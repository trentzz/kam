# Gradient Boosting for Variant Call Filtering in kam

## Overview

This document explores using gradient boosting (XGBoost and LightGBM) to improve
variant calling in the kam pipeline. The goal is to replace or augment the current
rule-based filters with a trained classifier that better separates true somatic
variants from noise, particularly at the 0.05–0.5% VAF range where the current
filters struggle most.

---

## 1. Where ML Fits in This Pipeline

### The current filter landscape

The pipeline emits calls with a `filter` column that takes values like `PASS`,
`LowConfidence`, `NotTargeted`, and `StrandBias`. These are rule-based decisions.
Looking at `snv_vaf0050_a/discovery.tsv`, the pattern is clear: any call with
`nalt == 1` gets `LowConfidence` and `conf ~= 0.859`. Any call with `nalt == 2`
gets `PASS` and `conf ~= 0.993`. The thresholds are fixed across all samples.

This works well at high VAF. It breaks down at low VAF (0.05–0.5%) because:

- True positives at 0.05% VAF legitimately have `nalt == 1`.
- Sequencing errors and UMI collisions also have `nalt == 1`.
- The current filter cannot distinguish these cases.

### Decisions a classifier could improve

**False positive reduction.** The largest gain. At low VAF, the output contains
many `LowConfidence` calls that are noise. A classifier trained on position,
allele context, duplex confirmation, strand balance, and coverage can score each
call and suppress the noise while preserving true variants.

**Threshold optimisation by sample.** The optimal `conf` cutoff and minimum `nalt`
threshold vary with total coverage (`nref + nalt`). A classifier naturally learns
this relationship from the data rather than requiring manual tuning per sample.

**Variant type disambiguation.** The discovery TSV mixes SNVs, short indels (MNVs
and multi-base substitutions with long ref/alt strings), and SV-adjacent events
like `<DUP>` in a single output. A multi-class classifier can label each call as
SNV, indel, SV, or artefact, enabling downstream tools to handle each type
correctly.

**VAF refinement at low depth.** The current `vaf_lo`/`vaf_hi` bounds from a
Beta-binomial model are wide at low counts. A classifier can use those bounds as
features alongside position-level context to produce a calibrated posterior
probability of the call being real, which is more useful clinically than a
binary PASS/FAIL.

### What ML should not replace

The `NotTargeted` filter in tumour-informed mode is deterministic and correct. ML
should not override it. Similarly, `StrandBias` failures backed by a significant
`sbp` value (e.g. `sbp < 0.05`) reflect a real data quality problem. ML can
learn to be cautious around strand-biased calls but should not override a hard
biological signal.

---

## 2. Feature Engineering

### Raw features available per call

From the TSV columns:

| Feature | Notes |
|---------|-------|
| `vaf` | Observed frequency |
| `vaf_lo`, `vaf_hi` | Beta-binomial CI bounds |
| `nref` | Reference molecule count |
| `nalt` | Alt molecule count |
| `ndupalt` | Duplex-confirmed alt count |
| `nsimalt` | Simplex-only alt count |
| `sbp` | Strand bias p-value |
| `conf` | Existing confidence score |

### Derived features to compute before training

**Coverage and depth:**
- `depth = nref + nalt` — total molecule depth at the site
- `log_depth = log2(depth)` — tree models split on log scale more naturally at this range

**Duplex evidence:**
- `duplex_frac = ndupalt / nalt` — fraction of alt molecules with duplex confirmation. Ranges 0–1. A true somatic variant at even 0.1% VAF in a well-covered sample should accumulate at least some duplex support. Pure noise rarely does.
- `has_duplex = (ndupalt > 0)` — binary flag, often more useful than the fraction when `nalt` is small

**Confidence interval:**
- `ci_width = vaf_hi - vaf_lo` — wide intervals indicate uncertain calls
- `ci_rel_width = ci_width / vaf` — normalised uncertainty

**Allele balance:**
- `alt_frac_of_depth = nalt / depth` — equivalent to `vaf` but keeps integer-denominator intuition
- `simplex_frac = nsimalt / nalt` — complement of `duplex_frac`

**Strand bias:**
- `sbp` directly (already a p-value, no transformation needed)
- `strand_biased = (sbp < 0.05)` — binary flag

**Variant type:**
- `ref_len = len(ref)` — 1 for SNV, >1 for indel or MNV
- `alt_len = len(alt)` — same
- `is_symbolic = (alt starts with "<")` — flags SV-type calls like `<DUP>`
- `indel_len = abs(len(alt) - len(ref))` — 0 for SNV/MNV, >0 for indel

**Position-level clustering:**
- `nearby_alt_count` — number of other alt calls within ±50 bp in the same sample. A cluster of low-confidence calls near a single high-confidence one is a common artefact pattern. This requires a pass over the sample's call set before training.

### Features not to include

Do not include `chrom` and `pos` as raw integers. These would cause the model to
memorise genomic coordinates from the simulation, which will not generalise to real
samples. The `nearby_alt_count` feature captures positional context without leaking
coordinates.

Do not include the existing `filter` label as a feature. It is the quantity you are
trying to improve on.

---

## 3. LightGBM vs XGBoost

### Dataset size

With 100 SNV/indel samples and 150 SV samples, each producing 20–200 calls, the
full labelled dataset is roughly 5,000–35,000 rows. This is small by ML standards.
Both XGBoost and LightGBM work at this scale without special treatment, but the
small size means overfitting is the primary risk.

### Class imbalance

True positives are rare. In `snv_vaf0050_a/discovery.tsv`, one call at position 149
is the planted true positive (VAF 0.046). The remaining 100+ calls are noise. The
positive:negative ratio is roughly 1:50 to 1:200, depending on VAF tier.

Both libraries handle imbalance via `scale_pos_weight` (XGBoost) or `is_unbalance`
/ `class_weight` (LightGBM). LightGBM also supports `pos_bagging_fraction` and
`neg_bagging_fraction` for more granular control.

### Speed

LightGBM trains faster than XGBoost on tabular data at this scale. At 30,000 rows
and ~15 features, both finish in under a second. Speed is not a differentiating
factor here.

### Interpretability

Both produce feature importance scores and support SHAP values. In a clinical
context, being able to say "this call was downscored because `ndupalt == 0` and
`nalt == 1`" is important. LightGBM integrates with `shap` as well as XGBoost
does.

### Recommendation for this project

Use LightGBM as the primary model. Reasons:

1. LightGBM's leaf-wise growth strategy finds deep, narrow splits efficiently.
   This suits the call set structure: most noise is distinguishable from true
   variants by a conjunction of multiple weak signals (low `nalt`, no duplex,
   wide CI, strand bias).
2. LightGBM's native categorical feature support handles the `filter` column (if
   used as a feature in future) and the `is_symbolic` flag without one-hot encoding.
3. LightGBM produces well-calibrated probabilities with isotonic regression or
   Platt scaling, which are preferable to binary PASS/FAIL in a clinical context.

XGBoost is a solid alternative and should be trained in parallel as a cross-check.
If the two models disagree on a call, flag it for manual review rather than
emitting a confident label.

---

## 4. Training Strategy

### Label construction

For each sample directory under `docs/benchmarking/snvindel/samples/` and
`docs/benchmarking/sv/samples/`, match calls from `discovery.tsv` against
`truth.tsv` by `(chrom, pos, ref, alt)` exact match. Assign label `1` to matched
calls and label `0` to all others.

For multi-allelic positions, a call is a true positive only if both `ref` and `alt`
match the truth entry exactly. Partial allele matches are false positives.

The SV samples include `<DUP>` and other symbolic alleles (visible in
`sv_vaf0100_a/discovery.tsv`). These should match on the symbolic allele string
exactly.

### Sample weighting

Weight samples by VAF tier. Calls from the 0.05% VAF samples are more clinically
important than calls from the 5% VAF samples. Assign per-sample weight inversely
proportional to the planted VAF. This prevents the model from optimising only for
easy high-VAF cases.

### Cross-validation strategy

With 250 samples, use group k-fold cross-validation with `k = 5`. Each fold holds
out 50 samples entirely. Never split a single sample across train and test folds.
This is critical: the noise calls within a sample are correlated by coverage, UMI
composition, and panel design. Splitting at the call level would leak information
and inflate performance estimates.

Use stratified grouping: ensure each fold contains a proportional mix of VAF tiers
and variant types.

### Evaluation metrics

Use area under the precision-recall curve (AUPRC) as the primary metric.
AUROC is misleading at high class imbalance. Report sensitivity and FDR at the
operating threshold that achieves 95% sensitivity on the validation fold. This
mirrors the clinical requirement.

---

## 5. Integration Plan

### CLI design

Add an optional `--ml-model <path>` flag to the `kam call` subcommand:

```
kam call --targets targets.bed --ml-model /path/to/kam_call_filter.lgbm [options]
```

When the flag is absent, the pipeline runs exactly as today. When present, each
call passes through the model after the rule-based filters. The model produces a
probability `ml_prob` (0–1). Calls with `ml_prob < threshold` receive an
additional `MLFiltered` filter tag. The threshold defaults to 0.5 and is
configurable with `--ml-threshold`.

The output TSV gains two columns: `ml_prob` and `ml_filter` (PASS or MLFiltered).

### Model distribution

The model file is a LightGBM binary model (`.lgbm` or `.bin`, produced by
`booster.save_model()`). It is not embedded in the Rust binary. It is distributed
as a separate file alongside the binary in the Nextflow container image, placed at
a standard path such as `/opt/kam/models/call_filter_v1.lgbm`.

The `--ml-model` flag accepts an absolute path. The Nextflow module passes the
container path automatically when ML filtering is enabled.

### Rust inference options

Three options exist for running the model from Rust:

**Option A: Python subprocess.** The Rust binary writes the call TSV to a temp
file, invokes a Python script (`kam-ml-filter.py`) as a child process, reads the
scored output. Simple to implement and debug. Adds a Python dependency to the
container. Acceptable for a pipeline tool where latency is not critical.

**Option B: `lightgbm` Rust crate.** The `lightgbm` crate (version 0.3+) wraps
the LightGBM C API. Load the model with `Booster::from_file()`, score each call
row with `predict()`. No Python required. The crate is maintained and builds
cleanly on Linux. This is the preferred option for production.

**Option C: ONNX export + `tract`.** Export the LightGBM model to ONNX format,
then use the `tract-onnx` crate for inference. ONNX export from LightGBM is
supported but adds a conversion step. Overkill for a gradient boosted tree.

The recommended path is to prototype with Option A (subprocess), validate the
model, then migrate to Option B (`lightgbm` crate) for the production release.

### Training pipeline (separate from Rust)

Training runs as a Python script outside the Rust workspace:

```
tools/ml/train_call_filter.py
  --snvindel-samples docs/benchmarking/snvindel/samples/
  --sv-samples       docs/benchmarking/sv/samples/
  --output           models/call_filter_v1.lgbm
```

The script reads all `discovery.tsv` and `truth.tsv` files, constructs labels,
computes derived features, runs group k-fold CV, trains the final model on all
data, and saves it. A companion script `tools/ml/evaluate_call_filter.py` produces
the AUPRC curve and per-VAF-tier sensitivity tables.

---

## 6. Risks and Limitations

### Overfitting on small datasets

With 250 samples and sparse positives (1–5 true calls per sample), the training
set contains 250–1,250 positive examples. This is small. LightGBM can memorise
these examples, particularly if any positives share unusual feature combinations
that are artefacts of the varforge simulation rather than real biology.

Mitigations:
- Limit tree depth (`max_depth = 4` or `num_leaves = 15`).
- Use early stopping on a held-out validation fold.
- Require minimum 20 samples per leaf.
- Inspect feature importances: if `pos` or any coordinate-derived feature ranks
  highly, the model is memorising simulation structure.

### Distributional shift

The benchmark data is simulated by varforge at fixed VAF tiers. Real ctDNA samples
have: variable coverage across the panel, non-uniform error profiles driven by
GC content and polymerase error, real UMI collisions at high depth, and multiple
true variants at variable VAFs. The model trained on simulated data may perform
differently on real data.

Mitigations:
- Include `depth` and `log_depth` as features so the model learns coverage
  dependence explicitly.
- When real data becomes available, retrain or fine-tune on a mix of simulated
  and real samples.
- Report confidence intervals on all benchmark metrics, not point estimates.

### Interpretability in clinical settings

A black-box ML filter on a diagnostic assay raises auditability questions. If a
true positive is suppressed by `MLFiltered`, a clinical scientist needs to
understand why.

Mitigations:
- Use SHAP values to generate a per-call explanation: "This call was filtered
  because duplex_frac = 0.0 and nalt = 1 and nearby_alt_count = 12."
- Emit the `ml_prob` score in the output TSV so lab staff can inspect borderline
  calls.
- Never fully suppress calls below threshold. Emit them with `MLFiltered` tag,
  allowing downstream review. The model is a ranking tool, not a hard gate.
- Maintain a changelog of model versions alongside the call filter binary. Each
  clinical result must record the model version used.

### Model version drift

As the pipeline evolves, new filters and features will be added. A model trained
on today's TSV schema will fail or behave unpredictably on a TSV with different
columns or altered `conf` semantics.

Mitigation: version the model file name (e.g. `call_filter_v1.lgbm`) and tie it
to a specific kam release version. The `--ml-model` flag should fail loudly if the
model was trained on a different schema version. Store the schema version and
training date in the model file metadata.

---

## Summary of Recommendations

1. Train a LightGBM binary classifier on all 250 benchmark samples using group
   k-fold CV, with `(chrom, pos, ref, alt)` matching against truth VCFs.
2. Use 15 features: the 9 raw TSV columns plus `duplex_frac`, `has_duplex`,
   `ci_width`, `ref_len`, `alt_len`, and `nearby_alt_count`.
3. Constrain tree depth (`max_depth = 4`, `num_leaves = 15`) and use early
   stopping to control overfitting.
4. Evaluate on AUPRC and per-VAF-tier sensitivity. The target is 95% sensitivity
   at the 0.1% VAF tier with a false discovery rate below 20%.
5. Integrate via `--ml-model` flag, using Python subprocess for the prototype and
   the `lightgbm` Rust crate for production.
6. Emit `ml_prob` and `ml_filter` in the output TSV. Never hard-suppress calls.
7. Retrain when real clinical data becomes available. Treat the simulated-data
   model as a baseline, not a production filter.
