# Experiment 03: twist-duplex-v2 — Real-data Retraining

**Status:** Complete
**Hypothesis:** A ML model trained on real confirmed TP/FP calls from the titration dataset, with new sequence-context features, will produce discriminating ml_prob distributions in discovery mode (unlike twist-duplex-v1 which assigns ≥ 0.9 to all PASS calls).
**Outcome:** Confirmed. See `results/README.md`. AUPRC 0.973, FP reduction 79–93%, precision 0.969–1.000 at 1–2% VAF, sensitivity retention 97–100%.
**Motivation:** `docs/benchmarking/07-snvindel-ml-boost-v1/discovery-precision-analysis.md`

---

## Problem Statement

`twist-duplex-v1` was trained on varforge-simulated data where all statistical PASS calls are real variants. At inference on real data, 100% of PASS calls — both true positives and false positives — receive ml_prob in [0.9, 1.0]. No threshold separates them. The model learned "passes the confidence ≥ 0.99 filter" as the positive class rather than "is actually somatic."

The titration dataset provides the first real labelled training data: ~200–250 confirmed FPs from 0% VAF negative controls, and ~700+ confirmed TPs from high-VAF spiked-in samples, all from real Twist duplex sequencing. The key new signal is that real background errors have systematic sequence signatures (deamination context, substitution spectrum) that synthetic data does not model.

---

## New Features

Sixteen features are added to `extract_features()` in `kam-call/src/ml.rs`, growing the vector from 33 to 49. Backward compatibility is preserved: `extract_features()` uses the `feature_names` list from the model JSON to build the lookup. Old models (`twist-duplex-v1`, `single-strand-v1`) list only their original 33 features and are unaffected.

### Category A — Sequence context (5 features, computed from `ref_sequence`/`alt_sequence`)

| Feature | Range | Rationale |
|---------|-------|-----------|
| `subst_type` | 0–12 | 12 canonical SNV classes (C>A=0 … A>C=11) + non-SNV sentinel 12. Captures damage signature. |
| `trinuc_context` | 0–64 | Base-4 encoded 3-mer at variant site: `ref[i-1]*16 + ref[i]*4 + ref[i+1]`. Sentinel 64 for non-SNV/edge. |
| `is_cpg` | 0/1 | 1 if C>T or G>A change at a CpG dinucleotide. Key deamination indicator. |
| `gc_content_ref` | 0–1 | GC fraction of `ref_sequence`. GC-rich regions have elevated oxidative damage rates. |
| `homopolymer_run` | 0+ | Longest run of identical bases adjacent to the variant in `ref_sequence`. PCR slippage indicator. |

### Category B — Existing VariantCall fields not previously exposed (7 features)

| Feature | Field | Notes |
|---------|-------|-------|
| `n_simplex_fwd_alt` | `call.n_simplex_fwd_alt` | Simplex forward-strand alt molecules |
| `n_simplex_rev_alt` | `call.n_simplex_rev_alt` | Simplex reverse-strand alt molecules |
| `n_duplex_ref` | `call.n_duplex_ref` | Duplex-confirmed reference molecules |
| `n_simplex_ref` | `call.n_simplex_ref` | Simplex-only reference molecules |
| `mean_alt_error_prob` | `call.mean_alt_error_prob` | Mean consensus error probability of alt-supporting reads |
| `min_variant_specific_duplex` | `call.min_variant_specific_duplex` | Minimum duplex count per variant-specific k-mer |
| `mean_variant_specific_molecules` | `call.mean_variant_specific_molecules` | Mean molecule count per variant-specific k-mer |

### Category C — Derived strand/duplex features (4 features)

| Feature | Formula | Rationale |
|---------|---------|-----------|
| `strand_asymmetry_alt` | `(fwd - rev) / (fwd + rev + ε)` | Alt-specific strand imbalance. Damage artefacts are often single-strand. |
| `duplex_vaf` | `ndupalt / (ndupalt + n_duplex_ref + ε)` | VAF estimated from duplex molecules only. |
| `simplex_vaf` | `nsimalt / (nsimalt + n_simplex_ref + ε)` | VAF estimated from simplex molecules only. |
| `duplex_simplex_vaf_delta` | `duplex_vaf - simplex_vaf` | Real variants: delta ≈ 0. Damage artefacts: large positive or negative delta. |

### Feature vector size

Total: 33 + 5 + 7 + 4 = **49 features**.

---

## Training Data Construction

### Source

`docs/benchmarking/07-snvindel-ml-boost-v1/results/tsvs_disc_ml/` — 24 per-sample TSVs from the discovery + ML run. Each file contains all PASS calls for one sample, including `ml_prob` and `ml_filter` columns, and the `ref_seq`/`alt_seq` columns needed for the new features.

### Labelling

| Sample type | Label assignment |
|-------------|-----------------|
| VAF = 0% (3 samples) | All PASS calls → FP (label = 0). Ground truth: zero spiked-in variants. |
| VAF > 0% (21 samples) | PASS call matches truth VCF → TP (label = 1). Otherwise → FP (label = 0). |

For VAF = 0% samples, labels are assigned unconditionally without coordinate matching. This avoids false TPs from incidental background errors that happen to match panel coordinates.

Truth matching uses the coordinate-recovery algorithm from `run_titration_batch.py` lines 110–190 (the diff_pos → VCF-position pipeline).

### Expected class distribution

Based on the precision study results:

| Source | Approx TPs | Approx FPs |
|--------|-----------|-----------|
| 0% VAF negative controls (3 samples) | 0 | ~200–250 |
| Positive VAF samples (21 samples) | ~1,800 | ~1,300 |
| **Total** | **~1,800** | **~1,500** |

Roughly balanced. Do not apply downsampling or class weights initially — verify actual counts from the label summary output before deciding.

### Train/test split

Hold out all 30ng samples (8 samples covering all VAF levels) as the test set. Train on 5ng (8 samples) and 15ng (8 samples). This ensures the model sees both the weakest and middle input conditions during training, and the test set evaluates generalisation to higher molecule counts and duplex rates. Holding out 30ng avoids the original concern of testing on the weakest condition (5ng) when the model has only seen stronger inputs.

---

## Scripts

| Script | Location | Purpose |
|--------|----------|---------|
| `build_real_training_data.py` | `scripts/ml/` | Extract features + labels from per-sample TSVs |
| `train_twist_duplex_v2.py` | `scripts/ml/` | HPO + training + calibration + ONNX export |
| `export_model_twist_duplex_v2.py` | `scripts/ml/` | Re-export from saved model (if needed) |

### `build_real_training_data.py` design

```
Inputs:
  --tsv-dir        tsvs_disc_ml/ directory (24 TSVs)
  --truth-vcf      docs/benchmarking/01-snvindel/scripts/truth_variants.vcf
  --output-dir     bigdata/experiments/03-ml-twist-duplex-v2/

Steps:
  1. Load truth VCF to (chrom, pos, ref, alt) set.
  2. For each TSV:
     a. Parse sample name → ng, vaf_nominal.
     b. For each PASS row: recover genomic key via diff_pos algorithm.
     c. Assign label: 1 if key in truth, 0 otherwise.
        For vaf_nominal == 0.0: force label = 0 unconditionally.
     d. Derive all 49 features (33 base + 7 category B + 4 category C + 5 category A).
  3. Emit train_features.csv.gz (5ng + 15ng samples).
  4. Emit test_features.csv.gz (30ng samples).
  5. Emit label_summary.csv (TP/FP counts per sample).

Output columns include: all 49 features, label, sample_id, ng_condition, vaf_nominal.
```

Important: do NOT include `ng_condition`, `vaf_nominal`, `sample_id` as model features. These are identifiers. Drop them from the feature matrix before training.

### `train_twist_duplex_v2.py` key differences from v1

- No `params.json` files — drop all `add_param_features()` logic.
- `GroupKFold(n_splits=4)` grouped by `sample_id`. With 16 training samples this gives 4 samples per fold.
- 50 Optuna trials per model (small dataset, fast evaluation).
- Feature list: 49 features (33 base + 16 new). Treat all new features as numeric, NOT categorical.
- Export metadata JSON with `"version": "twist-duplex-2"` and 49-element `feature_names` list.
- Output ONNX to `bigdata/experiments/03-ml-twist-duplex-v2/models/twist-duplex-v2.onnx`.

---

## Rust Changes

All changes confined to two files.

### `kam-call/src/ml.rs`

Five helper functions added after `variant_type_str()`, all gated by `#[cfg(feature = "ml")]`:

```rust
fn snv_subst_type(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32
fn snv_trinuc_context(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32
fn snv_is_cpg(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32
fn ref_gc_content(ref_seq: &[u8]) -> f32
fn ref_homopolymer_run(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32
```

16 new entries added to the `lookup` HashMap in `extract_features()` (categories A, B, C as above).

Feature vector grows from 33 to 49 elements. Old models unaffected (metadata-driven lookup).

### `kam/src/models.rs`

Add two static byte slices for the bundled v2 model, extend `builtin_names()` to include `"twist-duplex-v2"`, and add a match arm in `resolve()`. Total: ~6 lines added.

### No other Rust files need changes.

---

## Evaluation

### Minimum bar

AUPRC > 0.5 on the held-out 5ng test set. A model at AUPRC = 0.5 is no better than random — worse than v1 only in the sense that v1 at least does not reduce sensitivity. The model must improve over 0.5 to justify deployment.

### End-to-end benchmark

Re-run the discovery precision study (condition B: `--ml-model twist-duplex-v2`) using `run_titration_batch.py`. Compare against the v1 results:

| Metric | v1 target | v2 target |
|--------|-----------|-----------|
| FP calls with ML_PASS | 100% | < 50% |
| TP calls with ML_PASS | 100% | > 85% |
| Precision at 0.5 threshold (15ng 1% VAF) | 0.755 | > 0.85 |
| Sensitivity at 0.5 threshold | same as no-model | > 0.90 × no-model |

Also produce the TP vs FP ml_prob histogram from `tsvs_disc_ml/` (following the same analysis as `discovery-precision-analysis.md`). The v2 histogram must show separated distributions.

### Threshold sweep

Report precision, recall, and FP count at thresholds: 0.3, 0.5, 0.7, 0.9. Set `ml_pass_threshold` in the JSON metadata to the value that maximises F1 on the 5ng test set.

---

## Known Limitations

**Small training dataset.** 16 training samples, ~3,300 rows. CV variance will be high. Report confidence intervals across folds. The model may not generalise well to substantially different libraries or panel designs.

**Label noise in positive-VAF samples.** PASS calls that don't match the 375-variant truth panel are labelled FP, but some may be real variants outside the panel. This is unavoidable without deeper characterisation. Negative control samples provide clean labels and should be considered the most reliable FP examples.

**Indels.** The `trinuc_context` and `subst_type` features return sentinel values for all non-SNVs. The training data may have few indel TPs (the titration panel has 170 truth indels but many are missed at 2% VAF). The model's ability to discriminate indel artefacts is limited.

**Distribution shift.** The training data comes from one panel (375 loci, 100 bp targets) and one chemistry condition (Twist duplex, 2M reads). The model may not transfer to different panels with different background error profiles.

---

## Bigdata Layout

```
bigdata/experiments/03-ml-twist-duplex-v2/
  train_features.csv.gz
  test_features.csv.gz
  label_summary.csv
  optuna_study_lgb.pkl
  optuna_study_xgb.pkl
  models/
    twist-duplex-v2.onnx
    twist-duplex-v2.json
  results/
    cv_results_v2.csv
    test_results_v2.json
    feature_importance_lgb_v2.csv
    feature_importance_xgb_v2.csv
```

Trained artifacts go to `bigdata/` (gitignored), then copied to `kam/models/` and committed.

---

## Steps to Execute

1. Add `snv_subst_type()` and `snv_trinuc_context()` to `kam-call/src/ml.rs`. Run `cargo test -p kam-call`.
2. Write `scripts/ml/build_real_training_data.py`. Run it, inspect `label_summary.csv`.
3. Write `scripts/ml/train_twist_duplex_v2.py`. Train model, check AUPRC on test set.
4. Copy ONNX + JSON to `kam/models/twist-duplex-v2.{onnx,json}`.
5. Update `kam/src/models.rs`. Build with `--features ml`.
6. Re-run condition B benchmark with `--ml-model twist-duplex-v2`. Compare against v1 results.
7. Write results in `docs/project/experiments/03-ml-twist-duplex-v2/results/`.
8. Commit model artifacts and Rust changes.
