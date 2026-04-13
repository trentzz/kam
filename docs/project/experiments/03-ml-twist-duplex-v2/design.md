# Experiment 03: twist-duplex-v2 — Real-data Retraining

**Status:** Planned
**Hypothesis:** A ML model trained on real confirmed TP/FP calls from the titration dataset, with new sequence-context features, will produce discriminating ml_prob distributions in discovery mode (unlike twist-duplex-v1 which assigns ≥ 0.9 to all PASS calls).
**Motivation:** `docs/benchmarking/07-snvindel-ml-boost-v1/discovery-precision-analysis.md`

---

## Problem Statement

`twist-duplex-v1` was trained on varforge-simulated data where all statistical PASS calls are real variants. At inference on real data, 100% of PASS calls — both true positives and false positives — receive ml_prob in [0.9, 1.0]. No threshold separates them. The model learned "passes the confidence ≥ 0.99 filter" as the positive class rather than "is actually somatic."

The titration dataset provides the first real labelled training data: ~200–250 confirmed FPs from 0% VAF negative controls, and ~700+ confirmed TPs from high-VAF spiked-in samples, all from real Twist duplex sequencing. The key new signal is that real background errors have systematic sequence signatures (deamination context, substitution spectrum) that synthetic data does not model.

---

## New Features

Two features are added to `extract_features()` in `kam-call/src/ml.rs`. Both are computable from existing `VariantCall` fields (`ref_sequence`, `alt_sequence`, `variant_type`) with no struct changes.

### `subst_type` (integer 0–12)

Encodes the substitution class for SNVs using 12 canonical classes, plus a non-SNV sentinel.

| Value | Class | Note |
|-------|-------|------|
| 0 | C>A | |
| 1 | C>G | |
| 2 | C>T | Common deamination artefact |
| 3 | T>A | |
| 4 | T>C | |
| 5 | T>G | |
| 6 | G>T | |
| 7 | G>C | |
| 8 | G>A | Complement of C>T |
| 9 | A>T | |
| 10 | A>G | |
| 11 | A>C | |
| 12 | non-SNV | Sentinel for indels, SVs |

Implementation: find the first differing position in `ref_sequence`/`alt_sequence`, look up the (ref_base, alt_base) pair. Non-SNV variants and any unrecognised pair return 12.

### `trinuc_context` (integer 0–64)

Encodes the 3-mer reference context centred on the SNV site. Using base-4 encoding (A=0, C=1, G=2, T=3): `context[i-1] * 16 + context[i] * 4 + context[i+1]`. Returns 64 for non-SNVs, edge-position SNVs, or unrecognised bases.

Rationale: C>T deamination preferentially occurs at CpG sites (5'-CG-3'). Encoding the trinucleotide lets the model learn that C>T at `_CG` context (trinuc values 9, 11, 13, 15 for `ACG`, `CCG`, `GCG`, `TCG`) are high-risk background errors compared to C>T in non-CpG contexts.

### Feature vector size

The full vector grows from 33 to 35 elements. Backward compatibility is preserved: `extract_features()` uses the `feature_names` list from the model JSON to build the lookup. Models that don't list `subst_type` or `trinuc_context` in `feature_names` are unaffected. Existing `twist-duplex-v1` and `single-strand-v1` continue to work with their 33-element vectors.

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

Hold out all 5ng samples (8 samples covering all VAF levels) as the test set. Train on 15ng (8 samples) and 30ng (8 samples). This ensures the test set is from a genuinely different DNA input condition — the model must generalise across input amounts, not just across VAF levels.

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
     d. Derive 33 base features (same as build_twist_duplex_features.py).
     e. Compute subst_type and trinuc_context from ref_seq/alt_seq/diff_pos.
  3. Emit train_features.csv.gz (15ng + 30ng samples).
  4. Emit test_features.csv.gz (5ng samples).
  5. Emit label_summary.csv (TP/FP counts per sample).

Output columns include: all 35 features, label, sample_id, ng_condition, vaf_nominal.
```

Important: do NOT include `ng_condition`, `vaf_nominal`, `sample_id` as model features. These are identifiers. Drop them from the feature matrix before training.

### `train_twist_duplex_v2.py` key differences from v1

- No `params.json` files — drop all `add_param_features()` logic.
- `GroupKFold(n_splits=3)` grouped by `sample_id`. With 16 training samples this gives ~5–6 samples per fold.
- 50 Optuna trials per model (small dataset, fast evaluation).
- Feature list: base 33 + `subst_type` + `trinuc_context` = 35 features. Treat both new features as numeric integers, NOT categorical.
- Export metadata JSON with `"version": "twist-duplex-2"` and 35-element `feature_names` list.
- Output ONNX to `bigdata/experiments/03-ml-twist-duplex-v2/models/twist-duplex-v2.onnx`.

---

## Rust Changes

All changes confined to two files.

### `kam-call/src/ml.rs`

Add two helper functions (gated by `#[cfg(feature = "ml")]`) after the existing `variant_type_str()` function:

```rust
fn snv_subst_type(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32
fn snv_trinuc_context(ref_seq: &[u8], alt_seq: &[u8], vt: VariantType) -> f32
```

Add two entries to the `lookup` HashMap in `extract_features()`:
```rust
("subst_type", snv_subst_type(&call.ref_sequence, &call.alt_sequence, call.variant_type)),
("trinuc_context", snv_trinuc_context(&call.ref_sequence, &call.alt_sequence, call.variant_type)),
```

Feature vector grows from 33 to 35 elements. Old models unaffected (metadata-driven lookup).

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
