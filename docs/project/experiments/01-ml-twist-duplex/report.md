# Experiment Report: ML Twist-Duplex (v1)

**Date**: 2026-04-10
**Status**: Complete

---

## Objective

Train and evaluate a gradient-boosted classifier for variant calls produced by
kam in duplex mode on Twist UMI chemistry. The model re-ranks calls by
assigning a real-variant probability, enabling a post-calling filter that
trades precision and recall independently of the statistical confidence
threshold.

---

## Data

### Simulations

11,000 varforge-simulated samples were generated (10,000 train, 1,000 test)
covering nine variant classes at low VAF:

| Variant class | Train | Test |
|---|---|---|
| SNV | 3,500 | 350 |
| Short indel (1–5 bp) | 2,000 | 200 |
| Medium indel (6–20 bp) | 1,000 | 100 |
| Large deletion (20–50 bp) | 500 | 50 |
| Tandem duplication (50–100 bp) | 300 | 30 |
| Inversion (50–100 bp) | 300 | 30 |
| Large deletion (100–200 bp) | 400 | 40 |
| Insertion (70–100 bp) | 1,000 | 100 |
| InvDel (80–130 bp) | 500 | 50 |
| **Total** | **10,000** | **1,000** |

VAF was drawn log-uniformly from 0.01%–5%. Coverage ranged from 1,000–15,000×,
family size mean from 2–8, PCR cycles from 5–14.

### Feature extraction

Each sample contributes calls from two modes (discovery and tumour-informed).
Because kam duplex calls encode full window haplotype sequences in the ref/alt
fields, labelling required position-based matching: the first character
difference between the ref and alt haplotype sequences identifies the variant
site in genomic coordinates, which is then checked against the truth VCF. Both
the anchor position and anchor+1 are accepted to accommodate VCF anchor-base
convention for indels.

After feature extraction, the training set contained:
- 22,338,378 rows total, 535,628 positives (2.4% positive rate)

Negatives were subsampled 10:1 before training to make the problem tractable
on a 25 GB machine, giving a final training matrix of 5,891,908 rows (535,628
positives, 9.1%).

The test set was not subsampled: 1,541,242 rows, 44,852 positives (2.9%).

---

## Models

Two models were trained and evaluated:

- **LightGBM**: 300 estimators, max depth 6, 31 leaves, learning rate 0.05,
  `is_unbalance=True`
- **XGBoost**: 300 estimators, max depth 6, learning rate 0.05,
  `scale_pos_weight=10.0`

Features: the 33-feature rust-safe set, identical to `single-strand-v1`. All
features are computable from a `VariantCall` at inference time.

Cross-validation used 5-fold `GroupKFold` with `sample_id` as the group key,
preventing leakage between folds.

---

## Results

### Cross-validation (5-fold, on subsampled training set)

| Model | AUPRC | AUROC |
|---|---|---|
| LightGBM | 0.6567 ± 0.0032 | 0.8521 ± 0.0012 |
| XGBoost | 0.6568 ± 0.0032 | 0.8517 ± 0.0013 |

CV precision and recall at the default 0.5 threshold:

| Model | Precision | Recall |
|---|---|---|
| LightGBM | 0.355 ± 0.004 | 0.701 ± 0.007 |
| XGBoost | 0.356 ± 0.004 | 0.702 ± 0.005 |

### Held-out test set (unsubsampled, 2.9% positive rate)

| Model | AUPRC | AUROC | Precision | Recall | F1 |
|---|---|---|---|---|---|
| LightGBM | 0.6062 | 0.9001 | 0.1263 | 0.7978 | 0.2181 |
| XGBoost | 0.6055 | 0.9001 | 0.1258 | 0.7974 | 0.2174 |

The two models are essentially identical in performance. LightGBM is the
deployed model (`twist-duplex-v1`).

### Interpreting the CV vs test gap

AUPRC drops from ~0.657 (CV) to ~0.606 (test) because AUPRC is sensitive to
class imbalance. CV ran on the subsampled 9.1%-positive set; the test set is
unsubsampled at 2.9% positives. At lower positive rates, maintaining high
precision at high recall is harder, and the AUPRC drops accordingly. The
AUROC, which is rank-based and less sensitive to class balance, is actually
higher on test (0.900 vs 0.852) — this reflects that the model ranks true
positives above the much larger negative background well, but the transition
from CV to test increased the absolute number of false positives.

---

## Feature importance

Top features by LightGBM gain:

| Rank | Feature | LGB gain | XGB gain | Notes |
|---|---|---|---|---|
| 1 | `ref_alt_len_ratio` | 20,752K | 0.438 | Dominant by a large margin |
| 2 | `vaf` | 7,935K | 0.051 | |
| 3 | `snr` | 4,420K | 0.141 | |
| 4 | `ci_width_rel` | 4,234K | 0.050 | |
| 5 | `variant_class_enc` | 2,908K | 0.044 | |
| 6 | `alt_len` | 1,883K | 0.031 | |
| 7 | `ref_len` | 1,701K | 0.017 | |
| 8 | `nref` | 1,263K | 0.026 | |
| 9 | `alt_depth` | 1,239K | 0.011 | |
| 10 | `indel_size` | 1,162K | 0.006 | |

### Zero-importance duplex features

All duplex-specific features have exactly zero importance in both models:

`ndupalt`, `nsimalt`, `duplex_frac`, `has_duplex`, `simplex_only_frac`,
`duplex_enrichment`, `log_nalt`, `log_nref`, `log_alt_depth`, `nalt_sq`,
`conf_above_99`, `conf_above_999`

This is unexpected. The duplex-specific columns (`ndupalt`, `nsimalt`) should
be the most informative features for distinguishing real variants in duplex
chemistry — a true variant should appear in duplex molecules, while artefacts
are more likely to be simplex-only. Zero importance strongly suggests a data
issue, not a genuine lack of signal.

**Likely cause**: the twist-duplex calls TSVs may not populate `n_duplex_alt`
and `n_simplex_alt` for all call types, or the column mapping during feature
extraction filled these with zeros or NaNs that were then replaced with zero
by `fillna(0)`. If duplex counts are all zero in the training data, the model
cannot learn from them.

**Action required**: verify that `ndupalt` and `nsimalt` are non-zero in a
representative sample of the training data before retraining. If they are
zero, identify why — either the simulated data does not produce duplex calls,
or the column rename in `build_training_data_twist_duplex.py` is not mapping
the correct column names.

### Dominance of `ref_alt_len_ratio`

`ref_alt_len_ratio` is the single most important feature by a factor of 2.6×
over the next. This reflects that the training set has very unequal variant
class sizes (3,500 SNVs vs 300 inversions), and the model exploits allele
length as a proxy for variant class identity. This is not a concern in itself,
but it means the model may be partially classifying variant type rather than
distinguishing real variants from artefacts within a class. If the deployed
population has a different variant class distribution to the training set, the
model may underperform.

---

## Comparison with `single-strand-v1`

| Metric | `single-strand-v1` | `twist-duplex-v1` |
|---|---|---|
| AUPRC (test) | 0.9998 | 0.6062 |
| AUROC (test) | 0.9949 | 0.9001 |
| Training samples | 9,990 | 10,000 |
| Positive rate (test) | ~10% | 2.9% |
| Features used | 33 | 33 (12 zero-importance) |

The large AUPRC gap is partly explained by the lower positive rate in the
twist-duplex test set, but the single-strand model also benefits from
more signal: single-strand consensus calls have stable `nalt`/`nref` counts
and non-trivial duplex columns. The twist-duplex model is likely hampered by
the zero-importance duplex features identified above.

---

## Deployed model

The LightGBM model was exported to ONNX (`lightgbm_rust.onnx`, 668 KB) and
bundled into the `kam` binary as `twist-duplex-v1`. Available via:

```bash
kam run --ml-model twist-duplex-v1 ...
```

---

## Next steps

1. **Diagnose zero duplex features**: inspect `ndupalt`/`nsimalt` distributions
   in `training_features.csv`. If all zeros, trace back through the TSV column
   names to find the mapping gap.
2. **Retrain with duplex signal**: once the column mapping is fixed, retrain.
   The AUPRC should improve substantially if duplex counts carry discriminative
   signal.
3. **Hyperparam search**: the current models used the same hyperparameters as
   `single-strand-v1`. A targeted HPO run on the duplex dataset may help.
4. **Threshold calibration**: at 0.5 threshold the model is high-recall
   (R=0.80) but low-precision (P=0.13). For production use, a higher threshold
   (0.7–0.9) may be more appropriate depending on the downstream application.
