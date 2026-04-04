# Experiment 02 — Single-Strand ML: v2b Model Results

## Note on naming

This model was trained on the v2 dataset (training_data_v2.csv, 888k rows, legacy
samples plus ML2 parameter-sweep samples) using the 33 rust-safe features. The
training script fell back to v2 data because training_data_v3.csv did not exist at
training time. The model was mislabelled "v3" at training time. It has been renamed
"v2b" to distinguish it from the existing v2 files and to reserve the v3 name for
future training on the ML3 train/test split.

## Objective

This experiment tests whether a gradient boosting classifier trained on simulated
single-strand sequencing data can distinguish true somatic variants from spurious
calls. The v2b model extends earlier versions by training on a larger and more diverse
dataset (888k rows, legacy + ML2 parameter-sweep samples) and evaluates two model
families (LightGBM and XGBoost) on a 33-feature set that is directly computable at
Rust inference time.

---

## Data

### Dataset composition

The v2 training dataset combines earlier generations of varforge simulations:

| Source | Approximate configs |
|---|---|
| Original SNV/indel sweep | 250 |
| ML1 dataset | 433 |
| ML2 dataset (indel_ml2) | 3,433 |

Each config produces one simulated sample. Across all training configs, the
combined CSV contains approximately 888k rows. Each row is one variant call from
`kam run`.

Variant types covered: SNV, Insertion, Deletion, MNV, Complex, LargeDeletion,
TandemDuplication, Inversion, Fusion, InvDel, NovelInsertion.

VAF range: 0.0003–0.15 (log-uniform sampling).

Coverage range: 500–20,000x (log-uniform sampling).

Family size mean: 1.5–10.0 (log-uniform).

### Positive/negative split and class imbalance

Positives are calls that match a spiked-in truth variant by (chrom, pos, ref, alt).
Each simulated sample contributes a small number of truth variants against a large
background of spurious low-evidence calls. Across the five CV folds, each fold
validation set contains roughly 5,000–5,200 positives out of approximately 970,000
rows — a positive rate of around 0.52%. Random-chance AUPRC is therefore
approximately 0.005.

### Train/test split

training_data_v2.csv has no explicit split column. All data is used as training.
There is no held-out test split for this model. All results below are from
5-fold cross-validation on the full dataset.

### Limitations

This model was trained on training_data_v2.csv. The training script (`train_eval_v3.py`)
looked for training_data_v3.csv first and silently fell back to v2 when it was absent.
That fallback has since been removed: the script now exits with an error if
training_data_v3.csv is missing.

All data is simulated. The models have not been validated against real liquid biopsy
sequencing.

---

## Model Training

### Models

Two gradient boosting frameworks were trained and compared: LightGBM and XGBoost.
Both used 300 estimators, max depth 6, and learning rate 0.05. LightGBM used
`is_unbalance=True` to handle class imbalance. XGBoost used
`scale_pos_weight = (1 - pos_rate) / pos_rate`, which has the same effect.

### Features

33 features were used, all directly computable from a `VariantCall` struct at Rust
inference time. They fall into five groups:

**Raw call fields (9):** `vaf`, `nref`, `nalt`, `ndupalt`, `nsimalt`, `sbp`, `conf`,
`ref_len`, `alt_len`

**Depth and duplex summaries (4):** `duplex_frac`, `has_duplex`, `ci_width`, `alt_depth`

**Log transforms (4):** `log_nalt`, `log_nref`, `log_alt_depth`, `log_vaf`

**Interaction and ratio terms (7):** `vaf_times_conf`, `vaf_times_nalt`,
`nalt_over_conf`, `ci_width_rel`, `snr`, `ref_alt_len_ratio`, `indel_size`

**Squared and binary features (9):** `conf_sq`, `nalt_sq`, `vaf_sq`,
`duplex_enrichment`, `simplex_only_frac`, `conf_above_99`, `conf_above_999`,
`sbp_above_05`, `variant_class_enc`

### Cross-validation

5-fold GroupKFold was used, with `sample_id` as the group key. GroupKFold ensures
that all rows from one simulated sample stay in the same fold. Without grouping,
rows from the same sample can appear in both train and validation sets for the same
fold, which inflates performance estimates. Using sample-level groups gives a more
realistic estimate of how the model generalises to unseen samples.

---

## Results

### CV metrics

All metrics are mean ± standard deviation across 5 folds. Precision and recall are
computed at a 0.5 probability threshold. AUPRC and AUROC are threshold-independent.

| Model | AUPRC | AUROC | Recall | Precision |
|---|---|---|---|---|
| LightGBM v2b | 0.919 ± 0.004 | 0.996 ± 0.000 | 0.947 ± 0.005 | 0.236 ± 0.022 |
| XGBoost v2b  | 0.920 ± 0.005 | 0.996 ± 0.000 | 0.946 ± 0.006 | 0.248 ± 0.030 |

The two models perform almost identically. AUPRC differs by less than 0.001.

These are CV results on the training set. They are not held-out test results.

### Feature importance

Top 10 features by LightGBM gain importance:

| Rank | Feature | LightGBM gain | XGBoost gain | Interpretation |
|---|---|---|---|---|
| 1 | `conf` | 363,387,000 | 0.0392 | Posterior confidence score from kam's binomial model |
| 2 | `indel_size` | 24,988,000 | 0.0037 | Absolute difference in allele length; separates indels from SNVs |
| 3 | `snr` | 9,572,000 | 0.0013 | Alt molecule count relative to ref + 1; signal-to-noise ratio |
| 4 | `vaf_times_conf` | 8,810,000 | 0.0026 | Joint signal: high VAF and high confidence together |
| 5 | `ref_len` | 7,800,000 | 0.0028 | Reference allele length; distinguishes SNVs (1) from longer alleles |
| 6 | `alt_len` | 3,635,000 | 0.0010 | Alternative allele length |
| 7 | `vaf_times_nalt` | 2,605,000 | 0.0007 | Joint signal: VAF weighted by absolute alt molecule count |
| 8 | `sbp` | 1,141,000 | 0.0002 | Strand bias p-value; low values flag strand-asymmetric calls |
| 9 | `ci_width_rel` | 658,000 | 0.0010 | Relative CI width (ci_width / vaf); normalised uncertainty |
| 10 | `ci_width` | 601,000 | 0.0003 | Absolute width of the VAF confidence interval |

`conf` dominates by a wide margin in LightGBM. This is expected: it already
integrates depth and allele count into a single posterior score. The model learns
that high-confidence calls are almost always true and low-confidence calls are almost
always noise. `indel_size`, `snr`, and `vaf_times_conf` add discriminating signal
beyond what `conf` alone provides.

The four duplex features (`ndupalt`, `has_duplex`, `duplex_frac`,
`duplex_enrichment`) have zero importance in both models. All simulated data uses
single-strand reads. None of the training samples contain duplex molecules, so these
features are identically zero or near-zero for all rows. The model cannot learn from
them. On real Twist duplex data, duplex confirmation is one of the strongest signals
for a true somatic variant. Duplex features are expected to become highly important
once the model is trained or fine-tuned on real data.

---

## Discussion

### What AUPRC 0.92 means

Random-chance AUPRC at the training set's positive rate (~0.52%) is approximately
0.005. An AUPRC of 0.92 represents a roughly 180-fold improvement over chance. The
v2 model on a much smaller dataset reached approximately 0.98 (on a dataset with
different class balance and fewer rows). The v2b result is harder to compare directly
because the class balance and dataset size differ.

AUPRC is the right primary metric here. AUROC at 0.996 sounds impressive, but at
low positive rates it is easier to achieve a high AUROC than a high AUPRC. AUPRC
penalises heavily for false positives, which is the main operational concern in
liquid biopsy: reporting a spurious somatic variant has clinical consequences.

### Why precision is low

At a 0.5 threshold, precision is approximately 24%. This means roughly three out of
four calls above threshold are still false positives. This is a consequence of the
class imbalance: even a well-calibrated model scoring 0.5 probability on a case
where the prior is 0.52% will produce many false positives. Operationally, the
threshold would not be set at 0.5. A higher threshold trades recall for precision.
The AUPRC curve characterises this trade-off fully.

Additionally, these are simulated data. The false positive rate in simulation may
differ from real data in ways that inflate or deflate the apparent precision.

### Duplex features and implications for real data

The zero importance of all four duplex features is not a bug or a training failure.
It is a direct consequence of training exclusively on single-strand simulations.
This means the deployed model cannot use duplex confirmation as a signal. On real
Twist duplex samples, where a call supported by duplex molecules is far more likely
to be a true variant, the model will underweight this information. The model must be
retrained or fine-tuned on real duplex data before it is suitable for production use
on Twist duplex chemistry.

### No held-out test evaluation

This model has no held-out test split. training_data_v2.csv does not contain a
`split` column. All published numbers are CV results on the full dataset. A proper
test evaluation requires training on the ML3 dataset (training_data_v3.csv), which
has an explicit 10,000/1,000 train/test split.

### Simulated data vs real sequencing

All training data is from varforge simulations. Real sequencing introduces noise
sources not present in simulation: systematic base-calling errors, polymerase
slippage artefacts, chimeric reads, index hopping, and chemistry-specific artefacts
from the Twist UMI design. The model will see real data as out-of-distribution. Its
practical precision and recall on clinical samples are unknown.

---

## Next Steps

1. **Train on ML3 data.** Generate training_data_v3.csv using the ML3 train/test
   pipeline. Run `train_eval_v3.py` to produce genuine v3 models with a held-out
   test evaluation. This is the minimum required before claiming a reliable
   performance estimate.

2. **Validate on real sequencing data.** Apply the ONNX model to a real Twist
   duplex liquid biopsy run with known positives (orthogonally validated by ddPCR or
   clinical report). Measure sensitivity and specificity against truth. This will
   reveal whether the model generalises and will activate the duplex features.

3. **Retrain with real duplex data.** Once real labelled data is available, include
   it in training (either as additional rows or as a fine-tuning dataset). Duplex
   features should gain substantial importance.

4. **Calibrate model outputs.** The raw probability scores are not calibrated.
   Apply isotonic regression or Platt scaling to convert scores to well-calibrated
   probabilities. This is needed before using the score as a VAF filter or
   confidence adjustment.

5. **Per-VAF and per-type breakdowns.** The current CV results aggregate all
   variant types and VAF levels. Stratified evaluation (by VAF bin and by variant
   class) would reveal where the model struggles, particularly at the lowest VAF
   levels (0.03%) where signal is weakest.
