# ML Experiment Results: Variant Ranking with Simulated Data

## Overview

This document summarises the first ML training run on kam's simulated benchmarking data. The goal is to learn a ranking function over variant calls that outperforms the raw PASS filter. All training used five-fold cross-validation on a combined SNV/indel dataset.

---

## Dataset Statistics

| Property | Value |
|---|---|
| Total rows | 240,134 |
| Positive labels (true variants) | 62 |
| Negative labels (false calls) | 240,072 |
| Positive rate | 0.026% |

The extreme imbalance is inherent to the task. Each sample contributes a small number of spiked-in variants (typically 1–5 per sample, across SNVs and indels) against a background of many thousands of spurious low-evidence calls. This mirrors real liquid biopsy conditions: somatic variants at 0.1–5% VAF are genuinely rare events in a high-depth panel context.

The SV samples contributed zero positives. The truth variants in those samples fell in size classes or at positions that were not recovered by the discovery or tumour-informed modes at the VAFs tested. This will be resolved by expanding the VAF sweep and including more SV types in future benchmarking runs.

---

## Model Comparison

The table below shows fold-mean metrics across five-fold CV.

| Model | AUROC | AUPRC | Recall | Precision | F1 |
|---|---|---|---|---|---|
| LightGBM | 0.7689 | 0.00096 | 0.780 | 0.00094 | 0.00188 |
| XGBoost | 0.8389 | 0.00209 | 0.230 | 0.00133 | 0.00261 |
| Baseline (PASS) | 0.5258 | 0.00037 | 0.075 | 0.00076 | 0.00149 |

Both ML models substantially outrank the PASS filter on AUROC. XGBoost achieves 0.84 vs 0.53 for the baseline. LightGBM recovers more true variants (recall 0.78 vs 0.08) but at lower precision. The baseline PASS filter performs near chance on AUROC because it assigns the same rank to all PASS calls and does not distinguish between high- and low-confidence events.

Figures: `docs/benchmarking/ml/results/figures/`.

---

## Why AUPRC, Not AUROC

AUROC measures ranking ability across all thresholds and all class ratios. At 0.026% positive rate, a model that perfectly ranks all negatives above all positives would score an AUROC close to 1.0 trivially. AUROC therefore overstates performance in highly imbalanced settings.

AUPRC (area under the precision-recall curve) directly measures the trade-off between finding true variants and avoiding false positives. At 0.026% positive rate, random chance gives an AUPRC of roughly 0.00026. The ML models reach AUPRC of 0.001–0.002, a three- to eight-fold improvement over chance and a two- to six-fold improvement over the PASS filter.

The absolute AUPRC values remain low. This reflects the dataset, not the models. With 62 positives spread across 240k rows and no variation in coverage or family size across samples, there is limited signal for the model to exploit. The ceiling will rise significantly with richer training data.

---

## Feature Importance

LightGBM gain importance (top features):

1. `alt_depth` — total depth at the variant site
2. `ci_width` — width of the VAF confidence interval
3. `nref` — reference allele count
4. `ref_len` — reference allele length (separates SNVs from indels/SVs)
5. `vaf` — variant allele fraction
6. `sbp` — strand bias p-value
7. `conf` — posterior confidence score

`alt_depth`, `ci_width`, and `nref` all capture depth and confidence in the call. The model learns that low-depth, wide-interval calls are more likely to be noise. `ref_len` separates variant classes without requiring an explicit class feature. `vaf`, `sbp`, and `conf` contribute expected biological and statistical signal.

The four duplex features (`ndupalt`, `nsimalt`, `has_duplex`, `duplex_frac`) have zero gain importance. This is expected: all simulated data uses a single molecule family size with no duplex reads, so these features carry no discriminating information. On real duplex sequencing data, where duplex confirmation is a strong true-variant signal, these features will become among the most important.

---

## What the Results Mean for the ML Integration Plan

The ML models demonstrate that it is possible to rank calls better than the PASS filter using the existing feature set, even on simulated data with limited diversity. The signal is real.

The limiting factor is data richness, not model quality. The current training set:
- Has uniform coverage across samples
- Has a single family size (no variation in sequencing depth or UMI grouping quality)
- Has no duplex reads (the strongest true-variant signal is absent)
- Has positives only in SNV/indel samples (SV positives are absent)

Models trained on this data will not generalise well to real samples.

---

## Recommended Next Steps

1. **Retrain on richer varforge data.** Run VAF sweeps at 0.5%, 1%, 2%, 5% with coverage variation (50x, 100x, 200x) and family size variation (2–8). This adds real variation across the depth-related features.

2. **Add SV positives.** Expand the SV benchmarking samples (INS, large DEL, INVDEL) so the model sees diverse variant classes with confirmed labels.

3. **Test on real samples.** Apply the model trained on simulated data to a real liquid biopsy run (known positive samples). Measure sensitivity and specificity against orthogonal validation (ddPCR, clinical report). The duplex features will activate here.

4. **Build training_data_v2.** The enhanced feature set in `scripts/ml/build_training_data_v2.py` adds 26 derived features including log transforms, ratios, squared terms, binned depth/strand-bias categories, and duplex enrichment scores. Retrain on this feature set.

5. **Calibrate outputs.** The model score should be calibrated to a probability before use as a filter or ranking signal. Use isotonic regression or Platt scaling post-hoc.
