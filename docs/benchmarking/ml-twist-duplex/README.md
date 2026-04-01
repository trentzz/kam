# Twist Duplex ML Dataset

This directory contains the configuration, truth VCFs, and outputs for the
Twist duplex ML training dataset. The dataset is designed to train a
chemistry-specific variant classifier for Twist Biosciences UMI duplex
sequencing.

All chemistry parameters are fixed (UMI length=5, duplex=true, inline=true).
Only biological and sequencing noise parameters vary.

## Dataset Statistics

| Variant type | Train | Test |
|---|---|---|
| snv | 3,500 | 350 |
| indel_short (1-5bp) | 2,000 | 200 |
| indel_medium (6-20bp) | 1,000 | 100 |
| sv_del (20-50bp) | 500 | 50 |
| sv_dup (50-100bp) | 300 | 30 |
| sv_inv (50-100bp) | 300 | 30 |
| sv_large_del (100-200bp) | 400 | 40 |
| ins (70-100bp) | 1,000 | 100 |
| invdel (80-130bp) | 500 | 50 |
| **Total** | **10,000** | **1,000** |

### VAF distribution

Log-uniform between 0.01% and 5%. Most samples fall at low VAF (<1%) where
classification is hardest.

### Parameter ranges

| Parameter | Distribution | Range |
|---|---|---|
| coverage | log-uniform | 1,000–15,000× |
| family_size_mean | uniform | 2.0–8.0 |
| pcr_cycles | uniform int | 5–14 |
| fragment_mean | uniform | 130–280 bp |
| fragment_sd | uniform | 20–50 bp |
| mean_quality | uniform int | 30–40 |

## Directory Structure

```
ml-twist-duplex/
├── configs/
│   ├── train/          # 10,000 varforge YAML configs
│   └── test/           # 1,000 varforge YAML configs
├── truth_vcfs/         # 11,000 per-sample truth VCFs
├── manifest.json       # Full parameter manifest for all 11,000 samples
├── simulations/        # (Nextcloud) varforge outputs + kam results
│   ├── train/
│   │   └── {sample_name}/
│   │       ├── *_R1.fastq.gz
│   │       ├── *_R2.fastq.gz
│   │       ├── *.truth.vcf
│   │       ├── discovery.tsv
│   │       ├── tumour_informed.tsv
│   │       └── params.json
│   └── test/
├── train_features.csv.gz   # (Nextcloud) extracted features, training set
├── test_features.csv.gz    # (Nextcloud) extracted features, test set
├── checkpoint.json         # (Nextcloud) pipeline checkpoint state
├── models/
│   ├── lightgbm_twist.txt  # (Nextcloud) trained LightGBM model
│   ├── xgboost_twist.json  # (Nextcloud) trained XGBoost model
│   ├── ensemble_twist.pkl  # (Nextcloud) calibrated ensemble
│   └── xgboost_twist.onnx  # (Nextcloud) ONNX export for deployment
└── results/
    ├── cv_results.csv
    ├── test_results.json
    ├── feature_importance_lgb.csv
    ├── feature_importance_xgb.csv
    └── figures/
        ├── pr_curves.png
        ├── roc_curves.png
        └── auprc_by_vaf.png
```

Items marked `(Nextcloud)` are not tracked in git. See [NEXTCLOUD.md](../../NEXTCLOUD.md).

## Reproducing the Dataset

### Prerequisites

- `varforge` on PATH
- `kam` built: `cargo build --release`
- Python packages: `numpy pandas lightgbm xgboost optuna scikit-learn`
- Optional: `skl2onnx` for ONNX export, `matplotlib` for figures

### Step 1: Generate configs and truth VCFs

```bash
python3 scripts/ml/generate_twist_duplex_configs.py
```

This writes 11,000 YAML configs to `docs/benchmarking/ml-twist-duplex/configs/`
and 11,000 truth VCFs to `docs/benchmarking/ml-twist-duplex/truth_vcfs/`.
Also writes `manifest.json`. Runtime: ~2 minutes.

### Step 2: Run the pipeline (varforge + kam)

```bash
python3 scripts/ml/run_twist_duplex_pipeline.py --workers 12
```

This runs varforge and kam for all 11,000 samples. Checkpointing allows
resuming if interrupted. Progress uploads to Nextcloud after each batch
of 500 samples.

Estimated runtime on a 12-core machine: 24–48 hours for all 11,000 samples,
depending on coverage settings.

To run only the training split:
```bash
python3 scripts/ml/run_twist_duplex_pipeline.py --split train --workers 12
```

### Step 3: Extract features

```bash
python3 scripts/ml/build_twist_duplex_features.py
```

Outputs `train_features.csv.gz` and `test_features.csv.gz`.

### Step 4: Train models

```bash
python3 scripts/ml/train_twist_duplex.py --hpo-trials 200
```

Runs LightGBM and XGBoost HPO (200 trials each), 5-fold CV, final training,
calibration, and ONNX export. Results go to `results/` and `models/`.

### Step 5: Upload to Nextcloud

```bash
bash scripts/ml/upload_twist_duplex.sh all
```

## Feature List

The feature matrix (~55 features) is built by `build_twist_duplex_features.py`.

### Raw features (from kam TSV)

| Feature | Description |
|---|---|
| `vaf` | Estimated variant allele frequency |
| `vaf_ci_low` | Lower 95% CI bound on VAF |
| `vaf_ci_high` | Upper 95% CI bound on VAF |
| `n_molecules_ref` | Molecules supporting reference |
| `n_molecules_alt` | Molecules supporting alt |
| `n_duplex_alt` | Duplex molecules supporting alt |
| `n_simplex_alt` | Simplex molecules supporting alt |
| `strand_bias_p` | Strand bias p-value |
| `confidence` | Posterior variant confidence |
| `n_simplex_fwd_alt` | Simplex alt on forward strand |
| `n_simplex_rev_alt` | Simplex alt on reverse strand |
| `n_duplex_ref` | Duplex molecules supporting ref |
| `n_simplex_ref` | Simplex molecules supporting ref |
| `mean_alt_error_prob` | Mean base error probability at alt position |
| `min_variant_specific_duplex` | Minimum variant-specific duplex count |
| `mean_variant_specific_molecules` | Mean variant-specific molecule count |
| `ref_len` | Length of REF allele |
| `alt_len` | Length of ALT allele |

### Duplex-specific derived features

| Feature | Formula | Description |
|---|---|---|
| `duplex_vaf` | `n_duplex_alt / (n_duplex_alt + n_duplex_ref)` | Duplex-only VAF estimate |
| `simplex_vaf` | `n_simplex_alt / (n_simplex_alt + n_simplex_ref)` | Simplex-only VAF estimate |
| `duplex_simplex_vaf_delta` | `duplex_vaf - simplex_vaf` | Duplex-simplex disagreement |
| `vaf_duplex_agreement` | `|vaf - duplex_vaf|` | Calibration quality |
| `duplex_enrichment` | `n_duplex_alt / n_molecules_alt` | Fraction of alt that is duplex |
| `simplex_frac` | `n_simplex_alt / n_molecules_alt` | Fraction of alt that is simplex |
| `strand_balance` | `min(fwd, rev) / (fwd + rev)` | Forward-reverse balance |
| `strand_asymmetry` | `(fwd - rev) / (fwd + rev)` | Signed strand imbalance |
| `duplex_depth` | `n_duplex_alt + n_duplex_ref` | Total duplex coverage |
| `simplex_depth` | `n_simplex_alt + n_simplex_ref` | Total simplex coverage |
| `duplex_ref_frac` | `n_duplex_ref / (n_duplex_ref + n_duplex_alt)` | Ref fraction in duplex |
| `simplex_only_frac` | `n_simplex_alt / n_molecules_alt` | Simplex-only alt fraction |
| `variant_specific_duplex_frac` | `min_vs_dup / n_duplex_alt` | Variant-specific duplex purity |
| `qual_confidence_interaction` | `(1 - mean_err) * confidence` | Quality-confidence product |
| `qual_vaf_interaction` | `(1 - mean_err) * vaf` | Quality-VAF product |
| `has_duplex` | `n_duplex_alt > 0` | Binary duplex presence flag |

### Log transforms

`log_vaf`, `log_nalt`, `log_nref`, `log_nduplex_alt`, `log_nsimplex_alt`,
`log_alt_depth`, `log_duplex_depth`, `log_simplex_depth`

### Ratio and interaction terms

`ci_width`, `ci_width_rel`, `vaf_times_conf`, `vaf_times_nalt`,
`nalt_over_conf`, `conf_sq`, `nalt_sq`, `vaf_sq`, `duplex_frac`, `snr`,
`ref_alt_len_ratio`, `indel_size`

### Binned features

| Feature | Values | Description |
|---|---|---|
| `depth_bucket` | 0–4 | Quintile bin of total depth |
| `sbp_cat` | low/medium/high | Strand bias p-value bin |
| `conf_above_99` | 0/1 | Confidence ≥ 0.99 |
| `conf_above_999` | 0/1 | Confidence ≥ 0.999 |
| `sbp_above_05` | 0/1 | Strand bias p ≥ 0.05 |

### Parameter features (from params.json)

`coverage`, `family_size_mean`, `pcr_cycles`, `fragment_mean`,
`vaf_target`, `vaf_log_target`

### Variant classification

`variant_class`: SNV / short_indel / medium_indel / SV (derived from allele lengths)

## Model Training Methodology

1. Load `train_features.csv.gz` and `test_features.csv.gz`.
2. Drop rows where `filter = 'NotTargeted'`.
3. Label: 1 if the call (ref_seq, alt_seq) matches a truth variant, 0 otherwise.
4. Optuna HPO: 200 trials of 3-fold GroupKFold CV per model, optimising AUPRC.
5. Final training on full training set with best HPO parameters.
6. Calibration with `CalibratedClassifierCV(method='isotonic', cv='prefit')`.
7. Ensemble: `0.5 * P(lgb) + 0.5 * P(xgb)`.
8. Evaluation: AUROC, AUPRC, sensitivity at 1% and 5% FPR, best F1.
9. Breakdown by VAF bin and variant type.

## Nextcloud Upload Locations

| Content | Nextcloud path |
|---|---|
| Train batches | `benchmarking/ml-twist-duplex/train/batch_NNN.tar.gz` |
| Test batches | `benchmarking/ml-twist-duplex/test/batch_NNN.tar.gz` |
| Models | `benchmarking/ml-twist-duplex/models/` |
| Feature CSVs | `benchmarking/ml-twist-duplex/` |

## Expected Runtime Estimates

| Step | Estimate |
|---|---|
| Config + VCF generation (11k samples) | ~2 minutes |
| varforge simulation (11k samples, 12 workers) | 8–16 hours |
| kam discovery + TI (11k samples, 12 workers) | 8–16 hours |
| Feature extraction | ~10 minutes |
| Model training (200 HPO trials each) | 4–8 hours |

Total wall clock on a 12-core machine: approximately 1–2 days.
