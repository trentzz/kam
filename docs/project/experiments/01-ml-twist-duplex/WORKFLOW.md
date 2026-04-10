# Twist Duplex ML Workflow

## 1. Overview

This workflow trains a Twist duplex-specific variant calling classifier for kam. It produces two gradient-boosted models (LightGBM and XGBoost), a calibrated ensemble, and an ONNX export for deployment inside the kam binary.

The dataset fixes all Twist chemistry parameters (UMI length 5, duplex, inline) and varies only biological and sequencing noise parameters. This lets the models focus entirely on variant-signal structure and exploit duplex-specific features — duplex/simplex VAF agreement, strand balance, variant-specific duplex counts — that are meaningless for non-duplex chemistries.

```
generate_twist_duplex_configs.py
  IN:  benchmark VCF pools (docs/benchmarking/)
  OUT: configs/ (11,000 YAML), truth_vcfs/ (11,000 VCFs), manifest.json  [docs/, tracked in git]
       |
       v
run_twist_duplex_pipeline.py
  IN:  configs/, manifest.json
  OUT: bigdata/experiments/01-ml-twist-duplex/simulations/{train,test}/{sample}/
         *_R1.fastq.gz, *_R2.fastq.gz          (NOT uploaded)
         calls_discovery.tsv, calls_discovery.vcf
         calls_tumour_informed.tsv, calls_tumour_informed.vcf
         params.json
       checkpoint.json                          [docs/, tracked in git]
       batch tarballs uploaded to Nextcloud     (TSVs, VCFs, params.json only)
       |
       v
build_twist_duplex_features.py
  IN:  bigdata/simulations/, truth_vcfs/
  OUT: bigdata/experiments/01-ml-twist-duplex/train_features.csv.gz
       bigdata/experiments/01-ml-twist-duplex/test_features.csv.gz
       |
       v
train_twist_duplex.py
  IN:  train_features.csv.gz
  OUT: bigdata/.../models/lightgbm_twist.txt
       bigdata/.../models/xgboost_twist.json
       bigdata/.../models/ensemble_twist.pkl
       bigdata/.../models/xgboost_twist.onnx
       docs/.../results/cv_results.csv          [tracked in git]
       docs/.../results/test_results.json       [tracked in git]
       docs/.../results/feature_importance_*.csv [tracked in git]
       |
       v
evaluate_twist_duplex.py
  IN:  test_features.csv.gz, models/
  OUT: docs/.../results/evaluation_report.json  [tracked in git]
       docs/.../results/{roc,pr,confusion,vaf_bin,vtype}_*.csv [tracked in git]
       bigdata/.../results/figures/             (large plots, not tracked)
       |
       v
upload_twist_duplex.sh
  Pushes models, feature CSVs, and configs tarball to Nextcloud.
```

---

## 2. Prerequisites

**Binaries**

- `varforge` on `PATH`. This is a separate tool; it is not built by this repo.
- `kam` at `target/release/kam`. Build it with:
  ```bash
  cargo build --release
  ```

**Python packages**

```
lightgbm
xgboost
scikit-learn
pandas
numpy
optuna
skl2onnx          # for ONNX export; optional but recommended
onnxmltools       # required by skl2onnx for XGBoost conversion
matplotlib        # for figures; optional
```

Install with:
```bash
pip install lightgbm xgboost scikit-learn pandas numpy optuna skl2onnx onnxmltools matplotlib
```

**Environment**

- `.env` file at the repository root with `NEXTCLOUD_EDIT_TOKEN` set. This is only needed for upload steps. See `docs/project/devmanual/nextcloud.md`.

**Hardware**

The full 11,000-sample run is resource-intensive. 16 cores and ~64 GB RAM are recommended. On a 12-core machine, expect 1–2 days of wall-clock time. See runtime estimates at the end of each step.

---

## 3. Step 1: Generate varforge configs

Run from the repository root.

```bash
python3 scripts/ml/generate_twist_duplex_configs.py
```

This script samples parameters for 11,000 unique simulations (10,000 train + 1,000 test). For each sample it:

1. Selects a truth VCF from the benchmark VCF pools under `docs/benchmarking/`.
2. Draws random simulation parameters (coverage, family size, PCR cycles, fragment size, mean quality, VAF) from the distributions defined in `design.md`.
3. Writes a varforge YAML config.

**Outputs**

| Path | Description |
|---|---|
| `docs/project/experiments/01-ml-twist-duplex/configs/train/` | 10,000 YAML configs |
| `docs/project/experiments/01-ml-twist-duplex/configs/test/` | 1,000 YAML configs |
| `docs/project/experiments/01-ml-twist-duplex/manifest.json` | Full parameter manifest for all 11,000 samples |

These outputs are tracked in git. Run this step only once. If you need to regenerate, delete the existing configs and re-run; the script is deterministic given the same random seed.

**Estimated runtime:** ~2 minutes.

---

## 4. Step 2: Run the simulation pipeline

Run from the repository root.

```bash
python3 scripts/ml/run_twist_duplex_pipeline.py --split all --workers 16 --batch-size 500
```

For each sample the pipeline runs:

1. `varforge simulate` using the sample's YAML config. Produces paired FASTQs and a truth VCF.
2. `kam discovery` on the FASTQs. Produces `calls_discovery.tsv` and `calls_discovery.vcf`.
3. `kam tumour-informed` on the FASTQs with the truth VCF as the target list. Produces `calls_tumour_informed.tsv` and `calls_tumour_informed.vcf`.
4. Writes `params.json` to the sample directory.

After each batch of 500 samples, the script tars the results (excluding FASTQs) and uploads to Nextcloud.

**Options**

| Flag | Default | Description |
|---|---|---|
| `--split` | `all` | Which split to run: `train`, `test`, or `all` |
| `--workers` | 4 | Number of parallel workers |
| `--batch-size` | 500 | Samples per upload batch |
| `--dry-run` | off | Print commands without executing |
| `--skip-upload` | off | Skip Nextcloud upload after each batch |

**Outputs**

```
bigdata/experiments/01-ml-twist-duplex/simulations/{train,test}/{sample_name}/
  {sample_name}_R1.fastq.gz
  {sample_name}_R2.fastq.gz
  {sample_name}.truth.vcf.gz
  calls_discovery.tsv
  calls_discovery.vcf
  calls_tumour_informed.tsv
  calls_tumour_informed.vcf
  params.json
```

FASTQs are not uploaded to Nextcloud (too large). Only TSVs, VCFs, and `params.json` are included in the tarballs.

**Checkpointing**

`docs/project/experiments/01-ml-twist-duplex/checkpoint.json` records completed sample names. If the run is interrupted, re-run the same command and it will skip completed samples automatically. See Section 10 for how to force a re-run of specific samples.

**Estimated runtime:** 8–16 hours for varforge simulation + 8–16 hours for kam, with 12 workers. Total: 24–48 hours for all 11,000 samples, depending on coverage settings. Use `--workers 16` on a 16-core machine to reduce this.

---

## 5. Step 3: Extract features

Run from the repository root.

```bash
python3 scripts/ml/build_twist_duplex_features.py --mode tumour_informed
```

The script reads every sample directory under `bigdata/.../simulations/`, loads the kam TSV outputs, labels each called variant against the truth VCF, and builds the feature matrix.

**Mode options**

| Mode | Description |
|---|---|
| `discovery` | Features from `calls_discovery.tsv` only |
| `tumour_informed` | Features from `calls_tumour_informed.tsv` only |
| `both` | Both modes; produces separate feature files for each |

Use `tumour_informed` for the standard classifier. Use `both` if you want to compare modes.

**Features extracted (~55 total)**

The feature matrix includes raw TSV columns (VAF, CI bounds, molecule counts, duplex/simplex counts, strand bias, confidence), duplex-specific derived features (duplex VAF, simplex VAF, duplex enrichment, strand balance), log transforms, ratio and interaction terms, binned features, and per-sample simulation parameters from `params.json`. See `README.md` for the full feature list.

**Outputs**

| Path | Description |
|---|---|
| `bigdata/experiments/01-ml-twist-duplex/train_features.csv.gz` | Training feature matrix |
| `bigdata/experiments/01-ml-twist-duplex/test_features.csv.gz` | Test feature matrix |

**Estimated runtime:** ~10 minutes.

---

## 6. Step 4: Train models

Run from the repository root.

```bash
python3 scripts/ml/train_twist_duplex.py --hpo-trials 200
```

Training proceeds in order:

1. Loads `train_features.csv.gz`. Drops rows where `filter = 'NotTargeted'`. Labels each row: 1 if the call matches a truth variant, 0 otherwise.
2. Runs Optuna HPO for LightGBM (200 trials, 3-fold GroupKFold CV, optimising AUPRC).
3. Runs Optuna HPO for XGBoost (200 trials, same protocol).
4. Trains final LightGBM and XGBoost models on the full training set with best HPO parameters.
5. Calibrates each model with `CalibratedClassifierCV(method='isotonic', cv='prefit')`.
6. Builds a 50/50 ensemble: `0.5 * P(lgb) + 0.5 * P(xgb)`.
7. Exports the XGBoost model to ONNX (requires `skl2onnx`). Pass `--no-onnx` to skip.

**Options**

| Flag | Default | Description |
|---|---|---|
| `--hpo-trials` | 200 | Optuna trials per model |
| `--no-onnx` | off | Skip ONNX export |

**Outputs**

| Path | Description |
|---|---|
| `bigdata/.../models/lightgbm_twist.txt` | LightGBM model |
| `bigdata/.../models/xgboost_twist.json` | XGBoost model |
| `bigdata/.../models/ensemble_twist.pkl` | Calibrated ensemble (pickle) |
| `bigdata/.../models/xgboost_twist.onnx` | ONNX model for kam deployment |
| `docs/.../results/cv_results.csv` | Cross-validation results (tracked in git) |
| `docs/.../results/test_results.json` | Held-out test metrics (tracked in git) |
| `docs/.../results/feature_importance_lgb.csv` | LightGBM feature importances (tracked in git) |
| `docs/.../results/feature_importance_xgb.csv` | XGBoost feature importances (tracked in git) |
| `docs/.../results/figures/` | Training plots (tracked in git) |

**Estimated runtime:** 4–8 hours for 200 HPO trials per model.

---

## 7. Step 5: Evaluate

Run from the repository root.

```bash
python3 scripts/ml/evaluate_twist_duplex.py
```

This loads all three models (LightGBM, XGBoost, ensemble) and optionally the ONNX model, runs them on `test_features.csv.gz`, and produces a comprehensive report.

**Metrics produced**

- Overall: AUROC, AUPRC, best F1, sensitivity at 1% FPR, sensitivity at 5% FPR.
- Per-VAF bin: AUPRC, F1, and counts for each VAF range.
- Per-variant type: AUPRC, F1, and counts for SNV, short indel, medium indel, SV.
- Curve data: ROC and precision-recall curves for each model.

**Outputs**

| Path | Description |
|---|---|
| `docs/.../results/evaluation_report.json` | Full metrics dict (tracked in git) |
| `docs/.../results/roc_data.csv` | ROC curve data per model (tracked in git) |
| `docs/.../results/pr_data.csv` | Precision-recall curve data per model (tracked in git) |
| `docs/.../results/confusion_matrix.csv` | Confusion matrix at optimal F1 threshold (tracked in git) |
| `docs/.../results/vaf_bin_breakdown.csv` | Per-VAF-bin breakdown (tracked in git) |
| `docs/.../results/vtype_breakdown.csv` | Per-variant-type breakdown (tracked in git) |
| `bigdata/.../results/figures/` | Full-resolution figures (not tracked) |

---

## 8. Step 6: Upload to Nextcloud

Run from the repository root. Requires `NEXTCLOUD_EDIT_TOKEN` in `.env`.

```bash
# Upload trained models
bash scripts/ml/upload_twist_duplex.sh models

# Upload feature CSVs
bash scripts/ml/upload_twist_duplex.sh features

# Upload configs tarball
bash scripts/ml/upload_twist_duplex.sh configs

# Upload everything
bash scripts/ml/upload_twist_duplex.sh all
```

The `train` and `test` sub-commands upload simulation batch tarballs. These are normally handled automatically by Step 2 after each batch. Run them manually only if the automatic upload was skipped.

**Nextcloud layout**

| Description | Local path | Nextcloud remote path |
|---|---|---|
| Train simulation batches | `bigdata/.../simulations/train/` (tarballs) | `experiments/01-ml-twist-duplex/train/batch_NNN.tar.gz` |
| Test simulation batches | `bigdata/.../simulations/test/` (tarballs) | `experiments/01-ml-twist-duplex/test/batch_NNN.tar.gz` |
| Trained models | `bigdata/.../models/` | `experiments/01-ml-twist-duplex/models/` |
| Feature CSVs | `bigdata/experiments/01-ml-twist-duplex/*.csv.gz` | `experiments/01-ml-twist-duplex/` |
| Configs tarball | `docs/.../configs/` | `experiments/01-ml-twist-duplex/configs.tar.gz` |

---

## 9. Step 7: Deploy to the kam binary

The ONNX model (`xgboost_twist.onnx`) and its metadata file (`xgboost_twist.onnx.meta.json`) are the deployment artefacts. Copy them to a stable path accessible at runtime.

**Build with ML support**

```bash
cargo build --release --features ml
```

The default build (`cargo build --release`) omits ML support to keep the binary small. Passing `--ml-model` to a no-ML binary will print a warning and ignore the flag.

**Call with the model**

```bash
kam call --ml-model /path/to/xgboost_twist.onnx [other flags]
```

Or via the integrated runner:

```bash
kam run --ml-model /path/to/xgboost_twist.onnx [other flags]
```

The model file and its `.meta.json` sidecar must be in the same directory.

---

## 10. Resuming an interrupted run

`checkpoint.json` lists the name of every sample that has completed Step 2 successfully. Re-running the pipeline command skips any sample already in the checkpoint.

```bash
# Resume exactly where it stopped
python3 scripts/ml/run_twist_duplex_pipeline.py --split all --workers 16
```

To force a specific sample to re-run:

1. Open `docs/project/experiments/01-ml-twist-duplex/checkpoint.json`.
2. Remove the sample name from the list.
3. Delete the sample output directory:
   ```bash
   rm -rf bigdata/experiments/01-ml-twist-duplex/simulations/{split}/{sample_name}/
   ```
4. Re-run the pipeline command.

The pipeline is otherwise idempotent: re-running it for already-completed samples is safe (they will be skipped).

---

## 11. Nextcloud layout

| Description | Local path | Nextcloud remote path |
|---|---|---|
| Train simulation batches | `bigdata/experiments/01-ml-twist-duplex/simulations/train/` | `experiments/01-ml-twist-duplex/train/batch_NNN.tar.gz` |
| Test simulation batches | `bigdata/experiments/01-ml-twist-duplex/simulations/test/` | `experiments/01-ml-twist-duplex/test/batch_NNN.tar.gz` |
| Train feature matrix | `bigdata/experiments/01-ml-twist-duplex/train_features.csv.gz` | `experiments/01-ml-twist-duplex/train_features.csv.gz` |
| Test feature matrix | `bigdata/experiments/01-ml-twist-duplex/test_features.csv.gz` | `experiments/01-ml-twist-duplex/test_features.csv.gz` |
| LightGBM model | `bigdata/experiments/01-ml-twist-duplex/models/lightgbm_twist.txt` | `experiments/01-ml-twist-duplex/models/lightgbm_twist.txt` |
| XGBoost model | `bigdata/experiments/01-ml-twist-duplex/models/xgboost_twist.json` | `experiments/01-ml-twist-duplex/models/xgboost_twist.json` |
| Ensemble model | `bigdata/experiments/01-ml-twist-duplex/models/ensemble_twist.pkl` | `experiments/01-ml-twist-duplex/models/ensemble_twist.pkl` |
| ONNX model | `bigdata/experiments/01-ml-twist-duplex/models/xgboost_twist.onnx` | `experiments/01-ml-twist-duplex/models/xgboost_twist.onnx` |
| Configs tarball | `docs/project/experiments/01-ml-twist-duplex/configs/` | `experiments/01-ml-twist-duplex/configs.tar.gz` |

---

## Runtime summary

| Step | Estimated time (12–16 workers) |
|---|---|
| Step 1: Generate configs | ~2 minutes |
| Step 2: varforge simulation (11k samples) | 8–16 hours |
| Step 2: kam discovery + tumour-informed (11k samples) | 8–16 hours |
| Step 3: Feature extraction | ~10 minutes |
| Step 4: Model training (200 HPO trials each) | 4–8 hours |
| Step 5: Evaluation | ~5 minutes |
| **Total** | **~1–2 days** |
