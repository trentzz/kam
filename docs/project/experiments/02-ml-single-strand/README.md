# Single-Strand ML

This experiment trains a gradient boosting classifier to filter variant calls from
single-strand (non-duplex) sequencing data. It preceded the Twist duplex classifier
in experiment 01.

**Hypothesis:** a gradient boosting model trained on simulated single-strand indel
and SNV data can improve precision over the raw PASS filter at low VAF, where
spurious low-evidence calls are common.

## Status

v1 and v2 models trained and evaluated. v3 model trained on the ML3 dataset (indel_ml2,
4.85M rows, AUPRC 0.920). ONNX export complete. Rust integration added via `--ml-model`
flag in `kam run`.

## Models

| Version | File | AUPRC | Dataset |
|---|---|---|---|
| v1 | `models/lightgbm_model.txt` | — | v1 (240k rows) |
| v2 | `models/lightgbm_v2.txt`, `models/xgboost_v2.json` | ~0.85 | v2 |
| v3 | `models/lightgbm_v3.txt`, `models/xgboost_v3.json` | 0.920 | ML3 (4.85M rows) |
| v3 ONNX | `models/lightgbm_rust.onnx` | 0.920 | ML3 |

Model files are excluded from git (see `.gitignore`). They live in
`bigdata/experiments/02-ml-single-strand/models/` and on Nextcloud under
`experiments/02-ml-single-strand/`.

## Datasets

### v1/v2 — indel/SNV sweep (samples/)

Per-sample kam outputs from the original VAF sweep. Each sample directory contains:

| File | Description |
|---|---|
| `config.yaml` | varforge simulation config |
| `params.json` | simulation parameters |
| `truth.tsv` | spiked-in variant truth set |
| `discovery.tsv` | kam discovery mode output |
| `tumour_informed.tsv` | kam tumour-informed mode output |
| `varforge_cmd.txt` | exact varforge command used |

Naming convention: `indel_vaf{VAF}_{rep}` or `sv_vaf{VAF}_{rep}`.

### ML3 — indel_ml2 dataset (configs/)

10,000 train + 1,000 test varforge configs for indel variants at VAF 0.0003–0.05.
Configs are excluded from git (stored on Nextcloud). Large simulation outputs
(results/train/, results/test/) and training CSVs are in bigdata/.

## Results

Training metrics (committed small files):

| File | Description |
|---|---|
| `results/cv_results.csv` | v1 cross-validation metrics |
| `results/cv_results_v2.csv` | v2 cross-validation metrics |
| `results/cv_results_v3.csv` | v3 cross-validation metrics |
| `results/feature_importance*.csv` | Feature importance by version |
| `results/figures/` | AUPRC/AUROC curves and feature importance plots |
| `results_summary.md` | Narrative summary of v1 training run |

Large training/test feature CSVs and raw simulation outputs live in
`bigdata/experiments/02-ml-single-strand/` and on Nextcloud.

## Nextcloud

Stored under `experiments/02-ml-single-strand/`. See
`docs/project/devmanual/nextcloud.md` for upload/download instructions.

## Scripts

| Script | Description |
|---|---|
| `scripts/ml/generate_twist_duplex_configs.py` | Generate varforge configs |
| `scripts/ml/run_twist_duplex_pipeline.py` | Run the training pipeline |
| `scripts/ml/build_twist_duplex_features.py` | Build feature CSVs from kam output |
| `scripts/ml/build_training_data.py` | Aggregate training data |
| `scripts/ml/train_twist_duplex.py` | Train and export models |
| `scripts/ml/upload_twist_duplex.sh` | Upload outputs to Nextcloud |
