# Single-Strand ML

This experiment trains a gradient boosting classifier to filter variant calls from
single-strand (non-duplex) sequencing data. It preceded the Twist duplex classifier
in experiment 01.

**Hypothesis:** a gradient boosting model trained on simulated single-strand indel
and SNV data can improve precision over the raw PASS filter at low VAF, where
spurious low-evidence calls are common.

## Status

v1 and v2 models trained and evaluated. v2b model trained on the v2 dataset (training_data_v2.csv,
888k rows, AUPRC 0.920). ONNX export complete. Rust integration added via `--ml-model`
flag in `kam run`. Note: v2b was mislabelled "v3" at training time. The v3 name is reserved
for future training on the ML3 train/test split (training_data_v3.csv).

## Models

| Version | File | AUPRC | Dataset |
|---|---|---|---|
| v1 | `models/lightgbm_model.txt` | — | v1 (240k rows) |
| v2 | `models/lightgbm_v2.txt`, `models/xgboost_v2.json` | ~0.85 | v2 |
| v2b | `models/lightgbm_v2b.txt`, `models/xgboost_v2b.json` | 0.920 | v2 (888k rows, legacy + ML2 param sweep) |
| v2b ONNX | `models/lightgbm_v2b.onnx` | 0.920 | v2 |

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
| `results/cv_results_v2b.csv` | v2b cross-validation metrics |
| `results/feature_importance*.csv` | Feature importance by version |
| `results/figures/` | AUPRC/AUROC curves and feature importance plots |
| `results_summary.md` | Narrative summary of v1 training run |

Large training/test feature CSVs and raw simulation outputs live in
`bigdata/experiments/02-ml-single-strand/` and on Nextcloud.

## Nextcloud

Stored under `experiments/02-ml-single-strand/`. See
`docs/project/devmanual/nextcloud.md` for upload/download instructions.

## ML3 Train Pipeline — Monitoring

The ML3 train pipeline runs in the background via `nohup`. It survives closing the terminal or Claude session. Use these commands to check progress:

```bash
# How many kam runs are complete (target: 10,000)
ls bigdata/experiments/02-ml-single-strand/results/train/ | grep "^kam_" | wc -l

# Breakdown by variant type
ls bigdata/experiments/02-ml-single-strand/results/train/ | grep "^kam_" | \
  sed 's/_ml3_.*//' | sed 's/^kam_//' | sort | uniq -c

# How many sample dirs are built (target: 10,000 — happens after all kam done)
ls bigdata/experiments/02-ml-single-strand/samples/train/ | wc -l

# Live log tail
tail -f ~/tmp/ml3_train_pipeline_final.log

# Check if pipeline process is still running
ps aux | grep run_ml3_train_pipeline | grep -v grep
```

Once `samples/train/` reaches 10,000, run in order:

```bash
python3 scripts/ml/build_training_data_v3.py
python3 scripts/ml/train_eval_v3.py
```

## Scripts

| Script | Description |
|---|---|
| `scripts/ml/generate_twist_duplex_configs.py` | Generate varforge configs |
| `scripts/ml/run_twist_duplex_pipeline.py` | Run the training pipeline |
| `scripts/ml/build_twist_duplex_features.py` | Build feature CSVs from kam output |
| `scripts/ml/build_training_data.py` | Aggregate training data |
| `scripts/ml/train_twist_duplex.py` | Train and export models |
| `scripts/ml/upload_twist_duplex.sh` | Upload outputs to Nextcloud |
