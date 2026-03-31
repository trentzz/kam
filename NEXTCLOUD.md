# Nextcloud Storage

Large files that are not tracked in git are stored on Nextcloud.

Public share: `https://nextcloudlocal.trentz.me/s/pTizAiSAJQsPcDo`
WebDAV root: `https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo/`

## Contents

| Path on Nextcloud | Description | Reproduced by |
|---|---|---|
| `benchmarking/sv-results/` | varforge-generated FASTQs for SV benchmark (fusion, monitoring) | `docs/benchmarking/sv/` scripts |
| `benchmarking/sv-new-results/` | varforge-generated FASTQs for SV new-type benchmark (INS, large DEL, INVDEL, novel insertion) | `docs/benchmarking/sv_new/` scripts |
| `benchmarking/ml-fastqs/` | varforge-generated FASTQs for ML training dataset | `scripts/ml/generate_ml_dataset.py` |
| `benchmarking/per-sample/` | Per-sample TSV results for SNV/indel benchmark | `docs/benchmarking/snvindel/` scripts |
| `lightgbm_model.txt` | Trained LightGBM classifier (binary format) | `scripts/ml/train_model.py` |
| `training_data_v2.csv.gz` | ML training features extracted from benchmarking runs | `scripts/ml/extract_features.py` |
| `ml-boost/` | ML boost pipeline outputs | `scripts/ml/` |
| `paper/` | Paper PDFs and supplementary materials | `docs/paper/` |

## Uploading

To re-upload a set of FASTQs, tar each sim directory and PUT it via WebDAV:

```bash
tar -czf /tmp/name.tar.gz -C <results-parent-dir> <sim-dir-name>
curl -T /tmp/name.tar.gz \
    "https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo/<path>/name.tar.gz" \
    -u "pTizAiSAJQsPcDo:"
rm /tmp/name.tar.gz
```
