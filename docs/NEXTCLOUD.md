# Nextcloud Storage

Large output files that exceed GitHub's limits are stored in Nextcloud.

**Access**: https://nextcloudlocal.trentz.me/s/pTizAiSAJQsPcDo

## Contents

| Folder | Contents | Reproduced by |
|--------|----------|---------------|
| `ml-boost/` | ML training data, LightGBM model, CV results | `scripts/ml/build_training_data_v2.py`, `scripts/ml/train_eval_v2.py` |
| `benchmarking/sv-new-results/` | sv_new simulation FASTQs, VCFs, QC JSON | `docs/benchmarking/sv_new/scripts/run_sv_suite.sh` |
| `benchmarking/sv-results/` | sv simulation FASTQs and truth VCFs | `docs/benchmarking/sv/scripts/run_sv_suite.sh` |
| `benchmarking/per-sample/` | Real clinical sample call TSVs | (clinical data, not regenerable) |
| `paper/` | Compiled paper PDFs | `docs/paper/` LaTeX source |

## Download

Use the WebDAV URL to download programmatically:

```bash
curl -O "https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo/<folder>/<file>" \
     -u "pTizAiSAJQsPcDo:"
```
