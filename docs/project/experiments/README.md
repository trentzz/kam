# experiments

Research experiments that go beyond routine benchmarking. Each experiment has a clear hypothesis, a defined dataset, and a recorded outcome. Configs, scripts, and metadata live here. Large generated files (FASTQs, models) live in `bigdata/experiments/` and on Nextcloud.

## Experiment log

| # | Name | Hypothesis | Status | Outcome |
|---|---|---|---|---|
| 01 | [ml-twist-duplex](01-ml-twist-duplex/README.md) | A chemistry-specific classifier trained only on Twist duplex data will outperform the generic ML3 model by exploiting duplex-vs-simplex evidence features | In progress | Pending training runs |
| 02 | [ml-single-strand](02-ml-single-strand/README.md) | A gradient boosting model trained on simulated single-strand indel/SNV data can improve precision over the raw PASS filter at low VAF | Complete | v3 LightGBM: AUPRC 0.920 on 4.85M rows; ONNX integrated into `kam run --ml-model` |
| 03 | [ml-twist-duplex-v2](03-ml-twist-duplex-v2/design.md) | A model retrained on real TP/FP calls from the titration dataset with new sequence-context features (subst_type, trinuc_context) will produce discriminating ml_prob distributions in discovery mode | Complete | Confirmed: AUPRC 0.973, FP reduction 79–93%, precision 0.969–1.000 at 1–2% VAF |

## Adding a new experiment

1. Create `docs/project/experiments/NN-name/` with a `README.md` describing hypothesis, dataset, and scripts.
2. Create a matching `bigdata/experiments/NN-name/` directory for large outputs.
3. Add a row to this table.
4. Upload outputs to Nextcloud under `experiments/NN-name/`.
