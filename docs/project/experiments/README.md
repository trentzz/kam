# experiments

Research experiments that go beyond routine benchmarking. Each experiment has a clear hypothesis, a defined dataset, and a recorded outcome. Configs, scripts, and metadata live here. Large generated files (FASTQs, models) live in `bigdata/experiments/` and on Nextcloud.

## Experiment log

| # | Name | Hypothesis | Status | Outcome |
|---|---|---|---|---|
| 01 | [ml-twist-duplex](01-ml-twist-duplex/README.md) | A chemistry-specific classifier trained only on Twist duplex data will outperform the generic ML3 model by exploiting duplex-vs-simplex evidence features | In progress | Pending training runs |

## Adding a new experiment

1. Create `docs/project/experiments/NN-name/` with a `README.md` describing hypothesis, dataset, and scripts.
2. Create a matching `bigdata/experiments/NN-name/` directory for large outputs.
3. Add a row to this table.
4. Upload outputs to Nextcloud under `experiments/NN-name/`.
