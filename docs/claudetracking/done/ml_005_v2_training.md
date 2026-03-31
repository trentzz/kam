# ML-005: Train v2 models on combined dataset

**Epic**: ML-BOOST (overallplans/ML-BOOST.md)
**Priority**: high
**Depends on**: ML-001 (done), ML data generation (in progress)
**Status**: todo

## Goal

Once `docs/benchmarking/ml/training_data_v2.csv` exists (produced by
`scripts/ml/build_training_data_v2.py`), run the full v2 training comparison:

```
python3 scripts/ml/train_eval_v2.py
```

This compares:
- LightGBM raw (9 features)
- LightGBM v2 (30+ features with log transforms, interactions, ratios)
- XGBoost raw
- XGBoost v2
- Baseline (PASS filter)

And saves:
- `docs/benchmarking/ml/results/cv_results_v2.csv`
- `docs/benchmarking/ml/results/feature_importance_v2.csv`
- `docs/benchmarking/ml/models/lightgbm_v2.txt`
- `docs/benchmarking/ml/models/xgboost_v2.json`

## Steps

1. Confirm `training_data_v2.csv` exists and has rows from both the original
   250 samples and the new ML-dataset samples (look for param columns).
2. Run `python3 scripts/ml/train_eval_v2.py`.
3. Check the summary table output. Note whether v2 features improve AUPRC
   over raw features.
4. If params columns (coverage, family_size_mean) have non-zero variance in
   the new data, check their feature importance.
5. Update this task as done.

## Success criteria

- [ ] `cv_results_v2.csv` exists with 5 configs × 5 folds = 25 rows
  (plus 5 baseline rows)
- [ ] `feature_importance_v2.csv` exists
- [ ] Both LightGBM and XGBoost models saved
- [ ] Summary table shows whether v2 features improve over raw features
