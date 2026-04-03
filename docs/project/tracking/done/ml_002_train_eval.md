# ML-002: Train and evaluate LightGBM and XGBoost classifiers

**Epic**: ML-BOOST (overallplans/ML-BOOST.md)
**Priority**: high
**Depends on**: ML-001
**Status**: todo

## Goal

Train LightGBM and XGBoost classifiers on the labelled call dataset. Evaluate
using group k-fold CV. Report precision, recall, F1, and AUPRC. Compare
against the current rule-based baseline (PASS filter).

## Steps

1. Load `docs/benchmarking/ml/training_data.csv`.
2. Define features (see ML-BOOST epic for feature list).
3. Encode `variant_class`, `mode`, `dataset_type` as categoricals.
4. Run 5-fold group k-fold CV (group = sample_id).
5. Train LightGBM with class weight balancing (scale_pos_weight or
   is_unbalance=True).
6. Train XGBoost with scale_pos_weight.
7. For each fold: compute precision, recall, F1 at default threshold (0.5),
   AUPRC, AUROC. Average across folds.
8. Compute baseline metrics: treat PASS filter as the classifier (label=1
   if filter==PASS, else 0).
9. Write a per-fold result table to
   `docs/benchmarking/ml/results/cv_results.csv`.
10. Save the best LightGBM model (by mean AUPRC) to
    `docs/benchmarking/ml/models/lightgbm_model.txt`.
11. Print feature importance (gain) for top 10 features.

## Notes

- Use `lightgbm` and `xgboost` Python packages.
- Use `scikit-learn` for GroupKFold, metrics.
- If the dataset is very small (<500 rows), note this in the output and
  interpret results accordingly.
- Do not tune hyperparameters in this task — use sensible defaults. Tuning
  is out of scope.

## Success criteria

- [ ] `docs/benchmarking/ml/results/cv_results.csv` exists with fold-level
  metrics for LightGBM, XGBoost, and baseline
- [ ] `docs/benchmarking/ml/models/lightgbm_model.txt` saved
- [ ] Feature importance table printed to stdout and written to
  `docs/benchmarking/ml/results/feature_importance.csv`
- [ ] Mean AUPRC reported for each model and the baseline
