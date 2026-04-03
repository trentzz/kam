# ML-004: Figures and summary report

**Epic**: ML-BOOST (overallplans/ML-BOOST.md)
**Priority**: medium
**Depends on**: ML-002, ML-003
**Status**: todo

## Goal

Produce figures and a summary report from the ML experiment results. The
report should answer: does gradient boosting improve on the rule-based filter?
If so, by how much, and at what VAF levels?

## Steps

1. Load `docs/benchmarking/ml/results/cv_results.csv`.
2. Plot: precision-recall curve for LightGBM, XGBoost, and baseline on the
   held-out fold data.
3. Plot: AUPRC by variant class (SNV, indel, SV) for each model.
4. Plot: feature importance bar chart (top 10 features by gain).
5. Plot: sensitivity at fixed FPR (0.01 calls per sample) vs VAF, comparing
   LightGBM to baseline. Group by variant class.
6. Write a 1-2 page summary to `docs/research/ml_experiment_results.md`
   covering: dataset stats, model comparison, key features, conclusions, and
   recommended next steps.
7. Write figures to `docs/benchmarking/ml/results/figures/`.

## Notes

- Use matplotlib or seaborn. Follow `docs/claudeguide/` graph style if a
  graph-style guide exists.
- The summary should be direct and factual. State what worked and what did not.
- If results are poor, that is a valid and useful finding.

## Success criteria

- [ ] At least 4 figures in `docs/benchmarking/ml/results/figures/`
- [ ] `docs/research/ml_experiment_results.md` written
- [ ] Summary includes a clear recommendation (pursue or defer Rust integration)
