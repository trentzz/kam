# Comparisons and Baselines

## Primary Baseline

**Thesis alignment-based pipeline** (RaSCALL with HUMID/Redux deduplication).

- Same samples as the kam titration experiments.
- Results are in `/mnt/tzeng-local/tzeng-thesis/`:
  - `titration.ref.v2.redux.detection.xlsx` — per-variant detection table.
  - `benchmark-dataframes/` — structured benchmark data.
- This is the comparison that matters most. All sensitivity and precision claims are against this baseline.

## Secondary Baseline

**Public tool benchmarks** for public datasets.

- Use whatever alignment-based results exist for the public datasets included in the paper.
- Cite published benchmarks from the original dataset papers where available.
- Do not run new alignment-based pipelines on public data unless necessary.

---

## What to Compare

For each baseline, report:

- Sensitivity at each VAF level.
- Precision.
- F1 score.
- Runtime on matched input.

For per-variant concordance, document each variant individually: detected by kam, detected by baseline, detected by both, detected by neither.
