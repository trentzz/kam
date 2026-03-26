# SV-BASE-008: Generate Baseline Sensitivity Curves

**Epic**: SV-BASELINE
**Priority**: high
**Depends on**: SV-BASE-007
**Status**: todo

## Goal

Plot sensitivity vs VAF curves for each SV type in both discovery and
tumour-informed modes. These are the baseline figures for all subsequent
experiments.

## Success Criteria

- [ ] PNG figures in `docs/benchmarking/sv_new/figures/` showing sensitivity curves per SV type.
- [ ] Each figure includes both discovery and tumour-informed mode curves.
- [ ] `/update` has been run after changes.

## Steps

1. Read sv_new_summary.tsv and sv_new_per_dataset.tsv.
2. Plot sensitivity vs VAF for each of the three SV types.
3. Include error bars or replicate scatter where applicable.
4. Save figures as PNG.

## Notes

Follow the project graph style guide. Use consistent colours across SV types.
