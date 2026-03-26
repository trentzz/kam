# SV-SWP-007: Aggregate Sweep Results

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-SWP-002, SV-SWP-004, SV-SWP-006
**Status**: todo

## Goal

Produce heatmaps and line plots aggregating all sweep results into actionable
figures for parameter selection.

## Success Criteria

- [ ] Heatmap: k-mer size vs VAF, coloured by sensitivity, per SV type in `docs/benchmarking/sv_new/figures/`.
- [ ] Heatmap: confidence threshold vs min_alt_molecules per SV type.
- [ ] Line plot: SV size vs sensitivity per type.
- [ ] `/update` has been run after changes.

## Steps

1. Collate scored results from all three sweeps.
2. Generate heatmaps for k-mer and threshold sweeps.
3. Generate line plots for size sweep.
4. Add recommended configuration per SV type based on results.

## Notes

Follow the project graph style guide for consistent figure formatting.
