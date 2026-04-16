# SV-SWP-007: Aggregate Sweep Results

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-SWP-002, SV-SWP-004, SV-SWP-006
**Status**: done

## Goal

Produce heatmaps and line plots aggregating all sweep results into actionable
figures for parameter selection.

## Success Criteria

- [x] Aggregation script at `scripts/benchmarking/aggregate_sweep.py`.
- [x] Produces summary CSV with columns: parameter, value, vaf, sv_type, sensitivity, precision, n_tp, n_fp.
- [x] Generates heatmap CSVs and matplotlib plots (when available).
- [x] Prints recommended optimal parameters per SV type.
- [ ] `/update` has been run after changes.

## Steps

1. Collate scored results from all three sweeps.
2. Generate heatmaps for k-mer and threshold sweeps.
3. Generate line plots for size sweep.
4. Add recommended configuration per SV type based on results.

## Notes

Follow the project graph style guide for consistent figure formatting.
