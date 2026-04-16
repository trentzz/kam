# SV-SWP-003: Write Threshold Sweep Script

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-BASE-007
**Status**: done

## Goal

Write a script that varies sv_min_confidence (0.80, 0.90, 0.95, 0.99) and
sv_min_alt_molecules (1, 2, 3) across the same 6 key VAFs x 3 types x 2
replicates.

## Success Criteria

- [x] Executable script at `scripts/benchmarking/sv_threshold_sweep.sh`.
- [x] Script parameterises confidence and min_alt_molecules values.
- [ ] `/update` has been run after changes.

## Steps

1. Design the parameter grid (4 confidence x 3 molecule thresholds = 12 combos).
2. Write the script, iterating over the grid and dataset selection.
3. Test on one dataset.

## Notes

Use the same 6 key VAFs as the k-mer sweep for consistency.
