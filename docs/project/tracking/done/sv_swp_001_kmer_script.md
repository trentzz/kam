# SV-SWP-001: Write K-mer Sweep Script

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-BASE-007
**Status**: done

## Goal

Write a script that runs kam with k=21,25,27,31,35,41 on 6 key VAFs x 3 SV
types x 2 replicates (216 combinations, 432 total with discovery + TI modes).

## Success Criteria

- [x] Executable script at `scripts/benchmarking/sv_kmer_sweep.sh`.
- [x] Script accepts `--data-dir`, `--output-dir`, `--reference` arguments.
- [ ] `/update` has been run after changes.

## Steps

1. Design the directory naming convention for sweep results.
2. Write the script, parameterising k-mer size and dataset selection.
3. Test on one dataset to verify correctness.

## Notes

Select 6 key VAFs that span the sensitivity curve: e.g. 0.1%, 0.25%, 0.5%,
1%, 2%, 5%. Use both replicates at each VAF.
