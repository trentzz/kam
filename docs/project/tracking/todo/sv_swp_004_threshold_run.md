# SV-SWP-004: Run Threshold Sweep

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-SWP-003
**Status**: todo

## Goal

Execute the threshold sweep script and collect all results.

## Success Criteria

- [ ] All threshold sweep result directories are populated with VCF and TSV output.
- [ ] No empty or truncated output files.
- [ ] `/update` has been run after changes.

## Steps

1. Run the threshold sweep script from SV-SWP-003.
2. Monitor for failures and re-run as needed.
3. Verify output completeness.

## Notes

The 12-combination grid x 36 datasets = 432 runs, similar in scale to the
k-mer sweep.
