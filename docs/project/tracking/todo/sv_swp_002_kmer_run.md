# SV-SWP-002: Run K-mer Sweep

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-SWP-001
**Status**: todo

## Goal

Execute the k-mer sweep script and collect all results.

## Success Criteria

- [ ] All 432 result directories contain VCF and TSV output.
- [ ] No empty or truncated output files.
- [ ] `/update` has been run after changes.

## Steps

1. Run the sweep script from SV-SWP-001.
2. Monitor for failures and re-run as needed.
3. Verify output completeness.

## Notes

This is a compute-intensive task. Consider running in parallel batches.
Each kam run takes approximately 1-2 minutes per dataset.
