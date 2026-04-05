# SV-BASE-004: Run kam on InvDel Datasets

**Epic**: SV-BASELINE
**Priority**: critical
**Depends on**: SV-BASE-001
**Status**: todo

## Goal

Run kam in both discovery and tumour-informed modes on all 50 InvDel simulated
datasets to produce variant calls.

## Success Criteria

- [ ] All 50 kam_invdel_vaf* directories contain calls_discovery.vcf, calls_discovery.tsv, calls_tumour_informed.vcf, and calls_tumour_informed.tsv.
- [ ] `/update` has been run after changes.

## Steps

1. Build kam with the latest code.
2. Run kam discovery mode on all 50 InvDel datasets.
3. Run kam tumour-informed mode on all 50 InvDel datasets.
4. Verify output completeness.

## Notes

Use the run script in `docs/benchmarking/sv_new/scripts/` if available,
otherwise adapt from the existing SNV/indel run scripts.
