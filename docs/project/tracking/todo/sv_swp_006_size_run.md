# SV-SWP-006: Run Size Sweep

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-SWP-005
**Status**: todo

## Goal

Run varforge simulations, kam, and scoring for all SV size sweep configs.

## Success Criteria

- [ ] All size sweep simulations produce FASTQ + truth VCF.
- [ ] All kam runs produce discovery and TI output.
- [ ] All results are scored.
- [ ] `/update` has been run after changes.

## Steps

1. Run varforge for all size sweep configs.
2. Run kam on all simulated datasets.
3. Score results using the scoring script.

## Notes

This sweep answers the question: at what SV size does detection break down
for each type?
