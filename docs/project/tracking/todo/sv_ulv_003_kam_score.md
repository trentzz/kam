# SV-ULV-003: Run kam and Score Ultra-Low VAF Results

**Epic**: SV-ULTRAVAF
**Priority**: medium-low
**Depends on**: SV-ULV-002
**Status**: todo

## Goal

Run kam in both modes on all 24 ultra-low VAF datasets and score the results.

## Success Criteria

- [ ] All 24 kam result directories contain discovery and TI output.
- [ ] Scored results exist with sensitivity per VAF per type.
- [ ] `/update` has been run after changes.

## Steps

1. Run kam discovery and tumour-informed modes on all 24 datasets.
2. Run the scoring script.
3. Record results.

## Notes

Expect very low or zero sensitivity at 0.01-0.02% VAF. The value is in
establishing the exact detection floor.
