# SV-ULV-004: Produce Detection Limit Curves

**Epic**: SV-ULTRAVAF
**Priority**: medium-low
**Depends on**: SV-ULV-003
**Status**: todo

## Goal

Plot detection limit curves showing sensitivity vs VAF at the ultra-low end,
combined with the baseline data from SV-BASE-008 for a complete picture.

## Success Criteria

- [ ] Figures in `docs/benchmarking/sv_new/figures/` showing detection floor per SV type.
- [ ] Figures combine ultra-low and baseline VAF ranges on the same plot.
- [ ] Minimum detectable VAF at 80% sensitivity noted per type.
- [ ] `/update` has been run after changes.

## Steps

1. Merge ultra-low scored results with baseline scored results.
2. Plot extended sensitivity curves from 0.01% to 10% VAF.
3. Annotate the 80% sensitivity threshold.
4. Save figures.

## Notes

Use a log scale on the x-axis to show the full VAF range clearly.
