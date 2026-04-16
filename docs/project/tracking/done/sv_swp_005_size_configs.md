# SV-SWP-005: Generate SV Size Sweep Configs

**Epic**: SV-SWEEP
**Priority**: medium
**Depends on**: SV-BASE-007
**Status**: done

## Goal

Create varforge configs for SV sizes 20, 50, 100, 200, 500 bp at 1% VAF for
all three SV types with 2 replicates each. Generate truth VCFs.

## Success Criteria

- [x] 28 varforge configs generated in `scripts/benchmarking/sv_size_sweep_configs/` (DEL: 5 sizes x 2 reps, DUP: 5 sizes x 2 reps, INV: 4 sizes x 2 reps).
- [x] Generator script at `scripts/benchmarking/generate_size_sweep_configs.py`.
- [ ] `/update` has been run after changes.

## Steps

1. Adapt existing config generation scripts to parameterise SV size.
2. Generate configs for all size x type x replicate combinations.
3. Generate truth VCFs.

## Notes

Not all SV types support all sizes equally. InvDel and NovelInsertion scale
naturally with size. Fusion junction size is fixed, but the flanking deletion
can vary.
