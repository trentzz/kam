# SV-ULV-001: Generate Ultra-Low VAF Configs

**Epic**: SV-ULTRAVAF
**Priority**: medium-low
**Depends on**: SV-BASE-007
**Status**: todo

## Goal

Create varforge configs for 0.01%, 0.02%, 0.03%, 0.04% VAF across all 3 SV
types with 2 replicates each (24 configs total). Generate truth VCFs.

## Success Criteria

- [ ] 24 varforge configs generated (4 VAFs x 3 types x 2 replicates).
- [ ] Truth VCFs generated for each config.
- [ ] `/update` has been run after changes.

## Steps

1. Adapt existing config generation scripts for sub-0.05% VAF levels.
2. Generate all 24 configs.
3. Generate truth VCFs.

## Notes

At 5000x coverage, 0.01% VAF means ~0.5 variant molecules per dataset.
Detection at this level is expected to be near zero for most types. The
purpose is to find the floor, not to achieve high sensitivity.
