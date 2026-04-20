# DEPS-002: Upgrade statrs to clear rand 0.8 and paste advisories

**Epic**: DEPS
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

Upgrade `statrs` (or its transitive deps) so `rand 0.8.5` (RUSTSEC-2026-0097) and `paste 1.0.15` (RUSTSEC-2024-0436) no longer appear in `cargo audit`.

## Success Criteria

- [ ] `cargo audit` reports zero advisories for `rand` and `paste`.
- [ ] All `call_variant`, `estimate_vaf`, and `strand_bias_test` tests still pass without numerical drift.
- [ ] All tests pass across the workspace.
- [ ] `/update` has been run after changes.

## Steps

1. Check statrs CHANGELOG for a release on rand 0.9. If none, pin via Cargo patch.
2. Bump `statrs` in `Cargo.toml`.
3. Run full test suite; spot-check numeric outputs against 01-snvindel golden summary.

## Notes

`statrs::distribution::Beta::inverse_cdf` is load-bearing for VAF credible intervals. Verify no precision regression.
