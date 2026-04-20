# DEPS-003: Wire cargo-audit into CI

**Epic**: DEPS
**Priority**: medium
**Depends on**: QUALITY-002
**Status**: todo

## Goal

Run `cargo audit` on every pull request and fail the build on new advisories. Maintain an `audit.toml` to document intentionally allowed advisories.

## Success Criteria

- [ ] `.github/workflows/ci.yml` runs `cargo audit --deny warnings` (or equivalent gating).
- [ ] `audit.toml` documents any currently-allowed advisories with rationale.
- [ ] Failing audit in a PR blocks merge.
- [ ] `/update` has been run after changes.

## Steps

1. Add a dedicated `audit` job to CI.
2. Seed `audit.toml` with DEPS-001 / DEPS-002 outcomes.

## Notes

Gate on QUALITY-002 landing the initial CI skeleton.
