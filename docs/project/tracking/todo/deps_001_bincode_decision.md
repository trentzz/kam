# DEPS-001: Resolve bincode 1.x Unmaintained Advisory

**Epic**: DEPS
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

Decide and implement a path for RUSTSEC-2025-0141 (bincode 1.3.3 unmaintained). Either migrate to bincode 2.x with on-disk version header bump and a converter, or pin the version and document the decision.

## Success Criteria

- [ ] `docs/project/vision/decisions/log.md` has a dated entry for the bincode choice.
- [ ] `cargo audit` either does not report bincode, or the advisory is documented as allowed in `deny.toml` / `audit.toml`.
- [ ] If migrating: all inter-stage file formats round-trip through the new version; existing v1 files are handled by a `kam convert` sub-command or clearly rejected.
- [ ] If pinning: rationale is in the devmanual release doc.
- [ ] All tests pass.
- [ ] `/update` has been run after changes.

## Steps

1. Inventory every call to bincode across all crates. Confirm the set of file formats affected.
2. Prototype a bincode 2.x migration in a branch to estimate work.
3. Make the call; document in decisions log.
4. Implement.

## Notes

See `needs-review/002_bincode_migration.md` for options.
