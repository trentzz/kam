# QUALITY-005: Remove dead else-if branch in targeting.rs

**Epic**: QUALITY
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

`kam-call/src/targeting.rs::apply_target_filter_with_seq_fallback` has an unreachable `else if key.is_none()` branch after an `if let Some(...) = key` arm. Remove it or replace with an explicit `else {}` with a clear comment.

## Success Criteria

- [ ] The dead branch is removed or replaced with a documented explicit `else {}`.
- [ ] All existing sequence-fallback tests still pass.
- [ ] `cargo clippy` reports no redundant-pattern-matching warnings on targeting.rs.
- [ ] `/update` has been run after changes.

## Steps

1. Read `apply_target_filter_with_seq_fallback` body in targeting.rs.
2. Remove or replace the dead branch.
