# QUALITY-005: Remove dead else-if branch in targeting.rs

**Epic**: QUALITY
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

`kam-call/src/targeting.rs::apply_target_filter_with_seq_fallback` has an unreachable `else if key.is_none()` branch after an `if let Some(...) = key` arm. Remove it or replace with an explicit `else {}` with a clear comment.

## Success Criteria

- [x] The dead branch is removed or replaced with a documented explicit `else {}`. (branch was removed; code falls through cleanly to sequence-level fallback)
- [x] All existing sequence-fallback tests still pass.
- [x] `cargo clippy` reports no redundant-pattern-matching warnings on targeting.rs.
- [x] `/update` has been run after changes.

## Steps

1. Read `apply_target_filter_with_seq_fallback` body in targeting.rs.
2. Remove or replace the dead branch.
