# BUG-001: FamilyType::from_family_size((0,0)) returns Duplex

**Epic**: none (standalone bug fix)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Fix `FamilyType::from_family_size` so that a family with zero reads on both
strands does not return `Duplex`. Done looks like: the catch-all arm is removed
or narrowed; `(0, 0)` returns `Singleton` or a dedicated `Empty` variant; a
test covers this case.

## Steps

1. Read `kam-core/src/molecule.rs`, focusing on `FamilyType` and
   `from_family_size`.

2. Determine whether a `(0, 0)` family can actually be constructed in practice
   (search for callers of `from_family_size` or constructors of the family
   size fields). If it cannot occur, add a `debug_assert!` and document why.

3. If it can occur, fix the match arm. Options:
   - Map `(0, 0)` to `Singleton` (most conservative — matches a single-read
     molecule that has been downgraded).
   - Add an `Empty` variant to `FamilyType` if downstream code needs to
     distinguish zero-read molecules.

4. Remove the stale doc comment referencing "task-002 / task-004" (finding
   7.1 from codebase review).

5. Add a unit test: `assert_eq!(FamilyType::from_family_size(0, 0), ...)`.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding 1.1: the catch-all `_ => FamilyType::Duplex`
  at line 153 of `molecule.rs` is the bug.
- Changing `FamilyType` is a `kam-core` type change. Confirm before committing
  that adding a variant (if needed) does not break any match expressions
  elsewhere.
