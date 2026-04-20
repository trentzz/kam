# QUALITY-001: Replace unwrap/expect in library code with typed errors

**Epic**: QUALITY
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Bring non-test code in kam-core, kam-assemble, kam-index, kam-pathfind, kam-call, and kam-ml in line with CLAUDE.md rule: "No `unwrap()` in library code — use proper error handling with `thiserror`". Binary `main.rs`/command entry points may keep `?`-friendly `expect()` with documented invariants.

## Success Criteria

- [x] `grep -rn "\.unwrap()" kam-*/src/` returns only test code, const expressions, or calls annotated with a `// SAFETY:` comment explaining the infallible invariant.
- [x] `KamError` enum in `kam-core/src/error.rs` covers the new error paths.
- [x] `cargo clippy --all --all-targets -- -D warnings` passes.
- [x] All tests pass.
- [x] `/update` has been run after changes.

## Steps

1. Run the grep; classify each hit into: test code, infallible invariant, or genuine error path.
2. Add `KamError` variants for genuine error paths; propagate with `?`.
3. Annotate infallible call sites with `// SAFETY:` comments.

## Notes

Hottest files: `kam/src/commands/run.rs` (~35 hits), `kam-call/src/output.rs` (~58 hits), `kam-assemble/src/parser.rs` (~32 hits).
