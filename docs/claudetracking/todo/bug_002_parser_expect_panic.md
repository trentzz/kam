# BUG-002: parser.rs uses .expect() in library code

**Epic**: none (standalone bug fix)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Replace `.expect()` / `.unwrap()` calls in `kam-assemble/src/parser.rs` with
proper error propagation. Done looks like: no `.expect()` or `.unwrap()` in
library code; all fallible operations use `?`; the function signatures return
`Result`.

## Steps

1. Read `kam-assemble/src/parser.rs` in full.

2. Identify all `.expect()` calls (specifically the UMI/skip extraction at
   lines ~215-221 and any others).

3. Replace each with `?` or `.map_err(|e| ...)` propagating through the
   function's `Result` return type.

4. Update function signatures if they currently return `()` instead of
   `Result<_, Error>`.

5. Update callers as needed to handle the returned `Result`.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding 2.1.
- Project rule: "No `unwrap()` in library code — use proper error handling
  with `thiserror`."
