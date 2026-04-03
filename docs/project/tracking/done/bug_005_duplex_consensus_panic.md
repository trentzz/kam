# BUG-005: assert_eq! panic in duplex_consensus (library code)

**Epic**: none (standalone bug fix)
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Replace the `assert_eq!` in `duplex_consensus` in `kam-assemble/src/consensus.rs`
with proper error handling. Done looks like: the function returns `Option` or
`Result` for the length-mismatch case; no panicking assertions remain in library
code.

## Steps

1. Read `kam-assemble/src/consensus.rs`, focusing on `duplex_consensus` (around
   line 311).

2. Replace:
   ```rust
   assert_eq!(fwd_bases.len(), rev_bases.len());
   ```
   with a graceful error path. Options:
   - Change return type to `Option<ConsensusRead>` and return `None` with an
     `eprintln!` warning if lengths differ.
   - Change return type to `Result<ConsensusRead, ConsensusError>`.
   Prefer `Option` if the function currently returns `ConsensusRead` — it is
   a simpler change and a length mismatch is genuinely unexpected rather than
   a routine error.

3. Update all callers to handle the `Option`/`Result`.

4. Add a unit test that passes mismatched-length SSC inputs and confirms no
   panic occurs (the function returns `None` or an error).

5. Move this task file to done and update Status.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding NEW-003.
- Project rule: "No `unwrap()` in library code." Panicking asserts are
  equivalent to unwrap.
