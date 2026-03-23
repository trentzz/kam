# REFACTOR-001: Deduplicate CallerConfig construction in call.rs and run.rs

**Epic**: none (standalone refactor)
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

Extract the shared `CallerConfig` construction logic from
`kam/src/commands/call.rs` and `kam/src/commands/run.rs` into a single
function or builder. Done looks like: the construction appears once; both
commands use all fields including `sv_min_confidence` and
`sv_min_alt_molecules`; no fields are silently missing from either path.

## Steps

1. Read `kam/src/commands/call.rs` and `kam/src/commands/run.rs`, focusing
   on the `CallerConfig` construction blocks.

2. Identify which fields are present in `run.rs` but missing from `call.rs`
   (currently: `sv_min_confidence`, `sv_min_alt_molecules`).

3. Extract a shared function `caller_config_from_args(args: &impl SharedArgs)
   -> CallerConfig` or implement a trait for the shared CLI arg fields.

4. Use the shared function in both `call.rs` and `run.rs`.

5. Also deduplicate `parse_output_formats` and `format_extension` (codebase
   review finding 4.2) — move both to a shared module (e.g.
   `kam/src/output.rs`).

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review findings 4.2 and 4.4.
- The missing `sv_min_confidence` in `call.rs` is a correctness issue: the
  standalone `kam call` subcommand would silently use the default even if a
  custom value was intended.
