# CONFIG-005: Integration Tests and Documentation

**Epic**: CONFIG-TOML (docs/claudetracking/overallplans/CONFIG-TOML.md)
**Priority**: critical
**Depends on**: config_003_wire_pipeline.md, config_004_example_configs.md
**Status**: todo

## Goal

Write integration tests that confirm config file and CLI flag routes produce
identical output, and add doc comments to all public config types. Done looks
like: a `tests/config_integration.rs` test that runs the pipeline twice (once
with `--config`, once with equivalent CLI flags) and asserts the VCF output is
identical; all public config structs and fields have doc comments.

## Steps

1. Write `tests/config_integration.rs` (or add to an existing integration test
   file in `kam/tests/`):
   - Generate a minimal synthetic FASTQ pair (or reuse an existing fixture).
   - Run 1: `kam run --min-confidence 0.95 --r1 ... --r2 ...`
   - Run 2: write a temp TOML with `min_confidence = 0.95`, run
     `kam run --config temp.toml --r1 ... --r2 ...`
   - Assert output VCFs are identical (line by line, ignoring timestamp
     headers).

2. Add a second test for CLI override:
   - Config file sets `min_confidence = 0.90`.
   - CLI flag sets `--min-confidence 0.99`.
   - Assert the output matches a run with `--min-confidence 0.99` only.

3. Review all public items in `kam/src/config.rs` and add `///` doc comments
   to any that are missing them. Every struct and every field must have a
   comment.

4. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

5. Run `cargo doc --no-deps` and check that `KamConfig` and all sub-structs
   appear in the generated docs with readable descriptions.

## Notes

- Use `tempfile::tempdir()` to write the temp TOML in the integration test so
  the test is self-contained and leaves no artefacts.
- If the VCF includes a timestamp or run-ID header line, strip or ignore it
  before comparison.
- Doc comments on config fields should be written for end users (the person
  editing the TOML), not for Rust developers.
