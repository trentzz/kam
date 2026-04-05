# CONFIG-003: Wire KamConfig into run_pipeline

**Epic**: CONFIG-TOML (docs/claudetracking/overallplans/CONFIG-TOML.md)
**Priority**: critical
**Depends on**: config_002_loader_merge.md
**Status**: todo

## Goal

Replace the ad-hoc `ParserConfig`, `AssemblerConfig`, and `CallerConfig`
construction in `run.rs` with construction from a `KamConfig`. Done looks like:
`run_pipeline` accepts a `KamConfig`; the existing CLI-only path still works
(because `merge_cli_overrides` populates the config from flags); all tests pass.

## Steps

1. Read `kam/src/run.rs` in full to understand how the per-crate configs are
   currently constructed.

2. Change the signature of `run_pipeline` (or its internal logic) to accept
   `KamConfig` instead of (or in addition to) the raw `RunArgs`.

3. Construct each per-crate config from `KamConfig` fields:
   - `ParserConfig`: from `config.chemistry` and `config.assembly`
   - `AssemblerConfig`: from `config.assembly`
   - `IndexingConfig` (if separate): from `config.indexing`
   - `CallerConfig`: from `config.calling`

4. In `main.rs` (or wherever `run_pipeline` is called):
   - Call `load_config(args.config.as_deref())?`
   - Call `merge_cli_overrides(&mut config, &args)`
   - Pass `config` to `run_pipeline`

5. Verify that running without `--config` produces identical results to the
   previous behaviour. Do this by running the existing integration tests.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- The per-crate config structs (`ParserConfig`, `CallerConfig`, etc.) are
  unchanged. This task only changes how they are constructed.
- If `run_pipeline` currently takes `RunArgs` directly and uses many fields,
  consider creating a `RunConfig` struct in `run.rs` that is populated from
  `KamConfig + RunArgs`. This keeps `run_pipeline` testable without a full
  `RunArgs`.
- Do not change kam-core types.
