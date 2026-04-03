# CONFIG-002: Config File Loader and CLI Override Merge

**Epic**: CONFIG-TOML (docs/claudetracking/overallplans/CONFIG-TOML.md)
**Priority**: critical
**Depends on**: config_001_toml_schema.md
**Status**: todo

## Goal

Add a `--config` flag to `RunArgs` and implement a `load_config` function that
reads a TOML file and merges CLI overrides on top. Done looks like: `kam run
--config my.toml` loads and uses the file; any explicitly set CLI flag
overrides the file value; omitted CLI flags fall back to the file; the file is
optional (omitting `--config` uses defaults).

## Steps

1. Read `kam/src/cli.rs` (or wherever `RunArgs` is defined) and
   `kam/src/run.rs` in full.

2. Add `--config` to `RunArgs`:
   ```rust
   #[arg(long, value_name = "FILE")]
   pub config: Option<PathBuf>,
   ```

3. Write `load_config(path: Option<&Path>) -> Result<KamConfig>` in
   `kam/src/config.rs`:
   - If `path` is `None`, return `KamConfig::default()`.
   - Read the file to a string; parse with `toml::from_str`. Return an
     `Err` with a clear message if the file is missing or malformed.

4. Write `merge_cli_overrides(config: &mut KamConfig, args: &RunArgs)`:
   - For every CLI flag that can override a config field, check
     `args.value_source(field) == ValueSource::CommandLine` (or use
     `Option<T>` fields and check `is_some()`).
   - If the flag was explicitly set, write the CLI value into `config`.
   - If not set, leave the config file value unchanged.

5. Write unit tests in `config.rs`:
   - `test_default_fallback`: no config file, no CLI flags → defaults.
   - `test_file_override`: config file sets `min_confidence = 0.90`,
     no CLI flag → config uses 0.90.
   - `test_cli_override`: config file sets `min_confidence = 0.90`,
     CLI sets `--min-confidence 0.99` → config uses 0.99.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- CLI flags that are already `Option<T>` (e.g. `--target-variants`) are easy:
  `if let Some(v) = args.target_variants { config.input.target_variants = Some(v); }`.
- For flags with non-optional primitives, either change them to `Option<T>` or
  use `ArgMatches::value_source`. Prefer `Option<T>` for clarity.
- Do not wire into `run_pipeline` yet; that is CONFIG-003.
