# CONFIG-001: TOML Schema and Rust Structs

**Epic**: CONFIG-TOML (docs/claudetracking/overallplans/CONFIG-TOML.md)
**Priority**: critical
**Depends on**: none
**Status**: todo

## Goal

Create `kam/src/config.rs` with Rust structs for the TOML configuration schema
and add the `toml` dependency to the workspace. Done looks like: a new
`KamConfig` struct with sub-structs for each pipeline section, all fields using
`#[serde(default)]`, and the file compiling cleanly under `cargo build`.

## Steps

1. Read `kam/Cargo.toml` and `Cargo.toml` (workspace root) to understand
   current dependencies.

2. Add `toml = "1"` and confirm `serde` with `derive` feature is already
   present. Add to `[workspace.dependencies]` if using workspace inheritance,
   or to `kam/Cargo.toml` directly.

3. Create `kam/src/config.rs` with these structs:

   ```rust
   #[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
   #[serde(default)]
   pub struct KamConfig {
       pub chemistry: ChemistryConfig,
       pub assembly: AssemblyConfig,
       pub indexing: IndexingConfig,
       pub calling: CallingConfig,
       pub output: OutputConfig,
       pub input: InputConfig,
   }
   ```

   Sub-structs:
   - `ChemistryConfig`: `umi_length: usize`, `skip_length: usize`
   - `AssemblyConfig`: `max_hamming: u8`, `min_reads_per_molecule: u32`
   - `IndexingConfig`: `kmer_size: usize`, `allowlist_path: Option<PathBuf>`
   - `CallingConfig`: `min_confidence: f64`, `min_alt_molecules: u32`,
     `sv_min_confidence: f64`, `sv_min_alt_molecules: u32`,
     `strand_bias_threshold: f64`, `sv_strand_bias_threshold: f64`
   - `OutputConfig`: `tsv: bool`
   - `InputConfig`: `r1: Option<PathBuf>`, `r2: Option<PathBuf>`,
     `targets: Option<PathBuf>`, `target_variants: Option<PathBuf>`

4. Implement `Default` for each struct using the same values as the existing
   CLI defaults. Match `CallerConfig` defaults exactly.

5. Add a doc comment on every public field explaining what it controls.

6. Declare `mod config;` in `kam/src/main.rs` or `kam/src/lib.rs`.

7. Run `cargo build` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- Do not wire the config into `run_pipeline` yet; that is CONFIG-003.
- Do not add the `--config` CLI flag yet; that is CONFIG-002.
- Match default values to existing CLI arg defaults exactly to avoid silent
  behaviour changes when users migrate from flags to config files.
- `serde(default)` on the outer struct and each sub-struct ensures a partial
  TOML file (e.g. only `[calling]` section) is valid.
