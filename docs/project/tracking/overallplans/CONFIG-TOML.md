# CONFIG-TOML: TOML Configuration File Support

**Status**: todo
**Priority**: critical
**Branch**: epic/CONFIG-TOML

## Goal

Add a `--config` flag to `kam run` that loads pipeline parameters from a TOML
file. Implement a layered resolution strategy: built-in Rust defaults, then
config file, then CLI flag overrides (CLI always wins). This removes the need
to pass every parameter on the command line and enables reproducible runs by
committing a config file alongside data.

## Motivation

The `kam run` command currently requires all parameters to be passed as CLI
flags. A typical tumour-informed run with custom thresholds requires 10+ flags.
TOML config files allow a single `--config my-run.toml` to set all parameters,
with CLI flags reserved for per-run overrides. This is the standard pattern for
bioinformatics pipelines and makes Nextflow integration straightforward.

## Design

### TOML schema

Sections map to pipeline stages:

```toml
[chemistry]
umi_length = 5
skip_length = 2

[assembly]
max_hamming = 1
min_reads_per_molecule = 2

[indexing]
kmer_size = 31
allowlist_path = ""

[calling]
min_confidence = 0.99
min_alt_molecules = 2
sv_min_confidence = 0.95
sv_min_alt_molecules = 1
strand_bias_threshold = 0.01
sv_strand_bias_threshold = 1.0

[output]
format = "vcf"
tsv = false

[input]
r1 = ""
r2 = ""
targets = ""
target_variants = ""
```

All fields use `#[serde(default)]` so a partial config file is valid.

### Layered resolution

1. Start with Rust struct defaults (`Default::default()`).
2. If `--config path.toml` is given, deserialise and merge (overwrite fields
   that are present in the file).
3. For each CLI flag that is explicitly set, overwrite the corresponding field.

Step 3 requires distinguishing "user set this flag" from "this is the CLI
default". Use `Option<T>` wrappers for all overridable CLI args, or check
`ArgMatches::value_source`.

### File location

New file: `kam/src/config.rs`. The existing `ParserConfig`, `AssemblerConfig`,
`CallerConfig` types in downstream crates remain unchanged. The TOML config is
a front-end concern: `run.rs` reads it and constructs the per-crate configs
from it.

## Child tasks

| ID | File | Status |
|----|------|--------|
| CONFIG-001 | todo/config_001_toml_schema.md | todo |
| CONFIG-002 | todo/config_002_loader_merge.md | todo |
| CONFIG-003 | todo/config_003_wire_pipeline.md | todo |
| CONFIG-004 | todo/config_004_example_configs.md | todo |
| CONFIG-005 | todo/config_005_tests_docs.md | todo |

## Scope

- `kam/src/config.rs` — new file: TOML schema structs
- `kam/src/run.rs` — wire config into pipeline construction
- `kam/src/cli.rs` or `kam/src/main.rs` — add `--config` flag to `RunArgs`
- `Cargo.toml` (workspace or `kam/`) — add `toml` and `serde` dependencies
- `examples/` — two example config files

## Out of scope

- Changes to `kam-core` types
- Changes to per-crate config structs (`ParserConfig`, `CallerConfig`, etc.)
- Nextflow DSL2 integration (a separate epic)
