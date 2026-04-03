# Task 028: Wire Up `kam run` Full Pipeline

## Location
`kam/src/commands/run.rs`

## What to implement

Connect the `kam run` CLI subcommand: runs the full pipeline (assemble → index → pathfind → call) in a single process with zero-copy data passing between stages.

## Interface

```rust
use crate::cli::RunArgs;

pub fn run_pipeline(args: RunArgs) -> Result<(), Box<dyn std::error::Error>>;
```

## Logic

1. Create output directory
2. **Assemble:** read FASTQ pairs → parse → assemble molecules (in memory, no bincode needed)
3. **Index:** build allowlist from targets → extract k-mers from molecules (in memory)
4. **Pathfind:** for each target, build graph → walk paths → score (in memory)
5. **Call:** call variants from scored paths
6. Write final variant calls to output file in requested format
7. Write QC JSON files for each stage to output directory
8. Print summary to stderr

This is the zero-copy hot path — data flows as Rust structs between stages without serialization.

## Tests

Integration test: create small temp FASTQ and FASTA files, run full pipeline, verify variant output file and QC JSONs exist.

## Definition of done

- `cargo test -p kam` passes
- `cargo clippy -p kam -- -D warnings` passes
- `kam run --r1 R1.fq --r2 R2.fq --targets panel.fa --output-dir results/` runs end-to-end
