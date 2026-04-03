# Task 024: Wire Up `kam assemble` Subcommand

## Location
`kam/src/commands/assemble.rs`

## What to implement

Connect the `kam assemble` CLI subcommand to the actual assemble pipeline: read FASTQ → parse → assemble molecules → write output + QC.

## Interface

```rust
use crate::cli::AssembleArgs;

pub fn run_assemble(args: AssembleArgs) -> Result<(), Box<dyn std::error::Error>>;
```

## Logic

1. Build `ParserConfig` from CLI args (chemistry, min_umi_quality, min_template_length)
2. Build `AssemblerConfig` from CLI args (min_family_size, consensus config)
3. Call `read_fastq_pairs(r1, r2, parser_config)` to get parsed read pairs
4. Call `assemble_molecules(read_pairs, assembler_config)` to get molecules
5. Write molecules to output path using bincode serialization
6. Build `AssemblyQc` from `AssemblyStats`
7. Write QC JSON to `{output_dir}/assembly_qc.json` (or derive path from --output)
8. Print summary to stderr

## Structure

Create a `kam/src/commands/` module directory with `mod.rs` and `assemble.rs`. Update `main.rs` to call `commands::assemble::run_assemble(args)` instead of the placeholder print.

## Tests

Integration test: create small temp FASTQ files, run `run_assemble`, verify output file exists and QC JSON is valid. Use `tempfile` crate for temp directories.

Add `tempfile = "3"` to kam's dev-dependencies.

## Definition of done

- `cargo test -p kam` passes
- `cargo clippy -p kam -- -D warnings` passes
- `kam assemble --r1 test.R1.fq --r2 test.R2.fq --output out.bin` actually runs
