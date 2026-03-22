# Task 025: Wire Up `kam index` Subcommand

## Location
`kam/src/commands/index.rs`

## What to implement

Connect the `kam index` CLI subcommand: read molecules → read targets → build allowlist → extract k-mers → write index + QC.

## Interface

```rust
use crate::cli::IndexArgs;

pub fn run_index(args: IndexArgs) -> Result<(), Box<dyn std::error::Error>>;
```

## Logic

1. Read molecules from bincode input file
2. Read target sequences from FASTA file (use `needletail` to parse FASTA)
3. Build allowlist from targets with `build_allowlist`
4. Create `FilteredKmerIndex<HashKmerIndex>`
5. For each molecule, create `ConsensusReadInfo` and call `extract_and_index`
6. Write index to output using bincode
7. Build and write `IndexQc`

## Add needletail to kam's dependencies for FASTA reading.

## Tests

Integration test with temp files: write a small molecules bincode file and a small targets FASTA, run `run_index`, verify output exists.

## Definition of done

- `cargo test -p kam` passes
- `cargo clippy -p kam -- -D warnings` passes
