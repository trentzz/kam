# Task 026: Wire Up `kam pathfind` Subcommand

## Location
`kam/src/commands/pathfind.rs`

## What to implement

Connect the `kam pathfind` CLI subcommand: read k-mer index → read targets → build graph → walk paths → score → write output + QC.

## Interface

```rust
use crate::cli::PathfindArgs;

pub fn run_pathfind(args: PathfindArgs) -> Result<(), Box<dyn std::error::Error>>;
```

## Logic

1. Read k-mer index from bincode
2. Read target sequences from FASTA
3. For each target:
   a. Validate anchors
   b. Build de Bruijn graph from relevant k-mers
   c. Walk paths from start to end anchor
   d. Score paths
4. Collect all scored paths
5. Write to output using bincode
6. Build and write `PathfindQc`

## Tests

Integration test with temp files.

## Definition of done

- `cargo test -p kam` passes
- `cargo clippy -p kam -- -D warnings` passes
