# Task 021: FASTQ I/O with needletail

## Location
`kam-assemble/src/io.rs`

## What to implement

Read paired FASTQ files (R1 and R2) using the `needletail` crate, parse each read pair through the parser, and return a Vec of ParsedReadPairs plus ParseStats.

## Interface

```rust
use std::path::Path;
use crate::parser::{ParsedReadPair, ParserConfig, ParseStats};

/// Read paired FASTQ files and parse all read pairs.
/// R1 and R2 must have the same number of records in the same order.
pub fn read_fastq_pairs(
    r1_path: &Path,
    r2_path: &Path,
    config: &ParserConfig,
) -> Result<(Vec<ParsedReadPair>, ParseStats), Box<dyn std::error::Error>>;
```

## Dependencies

Add `needletail = "0.6"` to [workspace.dependencies] in root Cargo.toml and to kam-assemble's Cargo.toml.

## Implementation

1. Open R1 and R2 with `needletail::parse_fastq_file`
2. Iterate both in lockstep
3. For each pair, call `parse_read_pair` with the sequences and qualities
4. Collect Ok results into Vec, increment ParseStats for drops
5. Error if files have different record counts

## Tests required

1. Reading a valid small FASTQ pair produces correct ParsedReadPairs (write temp files in test)
2. Mismatched record counts returns error
3. Empty files return empty vec with zero stats
4. ParseStats correctly counts drops

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
