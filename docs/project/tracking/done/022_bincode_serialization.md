# Task 022: Bincode Serialization for Inter-Stage Data

## Location
`kam-core/src/serialize.rs`

## What to implement

Serialize and deserialize the main data structures to/from bincode files for Nextflow inter-stage handoff. Each file has a small header with type tag and version for `kam explore` compatibility.

## Interface

```rust
use std::path::Path;
use std::io::{Read, Write};
use serde::{Serialize, Deserialize};

/// File type tag for bincode files
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum FileType {
    Molecules,
    KmerIndex,
    ScoredPaths,
    VariantCalls,
}

/// Header written at the start of every kam bincode file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileHeader {
    pub magic: [u8; 4],          // b"KAM\0"
    pub version: u32,             // format version (1)
    pub file_type: FileType,
    pub record_count: u64,
}

/// Write a header + records to a bincode file
pub fn write_bincode<T: Serialize>(
    path: &Path,
    file_type: FileType,
    records: &[T],
) -> Result<(), Box<dyn std::error::Error>>;

/// Read header from a bincode file (without loading records)
pub fn read_header(path: &Path) -> Result<FileHeader, Box<dyn std::error::Error>>;

/// Read header + all records from a bincode file
pub fn read_bincode<T: for<'de> Deserialize<'de>>(
    path: &Path,
) -> Result<(FileHeader, Vec<T>), Box<dyn std::error::Error>>;
```

## Dependencies

Add `bincode = "1"` to [workspace.dependencies] and kam-core's Cargo.toml.

All serializable types (Molecule, ConsensusRead, MoleculeEvidence, etc.) need `#[derive(Serialize, Deserialize)]`. Add these derives to the existing types in kam-core — they already have `serde` as a dependency.

## Tests required

1. Write + read roundtrip produces identical data
2. Header magic bytes are correct
3. Header record_count matches input
4. read_header alone works without loading records
5. Wrong magic bytes returns error
6. Empty records list writes valid file with count 0

## Definition of done

- `cargo test -p kam-core` passes
- `cargo clippy -p kam-core -- -D warnings` passes
- All core types have Serialize/Deserialize derives
