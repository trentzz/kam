# Task 023: QC JSON Output Per Stage

## Location
`kam-core/src/qc.rs`

## What to implement

Structured QC JSON output for each pipeline stage, used by Nextflow for validation between stages.

## Interface

```rust
use std::path::Path;
use serde::Serialize;

/// QC output from the assembly stage
#[derive(Debug, Clone, Serialize)]
pub struct AssemblyQc {
    pub stage: String,       // "molecule_assembly"
    pub version: String,
    pub n_input_read_pairs: u64,
    pub n_molecules: u64,
    pub n_duplex: u64,
    pub n_simplex_fwd: u64,
    pub n_simplex_rev: u64,
    pub n_singletons: u64,
    pub duplex_fraction: f64,
    pub mean_family_size: f64,
    pub n_dropped_reads: u64,
    pub passed: bool,
}

/// QC output from the indexing stage
#[derive(Debug, Clone, Serialize)]
pub struct IndexQc {
    pub stage: String,       // "kmer_indexing"
    pub version: String,
    pub n_target_kmers: u64,
    pub n_kmers_observed: u64,
    pub mean_molecule_depth: f64,
    pub passed: bool,
}

/// QC output from the pathfinding stage
#[derive(Debug, Clone, Serialize)]
pub struct PathfindQc {
    pub stage: String,       // "graph_walking"
    pub version: String,
    pub n_targets_queried: u64,
    pub n_targets_with_variants: u64,
    pub n_anchors_non_unique: u64,
    pub passed: bool,
}

/// QC output from the calling stage
#[derive(Debug, Clone, Serialize)]
pub struct CallQc {
    pub stage: String,       // "variant_calling"
    pub version: String,
    pub n_variants_called: u64,
    pub n_pass: u64,
    pub n_filtered: u64,
    pub passed: bool,
}

/// Write any QC struct to a JSON file
pub fn write_qc<T: Serialize>(path: &Path, qc: &T) -> Result<(), Box<dyn std::error::Error>>;
```

## Dependencies

Add `serde_json` to kam-core's Cargo.toml.

## Tests required

1. write_qc produces valid JSON
2. AssemblyQc serializes with all fields
3. Deserialized JSON has correct stage name
4. duplex_fraction computed correctly

## Definition of done

- `cargo test -p kam-core` passes
- `cargo clippy -p kam-core -- -D warnings` passes
