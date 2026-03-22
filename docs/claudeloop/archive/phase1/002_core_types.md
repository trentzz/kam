# Task 002: kam-core Shared Types

## Location
`kam-core/src/lib.rs` and submodules

## What to implement

Define the core data types that all other crates depend on. These types are the foundation of the entire project. Implement exactly the types specified below.

## Interface (implement exactly this)

```rust
// kam-core/src/lib.rs
pub mod molecule;
pub mod kmer;
pub mod chemistry;
pub mod error;

// kam-core/src/molecule.rs
#[derive(Debug, Clone, PartialEq)]
pub struct ConsensusRead {
    pub sequence: Vec<u8>,
    pub per_base_error_prob: Vec<f32>,
    pub per_base_strand_support: Vec<(u8, u8)>,  // (fwd_count, rev_count)
    pub family_size: (u8, u8),                    // (fwd_reads, rev_reads)
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CanonicalUmiPair {
    pub umi_a: [u8; 5],  // lexicographically smaller
    pub umi_b: [u8; 5],  // lexicographically larger
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FamilyType {
    Duplex,         // fwd > 0 && rev > 0
    SimplexFwd,     // fwd > 0 && rev == 0
    SimplexRev,     // fwd == 0 && rev > 0
    Singleton,      // total == 1
}

// kam-core/src/kmer.rs
#[derive(Debug, Clone, Default)]
pub struct MoleculeEvidence {
    pub n_molecules: u32,
    pub n_duplex: u32,
    pub n_simplex_fwd: u32,
    pub n_simplex_rev: u32,
    pub min_base_error_prob: f32,
    pub mean_base_error_prob: f32,
}

pub trait KmerIndex {
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);
    fn contains(&self, kmer: u64) -> bool;
    fn molecule_count(&self, kmer: u64) -> u32;
}

// kam-core/src/chemistry.rs
#[derive(Debug, Clone)]
pub struct ReadStructure {
    pub umi_length: usize,
    pub skip_length: usize,
    // template is the remainder
}

impl ReadStructure {
    pub fn twist_umi_duplex() -> Self {
        Self { umi_length: 5, skip_length: 2 }
    }

    pub fn template_start(&self) -> usize {
        self.umi_length + self.skip_length
    }
}

// kam-core/src/error.rs
use thiserror::Error;

#[derive(Error, Debug)]
pub enum KamError {
    #[error("Read too short: expected at least {expected} bases, got {actual}")]
    ReadTooShort { expected: usize, actual: usize },
    #[error("Low quality UMI bases")]
    LowQualityUmi,
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}
```

## Tests required

1. `CanonicalUmiPair::new()` with umi_r1 < umi_r2 — umi_a should be umi_r1
2. `CanonicalUmiPair::new()` with umi_r1 > umi_r2 — umi_a should be umi_r2 (swapped)
3. `CanonicalUmiPair::new()` with umi_r1 == umi_r2 — should handle gracefully
4. `ReadStructure::twist_umi_duplex()` returns umi_length=5, skip_length=2
5. `ReadStructure::template_start()` returns 7 for Twist
6. `MoleculeEvidence::default()` has all zeros
7. `FamilyType` classification from family size tuples

## Definition of done

- `cargo test -p kam-core` passes
- `cargo clippy -p kam-core -- -D warnings` passes
- All public types have doc comments
- All types derive Debug, Clone at minimum
