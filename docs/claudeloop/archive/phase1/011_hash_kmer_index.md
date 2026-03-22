# Task 011: HashKmerIndex Implementation

## Location
`kam-index/src/hash_index.rs`

## What to implement

The primary in-memory k-mer index backed by HashMap, implementing the `KmerIndex` trait from kam-core.

## Interface

```rust
use std::collections::HashMap;
use kam_core::kmer::{KmerIndex, MoleculeEvidence};

/// In-memory k-mer index backed by HashMap.
/// Fast insertion and lookup, ~48 bytes per entry.
/// Best for targeted panel work.
#[derive(Debug, Clone, Default)]
pub struct HashKmerIndex {
    map: HashMap<u64, MoleculeEvidence>,
}

impl HashKmerIndex {
    pub fn new() -> Self;
    pub fn with_capacity(capacity: usize) -> Self;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;

    /// Get a mutable reference to evidence for a k-mer, inserting default if absent
    pub fn entry(&mut self, kmer: u64) -> &mut MoleculeEvidence;

    /// Iterate over all k-mers and their evidence
    pub fn iter(&self) -> impl Iterator<Item = (&u64, &MoleculeEvidence)>;
}

impl KmerIndex for HashKmerIndex {
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);
    fn contains(&self, kmer: u64) -> bool;
    fn molecule_count(&self, kmer: u64) -> u32;
}
```

## Tests required

1. New index is empty
2. Insert and get back same evidence
3. contains returns true after insert, false before
4. molecule_count returns n_molecules from evidence
5. molecule_count returns 0 for absent k-mer
6. entry() creates default if absent, returns existing if present
7. Multiple inserts to same k-mer: last write wins
8. len() and is_empty() correct
9. iter() yields all entries
10. with_capacity doesn't change behaviour, just pre-allocates

## Definition of done

- `cargo test -p kam-index` passes
- `cargo clippy -p kam-index -- -D warnings` passes
- Doc comments with examples
