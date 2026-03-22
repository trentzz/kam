# Task 012: Allowlist Builder and FilteredKmerIndex

## Location
`kam-index/src/allowlist.rs`

## What to implement

Build an allowlist of k-mers from target sequences, and a FilteredKmerIndex wrapper that only stores k-mers in the allowlist. This is the key memory optimization: memory scales with panel size, not sequencing depth.

## Interface

```rust
use std::collections::HashSet;
use kam_core::kmer::{KmerIndex, MoleculeEvidence};
use crate::encode::{encode_kmer, canonical, KmerIterator};

/// Build an allowlist of canonical k-mers from target sequences.
/// Pads each target by `k - 1` bases on each side to avoid boundary gaps.
pub fn build_allowlist(targets: &[&[u8]], k: usize) -> HashSet<u64>;

/// A k-mer index that only stores k-mers present in an allowlist.
/// Wraps any inner KmerIndex implementation.
#[derive(Debug, Clone)]
pub struct FilteredKmerIndex<T: KmerIndex> {
    inner: T,
    allowlist: HashSet<u64>,
}

impl<T: KmerIndex> FilteredKmerIndex<T> {
    pub fn new(inner: T, allowlist: HashSet<u64>) -> Self;

    /// Number of allowed k-mers (the allowlist size)
    pub fn allowlist_size(&self) -> usize;

    /// Number of allowed k-mers that have been observed (have evidence)
    pub fn observed_count(&self) -> usize;

    /// Check if a k-mer is in the allowlist (regardless of whether it's been observed)
    pub fn is_allowed(&self, kmer: u64) -> bool;
}

impl<T: KmerIndex> KmerIndex for FilteredKmerIndex<T> {
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);  // no-op if not in allowlist
    fn contains(&self, kmer: u64) -> bool;
    fn molecule_count(&self, kmer: u64) -> u32;
}
```

## Key behaviour

- `build_allowlist` extracts all k-mers from each target, canonicalises them, and collects into a HashSet. It should handle the k-1 padding naturally — the caller provides sequences that already include flanking context.
- `FilteredKmerIndex::insert` silently ignores k-mers not in the allowlist (no error, just no-op). This is the filtering mechanism.
- `FilteredKmerIndex::get`, `contains`, `molecule_count` delegate to the inner index.

## Tests required

1. build_allowlist from a simple target sequence produces expected k-mers
2. build_allowlist canonicalises k-mers (forward and revcomp map to same entry)
3. FilteredKmerIndex insert + get works for allowed k-mers
4. FilteredKmerIndex insert is no-op for non-allowed k-mers
5. allowlist_size returns correct count
6. is_allowed true for allowed, false for others
7. Multiple targets: allowlist is union of all target k-mers
8. observed_count starts at 0, increases with inserts

## Definition of done

- `cargo test -p kam-index` passes
- `cargo clippy -p kam-index -- -D warnings` passes
- Doc comments with examples
