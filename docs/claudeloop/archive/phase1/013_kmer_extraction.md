# Task 013: K-mer Extraction from Consensus Reads

## Location
`kam-index/src/extract.rs`

## What to implement

Extract k-mers from consensus molecules and accumulate MoleculeEvidence in a KmerIndex. This bridges kam-assemble output to kam-index.

## Interface

```rust
use kam_core::kmer::{KmerIndex, MoleculeEvidence};
use kam_core::molecule::{FamilyType};
use crate::encode::{canonical_with_strand, KmerIterator, KmerStrand};

/// Information about a consensus read needed for k-mer extraction
#[derive(Debug, Clone)]
pub struct ConsensusReadInfo {
    pub sequence: Vec<u8>,
    pub per_base_error_prob: Vec<f32>,
    pub family_type: FamilyType,
}

/// Extract k-mers from a consensus read and update the index.
/// For each k-mer position:
///   - Canonicalize the k-mer
///   - Determine strand
///   - Update MoleculeEvidence: increment molecule count, strand counts,
///     update min/mean error probability from the k-mer's base qualities
pub fn extract_and_index(
    read: &ConsensusReadInfo,
    k: usize,
    index: &mut dyn KmerIndex,
);

/// Extract k-mers from multiple consensus reads
pub fn extract_all(
    reads: &[ConsensusReadInfo],
    k: usize,
    index: &mut dyn KmerIndex,
);

/// Compute the mean error probability across bases in a k-mer
fn kmer_mean_error_prob(error_probs: &[f32], start: usize, k: usize) -> f32;
```

## Algorithm for each k-mer from a consensus read

1. Encode and canonicalize the k-mer
2. Get or create MoleculeEvidence entry in the index
3. Increment `n_molecules` by 1
4. Based on family_type:
   - Duplex → increment `n_duplex`
   - SimplexFwd → increment `n_simplex_fwd`
   - SimplexRev → increment `n_simplex_rev`
   - Singleton → increment `n_simplex_fwd` (treat as simplex)
5. Compute mean error prob across the k bases
6. Update `min_base_error_prob` if this is lower
7. Update `mean_base_error_prob` as running average

## Tests required

1. Single consensus read produces L-k+1 k-mers in index
2. Duplex read increments n_duplex
3. SimplexFwd read increments n_simplex_fwd
4. Two different reads with overlapping k-mers: n_molecules incremented for shared k-mers
5. min_base_error_prob tracks the minimum across all observations
6. K-mers with N bases are skipped
7. Empty sequence produces no k-mers
8. extract_all processes multiple reads correctly

## Definition of done

- `cargo test -p kam-index` passes
- `cargo clippy -p kam-index -- -D warnings` passes
- Doc comments with examples
