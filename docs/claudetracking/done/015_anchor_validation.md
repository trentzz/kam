# Task 015: Anchor Validation

## Location
`kam-pathfind/src/anchor.rs`

## What to implement

Validate that anchor k-mers (the start and end of a target sequence) are sufficiently unique in the index to support reliable path walking. Non-unique anchors in repeat regions produce spurious paths.

## Interface

```rust
use kam_core::kmer::{KmerIndex};
use crate::encode::{encode_kmer, canonical};

/// Result of anchor validation for a target
#[derive(Debug, Clone)]
pub struct AnchorValidation {
    pub start_kmer: u64,
    pub end_kmer: u64,
    pub start_count: u32,        // molecule count of start anchor in index
    pub end_count: u32,
    pub start_unique: bool,      // true if count < threshold
    pub end_unique: bool,
    pub warning: Option<String>,
}

/// Default threshold for anchor uniqueness
pub const DEFAULT_ANCHOR_THRESHOLD: u32 = 100;

/// Validate anchor k-mers for a target sequence.
/// Extracts the first and last k-mer from the target,
/// checks their molecule counts in the index.
pub fn validate_anchors(
    target_seq: &[u8],
    k: usize,
    index: &dyn KmerIndex,
    threshold: u32,
) -> Option<AnchorValidation>;
```

## Logic

1. Extract first k-mer from target sequence (positions 0..k)
2. Extract last k-mer from target sequence (positions len-k..len)
3. Canonicalize both
4. Look up molecule counts in the index
5. Flag as non-unique if count > threshold
6. Generate warning message if either anchor is non-unique

Returns None if target is shorter than k.

## Tests required

1. Unique anchors (low count) → both start_unique and end_unique true, no warning
2. Non-unique start anchor → start_unique false, warning present
3. Non-unique end anchor → end_unique false, warning present
4. Both non-unique → both false, warning mentions both
5. Target shorter than k → returns None
6. Threshold of 0 → everything is non-unique
7. Anchor at exactly threshold → unique (threshold is exclusive upper bound)

## Definition of done

- `cargo test -p kam-pathfind` passes
- `cargo clippy -p kam-pathfind -- -D warnings` passes
- Doc comments with examples
