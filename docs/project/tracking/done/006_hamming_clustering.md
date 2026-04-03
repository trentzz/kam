# Task 006: Hamming Distance UMI Clustering

## Location
`kam-assemble/src/clustering.rs`

## What to implement

Group `CanonicalUmiPair`s that are within Hamming distance 1 of each other. This handles sequencing errors in UMI bases — e.g., `ACGTA+TGCAT` and `ACGTT+TGCAT` should be grouped together if only 1 base differs in the combined UMI pair.

## Context

Twist uses 5bp random UMIs with no whitelist for error correction. The only way to correct UMI sequencing errors is Hamming distance clustering. With 1,024 possible 5bp UMIs, the expected number of UMIs within Hamming distance 1 of any given UMI is ~15 (5 positions × 3 alternatives). This is dense enough that clustering must be careful to avoid merging genuinely different molecules.

## Interface

```rust
use kam_core::molecule::CanonicalUmiPair;

/// Compute Hamming distance between two canonical UMI pairs.
/// Compares umi_a to umi_a and umi_b to umi_b, returns total mismatches.
pub fn umi_pair_hamming_distance(a: &CanonicalUmiPair, b: &CanonicalUmiPair) -> u32;

/// Group a set of canonical UMI pairs by Hamming distance.
/// Uses directional clustering: higher-count UMIs absorb lower-count neighbours
/// within the specified max_distance.
///
/// Returns: Vec of groups, where each group is a Vec of indices into the input.
pub fn cluster_umi_pairs(
    pairs: &[(CanonicalUmiPair, u32)],  // (umi_pair, read_count)
    max_distance: u32,                   // typically 1
) -> Vec<Vec<usize>>;
```

## Tests required

1. Hamming distance of identical pairs is 0
2. Hamming distance of pairs differing in 1 base of umi_a is 1
3. Hamming distance of pairs differing in 1 base of each UMI is 2
4. Clustering with max_distance=0 produces one group per unique pair
5. Clustering merges pairs within distance 1
6. Directional: high-count UMI absorbs low-count neighbour, not vice versa
7. Transitive clustering doesn't chain too far (A→B→C shouldn't merge A and C if A↔C distance > max)
8. Single-element input returns single group

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
- Doc comments with examples
