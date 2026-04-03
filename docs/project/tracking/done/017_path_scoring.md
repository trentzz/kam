# Task 017: Path Scoring with MoleculeEvidence

## Location
`kam-pathfind/src/score.rs`

## What to implement

Score each graph path using the MoleculeEvidence from the k-mer index. Determines how many independent molecules (especially duplex) support each path, and identifies the weakest point.

## Interface

```rust
use kam_core::kmer::{KmerIndex, MoleculeEvidence};
use crate::walk::GraphPath;

/// A scored path with molecule-level evidence
#[derive(Debug, Clone)]
pub struct ScoredPath {
    pub path: GraphPath,
    pub aggregate_evidence: PathEvidence,
    pub weakest_kmer: WeakestKmer,
    pub is_reference: bool,         // true if this is the reference/wildtype path
}

/// Aggregate evidence across all k-mers in a path
#[derive(Debug, Clone)]
pub struct PathEvidence {
    pub min_molecules: u32,         // minimum n_molecules across all k-mers in path
    pub mean_molecules: f32,        // mean n_molecules
    pub min_duplex: u32,
    pub mean_duplex: f32,
    pub total_simplex_fwd: u32,
    pub total_simplex_rev: u32,
    pub mean_error_prob: f32,       // mean of mean_base_error_prob across k-mers
}

/// The weakest point in the path (lowest molecule support)
#[derive(Debug, Clone)]
pub struct WeakestKmer {
    pub kmer: u64,
    pub position_in_path: usize,
    pub evidence: MoleculeEvidence,
}

/// Score a single path against the k-mer index
pub fn score_path(
    path: &GraphPath,
    index: &dyn KmerIndex,
) -> ScoredPath;

/// Score multiple paths and sort by evidence strength (strongest first)
pub fn score_and_rank_paths(
    paths: Vec<GraphPath>,
    index: &dyn KmerIndex,
    reference_seq: &[u8],
    k: usize,
) -> Vec<ScoredPath>;
```

## Algorithm

For each path:
1. Look up MoleculeEvidence for every k-mer in the path
2. Compute aggregate: min/mean of n_molecules, n_duplex across all k-mers
3. Find the weakest k-mer (lowest n_molecules) — this is the bottleneck
4. Determine if this path matches the reference sequence (is_reference = true)

For ranking:
1. Score all paths
2. Sort by `min_molecules` descending (strongest evidence first)
3. Mark the reference path

## Tests required

1. Single k-mer path → evidence equals that k-mer's evidence
2. Multi-kmer path → min_molecules is the minimum across k-mers
3. Weakest kmer identified correctly
4. Reference path detected by sequence comparison
5. Paths sorted by evidence strength
6. Path with k-mer missing from index → evidence has 0 molecules
7. Empty path → handle gracefully

## Definition of done

- `cargo test -p kam-pathfind` passes
- `cargo clippy -p kam-pathfind -- -D warnings` passes
- Doc comments with examples
