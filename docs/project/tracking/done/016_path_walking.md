# Task 016: Path Walking Between Anchors

## Location
`kam-pathfind/src/walk.rs`

## What to implement

Walk all paths through the de Bruijn graph from a start anchor k-mer to an end anchor k-mer. Each path represents a possible sequence (wildtype or variant) supported by the k-mer evidence.

## Interface

```rust
use kam_core::kmer::{KmerIndex, MoleculeEvidence};
use crate::graph::DeBruijnGraph;

/// A path through the graph from start to end anchor
#[derive(Debug, Clone)]
pub struct GraphPath {
    pub kmers: Vec<u64>,           // ordered k-mers in the path
    pub sequence: Vec<u8>,         // reconstructed sequence from the path
    pub length: usize,             // number of k-mers
}

/// Configuration for path walking
#[derive(Debug, Clone)]
pub struct WalkConfig {
    pub max_path_length: usize,    // maximum k-mers in a path (prevents infinite loops)
    pub max_paths: usize,          // maximum paths to enumerate (prevents combinatorial explosion)
}

impl Default for WalkConfig {
    fn default() -> Self {
        Self {
            max_path_length: 500,
            max_paths: 100,
        }
    }
}

/// Find all paths from start_kmer to end_kmer in the graph.
/// Uses BFS/DFS with depth limit.
pub fn walk_paths(
    graph: &DeBruijnGraph,
    start_kmer: u64,
    end_kmer: u64,
    config: &WalkConfig,
) -> Vec<GraphPath>;

/// Reconstruct the DNA sequence from an ordered list of k-mers.
/// Each successive k-mer overlaps the previous by k-1 bases.
pub fn reconstruct_sequence(kmers: &[u64], k: usize) -> Vec<u8>;
```

## Algorithm

BFS from start_kmer:
1. Initialize queue with `[start_kmer]`
2. For each partial path, extend by each successor of the last k-mer
3. If last k-mer == end_kmer, path is complete → collect it
4. If path length > max_path_length, abandon this branch
5. If total collected paths >= max_paths, stop
6. Track visited k-mers per path to avoid cycles

## Tests required

1. Linear path (no branches) → one path found
2. SNV creates two paths (ref and alt)
3. Deletion creates a shorter alternative path
4. Insertion creates a longer alternative path
5. No path exists → empty result
6. max_path_length limits long paths
7. max_paths limits total paths found
8. Cycle in graph doesn't cause infinite loop
9. reconstruct_sequence from k-mers matches original sequence
10. Start == end (zero-length variant) → single path of one k-mer

## Definition of done

- `cargo test -p kam-pathfind` passes
- `cargo clippy -p kam-pathfind -- -D warnings` passes
- Doc comments with examples
