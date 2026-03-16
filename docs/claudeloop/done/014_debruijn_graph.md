# Task 014: De Bruijn Graph Construction

## Location
`kam-pathfind/src/graph.rs`

## What to implement

Build a de Bruijn graph from a KmerIndex. Nodes are k-mers, edges connect k-mers that overlap by k-1 bases.

## Interface

```rust
use std::collections::HashMap;
use kam_core::kmer::{KmerIndex, MoleculeEvidence};

/// A de Bruijn graph built from a k-mer index.
/// Nodes are canonical k-mers. Edges connect k-mers with k-1 overlap.
#[derive(Debug)]
pub struct DeBruijnGraph {
    k: usize,
    /// Adjacency list: kmer → list of successor kmers
    successors: HashMap<u64, Vec<u64>>,
    /// Adjacency list: kmer → list of predecessor kmers
    predecessors: HashMap<u64, Vec<u64>>,
    /// Number of nodes
    n_nodes: usize,
    /// Number of edges
    n_edges: usize,
}

impl DeBruijnGraph {
    /// Build graph from a k-mer index. Only includes k-mers present in the index.
    pub fn from_index(index: &dyn KmerIndex, k: usize, all_kmers: &[u64]) -> Self;

    /// Get successors of a k-mer (k-mers that can follow it)
    pub fn successors(&self, kmer: u64) -> &[u64];

    /// Get predecessors of a k-mer
    pub fn predecessors(&self, kmer: u64) -> &[u64];

    /// Check if a k-mer is in the graph
    pub fn contains(&self, kmer: u64) -> bool;

    pub fn n_nodes(&self) -> usize;
    pub fn n_edges(&self) -> usize;
    pub fn k(&self) -> usize;
}
```

## How edges work

Two k-mers A and B are connected (A → B) if the last k-1 bases of A equal the first k-1 bases of B. In encoded form:
- suffix of A (k-1 bases) = `A & mask(k-1)` (lower 2*(k-1) bits)
- prefix of B (k-1 bases) = `B >> 2` (shift right by 2 to drop last base)
- Edge exists if suffix(A) == prefix(B)

Build the graph by:
1. For each k-mer, compute its (k-1)-prefix and (k-1)-suffix
2. Group k-mers by prefix → these are potential successors of anything with matching suffix
3. For each k-mer, look up its suffix in the prefix groups to find successors

## Tests required

1. Two overlapping k-mers form an edge
2. Non-overlapping k-mers have no edge
3. Linear sequence of k-mers forms a chain
4. Branch point (SNV) creates two successors
5. n_nodes and n_edges correct
6. Empty index → empty graph
7. Self-loop k-mer (e.g., AAAA with k=4) handled
8. successors and predecessors are consistent (A→B implies B has A as predecessor)

## Definition of done

- `cargo test -p kam-pathfind` passes
- `cargo clippy -p kam-pathfind -- -D warnings` passes
- Doc comments with examples
