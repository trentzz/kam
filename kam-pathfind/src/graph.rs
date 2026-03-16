//! De Bruijn graph construction from a k-mer index.
//!
//! Nodes are canonical k-mers. Edges connect k-mers that overlap by k-1 bases.
//! An edge A → B exists when the last k-1 bases of A equal the first k-1 bases of B.
//!
//! In 2-bit encoded form:
//! - suffix of A (k-1 bases) = `A & mask(k-1)` (lower 2*(k-1) bits)
//! - prefix of B (k-1 bases) = `B >> 2` (shift right by 2 to drop last base)
//! - Edge exists when suffix(A) == prefix(B)

use std::collections::HashMap;

use kam_core::kmer::KmerIndex;

/// A de Bruijn graph built from a k-mer index.
///
/// Nodes are canonical k-mers encoded as `u64`. An edge A → B exists when
/// `suffix(A, k-1) == prefix(B, k-1)` in 2-bit encoded form.
///
/// # Example
///
/// ```
/// use kam_pathfind::graph::DeBruijnGraph;
/// use kam_index::HashKmerIndex;
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
/// use kam_index::encode::encode_kmer;
///
/// let mut index = HashKmerIndex::new();
/// let acg = encode_kmer(b"ACG").unwrap();
/// let cgt = encode_kmer(b"CGT").unwrap();
/// index.insert(acg, MoleculeEvidence { n_molecules: 1, ..Default::default() });
/// index.insert(cgt, MoleculeEvidence { n_molecules: 1, ..Default::default() });
/// let all_kmers = vec![acg, cgt];
/// let graph = DeBruijnGraph::from_index(&index, 3, &all_kmers);
/// assert!(graph.successors(acg).contains(&cgt));
/// assert_eq!(graph.n_nodes(), 2);
/// assert_eq!(graph.n_edges(), 1);
/// ```
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
    /// Build a de Bruijn graph from a k-mer index.
    ///
    /// Only includes k-mers present in the index (filtered from `all_kmers`).
    /// Edges are derived from 2*(k-1) overlap: suffix(A) == prefix(B).
    ///
    /// # Arguments
    ///
    /// - `index`: the k-mer index (used to verify membership)
    /// - `k`: k-mer length
    /// - `all_kmers`: full set of encoded k-mers to consider as candidate nodes
    ///
    /// # Example
    ///
    /// ```
    /// use kam_pathfind::graph::DeBruijnGraph;
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    /// use kam_index::encode::encode_kmer;
    ///
    /// let mut index = HashKmerIndex::new();
    /// let acg = encode_kmer(b"ACG").unwrap();
    /// let cgt = encode_kmer(b"CGT").unwrap();
    /// index.insert(acg, MoleculeEvidence::default());
    /// index.insert(cgt, MoleculeEvidence::default());
    /// let kmers = vec![acg, cgt];
    /// let graph = DeBruijnGraph::from_index(&index, 3, &kmers);
    /// assert_eq!(graph.n_nodes(), 2);
    /// assert_eq!(graph.n_edges(), 1);
    /// ```
    pub fn from_index(index: &dyn KmerIndex, k: usize, all_kmers: &[u64]) -> Self {
        // Filter to nodes that are actually in the index.
        let nodes: Vec<u64> = all_kmers
            .iter()
            .copied()
            .filter(|&km| index.contains(km))
            .collect();

        // Initialise adjacency maps with empty vecs for every node.
        let mut successors: HashMap<u64, Vec<u64>> =
            nodes.iter().map(|&km| (km, Vec::new())).collect();
        let mut predecessors: HashMap<u64, Vec<u64>> =
            nodes.iter().map(|&km| (km, Vec::new())).collect();

        // Group nodes by their (k-1)-prefix for efficient edge lookup.
        // prefix(B) = B >> 2  (drop the last base, keeps first k-1 bases)
        let mut by_prefix: HashMap<u64, Vec<u64>> = HashMap::new();
        for &km in &nodes {
            let prefix = km >> 2;
            by_prefix.entry(prefix).or_default().push(km);
        }

        // Bitmask for the lower 2*(k-1) bits: this is the (k-1)-suffix of A.
        let suffix_mask = if k > 1 {
            (1u64 << (2 * (k - 1))) - 1
        } else {
            0
        };

        let mut n_edges = 0usize;

        for &a in &nodes {
            let a_suffix = a & suffix_mask;
            // All B whose (k-1)-prefix matches a_suffix are successors of A.
            if let Some(candidates) = by_prefix.get(&a_suffix) {
                for &b in candidates {
                    successors.entry(a).or_default().push(b);
                    predecessors.entry(b).or_default().push(a);
                    n_edges += 1;
                }
            }
        }

        DeBruijnGraph {
            k,
            successors,
            predecessors,
            n_nodes: nodes.len(),
            n_edges,
        }
    }

    /// Get successors of a k-mer (k-mers that can directly follow it).
    ///
    /// Returns an empty slice if the k-mer is not in the graph.
    ///
    /// # Example
    ///
    /// ```
    /// use kam_pathfind::graph::DeBruijnGraph;
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    /// use kam_index::encode::encode_kmer;
    ///
    /// let mut index = HashKmerIndex::new();
    /// let acg = encode_kmer(b"ACG").unwrap();
    /// let cgt = encode_kmer(b"CGT").unwrap();
    /// index.insert(acg, MoleculeEvidence::default());
    /// index.insert(cgt, MoleculeEvidence::default());
    /// let graph = DeBruijnGraph::from_index(&index, 3, &[acg, cgt]);
    /// assert!(graph.successors(acg).contains(&cgt));
    /// let empty: &[u64] = &[];
    /// assert_eq!(graph.successors(999), empty);
    /// ```
    pub fn successors(&self, kmer: u64) -> &[u64] {
        self.successors.get(&kmer).map(Vec::as_slice).unwrap_or(&[])
    }

    /// Get predecessors of a k-mer (k-mers that can directly precede it).
    ///
    /// Returns an empty slice if the k-mer is not in the graph.
    ///
    /// # Example
    ///
    /// ```
    /// use kam_pathfind::graph::DeBruijnGraph;
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    /// use kam_index::encode::encode_kmer;
    ///
    /// let mut index = HashKmerIndex::new();
    /// let acg = encode_kmer(b"ACG").unwrap();
    /// let cgt = encode_kmer(b"CGT").unwrap();
    /// index.insert(acg, MoleculeEvidence::default());
    /// index.insert(cgt, MoleculeEvidence::default());
    /// let graph = DeBruijnGraph::from_index(&index, 3, &[acg, cgt]);
    /// assert!(graph.predecessors(cgt).contains(&acg));
    /// let empty: &[u64] = &[];
    /// assert_eq!(graph.predecessors(999), empty);
    /// ```
    pub fn predecessors(&self, kmer: u64) -> &[u64] {
        self.predecessors
            .get(&kmer)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Return `true` if `kmer` is a node in the graph.
    ///
    /// # Example
    ///
    /// ```
    /// use kam_pathfind::graph::DeBruijnGraph;
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    /// use kam_index::encode::encode_kmer;
    ///
    /// let mut index = HashKmerIndex::new();
    /// let acg = encode_kmer(b"ACG").unwrap();
    /// index.insert(acg, MoleculeEvidence::default());
    /// let graph = DeBruijnGraph::from_index(&index, 3, &[acg]);
    /// assert!(graph.contains(acg));
    /// assert!(!graph.contains(999));
    /// ```
    pub fn contains(&self, kmer: u64) -> bool {
        self.successors.contains_key(&kmer)
    }

    /// Return the number of nodes in the graph.
    pub fn n_nodes(&self) -> usize {
        self.n_nodes
    }

    /// Return the number of directed edges in the graph.
    pub fn n_edges(&self) -> usize {
        self.n_edges
    }

    /// Return the k-mer length used to build this graph.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Create a new, empty graph for k-mers of length `k`.
    ///
    /// Intended for testing and incremental construction with [`add_edge`].
    /// For production use, prefer [`from_index`].
    ///
    /// # Example
    /// ```
    /// use kam_pathfind::graph::DeBruijnGraph;
    ///
    /// let g = DeBruijnGraph::new(5);
    /// assert_eq!(g.k(), 5);
    /// assert_eq!(g.n_nodes(), 0);
    /// ```
    ///
    /// [`add_edge`]: DeBruijnGraph::add_edge
    /// [`from_index`]: DeBruijnGraph::from_index
    pub fn new(k: usize) -> Self {
        DeBruijnGraph {
            k,
            successors: HashMap::new(),
            predecessors: HashMap::new(),
            n_nodes: 0,
            n_edges: 0,
        }
    }

    /// Add a directed edge from `from` to `to`, inserting both nodes if absent.
    ///
    /// This method is primarily intended for test graph construction. In
    /// production, edges are derived automatically by [`from_index`] based on
    /// k-mer overlap.
    ///
    /// # Example
    /// ```
    /// use kam_pathfind::graph::DeBruijnGraph;
    ///
    /// let mut g = DeBruijnGraph::new(3);
    /// g.add_edge(10, 20);
    /// assert!(g.contains(10));
    /// assert!(g.contains(20));
    /// assert_eq!(g.successors(10), &[20]);
    /// ```
    ///
    /// [`from_index`]: DeBruijnGraph::from_index
    pub fn add_edge(&mut self, from: u64, to: u64) {
        // Insert nodes if absent, tracking n_nodes separately.
        use std::collections::hash_map::Entry;
        if let Entry::Vacant(e) = self.successors.entry(from) {
            e.insert(Vec::new());
            self.predecessors.insert(from, Vec::new());
            self.n_nodes += 1;
        }
        if let Entry::Vacant(e) = self.successors.entry(to) {
            e.insert(Vec::new());
            self.predecessors.insert(to, Vec::new());
            self.n_nodes += 1;
        }
        self.successors.entry(from).or_default().push(to);
        self.predecessors.entry(to).or_default().push(from);
        self.n_edges += 1;
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use kam_core::kmer::MoleculeEvidence;
    use kam_index::encode::encode_kmer;
    use kam_index::HashKmerIndex;

    fn ev(n: u32) -> MoleculeEvidence {
        MoleculeEvidence {
            n_molecules: n,
            ..Default::default()
        }
    }

    fn build(kmers: &[u64], k: usize) -> DeBruijnGraph {
        let mut index = HashKmerIndex::new();
        for &km in kmers {
            index.insert(km, ev(1));
        }
        DeBruijnGraph::from_index(&index, k, kmers)
    }

    // Test 1: Two overlapping k-mers form an edge.
    #[test]
    fn two_overlapping_kmers_form_edge() {
        // ACG → CGT: suffix(ACG, k-1=2)=CG == prefix(CGT, k-1=2)=CG
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let graph = build(&[acg, cgt], 3);

        assert!(graph.successors(acg).contains(&cgt));
        assert!(graph.predecessors(cgt).contains(&acg));
    }

    // Test 2: Non-overlapping k-mers have no edge.
    #[test]
    fn non_overlapping_kmers_no_edge() {
        // ACG suffix=CG; TAC prefix=TA — no match.
        let acg = encode_kmer(b"ACG").unwrap();
        let tac = encode_kmer(b"TAC").unwrap();
        let graph = build(&[acg, tac], 3);

        assert!(graph.successors(acg).is_empty());
        assert!(graph.predecessors(tac).is_empty());
    }

    // Test 3: Linear sequence of k-mers forms a chain.
    #[test]
    fn linear_chain_of_kmers() {
        // ACGTA with k=3 → ACG → CGT → GTA
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let gta = encode_kmer(b"GTA").unwrap();
        let graph = build(&[acg, cgt, gta], 3);

        assert!(graph.successors(acg).contains(&cgt));
        assert!(graph.successors(cgt).contains(&gta));
        assert!(graph.predecessors(gta).contains(&cgt));
        assert!(graph.successors(gta).is_empty());
        assert!(graph.predecessors(acg).is_empty());
    }

    // Test 4: Branch point (SNV) creates two successors.
    #[test]
    fn branch_point_creates_two_successors() {
        // ACG → CGT and ACG → CGC (SNV at last base of successor)
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let cgc = encode_kmer(b"CGC").unwrap();
        let graph = build(&[acg, cgt, cgc], 3);

        let succs = graph.successors(acg);
        assert_eq!(succs.len(), 2);
        assert!(succs.contains(&cgt));
        assert!(succs.contains(&cgc));
    }

    // Test 5: n_nodes and n_edges correct.
    #[test]
    fn n_nodes_and_n_edges_correct() {
        // ACG → CGT → GTA: 3 nodes, 2 edges
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let gta = encode_kmer(b"GTA").unwrap();
        let graph = build(&[acg, cgt, gta], 3);

        assert_eq!(graph.n_nodes(), 3);
        assert_eq!(graph.n_edges(), 2);
    }

    // Test 6: Empty index → empty graph.
    #[test]
    fn empty_index_produces_empty_graph() {
        let index = HashKmerIndex::new();
        let graph = DeBruijnGraph::from_index(&index, 3, &[]);
        assert_eq!(graph.n_nodes(), 0);
        assert_eq!(graph.n_edges(), 0);
    }

    // Test 7: Self-loop k-mer (AAAA with k=4) handled.
    // AAAA encodes to 0. suffix(AAAA, k-1=3)=0; prefix(AAAA)=0 → self-loop.
    #[test]
    fn self_loop_kmer_handled() {
        let aaaa = encode_kmer(b"AAAA").unwrap(); // 0
        let graph = build(&[aaaa], 4);

        assert!(graph.successors(aaaa).contains(&aaaa));
        assert!(graph.predecessors(aaaa).contains(&aaaa));
        assert_eq!(graph.n_nodes(), 1);
        assert_eq!(graph.n_edges(), 1);
    }

    // Test 8: successors and predecessors are consistent (A→B ⟹ B has A as predecessor).
    #[test]
    fn successors_and_predecessors_consistent() {
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let cgc = encode_kmer(b"CGC").unwrap();
        let graph = build(&[acg, cgt, cgc], 3);

        // For every A→B edge, B must list A as predecessor.
        for (&a, succs) in &graph.successors {
            for &b in succs {
                assert!(
                    graph.predecessors(b).contains(&a),
                    "A→B edge exists but B does not list A as predecessor"
                );
            }
        }

        // And vice versa: for every B with predecessor A, A must list B as successor.
        for (&b, preds) in &graph.predecessors {
            for &a in preds {
                assert!(
                    graph.successors(a).contains(&b),
                    "B has A as predecessor but A does not list B as successor"
                );
            }
        }
    }
}
