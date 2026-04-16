//! De Bruijn graph construction from a k-mer index.
//!
//! Nodes are raw (non-canonical) k-mers. Edges connect k-mers that overlap by k-1 bases.
//! An edge A → B exists when the last k-1 bases of A equal the first k-1 bases of B.
//!
//! In 2-bit encoded form:
//! - suffix of A (k-1 bases) = `A & mask(k-1)` (lower 2*(k-1) bits)
//! - prefix of B (k-1 bases) = `B >> 2` (shift right by 2 to drop last base)
//! - Edge exists when suffix(A) == prefix(B)

use std::collections::{HashMap, HashSet, VecDeque};

use kam_core::kmer::KmerIndex;

/// A de Bruijn graph built from a k-mer index.
///
/// Nodes are raw (non-canonical) k-mers encoded as `u64`. An edge A → B exists when
/// `suffix(A, k-1) == prefix(B, k-1)` in 2-bit encoded form.
///
/// Only successor adjacency is stored. Predecessor adjacency was removed
/// because path walking only traverses forward edges.
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
/// let graph = DeBruijnGraph::from_index(&index, 3, &all_kmers, 1);
/// assert!(graph.successors(acg).contains(&cgt));
/// assert_eq!(graph.n_nodes(), 2);
/// assert_eq!(graph.n_edges(), 1);
/// ```
#[derive(Debug)]
pub struct DeBruijnGraph {
    k: usize,
    /// Forward adjacency list: kmer → list of successor kmers
    successors: HashMap<u64, Vec<u64>>,
    /// Reverse adjacency list: kmer → list of predecessor kmers.
    /// Built once at construction time so backward BFS does not need to
    /// rebuild it on every `backward_reachable` call.
    predecessors: HashMap<u64, Vec<u64>>,
    /// Number of nodes
    n_nodes: usize,
    /// Number of edges
    n_edges: usize,
}

impl DeBruijnGraph {
    /// Build a de Bruijn graph from a k-mer index.
    ///
    /// Only includes k-mers present in the index with at least `min_molecules`
    /// supporting molecules. Edges are derived from 2*(k-1) overlap:
    /// suffix(A) == prefix(B).
    ///
    /// The `min_molecules` threshold filters out low-evidence k-mers (PCR and
    /// sequencing errors) before graph construction. This reduces spurious
    /// branching from error k-mers, which is the primary driver of memory
    /// explosion in the path-walking stage. Use `min_molecules = 1` to include
    /// all observed k-mers.
    ///
    /// # Arguments
    ///
    /// - `index`: the k-mer index (used to verify membership and check evidence)
    /// - `k`: k-mer length
    /// - `all_kmers`: full set of encoded k-mers to consider as candidate nodes
    /// - `min_molecules`: minimum molecule count to include a node in the graph
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
    /// let kmers = vec![acg, cgt];
    /// let graph = DeBruijnGraph::from_index(&index, 3, &kmers, 1);
    /// assert_eq!(graph.n_nodes(), 2);
    /// assert_eq!(graph.n_edges(), 1);
    /// ```
    pub fn from_index(
        index: &dyn KmerIndex,
        k: usize,
        all_kmers: &[u64],
        min_molecules: u32,
    ) -> Self {
        // Filter to nodes that are in the index with sufficient molecule support.
        let nodes: Vec<u64> = all_kmers
            .iter()
            .copied()
            .filter(|&km| index.molecule_count(km) >= min_molecules)
            .collect();

        // Initialise successor adjacency map with empty vecs for every node.
        let mut successors: HashMap<u64, Vec<u64>> =
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
                    n_edges += 1;
                }
            }
        }

        // Build reverse adjacency (predecessor) map once at construction time
        // so that backward_reachable can BFS without rebuilding it on every call.
        let mut predecessors: HashMap<u64, Vec<u64>> =
            nodes.iter().map(|&km| (km, Vec::new())).collect();
        for (&a, succs) in &successors {
            for &b in succs {
                predecessors.entry(b).or_default().push(a);
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
    /// index.insert(acg, MoleculeEvidence { n_molecules: 1, ..Default::default() });
    /// index.insert(cgt, MoleculeEvidence { n_molecules: 1, ..Default::default() });
    /// let graph = DeBruijnGraph::from_index(&index, 3, &[acg, cgt], 1);
    /// assert!(graph.successors(acg).contains(&cgt));
    /// let empty: &[u64] = &[];
    /// assert_eq!(graph.successors(999), empty);
    /// ```
    pub fn successors(&self, kmer: u64) -> &[u64] {
        self.successors.get(&kmer).map(Vec::as_slice).unwrap_or(&[])
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
    /// index.insert(acg, MoleculeEvidence { n_molecules: 1, ..Default::default() });
    /// let graph = DeBruijnGraph::from_index(&index, 3, &[acg], 1);
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
    /// Compute the set of nodes reachable from `start` via forward edges within
    /// `max_hops` steps.
    ///
    /// Used together with [`backward_reachable`] to identify nodes that lie on
    /// at least one valid start→end path, enabling aggressive DFS pruning.
    ///
    /// [`backward_reachable`]: DeBruijnGraph::backward_reachable
    pub fn forward_reachable(&self, start: u64, max_hops: usize) -> HashSet<u64> {
        let mut reachable: HashSet<u64> = HashSet::new();
        let mut queue: VecDeque<(u64, usize)> = VecDeque::new();
        reachable.insert(start);
        queue.push_back((start, 0));

        while let Some((node, depth)) = queue.pop_front() {
            if depth >= max_hops {
                continue;
            }
            for &succ in self.successors(node) {
                if reachable.insert(succ) {
                    queue.push_back((succ, depth + 1));
                }
            }
        }

        reachable
    }

    /// Compute the set of nodes from which `end` is reachable within `max_hops`
    /// forward steps.
    ///
    /// This is a backward BFS: it traverses predecessor edges (computed on the
    /// fly from the forward adjacency list) starting from `end`. The resulting
    /// set is used together with [`forward_reachable`] to identify nodes that
    /// lie on at least one valid start→end path.
    ///
    /// Runs in O(n_edges) time. Building the predecessor map is O(n_edges);
    /// BFS is O(n_nodes + n_edges).
    ///
    /// [`forward_reachable`]: DeBruijnGraph::forward_reachable
    pub fn backward_reachable(&self, end: u64, max_hops: usize) -> HashSet<u64> {
        // BFS backward from `end` using the pre-built predecessor map.
        let mut reachable: HashSet<u64> = HashSet::new();
        let mut queue: VecDeque<(u64, usize)> = VecDeque::new();
        reachable.insert(end);
        queue.push_back((end, 0));

        while let Some((node, depth)) = queue.pop_front() {
            if depth >= max_hops {
                continue;
            }
            if let Some(preds) = self.predecessors.get(&node) {
                for &pred in preds {
                    if reachable.insert(pred) {
                        queue.push_back((pred, depth + 1));
                    }
                }
            }
        }

        reachable
    }

    pub fn add_edge(&mut self, from: u64, to: u64) {
        use std::collections::hash_map::Entry;
        if let Entry::Vacant(e) = self.successors.entry(from) {
            e.insert(Vec::new());
            self.n_nodes += 1;
        }
        if let Entry::Vacant(e) = self.successors.entry(to) {
            e.insert(Vec::new());
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
        DeBruijnGraph::from_index(&index, k, kmers, 1)
    }

    // Test 1: Two overlapping k-mers form an edge.
    #[test]
    fn two_overlapping_kmers_form_edge() {
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let graph = build(&[acg, cgt], 3);

        assert!(graph.successors(acg).contains(&cgt));
    }

    // Test 2: Non-overlapping k-mers have no edge.
    #[test]
    fn non_overlapping_kmers_no_edge() {
        let acg = encode_kmer(b"ACG").unwrap();
        let tac = encode_kmer(b"TAC").unwrap();
        let graph = build(&[acg, tac], 3);

        assert!(graph.successors(acg).is_empty());
    }

    // Test 3: Linear sequence of k-mers forms a chain.
    #[test]
    fn linear_chain_of_kmers() {
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let gta = encode_kmer(b"GTA").unwrap();
        let graph = build(&[acg, cgt, gta], 3);

        assert!(graph.successors(acg).contains(&cgt));
        assert!(graph.successors(cgt).contains(&gta));
        assert!(graph.successors(gta).is_empty());
        assert!(!graph.contains(999));
    }

    // Test 4: Branch point (SNV) creates two successors.
    #[test]
    fn branch_point_creates_two_successors() {
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
        let graph = DeBruijnGraph::from_index(&index, 3, &[], 1);
        assert_eq!(graph.n_nodes(), 0);
        assert_eq!(graph.n_edges(), 0);
    }

    // Test 7: Self-loop k-mer (AAAA with k=4) handled.
    #[test]
    fn self_loop_kmer_handled() {
        let aaaa = encode_kmer(b"AAAA").unwrap();
        let graph = build(&[aaaa], 4);

        assert!(graph.successors(aaaa).contains(&aaaa));
        assert_eq!(graph.n_nodes(), 1);
        assert_eq!(graph.n_edges(), 1);
    }

    // Test 8: min_molecules filters low-evidence k-mers.
    #[test]
    fn min_molecules_filters_low_evidence() {
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let mut index = HashKmerIndex::new();
        index.insert(acg, ev(3));
        index.insert(cgt, ev(1)); // below threshold

        // With min_molecules=2, cgt is excluded.
        let graph = DeBruijnGraph::from_index(&index, 3, &[acg, cgt], 2);
        assert!(graph.contains(acg));
        assert!(!graph.contains(cgt));
        assert_eq!(graph.n_nodes(), 1);
        assert_eq!(graph.n_edges(), 0);
    }

    // Test 9: add_edge builds graph correctly (no predecessors map required).
    #[test]
    fn add_edge_builds_successors() {
        let mut g = DeBruijnGraph::new(3);
        g.add_edge(10, 20);
        assert!(g.contains(10));
        assert!(g.contains(20));
        assert_eq!(g.successors(10), &[20]);
        assert_eq!(g.n_nodes(), 2);
        assert_eq!(g.n_edges(), 1);
    }

    // ── New edge-case tests ──────────────────────────────────────────────────

    // Test 10: Single k-mer in the graph produces one node and no edges
    // (unless it self-loops, but a non-self-looping k-mer should have 0 edges).
    #[test]
    fn single_non_looping_kmer_one_node_no_edges() {
        // ACG has suffix CG and prefix AC. They do not match, so no self-loop.
        let acg = encode_kmer(b"ACG").expect("valid k-mer");
        let graph = build(&[acg], 3);
        assert_eq!(graph.n_nodes(), 1);
        assert_eq!(graph.n_edges(), 0);
        assert!(graph.successors(acg).is_empty());
    }

    // Test 11: Graph from a linear sequence of length 2k.
    // For k=3 and sequence "ACGTAC" (length 6 = 2*3), there are 4 k-mers:
    // ACG → CGT → GTA → TAC. These form 3 linear edges. However, TAC has
    // suffix AC, and ACG has prefix AC (ACG >> 2), so an additional edge
    // TAC → ACG exists, creating a cycle. Total edges = 4.
    #[test]
    fn linear_sequence_of_length_2k() {
        let k = 3;
        let seq = b"ACGTAC";
        assert_eq!(seq.len(), 2 * k);

        let kmers: Vec<u64> = (0..=seq.len() - k)
            .map(|i| encode_kmer(&seq[i..i + k]).expect("valid k-mer"))
            .collect();
        assert_eq!(kmers.len(), k + 1, "expected k+1 = {} k-mers", k + 1);

        let graph = build(&kmers, k);
        assert_eq!(graph.n_nodes(), k + 1);
        // 3 linear edges plus a wrap-around edge TAC → ACG = 4 edges total.
        assert_eq!(
            graph.n_edges(),
            k + 1,
            "linear chain plus wrap-around edge TAC → ACG"
        );

        // Each consecutive pair has a forward edge.
        for w in kmers.windows(2) {
            assert!(
                graph.successors(w[0]).contains(&w[1]),
                "expected edge {} → {}",
                w[0],
                w[1]
            );
        }
        // TAC → ACG wrap-around edge.
        let tac = *kmers.last().expect("non-empty");
        let acg = kmers[0];
        assert!(
            graph.successors(tac).contains(&acg),
            "TAC should have a wrap-around edge to ACG"
        );
    }

    // Test 12: forward_reachable from the start of a linear graph returns
    // only downstream nodes within the hop limit.
    #[test]
    fn forward_reachable_linear_graph() {
        let mut g = DeBruijnGraph::new(3);
        // Chain: 0 → 1 → 2 → 3 → 4
        for i in 0u64..4 {
            g.add_edge(i, i + 1);
        }

        // From node 0 with max_hops=2, reachable nodes are {0, 1, 2}.
        let reach = g.forward_reachable(0, 2);
        assert!(reach.contains(&0));
        assert!(reach.contains(&1));
        assert!(reach.contains(&2));
        assert!(!reach.contains(&3), "node 3 is 3 hops away, beyond limit");
        assert!(!reach.contains(&4));

        // With unlimited hops, all nodes are reachable.
        let reach_all = g.forward_reachable(0, 100);
        for i in 0u64..=4 {
            assert!(reach_all.contains(&i), "node {i} should be reachable");
        }
    }

    // Test 13: backward_reachable from the end of a linear graph returns
    // only upstream nodes within the hop limit.
    #[test]
    fn backward_reachable_linear_graph() {
        let mut g = DeBruijnGraph::new(3);
        // Chain: 0 → 1 → 2 → 3 → 4
        for i in 0u64..4 {
            g.add_edge(i, i + 1);
        }

        // From node 4 with max_hops=2, backward-reachable nodes are {4, 3, 2}.
        let reach = g.backward_reachable(4, 2);
        assert!(reach.contains(&4));
        assert!(reach.contains(&3));
        assert!(reach.contains(&2));
        assert!(!reach.contains(&1), "node 1 is 3 hops upstream, beyond limit");
        assert!(!reach.contains(&0));
    }

    // Test 14: Successors of a node not in the graph returns an empty slice.
    #[test]
    fn successors_of_absent_node_empty() {
        let graph = build(&[], 3);
        let empty: &[u64] = &[];
        assert_eq!(graph.successors(12345), empty);
    }

    // Test 15: add_edge with the same pair twice creates duplicate edges.
    // This mirrors real graph construction where a k-mer can overlap itself.
    #[test]
    fn add_edge_allows_duplicates() {
        let mut g = DeBruijnGraph::new(3);
        g.add_edge(10, 20);
        g.add_edge(10, 20);
        assert_eq!(g.successors(10).len(), 2, "duplicate edges are preserved");
        assert_eq!(g.n_edges(), 2);
        // Node count should still be 2 (not 4).
        assert_eq!(g.n_nodes(), 2);
    }

    // Test 16: With k=1 and suffix_mask=0, every node's suffix is 0 and every
    // node's prefix is also 0 (B >> 2 for a 2-bit value is 0). So every node
    // is a successor of every other node (and itself). A full graph on 4
    // single-base k-mers has 4 nodes and 16 edges (including self-loops).
    #[test]
    fn k_equals_one_fully_connected_graph() {
        let a = encode_kmer(b"A").expect("valid 1-mer");
        let c = encode_kmer(b"C").expect("valid 1-mer");
        let g_base = encode_kmer(b"G").expect("valid 1-mer");
        let t = encode_kmer(b"T").expect("valid 1-mer");
        let all = vec![a, c, g_base, t];

        let graph = build(&all, 1);
        assert_eq!(graph.n_nodes(), 4);
        // Every pair (including self-loops) is an edge: 4 * 4 = 16.
        assert_eq!(graph.n_edges(), 16);

        for &src in &all {
            let succs = graph.successors(src);
            assert_eq!(succs.len(), 4, "every node should connect to all 4 nodes");
            for &dst in &all {
                assert!(
                    succs.contains(&dst),
                    "edge {} -> {} should exist",
                    src,
                    dst
                );
            }
        }
    }

    // Test 17: Diamond graph (branch then merge). Two paths from start to end.
    // ACG branches to CGA and CGT, both of which converge on a common
    // successor. Verify correct node and edge counts.
    #[test]
    fn diamond_graph_branch_and_merge() {
        let k = 3;
        let acg = encode_kmer(b"ACG").expect("valid k-mer");
        let cga = encode_kmer(b"CGA").expect("valid k-mer");
        let cgt = encode_kmer(b"CGT").expect("valid k-mer");
        // Both CGA and CGT share suffix GA/GT? No. We need them to converge.
        // CGA has suffix GA. CGT has suffix GT. They converge on different
        // targets. Instead, build explicitly with add_edge.
        let end = encode_kmer(b"GTA").expect("valid k-mer");

        let mut g = DeBruijnGraph::new(k);
        g.add_edge(acg, cga);
        g.add_edge(acg, cgt);
        g.add_edge(cga, end);
        g.add_edge(cgt, end);

        assert_eq!(g.n_nodes(), 4);
        assert_eq!(g.n_edges(), 4);
        assert_eq!(g.successors(acg).len(), 2);
        assert!(g.successors(cga).contains(&end));
        assert!(g.successors(cgt).contains(&end));
    }

    // Test 18: forward_reachable with max_hops=0 returns only the start node.
    #[test]
    fn forward_reachable_zero_hops_returns_only_start() {
        let mut g = DeBruijnGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);

        let reach = g.forward_reachable(0, 0);
        assert_eq!(reach.len(), 1);
        assert!(reach.contains(&0));
    }

    // Test 19: backward_reachable with max_hops=0 returns only the end node.
    #[test]
    fn backward_reachable_zero_hops_returns_only_end() {
        let mut g = DeBruijnGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);

        let reach = g.backward_reachable(2, 0);
        assert_eq!(reach.len(), 1);
        assert!(reach.contains(&2));
    }

    // Test 20: from_index with all k-mers below min_molecules produces an
    // empty graph. No nodes survive the filter.
    #[test]
    fn all_below_min_molecules_produces_empty_graph() {
        let mut index = HashKmerIndex::new();
        let acg = encode_kmer(b"ACG").expect("valid k-mer");
        let cgt = encode_kmer(b"CGT").expect("valid k-mer");
        index.insert(acg, ev(2));
        index.insert(cgt, ev(3));

        // min_molecules=10 filters out both k-mers.
        let graph = DeBruijnGraph::from_index(&index, 3, &[acg, cgt], 10);
        assert_eq!(graph.n_nodes(), 0);
        assert_eq!(graph.n_edges(), 0);
    }

    // Test 21: Disconnected components. Two separate linear chains share no
    // reachable nodes.
    #[test]
    fn disconnected_components_independent_reachability() {
        let mut g = DeBruijnGraph::new(3);
        // Component A: 0 → 1 → 2
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        // Component B: 10 → 11 → 12
        g.add_edge(10, 11);
        g.add_edge(11, 12);

        let reach_a = g.forward_reachable(0, 100);
        assert!(reach_a.contains(&0));
        assert!(reach_a.contains(&1));
        assert!(reach_a.contains(&2));
        assert!(
            !reach_a.contains(&10),
            "component B nodes should not be reachable from component A"
        );
        assert!(!reach_a.contains(&11));
        assert!(!reach_a.contains(&12));

        let reach_b = g.backward_reachable(12, 100);
        assert!(reach_b.contains(&10));
        assert!(reach_b.contains(&11));
        assert!(reach_b.contains(&12));
        assert!(!reach_b.contains(&0));
    }

    // Test 22: forward_reachable on a node not in the graph returns a set
    // containing only that node (the start is always included).
    #[test]
    fn forward_reachable_absent_node() {
        let g = DeBruijnGraph::new(3);
        let reach = g.forward_reachable(999, 10);
        assert_eq!(reach.len(), 1);
        assert!(reach.contains(&999));
    }
}
