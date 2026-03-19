//! Path walking through a de Bruijn graph between anchor k-mers.
//!
//! Given a start anchor k-mer and an end anchor k-mer, [`walk_paths`] enumerates
//! paths through the graph connecting them. Each path corresponds to a possible
//! sequence (wildtype or variant).
//!
//! The walker uses iterative DFS with a single path buffer and a `HashSet` for
//! O(1) cycle detection. This uses O(B × L) memory regardless of branching
//! factor B and path length L, compared to O(B^(L/2) × L) for BFS with full
//! path cloning.
//!
//! # Example
//! ```
//! use kam_pathfind::graph::DeBruijnGraph;
//! use kam_pathfind::walk::{walk_paths, WalkConfig};
//!
//! // Build a simple linear graph: 0 → 1 → 2
//! let mut g = DeBruijnGraph::new(4);
//! g.add_edge(0, 1);
//! g.add_edge(1, 2);
//!
//! let paths = walk_paths(&g, 0, 2, &WalkConfig::default());
//! assert_eq!(paths.len(), 1);
//! assert_eq!(paths[0].kmers, vec![0, 1, 2]);
//! ```

use std::collections::HashSet;

use crate::graph::DeBruijnGraph;
use kam_index::encode::decode_kmer;

/// A path through the de Bruijn graph from a start to an end anchor k-mer.
#[derive(Debug, Clone)]
pub struct GraphPath {
    /// Ordered k-mers in the path (start k-mer first, end k-mer last).
    pub kmers: Vec<u64>,
    /// Reconstructed DNA sequence from the path.
    pub sequence: Vec<u8>,
    /// Number of k-mers in the path.
    pub length: usize,
}

/// Configuration controlling DFS path enumeration.
///
/// Set `max_path_length` based on the target sequence length rather than a
/// fixed bound. For a 100bp target with k=31, the expected path is ~70 k-mers;
/// a value of 150 gives generous headroom for large indels while preventing
/// the walker from exploring far outside the target window.
#[derive(Debug, Clone)]
pub struct WalkConfig {
    /// Maximum number of k-mers in a single path. Paths exceeding this length
    /// are abandoned. Set to approximately `target_len / (k - 1) + 50`.
    pub max_path_length: usize,
    /// Maximum total paths to return. Enumeration stops once this many
    /// complete paths have been found.
    pub max_paths: usize,
}

impl Default for WalkConfig {
    /// Conservative defaults suitable when the target length is not known.
    /// In the pathfind command, `max_path_length` is overridden per-target.
    fn default() -> Self {
        Self {
            max_path_length: 150,
            max_paths: 100,
        }
    }
}

/// Find all paths from `start_kmer` to `end_kmer` in the graph.
///
/// Uses iterative DFS with a single path buffer and a `HashSet` for O(1) cycle
/// detection. Memory use is O(B × L) where B is the maximum branching factor
/// and L is `max_path_length` — constant regardless of the number of paths
/// enumerated.
///
/// Enumeration stops once `config.max_paths` complete paths have been found.
/// Returns an empty `Vec` if no path exists or either k-mer is not in the graph.
///
/// # Example
/// ```
/// use kam_pathfind::graph::DeBruijnGraph;
/// use kam_pathfind::walk::{walk_paths, WalkConfig};
///
/// let mut g = DeBruijnGraph::new(3);
/// g.add_edge(1, 2);
/// g.add_edge(2, 3);
/// let paths = walk_paths(&g, 1, 3, &WalkConfig::default());
/// assert_eq!(paths.len(), 1);
/// assert_eq!(paths[0].length, 3);
/// ```
pub fn walk_paths(
    graph: &DeBruijnGraph,
    start_kmer: u64,
    end_kmer: u64,
    config: &WalkConfig,
) -> Vec<GraphPath> {
    let k = graph.k();
    let mut completed: Vec<GraphPath> = Vec::new();

    // Special case: start == end (single k-mer path).
    if start_kmer == end_kmer {
        if graph.contains(start_kmer) {
            let kmers = vec![start_kmer];
            let sequence = decode_kmer(start_kmer, k);
            completed.push(GraphPath {
                length: kmers.len(),
                kmers,
                sequence,
            });
        }
        return completed;
    }

    if !graph.contains(start_kmer) {
        return completed;
    }

    // DFS state: stack of (node, next_successor_index_to_try).
    // current_path holds the single in-progress path; visited tracks nodes in it.
    let mut stack: Vec<(u64, usize)> = vec![(start_kmer, 0)];
    let mut current_path: Vec<u64> = vec![start_kmer];
    let mut visited: HashSet<u64> = HashSet::from([start_kmer]);

    while let Some((node, succ_idx)) = stack.last_mut() {
        if completed.len() >= config.max_paths {
            break;
        }

        let succs = graph.successors(*node);

        if *succ_idx >= succs.len() {
            // All successors of this node tried — backtrack.
            let popped = current_path.pop().expect("path never empty during DFS");
            visited.remove(&popped);
            stack.pop();
            continue;
        }

        let next = succs[*succ_idx];
        *succ_idx += 1;

        // Skip cycles and paths that are already too long.
        if visited.contains(&next) {
            continue;
        }
        if current_path.len() >= config.max_path_length {
            continue;
        }

        if next == end_kmer {
            // Complete path found — record it without pushing to the DFS stack.
            let mut path_kmers = current_path.clone();
            path_kmers.push(next);
            let sequence = reconstruct_sequence(&path_kmers, k);
            completed.push(GraphPath {
                length: path_kmers.len(),
                kmers: path_kmers,
                sequence,
            });
        } else {
            // Extend the current path and push a new DFS frame.
            visited.insert(next);
            current_path.push(next);
            stack.push((next, 0));
        }
    }

    completed
}

/// Reconstruct the DNA sequence from an ordered list of encoded k-mers.
///
/// Successive k-mers overlap by k−1 bases. The first k-mer contributes all `k`
/// bases; each subsequent k-mer contributes only its last base.
///
/// Returns an empty `Vec` when `kmers` is empty.
///
/// # Example
/// ```
/// use kam_pathfind::walk::reconstruct_sequence;
/// use kam_index::encode::encode_kmer;
///
/// // Two 4-mers: ACGT and CGTA overlap in CGT.
/// // Sequence should be ACGTA.
/// let k1 = encode_kmer(b"ACGT").unwrap();
/// let k2 = encode_kmer(b"CGTA").unwrap();
/// let seq = reconstruct_sequence(&[k1, k2], 4);
/// assert_eq!(seq, b"ACGTA");
/// ```
pub fn reconstruct_sequence(kmers: &[u64], k: usize) -> Vec<u8> {
    if kmers.is_empty() {
        return Vec::new();
    }

    let mut seq = decode_kmer(kmers[0], k);

    for &kmer in &kmers[1..] {
        let last_base_bits = (kmer & 0b11) as u8;
        let last_base = match last_base_bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        };
        seq.push(last_base);
    }

    seq
}

#[cfg(test)]
mod tests {
    use super::*;
    use kam_index::encode::encode_kmer;

    fn linear_graph(seq: &[u8], k: usize) -> (DeBruijnGraph, Vec<u64>) {
        let mut g = DeBruijnGraph::new(k);
        let mut kmers: Vec<u64> = Vec::new();
        for i in 0..=seq.len().saturating_sub(k) {
            kmers.push(encode_kmer(&seq[i..i + k]).unwrap());
        }
        for w in kmers.windows(2) {
            g.add_edge(w[0], w[1]);
        }
        (g, kmers)
    }

    // Test 1: Linear path (no branches) → one path found.
    #[test]
    fn linear_path_single_result() {
        let seq = b"ACGTCCAG";
        let k = 4;
        let (g, kmers) = linear_graph(seq, k);
        let unique: std::collections::HashSet<u64> = kmers.iter().copied().collect();
        assert_eq!(
            unique.len(),
            kmers.len(),
            "test sequence must have unique k-mers"
        );
        let paths = walk_paths(&g, kmers[0], *kmers.last().unwrap(), &WalkConfig::default());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].kmers, kmers);
    }

    // Test 2: SNV creates two paths (ref and alt).
    #[test]
    fn snv_two_paths() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);

        let aaa = encode_kmer(b"AAA").unwrap();
        let aac = encode_kmer(b"AAC").unwrap();
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let gtt = encode_kmer(b"GTT").unwrap();
        let ttt = encode_kmer(b"TTT").unwrap();
        let aat = encode_kmer(b"AAT").unwrap();
        let att = encode_kmer(b"ATT").unwrap();

        g.add_edge(aaa, aac);
        g.add_edge(aac, acg);
        g.add_edge(acg, cgt);
        g.add_edge(cgt, gtt);
        g.add_edge(gtt, ttt);
        g.add_edge(aaa, aat);
        g.add_edge(aat, att);
        g.add_edge(att, ttt);

        let paths = walk_paths(&g, aaa, ttt, &WalkConfig::default());
        assert_eq!(paths.len(), 2);
    }

    // Test 3: Deletion creates a shorter alternative path.
    #[test]
    fn deletion_shorter_path() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);

        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let gta = encode_kmer(b"GTA").unwrap();

        g.add_edge(acg, cgt);
        g.add_edge(cgt, gta);
        g.add_edge(acg, gta);

        let paths = walk_paths(&g, acg, gta, &WalkConfig::default());
        assert_eq!(paths.len(), 2);

        let lengths: std::collections::BTreeSet<usize> = paths.iter().map(|p| p.length).collect();
        assert!(lengths.contains(&2));
        assert!(lengths.contains(&3));
    }

    // Test 4: Insertion creates a longer alternative path.
    #[test]
    fn insertion_longer_path() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);

        let start = encode_kmer(b"ACG").unwrap();
        let mid1 = encode_kmer(b"CGT").unwrap();
        let mid2 = encode_kmer(b"GTA").unwrap();
        let mid3 = encode_kmer(b"TAC").unwrap();
        let end = encode_kmer(b"ACT").unwrap();

        g.add_edge(start, mid1);
        g.add_edge(mid1, end);
        g.add_edge(start, mid2);
        g.add_edge(mid2, mid3);
        g.add_edge(mid3, end);

        let paths = walk_paths(&g, start, end, &WalkConfig::default());
        assert_eq!(paths.len(), 2);

        let lengths: std::collections::BTreeSet<usize> = paths.iter().map(|p| p.length).collect();
        assert!(lengths.contains(&3));
        assert!(lengths.contains(&4));
    }

    // Test 5: No path exists → empty result.
    #[test]
    fn no_path_returns_empty() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let a = encode_kmer(b"ACG").unwrap();
        let b = encode_kmer(b"CGT").unwrap();
        let c = encode_kmer(b"GTA").unwrap();
        g.add_edge(a, b);
        g.add_edge(c, a);

        let paths = walk_paths(&g, b, c, &WalkConfig::default());
        assert!(paths.is_empty());
    }

    // Test 6: max_path_length limits long paths.
    #[test]
    fn max_path_length_respected() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let kmers: Vec<u64> = (0u64..20).collect();
        for w in kmers.windows(2) {
            g.add_edge(w[0], w[1]);
        }

        let config = WalkConfig {
            max_path_length: 5,
            max_paths: 100,
        };
        let paths = walk_paths(&g, kmers[0], *kmers.last().unwrap(), &config);
        assert!(paths.is_empty(), "path should be abandoned as too long");
    }

    // Test 7: max_paths limits total paths found.
    #[test]
    fn max_paths_limits_results() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let start = 0u64;
        let end = 999u64;

        for mid in 1u64..=10 {
            g.add_edge(start, mid);
            g.add_edge(mid, end);
        }

        let config = WalkConfig {
            max_path_length: 150,
            max_paths: 3,
        };
        let paths = walk_paths(&g, start, end, &config);
        assert_eq!(paths.len(), 3);
    }

    // Test 8: Cycle in graph doesn't cause infinite loop.
    #[test]
    fn cycle_does_not_loop_forever() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);

        let a = 1u64;
        let b = 2u64;
        let c = 3u64;
        let end = 4u64;

        g.add_edge(a, b);
        g.add_edge(b, c);
        g.add_edge(c, a);
        g.add_edge(a, end);

        let paths = walk_paths(&g, a, end, &WalkConfig::default());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].kmers, vec![a, end]);
    }

    // Test 9: reconstruct_sequence from k-mers matches original sequence.
    #[test]
    fn reconstruct_sequence_matches_original() {
        let seq = b"ACGTACGT";
        let k = 4;
        let kmers: Vec<u64> = (0..=seq.len() - k)
            .map(|i| encode_kmer(&seq[i..i + k]).unwrap())
            .collect();
        let reconstructed = reconstruct_sequence(&kmers, k);
        assert_eq!(reconstructed, seq.as_ref());
    }

    // Test 10: Start == end (zero-length variant) → single path of one k-mer.
    #[test]
    fn start_equals_end_single_path() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let kmer = encode_kmer(b"ACG").unwrap();
        g.add_edge(kmer, kmer + 1);

        let paths = walk_paths(&g, kmer, kmer, &WalkConfig::default());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].kmers, vec![kmer]);
        assert_eq!(paths[0].length, 1);
    }
}
