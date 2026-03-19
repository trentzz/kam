//! Path walking through a de Bruijn graph between anchor k-mers.
//!
//! Given a start anchor k-mer and an end anchor k-mer, [`walk_paths`] enumerates
//! all paths through the graph connecting them. Each path corresponds to a
//! possible sequence (wildtype or variant).
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

use std::collections::VecDeque;

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

/// Configuration controlling BFS path enumeration.
///
/// The defaults are conservative upper bounds that work well for typical panel
/// amplicons (≤500 bp, diploid alleles only). Reduce them for performance if
/// you know the region is small.
#[derive(Debug, Clone)]
pub struct WalkConfig {
    /// Maximum number of k-mers in a single path (prevents infinite loops on
    /// very long regions or graphs with many long cycles).
    pub max_path_length: usize,
    /// Maximum total paths to return. Enumeration stops once this many
    /// complete paths have been found.
    pub max_paths: usize,
}

impl Default for WalkConfig {
    fn default() -> Self {
        Self {
            max_path_length: 500,
            max_paths: 100,
        }
    }
}

/// Find all paths from `start_kmer` to `end_kmer` in the graph.
///
/// Uses BFS with per-path visited-set tracking to avoid cycles. Each path is
/// bounded to at most `config.max_path_length` k-mers; exploration stops once
/// `config.max_paths` complete paths have been found.
///
/// Returns an empty `Vec` if no path exists (including if either k-mer is not
/// in the graph).
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

    // Special case: start == end (zero-length variant — single k-mer path).
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

    // BFS queue: each element is the partial path so far (as a list of kmers).
    // We store the path as a Vec to allow cycle detection per-branch.
    let mut queue: VecDeque<Vec<u64>> = VecDeque::new();
    queue.push_back(vec![start_kmer]);

    while let Some(partial) = queue.pop_front() {
        if completed.len() >= config.max_paths {
            break;
        }

        let last = *partial.last().expect("partial path is never empty");

        for &next in graph.successors(last) {
            if completed.len() >= config.max_paths {
                break;
            }

            // Cycle detection: do not revisit a k-mer already in this path.
            if partial.contains(&next) {
                continue;
            }

            let mut new_path = partial.clone();
            new_path.push(next);

            if new_path.len() > config.max_path_length {
                // This branch is too long — abandon it.
                continue;
            }

            if next == end_kmer {
                // Complete path found.
                let sequence = reconstruct_sequence(&new_path, k);
                completed.push(GraphPath {
                    length: new_path.len(),
                    kmers: new_path,
                    sequence,
                });
            } else {
                queue.push_back(new_path);
            }
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

    // Start with the first k-mer fully decoded.
    let mut seq = decode_kmer(kmers[0], k);

    // Each subsequent k-mer contributes only its last (rightmost) base.
    for &kmer in &kmers[1..] {
        // The last base is encoded in the 2 least-significant bits.
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

    /// Build a simple linear graph from a DNA sequence using k-mers of length k.
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
        // Use a sequence with all unique k-mers (no repeats).
        let seq = b"ACGTCCAG";
        let k = 4;
        let (g, kmers) = linear_graph(seq, k);
        // Ensure all k-mers are distinct (required for cycle-free linear path).
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
        // Reference: ACGTA, Alt: ACTTA (G→T SNV at position 2).
        // k=3: ACG-CGT-GTA (ref), ACT-CTT-TTA (alt)
        // Build a graph with shared start and end k-mers but two middle paths.
        //
        // Shared prefix k-mer: ACA (acts as start anchor)
        // Shared suffix k-mer: TAT (acts as end anchor)
        // Ref path: ACA → CAG → AGT → GTA → TAT
        // Alt path: ACA → CAT → ATT → TTA → TAT  (wait, overlapping TAT is hard to share)
        //
        // Simpler: use k=3
        // ref seq: ACGTA  → ACG, CGT, GTA
        // alt seq: ACATA  → ACA, CAT, ATA
        // shared: start=ACG is not shared...
        //
        // Design: start anchor = AAA, end anchor = TTT, two intermediate paths.
        // ref: AAA → AAC → ACG → CGT → GTT → TTT
        // alt: AAA → AAT → ATT → TTT
        // (both share start=AAA and end=TTT but have different middle k-mers)
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

        // Ref branch: AAA → AAC → ACG → CGT → GTT → TTT
        g.add_edge(aaa, aac);
        g.add_edge(aac, acg);
        g.add_edge(acg, cgt);
        g.add_edge(cgt, gtt);
        g.add_edge(gtt, ttt);

        // Alt branch: AAA → AAT → ATT → TTT
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

        // ref: AACGT → AAC → ACG → CGT
        // del: AAGT  → AAG → AGT  (shorter — one k-mer in common at boundary)
        // anchor start = AA*, end = *GT
        // Use: start=AAC, end=CGT for ref; and AAG→AGT for del where end=AGT
        // Actually let's use start/end as fixed anchors and have two paths.
        //
        // Start: ACG, End: GTA
        // Long (ref): ACG → CGT → GTA
        // Short (del): ACG → GTA  (direct edge simulates 1-base deletion at kmer level)
        let acg = encode_kmer(b"ACG").unwrap();
        let cgt = encode_kmer(b"CGT").unwrap();
        let gta = encode_kmer(b"GTA").unwrap();

        // Ref path (longer)
        g.add_edge(acg, cgt);
        g.add_edge(cgt, gta);

        // Del path (shorter, direct)
        g.add_edge(acg, gta);

        let paths = walk_paths(&g, acg, gta, &WalkConfig::default());
        assert_eq!(paths.len(), 2);

        // One path should be length 2 (direct), one length 3 (through CGT).
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
        let mid3 = encode_kmer(b"TAC").unwrap(); // extra k-mer for insertion path
        let end = encode_kmer(b"ACT").unwrap();

        // Short (ref) path: start → mid1 → end
        g.add_edge(start, mid1);
        g.add_edge(mid1, end);

        // Long (ins) path: start → mid2 → mid3 → end
        g.add_edge(start, mid2);
        g.add_edge(mid2, mid3);
        g.add_edge(mid3, end);

        let paths = walk_paths(&g, start, end, &WalkConfig::default());
        assert_eq!(paths.len(), 2);

        let lengths: std::collections::BTreeSet<usize> = paths.iter().map(|p| p.length).collect();
        assert!(lengths.contains(&3)); // short path
        assert!(lengths.contains(&4)); // long path
    }

    // Test 5: No path exists → empty result.
    #[test]
    fn no_path_returns_empty() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let a = encode_kmer(b"ACG").unwrap();
        let b = encode_kmer(b"CGT").unwrap();
        g.add_edge(a, b); // a → b, but not b → anything useful

        let c = encode_kmer(b"GTA").unwrap();
        g.add_edge(c, a); // c → a, but a is not start here

        // Try to find a path from b to c — no such path exists.
        let paths = walk_paths(&g, b, c, &WalkConfig::default());
        assert!(paths.is_empty());
    }

    // Test 6: max_path_length limits long paths.
    #[test]
    fn max_path_length_respected() {
        // Build a long linear chain of 20 k-mers.
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let kmers: Vec<u64> = (0u64..20).collect();
        for w in kmers.windows(2) {
            g.add_edge(w[0], w[1]);
        }

        let config = WalkConfig {
            max_path_length: 5, // too short to reach the end
            max_paths: 100,
        };

        let paths = walk_paths(&g, kmers[0], *kmers.last().unwrap(), &config);
        assert!(paths.is_empty(), "path should be abandoned as too long");
    }

    // Test 7: max_paths limits total paths found.
    #[test]
    fn max_paths_limits_results() {
        // Build a graph with many alternative routes by adding lots of branches.
        let k = 3;
        let mut g = DeBruijnGraph::new(k);

        let start = 0u64;
        let end = 999u64;

        // Add 10 independent single-hop paths from start to end.
        for mid in 1u64..=10 {
            g.add_edge(start, mid);
            g.add_edge(mid, end);
        }

        let config = WalkConfig {
            max_path_length: 500,
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

        // Cycle: a → b → c → a (and a also → end)
        g.add_edge(a, b);
        g.add_edge(b, c);
        g.add_edge(c, a);
        g.add_edge(a, end);

        let paths = walk_paths(&g, a, end, &WalkConfig::default());
        // The only path without revisiting is a → end (direct).
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
        g.add_edge(kmer, kmer + 1); // ensure the node exists

        let paths = walk_paths(&g, kmer, kmer, &WalkConfig::default());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].kmers, vec![kmer]);
        assert_eq!(paths[0].length, 1);
    }
}
