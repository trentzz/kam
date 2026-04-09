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
//! let (paths, exceeded) = walk_paths(&g, 0, 2, &WalkConfig::default());
//! assert_eq!(paths.len(), 1);
//! assert_eq!(paths[0].kmers, vec![0, 1, 2]);
//! assert!(!exceeded);
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
    /// Total DFS node expansions before aborting. 0 = unlimited.
    pub max_expansions: usize,
}

impl Default for WalkConfig {
    /// Conservative defaults suitable when the target length is not known.
    /// In the pathfind command, `max_path_length` is overridden per-target.
    fn default() -> Self {
        Self {
            max_path_length: 150,
            max_paths: 100,
            max_expansions: 0,
        }
    }
}

/// Shared iterative DFS implementation parameterised by successor ordering.
///
/// `sort_successors` is called on each node's successor list before the DFS
/// explores it. Pass a no-op closure for unordered traversal or a sort by
/// molecule count for evidence-biased traversal.
fn walk_paths_inner<F>(
    graph: &DeBruijnGraph,
    start_kmer: u64,
    end_kmer: u64,
    config: &WalkConfig,
    sort_successors: F,
) -> (Vec<GraphPath>, bool)
where
    F: Fn(&mut Vec<u64>),
{
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
        return (completed, false);
    }

    if !graph.contains(start_kmer) {
        return (completed, false);
    }

    // Stack stores (node, pre-ordered successors, next index to try).
    // Successors are sorted once on first visit so the ordering is stable for
    // the lifetime of that DFS frame.
    let start_succs = {
        let mut s = graph.successors(start_kmer).to_vec();
        sort_successors(&mut s);
        s
    };
    let mut stack: Vec<(u64, Vec<u64>, usize)> = vec![(start_kmer, start_succs, 0)];
    let mut current_path: Vec<u64> = vec![start_kmer];
    let mut visited: HashSet<u64> = HashSet::from([start_kmer]);
    let mut n_expansions: usize = 0;
    let mut budget_exceeded = false;

    while let Some((_node, succs, succ_idx)) = stack.last_mut() {
        if completed.len() >= config.max_paths {
            break;
        }

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
            let mut next_succs = graph.successors(next).to_vec();
            sort_successors(&mut next_succs);
            stack.push((next, next_succs, 0));
            n_expansions += 1;
            if config.max_expansions > 0 && n_expansions >= config.max_expansions {
                budget_exceeded = true;
                break;
            }
        }
    }

    (completed, budget_exceeded)
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
/// let (paths, exceeded) = walk_paths(&g, 1, 3, &WalkConfig::default());
/// assert_eq!(paths.len(), 1);
/// assert_eq!(paths[0].length, 3);
/// assert!(!exceeded);
/// ```
pub fn walk_paths(
    graph: &DeBruijnGraph,
    start_kmer: u64,
    end_kmer: u64,
    config: &WalkConfig,
) -> (Vec<GraphPath>, bool) {
    walk_paths_inner(graph, start_kmer, end_kmer, config, |_succs| {})
}

/// Find alt paths by branching off the reference path at each position.
///
/// Exhaustive DFS (even evidence-biased) may not find early-branching SV alt
/// paths within `max_paths`: single-base error bubbles form complete paths for
/// every position in the reference path, exhausting the budget before the DFS
/// backtracks to the SV branch point (which may be only 20–50 k-mers from the
/// start for a large deletion).
///
/// This function walks `reference_kmers` and at each position checks for
/// successors NOT on the reference path with molecule evidence >=
/// `min_alt_evidence`. For each such "high-evidence branch", a short DFS
/// (`max_alt_paths_per_branch` paths) is run from the branch k-mer to
/// `end_kmer`. The resulting paths are stitched with the reference prefix to
/// form complete alt paths.
///
/// # Arguments
///
/// - `reference_kmers`: ordered k-mers of the identified reference path.
/// - `end_kmer`: the end anchor k-mer (same as used in the original DFS).
/// - `max_path_length`: maximum k-mers in a sub-path (same as `WalkConfig`).
/// - `max_alt_paths_per_branch`: `max_paths` for each sub-DFS. 3–5 is enough
///   for most SV alt paths.
/// - `branch_filter`: closure that returns `true` for non-reference successors
///   that should be explored. Use this to require junction k-mer membership
///   (preventing error-bubble DFS runs) and/or a minimum evidence threshold.
/// - `molecule_count`: closure returning molecule count for a raw k-mer.
pub fn find_alt_paths_from_reference(
    graph: &DeBruijnGraph,
    reference_kmers: &[u64],
    end_kmer: u64,
    max_path_length: usize,
    max_alt_paths_per_branch: usize,
    branch_filter: impl Fn(u64) -> bool,
    molecule_count: impl Fn(u64) -> u32 + Copy,
) -> Vec<GraphPath> {
    let k = graph.k();

    // Build a set of raw k-mers on the reference path for fast exclusion.
    let ref_kmer_set: std::collections::HashSet<u64> = reference_kmers.iter().copied().collect();

    let mut alt_paths: Vec<GraphPath> = Vec::new();
    // Deduplicate by sequence bytes to handle multiple branch points that
    // produce the same alt sequence.
    let mut seen_sequences: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();

    for (i, &ref_kmer) in reference_kmers.iter().enumerate() {
        for &successor in graph.successors(ref_kmer) {
            // Skip k-mers already on the reference path.
            if ref_kmer_set.contains(&successor) {
                continue;
            }
            // Apply branch filter (junction membership / evidence threshold).
            if !branch_filter(successor) {
                continue;
            }

            // Sub-path budget: cap the sub-walk tightly so the DFS stays
            // bounded. The sub-walk from `successor` to `end_kmer` can use at
            // most (max_path_length - (i+1)) k-mers (the full alt path minus
            // the already-committed reference prefix).  Hard-cap at 300 to
            // prevent exponential search when a spurious junction k-mer match
            // triggers a DFS into a highly connected region far from end_kmer.
            let sub_budget = max_path_length.saturating_sub(i + 1).clamp(1, 300);
            let walk_config = WalkConfig {
                max_path_length: sub_budget,
                max_paths: max_alt_paths_per_branch,
                max_expansions: 300_000,
            };

            // Walk from `successor` to `end_kmer`. Budget exceeded flag is
            // discarded here — sub-walks are small and bounded by design.
            let (rest_paths, _) =
                walk_paths_biased(graph, successor, end_kmer, &walk_config, molecule_count);
            for rest in rest_paths {
                // Stitch: reference_kmers[0..=i] + rest.kmers.
                // rest.kmers[0] == successor, which is a graph-valid successor
                // of reference_kmers[i], so the concatenation is a valid path.
                let total_len = (i + 1) + rest.kmers.len();
                let mut full_kmers = Vec::with_capacity(total_len);
                full_kmers.extend_from_slice(&reference_kmers[..=i]);
                full_kmers.extend_from_slice(&rest.kmers);

                let sequence = reconstruct_sequence(&full_kmers, k);

                if seen_sequences.contains(&sequence) {
                    continue;
                }
                seen_sequences.insert(sequence.clone());

                alt_paths.push(GraphPath {
                    length: full_kmers.len(),
                    kmers: full_kmers,
                    sequence,
                });
            }
        }
    }

    alt_paths
}

/// Find all paths from `start_kmer` to `end_kmer`, exploring high-evidence
/// successors first.
///
/// Identical to [`walk_paths`] except that at each DFS node the successors are
/// sorted by `molecule_count(kmer)` descending before exploration. This ensures
/// the high-evidence reference path is encountered early, before `max_paths` is
/// exhausted by low-evidence error-k-mer branches that branch off the main path.
///
/// The `molecule_count` closure should canonicalize the raw k-mer and query the
/// evidence index. Successors with equal evidence are explored in their natural
/// order (stable sort).
///
/// # Example
/// ```
/// use kam_pathfind::graph::DeBruijnGraph;
/// use kam_pathfind::walk::{walk_paths_biased, WalkConfig};
///
/// let mut g = DeBruijnGraph::new(4);
/// g.add_edge(0, 1);
/// g.add_edge(1, 2);
/// let (paths, exceeded) = walk_paths_biased(&g, 0, 2, &WalkConfig::default(), |_| 1);
/// assert_eq!(paths.len(), 1);
/// assert_eq!(paths[0].kmers, vec![0, 1, 2]);
/// assert!(!exceeded);
/// ```
pub fn walk_paths_biased(
    graph: &DeBruijnGraph,
    start_kmer: u64,
    end_kmer: u64,
    config: &WalkConfig,
    molecule_count: impl Fn(u64) -> u32,
) -> (Vec<GraphPath>, bool) {
    walk_paths_inner(graph, start_kmer, end_kmer, config, |succs| {
        succs.sort_by_key(|&km| std::cmp::Reverse(molecule_count(km)));
    })
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
        let (paths, _) = walk_paths(&g, kmers[0], *kmers.last().unwrap(), &WalkConfig::default());
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

        let (paths, _) = walk_paths(&g, aaa, ttt, &WalkConfig::default());
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

        let (paths, _) = walk_paths(&g, acg, gta, &WalkConfig::default());
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

        let (paths, _) = walk_paths(&g, start, end, &WalkConfig::default());
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

        let (paths, _) = walk_paths(&g, b, c, &WalkConfig::default());
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
            max_expansions: 0,
        };
        let (paths, _) = walk_paths(&g, kmers[0], *kmers.last().unwrap(), &config);
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
            max_expansions: 0,
        };
        let (paths, _) = walk_paths(&g, start, end, &config);
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

        let (paths, _) = walk_paths(&g, a, end, &WalkConfig::default());
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

    // Test 9b: Budget exceeded returns partial result and flag.
    #[test]
    fn budget_exceeded_returns_partial_and_flag() {
        // Build a graph where start branches into dead-end paths that never
        // reach end. The graph has exactly 4 possible node expansions total
        // (b1, b3, b2, b4), so a budget of 3 fires before all branches are
        // explored.
        //
        // DFS order (biased, all molecule counts equal):
        //   expand b1 (n=1) → expand b3 (n=2) → b3 dead-ends, backtrack
        //   → expand b2 (n=3) → budget triggered
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let start = encode_kmer(b"AAA").unwrap();
        let end = encode_kmer(b"TTT").unwrap();
        let b1 = encode_kmer(b"AAC").unwrap();
        let b2 = encode_kmer(b"AAG").unwrap();
        let b3 = encode_kmer(b"ACG").unwrap();
        let b4 = encode_kmer(b"AGC").unwrap();
        g.add_edge(start, b1);
        g.add_edge(start, b2);
        g.add_edge(b1, b3);
        g.add_edge(b2, b4);
        // Add end as an isolated node with no path from start to it.
        g.add_edge(end, end);

        let cfg = WalkConfig {
            max_path_length: 50,
            max_paths: 100,
            max_expansions: 3,
        };
        let (paths, exceeded) = walk_paths_biased(&g, start, end, &cfg, |_| 1);
        assert!(
            exceeded,
            "budget of 3 should be exceeded in this 4-expansion graph"
        );
        assert!(paths.is_empty(), "no complete paths reach end");
    }

    // Test 10: Start == end (zero-length variant) → single path of one k-mer.
    #[test]
    fn start_equals_end_single_path() {
        let k = 3;
        let mut g = DeBruijnGraph::new(k);
        let kmer = encode_kmer(b"ACG").unwrap();
        g.add_edge(kmer, kmer + 1);

        let (paths, _) = walk_paths(&g, kmer, kmer, &WalkConfig::default());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].kmers, vec![kmer]);
        assert_eq!(paths[0].length, 1);
    }
}
