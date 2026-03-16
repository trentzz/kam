//! Path scoring using molecule-level evidence from a k-mer index.
//!
//! Each [`GraphPath`] is scored by looking up the [`MoleculeEvidence`] for
//! every k-mer in the path. The aggregate statistics (minimum/mean molecule
//! counts, duplex fractions) summarise path confidence, and the weakest k-mer
//! identifies the bottleneck supporting the call.
//!
//! # Example
//! ```
//! use kam_pathfind::walk::GraphPath;
//! use kam_pathfind::score::{score_path, ScoredPath};
//! use kam_index::HashKmerIndex;
//! use kam_core::kmer::{KmerIndex, MoleculeEvidence};
//!
//! let mut index = HashKmerIndex::new();
//! index.insert(1, MoleculeEvidence { n_molecules: 10, n_duplex: 5, ..Default::default() });
//! index.insert(2, MoleculeEvidence { n_molecules: 8,  n_duplex: 4, ..Default::default() });
//!
//! let path = GraphPath { kmers: vec![1, 2], sequence: vec![], length: 2 };
//! let scored = score_path(&path, &index);
//! assert_eq!(scored.aggregate_evidence.min_molecules, 8);
//! assert_eq!(scored.weakest_kmer.kmer, 2);
//! ```

use kam_core::kmer::{KmerIndex, MoleculeEvidence};

use crate::walk::GraphPath;

/// A graph path annotated with molecule-level evidence.
#[derive(Debug, Clone)]
pub struct ScoredPath {
    /// The underlying path through the graph.
    pub path: GraphPath,
    /// Aggregate evidence statistics across all k-mers in the path.
    pub aggregate_evidence: PathEvidence,
    /// The weakest point in the path (lowest molecule support).
    pub weakest_kmer: WeakestKmer,
    /// `true` if the reconstructed sequence of this path matches the
    /// reference sequence.
    pub is_reference: bool,
}

/// Aggregate molecule-level evidence across all k-mers in a path.
#[derive(Debug, Clone)]
pub struct PathEvidence {
    /// Minimum `n_molecules` across all k-mers in the path.
    pub min_molecules: u32,
    /// Mean `n_molecules` across all k-mers in the path.
    pub mean_molecules: f32,
    /// Minimum `n_duplex` across all k-mers in the path.
    pub min_duplex: u32,
    /// Mean `n_duplex` across all k-mers in the path.
    pub mean_duplex: f32,
    /// Sum of `n_simplex_fwd` across all k-mers.
    pub total_simplex_fwd: u32,
    /// Sum of `n_simplex_rev` across all k-mers.
    pub total_simplex_rev: u32,
    /// Mean of `mean_base_error_prob` across all k-mers.
    pub mean_error_prob: f32,
}

/// The weakest k-mer in the path, i.e., the one with the fewest supporting
/// molecules.
#[derive(Debug, Clone)]
pub struct WeakestKmer {
    /// Encoded k-mer value.
    pub kmer: u64,
    /// Zero-based index of this k-mer in the path's `kmers` vector.
    pub position_in_path: usize,
    /// Evidence for this k-mer in the index.
    pub evidence: MoleculeEvidence,
}

/// Score a single path against the k-mer index.
///
/// Every k-mer in the path is looked up in `index`. K-mers absent from the
/// index are treated as having zero evidence (default [`MoleculeEvidence`]).
/// An empty path returns a default [`ScoredPath`] with all-zero statistics and
/// `weakest_kmer.kmer == 0`.
///
/// # Example
/// ```
/// use kam_pathfind::walk::GraphPath;
/// use kam_pathfind::score::score_path;
/// use kam_index::HashKmerIndex;
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
///
/// let mut index = HashKmerIndex::new();
/// index.insert(10, MoleculeEvidence { n_molecules: 5, ..Default::default() });
///
/// let path = GraphPath { kmers: vec![10], sequence: b"ACGT".to_vec(), length: 1 };
/// let scored = score_path(&path, &index);
/// assert_eq!(scored.aggregate_evidence.min_molecules, 5);
/// assert_eq!(scored.aggregate_evidence.mean_molecules, 5.0);
/// ```
pub fn score_path(path: &GraphPath, index: &dyn KmerIndex) -> ScoredPath {
    if path.kmers.is_empty() {
        let default_ev = MoleculeEvidence::default();
        return ScoredPath {
            path: path.clone(),
            aggregate_evidence: PathEvidence {
                min_molecules: 0,
                mean_molecules: 0.0,
                min_duplex: 0,
                mean_duplex: 0.0,
                total_simplex_fwd: 0,
                total_simplex_rev: 0,
                mean_error_prob: 0.0,
            },
            weakest_kmer: WeakestKmer {
                kmer: 0,
                position_in_path: 0,
                evidence: default_ev,
            },
            is_reference: false,
        };
    }

    // Collect evidence for every k-mer.
    let evidences: Vec<MoleculeEvidence> = path
        .kmers
        .iter()
        .map(|&kmer| {
            index
                .get(kmer)
                .cloned()
                .unwrap_or_default()
        })
        .collect();

    let n = evidences.len() as f32;

    let min_molecules = evidences.iter().map(|e| e.n_molecules).min().unwrap_or(0);
    let mean_molecules = evidences.iter().map(|e| e.n_molecules as f32).sum::<f32>() / n;
    let min_duplex = evidences.iter().map(|e| e.n_duplex).min().unwrap_or(0);
    let mean_duplex = evidences.iter().map(|e| e.n_duplex as f32).sum::<f32>() / n;
    let total_simplex_fwd = evidences.iter().map(|e| e.n_simplex_fwd).sum();
    let total_simplex_rev = evidences.iter().map(|e| e.n_simplex_rev).sum();
    let mean_error_prob =
        evidences.iter().map(|e| e.mean_base_error_prob).sum::<f32>() / n;

    // Find the weakest k-mer (lowest n_molecules; ties broken by first occurrence).
    let (weakest_idx, weakest_ev) = evidences
        .iter()
        .enumerate()
        .min_by_key(|(_, e)| e.n_molecules)
        .expect("evidences is non-empty");

    ScoredPath {
        path: path.clone(),
        aggregate_evidence: PathEvidence {
            min_molecules,
            mean_molecules,
            min_duplex,
            mean_duplex,
            total_simplex_fwd,
            total_simplex_rev,
            mean_error_prob,
        },
        weakest_kmer: WeakestKmer {
            kmer: path.kmers[weakest_idx],
            position_in_path: weakest_idx,
            evidence: weakest_ev.clone(),
        },
        is_reference: false, // set by score_and_rank_paths
    }
}

/// Score all paths, mark the reference path, and sort by evidence strength.
///
/// The reference path is determined by comparing each path's reconstructed
/// `sequence` against k-mers derived from `reference_seq` using a sliding
/// window of length `k`. A path is considered the reference if its `sequence`
/// equals `reference_seq` (byte-for-byte).
///
/// Paths are sorted by `min_molecules` descending (strongest evidence first).
/// Ties are broken by `mean_molecules` descending, then by path index for
/// determinism.
///
/// # Example
/// ```
/// use kam_pathfind::walk::GraphPath;
/// use kam_pathfind::score::score_and_rank_paths;
/// use kam_index::HashKmerIndex;
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
///
/// let mut index = HashKmerIndex::new();
/// index.insert(1, MoleculeEvidence { n_molecules: 10, ..Default::default() });
/// index.insert(2, MoleculeEvidence { n_molecules: 2, ..Default::default() });
///
/// let strong = GraphPath { kmers: vec![1], sequence: b"ACGT".to_vec(), length: 1 };
/// let weak   = GraphPath { kmers: vec![2], sequence: b"ACTT".to_vec(), length: 1 };
///
/// let ranked = score_and_rank_paths(vec![weak, strong], &index, b"ACGT", 4);
/// assert_eq!(ranked[0].aggregate_evidence.min_molecules, 10); // strong first
/// ```
pub fn score_and_rank_paths(
    paths: Vec<GraphPath>,
    index: &dyn KmerIndex,
    reference_seq: &[u8],
    _k: usize,
) -> Vec<ScoredPath> {
    let mut scored: Vec<ScoredPath> = paths
        .into_iter()
        .map(|p| score_path(&p, index))
        .collect();

    // Mark the reference path.
    for sp in &mut scored {
        if sp.path.sequence == reference_seq {
            sp.is_reference = true;
        }
    }

    // Sort by min_molecules descending, then mean_molecules descending for ties.
    scored.sort_by(|a, b| {
        b.aggregate_evidence
            .min_molecules
            .cmp(&a.aggregate_evidence.min_molecules)
            .then_with(|| {
                b.aggregate_evidence
                    .mean_molecules
                    .partial_cmp(&a.aggregate_evidence.mean_molecules)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    scored
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use kam_index::HashKmerIndex;

    fn ev(n_molecules: u32, n_duplex: u32) -> MoleculeEvidence {
        MoleculeEvidence {
            n_molecules,
            n_duplex,
            ..Default::default()
        }
    }

    fn simple_path(kmers: Vec<u64>, sequence: Vec<u8>) -> GraphPath {
        let length = kmers.len();
        GraphPath { kmers, sequence, length }
    }

    // Test 1: Single k-mer path → evidence equals that k-mer's evidence.
    #[test]
    fn single_kmer_path_evidence_matches() {
        let mut index = HashKmerIndex::new();
        index.insert(42, MoleculeEvidence {
            n_molecules: 7,
            n_duplex: 3,
            n_simplex_fwd: 2,
            n_simplex_rev: 2,
            min_base_error_prob: 0.001,
            mean_base_error_prob: 0.002,
        });

        let path = simple_path(vec![42], b"ACGT".to_vec());
        let scored = score_path(&path, &index);

        assert_eq!(scored.aggregate_evidence.min_molecules, 7);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 7.0);
        assert_eq!(scored.aggregate_evidence.min_duplex, 3);
        assert_eq!(scored.aggregate_evidence.mean_duplex, 3.0);
        assert_eq!(scored.aggregate_evidence.total_simplex_fwd, 2);
        assert_eq!(scored.aggregate_evidence.total_simplex_rev, 2);
        assert!((scored.aggregate_evidence.mean_error_prob - 0.002).abs() < 1e-6);
    }

    // Test 2: Multi-kmer path → min_molecules is the minimum across k-mers.
    #[test]
    fn multi_kmer_min_is_minimum() {
        let mut index = HashKmerIndex::new();
        index.insert(1, ev(10, 5));
        index.insert(2, ev(3, 1));
        index.insert(3, ev(8, 4));

        let path = simple_path(vec![1, 2, 3], b"ACGTA".to_vec());
        let scored = score_path(&path, &index);

        assert_eq!(scored.aggregate_evidence.min_molecules, 3);
        // mean = (10 + 3 + 8) / 3 = 7.0
        assert!((scored.aggregate_evidence.mean_molecules - 7.0).abs() < 1e-5);
    }

    // Test 3: Weakest k-mer identified correctly.
    #[test]
    fn weakest_kmer_identified() {
        let mut index = HashKmerIndex::new();
        index.insert(10, ev(15, 7));
        index.insert(20, ev(2, 0)); // weakest
        index.insert(30, ev(12, 6));

        let path = simple_path(vec![10, 20, 30], b"ACGTA".to_vec());
        let scored = score_path(&path, &index);

        assert_eq!(scored.weakest_kmer.kmer, 20);
        assert_eq!(scored.weakest_kmer.position_in_path, 1);
        assert_eq!(scored.weakest_kmer.evidence.n_molecules, 2);
    }

    // Test 4: Reference path detected by sequence comparison.
    #[test]
    fn reference_path_detected() {
        let ref_seq = b"ACGTACGT";
        let k = 4;

        let mut index = HashKmerIndex::new();
        index.insert(1, ev(10, 5));
        index.insert(2, ev(10, 5));

        let ref_path = simple_path(vec![1], ref_seq.to_vec());
        let alt_path = simple_path(vec![2], b"ACTTACGT".to_vec());

        let ranked = score_and_rank_paths(vec![alt_path, ref_path], &index, ref_seq, k);

        let ref_count = ranked.iter().filter(|p| p.is_reference).count();
        assert_eq!(ref_count, 1);
        assert!(ranked.iter().any(|p| p.is_reference && p.path.sequence == ref_seq));
    }

    // Test 5: Paths sorted by evidence strength.
    #[test]
    fn paths_sorted_strongest_first() {
        let mut index = HashKmerIndex::new();
        index.insert(1, ev(1, 0));   // weak
        index.insert(2, ev(50, 25)); // strong
        index.insert(3, ev(10, 5));  // medium

        let weak   = simple_path(vec![1], b"TTTT".to_vec());
        let strong = simple_path(vec![2], b"ACGT".to_vec());
        let medium = simple_path(vec![3], b"CCCC".to_vec());

        let ranked = score_and_rank_paths(vec![weak, strong, medium], &index, b"ACGT", 4);

        assert_eq!(ranked[0].aggregate_evidence.min_molecules, 50);
        assert_eq!(ranked[1].aggregate_evidence.min_molecules, 10);
        assert_eq!(ranked[2].aggregate_evidence.min_molecules, 1);
    }

    // Test 6: Path with k-mer missing from index → evidence has 0 molecules.
    #[test]
    fn missing_kmer_has_zero_evidence() {
        let index = HashKmerIndex::new(); // empty index

        let path = simple_path(vec![99999], b"ACGT".to_vec());
        let scored = score_path(&path, &index);

        assert_eq!(scored.aggregate_evidence.min_molecules, 0);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 0.0);
        assert_eq!(scored.weakest_kmer.evidence.n_molecules, 0);
    }

    // Test 7: Empty path → handle gracefully (all-zero evidence).
    #[test]
    fn empty_path_graceful() {
        let index = HashKmerIndex::new();
        let path = simple_path(vec![], vec![]);
        let scored = score_path(&path, &index);

        assert_eq!(scored.aggregate_evidence.min_molecules, 0);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 0.0);
        assert_eq!(scored.weakest_kmer.kmer, 0);
    }
}
