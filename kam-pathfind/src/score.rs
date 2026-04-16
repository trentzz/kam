//! Path scoring using molecule-level evidence from a k-mer index.
//!
//! Each [`GraphPath`] is scored by looking up the [`MoleculeEvidence`] for
//! every k-mer in the path. The aggregate statistics (minimum/mean molecule
//! counts, duplex fractions) summarise path confidence, and the weakest k-mer
//! identifies the bottleneck supporting the call.
//!
//! Path k-mers are raw (non-canonical) because the de Bruijn graph is built
//! from raw k-mers to preserve suffix/prefix overlap relationships. Before
//! each index lookup, the raw k-mer is canonicalized so it matches the
//! canonical evidence index entries.
//!
//! # Example
//! ```
//! use kam_pathfind::walk::GraphPath;
//! use kam_pathfind::score::{score_path, ScoredPath};
//! use kam_index::HashKmerIndex;
//! use kam_index::encode::{encode_kmer, canonical};
//! use kam_core::kmer::{KmerIndex, MoleculeEvidence};
//!
//! let mut index = HashKmerIndex::new();
//! // Insert under canonical keys.
//! let k = 4;
//! let km1 = encode_kmer(b"ACGT").unwrap();
//! let km2 = encode_kmer(b"CGTA").unwrap();
//! index.insert(canonical(km1, k), MoleculeEvidence { n_molecules: 10, n_duplex: 5, ..Default::default() });
//! index.insert(canonical(km2, k), MoleculeEvidence { n_molecules: 8,  n_duplex: 4, ..Default::default() });
//!
//! // Path stores raw k-mers.
//! let path = GraphPath { kmers: vec![km1, km2], sequence: vec![], length: 2 };
//! let scored = score_path(&path, &index, k);
//! assert_eq!(scored.aggregate_evidence.min_molecules, 8);
//! ```

use std::collections::HashSet;

use kam_core::kmer::{KmerIndex, MoleculeEvidence};
use kam_index::encode::canonical;

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
    /// Minimum duplex count across variant-specific k-mers only.
    ///
    /// Variant-specific k-mers are those present in this (alt) path but absent
    /// from the reference path k-mer set.  For these k-mers, molecule evidence
    /// comes only from alt-supporting reads, not from the many reference
    /// molecules covering the shared anchor k-mers.  This gives a calibrated
    /// duplex count at the actual variant site rather than the inflated
    /// `mean_duplex` computed across all ~70 path k-mers.
    ///
    /// Set to 0 for paths where all k-mers are shared with the reference
    /// (including the reference path itself).  Computed by
    /// [`score_and_rank_paths`] after the reference path is identified; not
    /// available from [`score_path`] alone.
    pub min_variant_specific_duplex: u32,
    /// Mean `n_molecules` across variant-specific k-mers only.
    ///
    /// Variant-specific k-mers are those not shared with the reference path.
    /// For large SV paths, this is a better evidence estimate than
    /// `min_molecules` because the minimum over 70–150 k-mer positions
    /// bottlenecks at 1–3 even at moderate VAF.
    ///
    /// For paths where all k-mers are shared with the reference (including the
    /// reference path itself and fusion junction paths that match the fusion
    /// target sequence), this falls back to `mean_molecules`. Computed by
    /// [`score_and_rank_paths`].
    pub mean_variant_specific_molecules: f32,
    /// Minimum `n_simplex_fwd` across all k-mers in the path.
    ///
    /// Using the minimum rather than the sum gives a per-molecule count rather
    /// than a per-k-mer count.  The sum inflates the apparent sample size by
    /// the number of k-mers in the path (~70 for k=31 on a 100 bp target),
    /// making Fisher's exact test ~70× too sensitive and causing genuine low-VAF
    /// variants to be filtered as strand-biased.
    pub min_simplex_fwd: u32,
    /// Minimum `n_simplex_rev` across all k-mers in the path.  See `min_simplex_fwd`.
    pub min_simplex_rev: u32,
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
/// Every k-mer in the path is canonicalized before the index lookup because
/// the de Bruijn graph stores raw (non-canonical) k-mers while the evidence
/// index stores canonical k-mers.  K-mers absent from the index are treated
/// as having zero evidence (default [`MoleculeEvidence`]).  An empty path
/// returns a default [`ScoredPath`] with all-zero statistics and
/// `weakest_kmer.kmer == 0`.
///
/// # Arguments
///
/// - `path`: path through the raw de Bruijn graph
/// - `index`: canonical k-mer evidence index
/// - `k`: k-mer length (required to compute the canonical form)
///
/// # Example
/// ```
/// use kam_pathfind::walk::GraphPath;
/// use kam_pathfind::score::score_path;
/// use kam_index::HashKmerIndex;
/// use kam_index::encode::{encode_kmer, canonical};
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
///
/// let k = 4;
/// let raw_km = encode_kmer(b"TTTT").unwrap();
/// let can_km = canonical(raw_km, k); // AAAA (the smaller form)
///
/// let mut index = HashKmerIndex::new();
/// index.insert(can_km, MoleculeEvidence { n_molecules: 5, ..Default::default() });
///
/// // Path uses the raw k-mer; score_path canonicalizes before lookup.
/// let path = GraphPath { kmers: vec![raw_km], sequence: b"TTTT".to_vec(), length: 1 };
/// let scored = score_path(&path, &index, k);
/// assert_eq!(scored.aggregate_evidence.min_molecules, 5);
/// assert_eq!(scored.aggregate_evidence.mean_molecules, 5.0);
/// ```
pub fn score_path(path: &GraphPath, index: &dyn KmerIndex, k: usize) -> ScoredPath {
    if path.kmers.is_empty() {
        let default_ev = MoleculeEvidence::default();
        return ScoredPath {
            path: path.clone(),
            aggregate_evidence: PathEvidence {
                min_molecules: 0,
                mean_molecules: 0.0,
                min_duplex: 0,
                mean_duplex: 0.0,
                min_variant_specific_duplex: 0,
                mean_variant_specific_molecules: 0.0,
                min_simplex_fwd: 0,
                min_simplex_rev: 0,
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
    // Path k-mers are raw (non-canonical); canonicalize before the index lookup.
    let evidences: Vec<MoleculeEvidence> = path
        .kmers
        .iter()
        .map(|&kmer| {
            let canon = canonical(kmer, k);
            index.get(canon).cloned().unwrap_or_default()
        })
        .collect();

    let n = evidences.len() as f32;

    let min_molecules = evidences.iter().map(|e| e.n_molecules).min().unwrap_or(0);
    let mean_molecules = evidences.iter().map(|e| e.n_molecules as f32).sum::<f32>() / n;
    let min_duplex = evidences.iter().map(|e| e.n_duplex).min().unwrap_or(0);
    let mean_duplex = evidences.iter().map(|e| e.n_duplex as f32).sum::<f32>() / n;
    let min_simplex_fwd = evidences.iter().map(|e| e.n_simplex_fwd).min().unwrap_or(0);
    let min_simplex_rev = evidences.iter().map(|e| e.n_simplex_rev).min().unwrap_or(0);
    let mean_error_prob = evidences
        .iter()
        .map(|e| e.mean_base_error_prob)
        .sum::<f32>()
        / n;

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
            min_variant_specific_duplex: 0, // set by score_and_rank_paths
            mean_variant_specific_molecules: 0.0, // set by score_and_rank_paths
            min_simplex_fwd,
            min_simplex_rev,
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
/// `sequence` against `reference_seq` (byte-for-byte).
///
/// Paths are sorted by `min_molecules` descending (strongest evidence first).
/// Ties are broken by `mean_molecules` descending, then by path index for
/// determinism.
///
/// The `k` parameter is passed to [`score_path`] so it can canonicalize raw
/// path k-mers before looking them up in the canonical evidence index.
///
/// # Example
/// ```
/// use kam_pathfind::walk::GraphPath;
/// use kam_pathfind::score::score_and_rank_paths;
/// use kam_index::HashKmerIndex;
/// use kam_index::encode::{encode_kmer, canonical};
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
///
/// let k = 4;
/// let km1 = encode_kmer(b"ACGT").unwrap();
/// let km2 = encode_kmer(b"ACTT").unwrap();
/// let mut index = HashKmerIndex::new();
/// index.insert(canonical(km1, k), MoleculeEvidence { n_molecules: 10, ..Default::default() });
/// index.insert(canonical(km2, k), MoleculeEvidence { n_molecules: 2, ..Default::default() });
///
/// let strong = GraphPath { kmers: vec![km1], sequence: b"ACGT".to_vec(), length: 1 };
/// let weak   = GraphPath { kmers: vec![km2], sequence: b"ACTT".to_vec(), length: 1 };
///
/// let ranked = score_and_rank_paths(vec![weak, strong], &index, b"ACGT", k);
/// assert_eq!(ranked[0].aggregate_evidence.min_molecules, 10); // strong first
/// ```
pub fn score_and_rank_paths(
    paths: Vec<GraphPath>,
    index: &dyn KmerIndex,
    reference_seq: &[u8],
    k: usize,
) -> Vec<ScoredPath> {
    let mut scored: Vec<ScoredPath> = paths
        .into_iter()
        .map(|p| score_path(&p, index, k))
        .collect();

    // Mark the reference path.
    for sp in &mut scored {
        if sp.path.sequence == reference_seq {
            sp.is_reference = true;
        }
    }

    // Fallback: if no exact sequence match was found (e.g. DFS exhausted
    // max_paths before reaching the true reference path), identify the
    // reference by evidence. The true reference path has far higher molecule
    // support than error-k-mer paths (min_mol=1). Sort by evidence and mark
    // the strongest path as reference so alt calling can proceed.
    if !scored.iter().any(|sp| sp.is_reference) && !scored.is_empty() {
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
        scored[0].is_reference = true;
    }

    // Build a set of canonical k-mers from the reference path.
    // These are the anchor k-mers shared with every alt path.  Duplex evidence
    // at these k-mers reflects reference molecules, not alt-supporting molecules,
    // so they must be excluded when computing variant-specific duplex support.
    let ref_kmer_set: HashSet<u64> = scored
        .iter()
        .find(|sp| sp.is_reference)
        .map(|sp| {
            sp.path
                .kmers
                .iter()
                .map(|&kmer| canonical(kmer, k))
                .collect()
        })
        .unwrap_or_default();

    // For each path, compute the minimum duplex count and mean molecule count
    // across variant-specific k-mers (those not shared with the reference path).
    // This gives calibrated evidence at the actual variant site.
    // For the reference path (and fusion junction paths that match the reference),
    // all k-mers are in ref_kmer_set, so vs_evidences is empty and the fallback
    // sets mean_variant_specific_molecules = mean_molecules.
    for sp in &mut scored {
        let vs_evidences: Vec<MoleculeEvidence> = sp
            .path
            .kmers
            .iter()
            .filter_map(|&kmer| {
                let canon = canonical(kmer, k);
                if ref_kmer_set.contains(&canon) {
                    return None; // anchor k-mer: skip
                }
                Some(index.get(canon).cloned().unwrap_or_default())
            })
            .collect();

        sp.aggregate_evidence.min_variant_specific_duplex =
            vs_evidences.iter().map(|e| e.n_duplex).min().unwrap_or(0);

        // Mean molecule count across variant-specific k-mers only.
        // Falls back to the overall path mean when no variant-specific k-mers
        // exist (e.g. when all path k-mers are shared with the reference).
        sp.aggregate_evidence.mean_variant_specific_molecules = if vs_evidences.is_empty() {
            sp.aggregate_evidence.mean_molecules
        } else {
            vs_evidences
                .iter()
                .map(|e| e.n_molecules as f32)
                .sum::<f32>()
                / vs_evidences.len() as f32
        };
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
    use kam_index::encode::{canonical, encode_kmer};
    use kam_index::HashKmerIndex;

    // k used across most tests.
    const K: usize = 4;

    fn ev(n_molecules: u32, n_duplex: u32) -> MoleculeEvidence {
        MoleculeEvidence {
            n_molecules,
            n_duplex,
            ..Default::default()
        }
    }

    fn simple_path(kmers: Vec<u64>, sequence: Vec<u8>) -> GraphPath {
        let length = kmers.len();
        GraphPath {
            kmers,
            sequence,
            length,
        }
    }

    // Convenience: encode a 4-base k-mer sequence.
    fn enc(seq: &[u8]) -> u64 {
        encode_kmer(seq).unwrap()
    }

    // Convenience: canonical form of a 4-base k-mer sequence.
    fn can(seq: &[u8]) -> u64 {
        canonical(enc(seq), K)
    }

    // Test 1: Single k-mer path → evidence equals that k-mer's evidence.
    // The path stores the raw k-mer; the index stores its canonical form.
    // score_path must canonicalize before the lookup.
    #[test]
    fn single_kmer_path_evidence_matches() {
        // TTTT → canonical is AAAA (the smaller value).
        let raw = enc(b"TTTT");
        let canon = can(b"TTTT");
        assert_ne!(raw, canon, "TTTT should not be its own canonical form");

        let mut index = HashKmerIndex::new();
        index.insert(
            canon,
            MoleculeEvidence {
                n_molecules: 7,
                n_duplex: 3,
                n_simplex_fwd: 2,
                n_simplex_rev: 2,
                min_base_error_prob: 0.001,
                mean_base_error_prob: 0.002,
            },
        );

        // Path carries the raw k-mer (as the graph would).
        let path = simple_path(vec![raw], b"TTTT".to_vec());
        let scored = score_path(&path, &index, K);

        assert_eq!(scored.aggregate_evidence.min_molecules, 7);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 7.0);
        assert_eq!(scored.aggregate_evidence.min_duplex, 3);
        assert_eq!(scored.aggregate_evidence.mean_duplex, 3.0);
        assert_eq!(scored.aggregate_evidence.min_simplex_fwd, 2);
        assert_eq!(scored.aggregate_evidence.min_simplex_rev, 2);
        assert!((scored.aggregate_evidence.mean_error_prob - 0.002).abs() < 1e-6);
    }

    // Test 2: Multi-kmer path → min_molecules is the minimum across k-mers.
    // Use raw k-mers in the path and canonical k-mers in the index.
    #[test]
    fn multi_kmer_min_is_minimum() {
        // Three raw k-mers from a linear sequence ACGTACGTA (k=4).
        // ACGT, CGTA, GTAC
        let raw1 = enc(b"ACGT");
        let raw2 = enc(b"CGTA");
        let raw3 = enc(b"GTAC");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw1, K), ev(10, 5));
        index.insert(canonical(raw2, K), ev(3, 1));
        index.insert(canonical(raw3, K), ev(8, 4));

        let path = simple_path(vec![raw1, raw2, raw3], b"ACGTACGTA".to_vec());
        let scored = score_path(&path, &index, K);

        assert_eq!(scored.aggregate_evidence.min_molecules, 3);
        // mean = (10 + 3 + 8) / 3 = 7.0
        assert!((scored.aggregate_evidence.mean_molecules - 7.0).abs() < 1e-5);
    }

    // Test 3: Weakest k-mer identified correctly.
    #[test]
    fn weakest_kmer_identified() {
        let raw1 = enc(b"AAAC");
        let raw2 = enc(b"TTTT"); // weakest — lowest molecule count
        let raw3 = enc(b"CCCC");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw1, K), ev(15, 7));
        index.insert(canonical(raw2, K), ev(2, 0));
        index.insert(canonical(raw3, K), ev(12, 6));

        let path = simple_path(vec![raw1, raw2, raw3], b"AAACTTTTCCCC".to_vec());
        let scored = score_path(&path, &index, K);

        // Weakest k-mer is raw2 (TTTT), which should be at index 1.
        assert_eq!(scored.weakest_kmer.kmer, raw2);
        assert_eq!(scored.weakest_kmer.position_in_path, 1);
        assert_eq!(scored.weakest_kmer.evidence.n_molecules, 2);
    }

    // Test 4: Reference path detected by sequence comparison.
    #[test]
    fn reference_path_detected() {
        let ref_seq = b"ACGTACGT";
        let k = 4;

        let raw_ref = enc(b"ACGT");
        let raw_alt = enc(b"TTTT");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw_ref, k), ev(10, 5));
        index.insert(canonical(raw_alt, k), ev(10, 5));

        let ref_path = simple_path(vec![raw_ref], ref_seq.to_vec());
        let alt_path = simple_path(vec![raw_alt], b"TTTTACGT".to_vec());

        let ranked = score_and_rank_paths(vec![alt_path, ref_path], &index, ref_seq, k);

        let ref_count = ranked.iter().filter(|p| p.is_reference).count();
        assert_eq!(ref_count, 1);
        assert!(ranked
            .iter()
            .any(|p| p.is_reference && p.path.sequence == ref_seq));
    }

    // Test 5: Paths sorted by evidence strength.
    #[test]
    fn paths_sorted_strongest_first() {
        let raw_weak = enc(b"TTTT");
        let raw_strong = enc(b"AAAC");
        let raw_medium = enc(b"CCCC");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw_weak, K), ev(1, 0));
        index.insert(canonical(raw_strong, K), ev(50, 25));
        index.insert(canonical(raw_medium, K), ev(10, 5));

        let weak = simple_path(vec![raw_weak], b"TTTT".to_vec());
        let strong = simple_path(vec![raw_strong], b"AAAC".to_vec());
        let medium = simple_path(vec![raw_medium], b"CCCC".to_vec());

        let ranked = score_and_rank_paths(vec![weak, strong, medium], &index, b"AAAC", K);

        assert_eq!(ranked[0].aggregate_evidence.min_molecules, 50);
        assert_eq!(ranked[1].aggregate_evidence.min_molecules, 10);
        assert_eq!(ranked[2].aggregate_evidence.min_molecules, 1);
    }

    // Test 6: Path with k-mer missing from index → evidence has 0 molecules.
    #[test]
    fn missing_kmer_has_zero_evidence() {
        let index = HashKmerIndex::new(); // empty index
        let raw = enc(b"ACGT");
        let path = simple_path(vec![raw], b"ACGT".to_vec());
        let scored = score_path(&path, &index, K);

        assert_eq!(scored.aggregate_evidence.min_molecules, 0);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 0.0);
        assert_eq!(scored.weakest_kmer.evidence.n_molecules, 0);
    }

    // Test 7: Empty path → handle gracefully (all-zero evidence).
    #[test]
    fn empty_path_graceful() {
        let index = HashKmerIndex::new();
        let path = simple_path(vec![], vec![]);
        let scored = score_path(&path, &index, K);

        assert_eq!(scored.aggregate_evidence.min_molecules, 0);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 0.0);
        assert_eq!(scored.weakest_kmer.kmer, 0);
    }

    // Test 8: mean_variant_specific_molecules uses only non-reference k-mers.
    //
    // Reference path has k-mers [ref1, ref2].
    // Alt path has k-mers [ref1, alt1, alt2] where ref1 is shared with reference.
    // Variant-specific k-mers are alt1 (20 molecules) and alt2 (10 molecules).
    // Expected mean_variant_specific_molecules = (20 + 10) / 2 = 15.0.
    // min_molecules for the alt path is min(ref1=50, alt1=20, alt2=10) = 10.
    // The mean (15.0) is a better estimate than the minimum (10) at moderate VAF.
    #[test]
    fn mean_variant_specific_molecules_excludes_anchor_kmers() {
        let k = 4;
        let raw_ref1 = enc(b"AAAC");
        let raw_ref2 = enc(b"CCCC");
        let raw_alt1 = enc(b"TTTT");
        // raw_alt2 is redefined below to avoid GGGG (canonical = CCCC = raw_ref2).

        let mut index = HashKmerIndex::new();
        // ref1 has high molecule count — it is an anchor k-mer shared with the
        // alt path.  Its count must not inflate mean_variant_specific_molecules.
        index.insert(canonical(raw_ref1, k), ev(50, 25));
        index.insert(canonical(raw_ref2, k), ev(48, 24));
        // alt-specific k-mers have lower but consistent support.
        // Note: raw_alt2 must not share a canonical form with any reference k-mer.
        // TTTT: canonical = min(TTTT, AAAA) = AAAA ≠ ref1(AAAC) or ref2(CCCC) — OK.
        // ACGT: canonical = min(ACGT, ACGT) = ACGT ≠ AAAC or CCCC — OK.
        let raw_alt2 = enc(b"ACGT"); // replacing GGGG (canonical = CCCC = raw_ref2)
        index.insert(canonical(raw_alt1, k), ev(20, 10));
        index.insert(canonical(raw_alt2, k), ev(10, 5));

        // Reference path: ref1 + ref2.
        let ref_seq = b"AAACCCCC";
        let ref_path = simple_path(vec![raw_ref1, raw_ref2], ref_seq.to_vec());

        // Alt path: ref1 (anchor) + alt1 + alt2 (variant-specific).
        let alt_path = simple_path(vec![raw_ref1, raw_alt1, raw_alt2], b"AAACTTTTACGT".to_vec());

        let ranked = score_and_rank_paths(vec![alt_path, ref_path], &index, ref_seq, k);

        let alt = ranked
            .iter()
            .find(|p| !p.is_reference)
            .expect("alt path present");

        // mean of alt1 (20) and alt2 (10) = 15.0.
        assert!(
            (alt.aggregate_evidence.mean_variant_specific_molecules - 15.0).abs() < 1e-5,
            "expected 15.0, got {}",
            alt.aggregate_evidence.mean_variant_specific_molecules
        );

        // min_molecules is over all three k-mers including the anchor, so it
        // equals min(50, 20, 10) = 10. The variant-specific mean (15.0) is
        // higher and more representative of actual alt support.
        assert_eq!(alt.aggregate_evidence.min_molecules, 10);
    }

    // Test 9: reference path falls back to mean_molecules for
    // mean_variant_specific_molecules because all its k-mers are in ref_kmer_set
    // (vs_evidences is empty).
    //
    // This is correct behaviour for fusion calling: the fusion junction path IS
    // the reference path, so it must get a non-zero mean_variant_specific_molecules
    // for VAF computation in call_fusion.
    #[test]
    fn reference_path_mean_vs_molecules_falls_back_to_mean_molecules() {
        let k = 4;
        let raw_ref = enc(b"AAAC");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw_ref, k), ev(50, 25));

        let ref_seq = b"AAAC";
        let ref_path = simple_path(vec![raw_ref], ref_seq.to_vec());
        let ranked = score_and_rank_paths(vec![ref_path], &index, ref_seq, k);

        let refp = ranked.iter().find(|p| p.is_reference).unwrap();
        // vs_evidences is empty (the sole k-mer is in ref_kmer_set), so the
        // fallback sets mean_variant_specific_molecules = mean_molecules = 50.0.
        assert!(
            (refp.aggregate_evidence.mean_variant_specific_molecules - 50.0).abs() < 1e-5,
            "expected 50.0, got {}",
            refp.aggregate_evidence.mean_variant_specific_molecules
        );
    }

    // Test 10: Score path with all k-mers at 0 molecules (absent from index).
    // All aggregate values should be 0 or 0.0, with no panics or NaN.
    #[test]
    fn score_path_all_zero_molecules_no_nan() {
        let index = HashKmerIndex::new(); // empty index
        let raw1 = enc(b"AAAA");
        let raw2 = enc(b"CCCC");
        let path = simple_path(vec![raw1, raw2], b"AAAACCCC".to_vec());
        let scored = score_path(&path, &index, K);

        assert_eq!(scored.aggregate_evidence.min_molecules, 0);
        assert_eq!(scored.aggregate_evidence.mean_molecules, 0.0);
        assert!(!scored.aggregate_evidence.mean_molecules.is_nan());
        assert_eq!(scored.aggregate_evidence.min_duplex, 0);
        assert_eq!(scored.aggregate_evidence.mean_duplex, 0.0);
        assert!(!scored.aggregate_evidence.mean_duplex.is_nan());
        assert_eq!(scored.weakest_kmer.evidence.n_molecules, 0);
    }

    // Test 11: score_and_rank_paths with an empty paths list returns empty,
    // no panic.
    #[test]
    fn score_and_rank_empty_paths_list() {
        let index = HashKmerIndex::new();
        let ranked = score_and_rank_paths(vec![], &index, b"ACGT", K);
        assert!(ranked.is_empty(), "empty input must produce empty output");
    }

    // Test 12: Tie-breaking. Paths with the same min_molecules are sorted by
    // mean_molecules descending.
    #[test]
    fn tie_breaking_by_mean_molecules() {
        let raw_a = enc(b"AAAC");
        let raw_b = enc(b"CCCC");
        let raw_c = enc(b"TTTT");

        let mut index = HashKmerIndex::new();
        // Path A: k-mers with molecules [5, 5]. min=5, mean=5.0
        index.insert(canonical(raw_a, K), ev(5, 2));
        // Path B has two k-mers: one at 5, one at 10. min=5, mean=7.5
        index.insert(canonical(raw_b, K), ev(5, 2));
        index.insert(canonical(raw_c, K), ev(10, 5));

        let path_a = simple_path(vec![raw_a], b"AAAC".to_vec());
        let path_b = simple_path(vec![raw_b, raw_c], b"CCCCTTTT".to_vec());

        // Neither matches reference sequence "ZZZZ", so fallback marks
        // strongest as reference.
        let ranked = score_and_rank_paths(vec![path_a, path_b], &index, b"ZZZZ", K);
        assert_eq!(ranked.len(), 2);

        // Both have min_molecules=5. Path B has mean=7.5 > path A mean=5.0.
        // Path B should be ranked first.
        assert!(
            ranked[0].aggregate_evidence.mean_molecules > ranked[1].aggregate_evidence.mean_molecules,
            "tie on min_molecules should break by mean_molecules descending"
        );
    }

    // Test 13: Duplex fields (min_duplex, mean_duplex) are correctly aggregated.
    #[test]
    fn duplex_fields_aggregated_correctly() {
        let raw1 = enc(b"ACGT");
        let raw2 = enc(b"CGTA");
        let raw3 = enc(b"GTAC");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw1, K), ev(10, 8));
        index.insert(canonical(raw2, K), ev(10, 3));
        index.insert(canonical(raw3, K), ev(10, 6));

        let path = simple_path(vec![raw1, raw2, raw3], b"ACGTACGTAC".to_vec());
        let scored = score_path(&path, &index, K);

        assert_eq!(
            scored.aggregate_evidence.min_duplex, 3,
            "min_duplex should be the smallest duplex count across all k-mers"
        );
        // mean_duplex = (8 + 3 + 6) / 3 ≈ 5.667
        let expected_mean = (8.0 + 3.0 + 6.0) / 3.0;
        assert!(
            (scored.aggregate_evidence.mean_duplex - expected_mean).abs() < 1e-3,
            "expected mean_duplex ≈ {expected_mean}, got {}",
            scored.aggregate_evidence.mean_duplex
        );
    }

    // Test 14: Fallback reference detection when no path's sequence matches
    // the provided reference_seq. The strongest path (by min_molecules then
    // mean_molecules) should be marked as reference.
    #[test]
    fn fallback_reference_detection_marks_strongest() {
        let raw_strong = enc(b"AAAC");
        let raw_weak = enc(b"TTTT");

        let mut index = HashKmerIndex::new();
        index.insert(canonical(raw_strong, K), ev(50, 25));
        index.insert(canonical(raw_weak, K), ev(2, 1));

        let strong = simple_path(vec![raw_strong], b"AAAC".to_vec());
        let weak = simple_path(vec![raw_weak], b"TTTT".to_vec());

        // Reference sequence "XXXX" matches neither path.
        let ranked = score_and_rank_paths(vec![weak, strong], &index, b"XXXX", K);
        assert_eq!(ranked.len(), 2);

        // Exactly one path should be marked as reference (the strongest).
        let ref_count = ranked.iter().filter(|p| p.is_reference).count();
        assert_eq!(ref_count, 1, "exactly one path should be marked as reference");

        let ref_path = ranked.iter().find(|p| p.is_reference).expect("reference present");
        assert_eq!(
            ref_path.aggregate_evidence.min_molecules, 50,
            "the strongest path should be the fallback reference"
        );
    }

    // Test 15: score_path with a single k-mer whose canonical form differs
    // from the raw form. The lookup must canonicalise before querying the index.
    // This tests the canonicalisation path directly.
    #[test]
    fn score_path_canonicalises_before_lookup() {
        let raw = enc(b"TTTT"); // raw = TTTT
        let canon = canonical(raw, K); // canonical = AAAA (smaller)
        assert_ne!(raw, canon, "TTTT should not equal its canonical form AAAA");

        let mut index = HashKmerIndex::new();
        // Insert under the canonical key only.
        index.insert(canon, ev(42, 21));

        let path = simple_path(vec![raw], b"TTTT".to_vec());
        let scored = score_path(&path, &index, K);
        assert_eq!(
            scored.aggregate_evidence.min_molecules, 42,
            "score_path must canonicalise the raw k-mer to find evidence"
        );
    }

    // Test 16: Score path with simplex strand counts. Verify min_simplex_fwd
    // and min_simplex_rev are the minimums across all k-mers.
    #[test]
    fn simplex_strand_counts_aggregated() {
        let raw1 = enc(b"ACGT");
        let raw2 = enc(b"CGTA");

        let mut index = HashKmerIndex::new();
        index.insert(
            canonical(raw1, K),
            MoleculeEvidence {
                n_molecules: 10,
                n_duplex: 5,
                n_simplex_fwd: 3,
                n_simplex_rev: 2,
                min_base_error_prob: 0.001,
                mean_base_error_prob: 0.002,
            },
        );
        index.insert(
            canonical(raw2, K),
            MoleculeEvidence {
                n_molecules: 10,
                n_duplex: 5,
                n_simplex_fwd: 1,
                n_simplex_rev: 4,
                min_base_error_prob: 0.001,
                mean_base_error_prob: 0.003,
            },
        );

        let path = simple_path(vec![raw1, raw2], b"ACGTCGTA".to_vec());
        let scored = score_path(&path, &index, K);

        assert_eq!(
            scored.aggregate_evidence.min_simplex_fwd, 1,
            "min_simplex_fwd should be the minimum across k-mers"
        );
        assert_eq!(
            scored.aggregate_evidence.min_simplex_rev, 2,
            "min_simplex_rev should be the minimum across k-mers"
        );
        // mean_error_prob = (0.002 + 0.003) / 2 = 0.0025
        assert!(
            (scored.aggregate_evidence.mean_error_prob - 0.0025).abs() < 1e-6,
            "mean_error_prob should be the average across k-mers"
        );
    }
}
