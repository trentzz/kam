//! K-mer extraction from consensus reads into a [`KmerIndex`].
//!
//! This module bridges [`kam-assemble`] output to [`kam-index`]: for each
//! consensus read it extracts every k-mer, canonicalises it, and accumulates
//! [`MoleculeEvidence`] in the provided index.

use kam_core::kmer::KmerIndex;
use kam_core::molecule::FamilyType;

use crate::encode::{canonical_with_strand, KmerIterator};

/// Information about a consensus read needed for k-mer extraction.
///
/// This is a lightweight view type intended to be constructed by the assembly
/// stage before passing reads to the indexing stage.
///
/// # Example
///
/// ```
/// use kam_index::extract::ConsensusReadInfo;
/// use kam_core::molecule::FamilyType;
///
/// let read = ConsensusReadInfo {
///     sequence: b"ACGT".to_vec(),
///     per_base_error_prob: vec![0.01, 0.01, 0.01, 0.01],
///     family_type: FamilyType::Duplex,
/// };
/// assert_eq!(read.sequence.len(), 4);
/// ```
#[derive(Debug, Clone)]
pub struct ConsensusReadInfo {
    /// The consensus base sequence (ASCII bytes, e.g. `b"ACGT..."`).
    pub sequence: Vec<u8>,
    /// Per-base probability of error (one value per base, 0.0–1.0).
    pub per_base_error_prob: Vec<f32>,
    /// Classification of the read family that produced this consensus.
    pub family_type: FamilyType,
}

/// Extract k-mers from a single consensus read and update the index.
///
/// For each k-mer position in `read.sequence`:
/// 1. Encode and canonicalise the k-mer.
/// 2. Retrieve or create a [`MoleculeEvidence`] entry in the index.
/// 3. Increment `n_molecules`.
/// 4. Increment the appropriate strand counter based on `read.family_type`.
/// 5. Compute the mean error probability across the k bases of the window.
/// 6. Update `min_base_error_prob` and `mean_base_error_prob` (running average).
///
/// K-mers containing `N` or other unrecognised bases are silently skipped.
///
/// # Arguments
///
/// * `read` — the consensus read to process.
/// * `k` — k-mer length (must be in `1..=31`).
/// * `index` — the index to update.
///
/// # Example
///
/// ```
/// use kam_index::extract::{ConsensusReadInfo, extract_and_index};
/// use kam_index::HashKmerIndex;
/// use kam_core::molecule::FamilyType;
/// use kam_core::kmer::KmerIndex;
///
/// let read = ConsensusReadInfo {
///     sequence: b"ACGTACGT".to_vec(),
///     per_base_error_prob: vec![0.01; 8],
///     family_type: FamilyType::Duplex,
/// };
/// let mut index = HashKmerIndex::new();
/// extract_and_index(&read, 4, &mut index);
/// // 8 - 4 + 1 = 5 k-mer positions; some canonical forms may overlap.
/// assert!(!index.is_empty());
/// ```
pub fn extract_and_index(read: &ConsensusReadInfo, k: usize, index: &mut dyn KmerIndex) {
    let seq = &read.sequence;
    let probs = &read.per_base_error_prob;

    for (pos, raw_kmer) in KmerIterator::new(seq, k) {
        let (canon, _strand) = canonical_with_strand(raw_kmer, k);
        let mean_err = kmer_mean_error_prob(probs, pos, k);

        // Retrieve or create the evidence entry.
        // We must read and then write to satisfy the trait interface (no entry() API).
        let existing = index.get(canon).cloned();
        let mut ev = existing.unwrap_or_default();

        // 3. Increment molecule count.
        ev.n_molecules += 1;

        // 4. Increment strand-specific counter.
        match read.family_type {
            FamilyType::Duplex => ev.n_duplex += 1,
            FamilyType::SimplexFwd => ev.n_simplex_fwd += 1,
            FamilyType::SimplexRev => ev.n_simplex_rev += 1,
            // Singleton is treated as simplex forward.
            FamilyType::Singleton => ev.n_simplex_fwd += 1,
        }

        // 5-7. Update error probability statistics.
        let prev_count = ev.n_molecules - 1; // observations before this one
        if prev_count == 0 {
            // First observation: initialise both stats.
            ev.min_base_error_prob = mean_err;
            ev.mean_base_error_prob = mean_err;
        } else {
            if mean_err < ev.min_base_error_prob {
                ev.min_base_error_prob = mean_err;
            }
            // Running average: new_mean = (old_mean * n + new_val) / (n + 1)
            ev.mean_base_error_prob =
                (ev.mean_base_error_prob * prev_count as f32 + mean_err) / ev.n_molecules as f32;
        }

        index.insert(canon, ev);
    }
}

/// Extract k-mers from multiple consensus reads into the index.
///
/// Equivalent to calling [`extract_and_index`] for each element of `reads`.
///
/// # Arguments
///
/// * `reads` — the consensus reads to process.
/// * `k` — k-mer length (must be in `1..=31`).
/// * `index` — the index to update.
///
/// # Example
///
/// ```
/// use kam_index::extract::{ConsensusReadInfo, extract_all};
/// use kam_index::HashKmerIndex;
/// use kam_core::molecule::FamilyType;
///
/// let reads = vec![
///     ConsensusReadInfo {
///         sequence: b"ACGT".to_vec(),
///         per_base_error_prob: vec![0.01; 4],
///         family_type: FamilyType::SimplexFwd,
///     },
/// ];
/// let mut index = HashKmerIndex::new();
/// extract_all(&reads, 4, &mut index);
/// ```
pub fn extract_all(reads: &[ConsensusReadInfo], k: usize, index: &mut dyn KmerIndex) {
    for read in reads {
        extract_and_index(read, k, index);
    }
}

/// Compute the mean error probability across the k bases starting at `start`.
///
/// If `error_probs` is shorter than `start + k` the available values are used;
/// missing values are treated as `0.0`.
fn kmer_mean_error_prob(error_probs: &[f32], start: usize, k: usize) -> f32 {
    if k == 0 {
        return 0.0;
    }
    let sum: f32 = (start..start + k)
        .map(|i| error_probs.get(i).copied().unwrap_or(0.0))
        .sum();
    sum / k as f32
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::{canonical, encode_kmer};
    use crate::hash_index::HashKmerIndex;
    use kam_core::kmer::KmerIndex;

    fn read(seq: &[u8], err: f32, family_type: FamilyType) -> ConsensusReadInfo {
        ConsensusReadInfo {
            sequence: seq.to_vec(),
            per_base_error_prob: vec![err; seq.len()],
            family_type,
        }
    }

    // Test 1: Single consensus read produces L-k+1 k-mers in the index.
    #[test]
    fn single_read_produces_correct_kmer_count() {
        let seq = b"ACGTACGT"; // length 8
        let k = 4;
        let r = read(seq, 0.01, FamilyType::Duplex);
        let mut idx = HashKmerIndex::new();
        extract_and_index(&r, k, &mut idx);

        // 8 - 4 + 1 = 5 k-mer positions; canonical forms may merge some.
        // All positions must produce at least one entry.
        assert!(!idx.is_empty());

        // Count unique canonical k-mers from the sequence to verify exactly.
        use crate::encode::KmerIterator;
        use std::collections::HashSet;
        let unique_canonical: HashSet<u64> = KmerIterator::new(seq, k)
            .map(|(_, km)| canonical(km, k))
            .collect();
        assert_eq!(idx.len(), unique_canonical.len());
    }

    // Test 2: Duplex read increments n_duplex.
    #[test]
    fn duplex_read_increments_n_duplex() {
        let seq = b"AAAA";
        let k = 4;
        let r = read(seq, 0.01, FamilyType::Duplex);
        let mut idx = HashKmerIndex::new();
        extract_and_index(&r, k, &mut idx);

        let kmer = canonical(encode_kmer(seq).unwrap(), k);
        let ev = idx.get(kmer).unwrap();
        assert_eq!(ev.n_duplex, 1);
        assert_eq!(ev.n_simplex_fwd, 0);
        assert_eq!(ev.n_simplex_rev, 0);
    }

    // Test 3: SimplexFwd read increments n_simplex_fwd.
    #[test]
    fn simplex_fwd_read_increments_n_simplex_fwd() {
        let seq = b"CCCC";
        let k = 4;
        let r = read(seq, 0.02, FamilyType::SimplexFwd);
        let mut idx = HashKmerIndex::new();
        extract_and_index(&r, k, &mut idx);

        let kmer = canonical(encode_kmer(seq).unwrap(), k);
        let ev = idx.get(kmer).unwrap();
        assert_eq!(ev.n_simplex_fwd, 1);
        assert_eq!(ev.n_duplex, 0);
    }

    // Test 4: Two reads with overlapping k-mers: n_molecules incremented for shared k-mers.
    #[test]
    fn overlapping_kmers_increment_n_molecules() {
        let shared_seq = b"ACGT";
        let k = 4;
        // Both reads contain the same single 4-mer.
        let r1 = read(shared_seq, 0.01, FamilyType::SimplexFwd);
        let r2 = read(shared_seq, 0.02, FamilyType::SimplexRev);
        let mut idx = HashKmerIndex::new();
        extract_and_index(&r1, k, &mut idx);
        extract_and_index(&r2, k, &mut idx);

        let kmer = canonical(encode_kmer(shared_seq).unwrap(), k);
        let ev = idx.get(kmer).unwrap();
        assert_eq!(ev.n_molecules, 2);
    }

    // Test 5: min_base_error_prob tracks the minimum across all observations.
    #[test]
    fn min_base_error_prob_tracks_minimum() {
        let seq = b"AAAA";
        let k = 4;
        let kmer = canonical(encode_kmer(seq).unwrap(), k);

        let r1 = read(seq, 0.05, FamilyType::SimplexFwd);
        let r2 = read(seq, 0.001, FamilyType::SimplexFwd); // lower error prob
        let r3 = read(seq, 0.1, FamilyType::SimplexFwd);

        let mut idx = HashKmerIndex::new();
        extract_and_index(&r1, k, &mut idx);
        extract_and_index(&r2, k, &mut idx);
        extract_and_index(&r3, k, &mut idx);

        let ev = idx.get(kmer).unwrap();
        assert!(
            (ev.min_base_error_prob - 0.001).abs() < 1e-6,
            "min was {}",
            ev.min_base_error_prob
        );
    }

    // Test 6: K-mers with N bases are skipped.
    #[test]
    fn kmers_with_n_are_skipped() {
        // "ACNT" contains an N — KmerIterator will skip it.
        let seq = b"ACNT";
        let k = 3;
        let r = read(seq, 0.01, FamilyType::Duplex);
        let mut idx = HashKmerIndex::new();
        extract_and_index(&r, k, &mut idx);

        // Only "ACN" and "CNT" contain N; "ACN" fails (N at pos 2) and "CNT" fails.
        // In reality KmerIterator resets at N, so no k-mer of length 3 can complete.
        // ACN: N at pos 2 → window reset. CNT: pos 1=C, pos 2=N → reset. NT: only 2 chars.
        assert!(
            idx.is_empty(),
            "no k-mers should be indexed when N prevents full windows"
        );
    }

    // Test 7: Empty sequence produces no k-mers.
    #[test]
    fn empty_sequence_produces_no_kmers() {
        let r = ConsensusReadInfo {
            sequence: vec![],
            per_base_error_prob: vec![],
            family_type: FamilyType::Duplex,
        };
        let mut idx = HashKmerIndex::new();
        extract_and_index(&r, 4, &mut idx);
        assert!(idx.is_empty());
    }

    // Test 8: extract_all processes multiple reads correctly.
    #[test]
    fn extract_all_processes_multiple_reads() {
        let r1 = read(b"AAAA", 0.01, FamilyType::Duplex);
        let r2 = read(b"CCCC", 0.02, FamilyType::SimplexFwd);
        let r3 = read(b"AAAA", 0.03, FamilyType::SimplexRev); // overlaps with r1

        let k = 4;
        let mut idx = HashKmerIndex::new();
        extract_all(&[r1, r2, r3], k, &mut idx);

        let kmer_aaaa = canonical(encode_kmer(b"AAAA").unwrap(), k);
        let kmer_cccc = canonical(encode_kmer(b"CCCC").unwrap(), k);

        // AAAA was seen by r1 and r3.
        assert_eq!(idx.molecule_count(kmer_aaaa), 2);
        // CCCC was seen only by r2.
        assert_eq!(idx.molecule_count(kmer_cccc), 1);
    }
}
