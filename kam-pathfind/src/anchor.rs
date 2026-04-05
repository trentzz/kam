//! Anchor k-mer validation for targeted path walking.
//!
//! Validates that the start and end k-mers of a target sequence are sufficiently
//! unique in the index to support reliable path walking. Non-unique anchors in
//! repeat regions can produce spurious paths.

use kam_core::kmer::KmerIndex;
use kam_index::encode::{canonical, encode_kmer};

/// Default threshold for anchor uniqueness.
///
/// Anchors with molecule counts at or above this value are considered non-unique.
pub const DEFAULT_ANCHOR_THRESHOLD: u32 = 100;

/// Result of anchor validation for a target sequence.
///
/// Holds the encoded canonical k-mers at each end of the target, their molecule
/// counts from the index, uniqueness flags, and an optional warning message.
#[derive(Debug, Clone)]
pub struct AnchorValidation {
    /// Canonical encoding of the start (first k bases) anchor k-mer.
    pub start_kmer: u64,
    /// Canonical encoding of the end (last k bases) anchor k-mer.
    pub end_kmer: u64,
    /// Molecule count of the start anchor in the index (0 if absent).
    pub start_count: u32,
    /// Molecule count of the end anchor in the index (0 if absent).
    pub end_count: u32,
    /// `true` if `start_count < threshold` (exclusive upper bound).
    pub start_unique: bool,
    /// `true` if `end_count < threshold` (exclusive upper bound).
    pub end_unique: bool,
    /// Warning message if either anchor is non-unique; `None` otherwise.
    pub warning: Option<String>,
}

/// Validate anchor k-mers for a target sequence.
///
/// Extracts the first and last k-mer from `target_seq`, canonicalises them,
/// and looks up their molecule counts in `index`. An anchor is considered
/// non-unique when its count is **≥ threshold**.
///
/// Returns `None` if `target_seq` is shorter than `k`.
///
/// # Arguments
///
/// - `target_seq`: raw DNA bytes for the target region
/// - `k`: k-mer length
/// - `index`: the k-mer index to query
/// - `threshold`: molecule count above which an anchor is non-unique (exclusive)
///
/// # Example
///
/// ```
/// use kam_pathfind::anchor::{validate_anchors, DEFAULT_ANCHOR_THRESHOLD};
/// use kam_index::HashKmerIndex;
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
/// use kam_index::encode::{encode_kmer, canonical};
///
/// let mut index = HashKmerIndex::new();
/// // Insert low-count evidence for the start and end k-mers.
/// let start = canonical(encode_kmer(b"ACG").unwrap(), 3);
/// let end   = canonical(encode_kmer(b"GTA").unwrap(), 3);
/// index.insert(start, MoleculeEvidence { n_molecules: 1, ..Default::default() });
/// index.insert(end,   MoleculeEvidence { n_molecules: 1, ..Default::default() });
///
/// let result = validate_anchors(b"ACGGTA", 3, &index, DEFAULT_ANCHOR_THRESHOLD).unwrap();
/// assert!(result.start_unique);
/// assert!(result.end_unique);
/// assert!(result.warning.is_none());
/// ```
pub fn validate_anchors(
    target_seq: &[u8],
    k: usize,
    index: &dyn KmerIndex,
    threshold: u32,
) -> Option<AnchorValidation> {
    if target_seq.len() < k {
        return None;
    }

    // Extract and canonicalise start and end k-mers.
    let start_raw = encode_kmer(&target_seq[..k])?;
    let end_raw = encode_kmer(&target_seq[target_seq.len() - k..])?;

    let start_kmer = canonical(start_raw, k);
    let end_kmer = canonical(end_raw, k);

    let start_count = index.molecule_count(start_kmer);
    let end_count = index.molecule_count(end_kmer);

    let start_unique = start_count < threshold;
    let end_unique = end_count < threshold;

    let warning = match (start_unique, end_unique) {
        (true, true) => None,
        (false, true) => Some(format!(
            "Start anchor has high molecule count ({start_count} >= {threshold}); \
             path walking may be unreliable"
        )),
        (true, false) => Some(format!(
            "End anchor has high molecule count ({end_count} >= {threshold}); \
             path walking may be unreliable"
        )),
        (false, false) => Some(format!(
            "Both anchors have high molecule counts \
             (start={start_count}, end={end_count}, threshold={threshold}); \
             path walking may be unreliable"
        )),
    };

    Some(AnchorValidation {
        start_kmer,
        end_kmer,
        start_count,
        end_count,
        start_unique,
        end_unique,
        warning,
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use kam_core::kmer::MoleculeEvidence;
    use kam_index::encode::encode_kmer;
    use kam_index::HashKmerIndex;

    /// Build an index that maps each of the provided (kmer, count) pairs.
    fn make_index(entries: &[(u64, u32)]) -> HashKmerIndex {
        let mut index = HashKmerIndex::new();
        for &(km, n) in entries {
            index.insert(
                km,
                MoleculeEvidence {
                    n_molecules: n,
                    ..Default::default()
                },
            );
        }
        index
    }

    // Helper: canonical encode from raw bytes.
    fn ck(seq: &[u8]) -> u64 {
        canonical(encode_kmer(seq).unwrap(), seq.len())
    }

    // Test 1: Unique anchors (low count) → both start_unique and end_unique true, no warning.
    #[test]
    fn unique_anchors_no_warning() {
        let k = 3;
        let start = ck(b"ACG");
        let end = ck(b"GTA");
        let index = make_index(&[(start, 1), (end, 2)]);
        let result = validate_anchors(b"ACGGTA", k, &index, 100).unwrap();
        assert!(result.start_unique);
        assert!(result.end_unique);
        assert!(result.warning.is_none());
    }

    // Test 2: Non-unique start anchor → start_unique false, warning present.
    #[test]
    fn non_unique_start_anchor_warns() {
        let k = 3;
        let start = ck(b"ACG");
        let end = ck(b"GTA");
        let index = make_index(&[(start, 200), (end, 5)]);
        let result = validate_anchors(b"ACGGTA", k, &index, 100).unwrap();
        assert!(!result.start_unique);
        assert!(result.end_unique);
        assert!(result.warning.is_some());
        let w = result.warning.unwrap();
        assert!(
            w.contains("Start anchor"),
            "warning should mention Start anchor: {w}"
        );
    }

    // Test 3: Non-unique end anchor → end_unique false, warning present.
    #[test]
    fn non_unique_end_anchor_warns() {
        let k = 3;
        let start = ck(b"ACG");
        let end = ck(b"GTA");
        let index = make_index(&[(start, 5), (end, 150)]);
        let result = validate_anchors(b"ACGGTA", k, &index, 100).unwrap();
        assert!(result.start_unique);
        assert!(!result.end_unique);
        assert!(result.warning.is_some());
        let w = result.warning.unwrap();
        assert!(
            w.contains("End anchor"),
            "warning should mention End anchor: {w}"
        );
    }

    // Test 4: Both non-unique → both false, warning mentions both.
    #[test]
    fn both_non_unique_warns_both() {
        let k = 3;
        let start = ck(b"ACG");
        let end = ck(b"GTA");
        let index = make_index(&[(start, 200), (end, 300)]);
        let result = validate_anchors(b"ACGGTA", k, &index, 100).unwrap();
        assert!(!result.start_unique);
        assert!(!result.end_unique);
        let w = result.warning.unwrap();
        // The warning must reference both anchors.
        assert!(
            w.contains("start") && w.contains("end"),
            "warning should mention both anchors: {w}"
        );
    }

    // Test 5: Target shorter than k → returns None.
    #[test]
    fn target_shorter_than_k_returns_none() {
        let index = make_index(&[]);
        assert!(validate_anchors(b"AC", 3, &index, 100).is_none());
        assert!(validate_anchors(b"", 3, &index, 100).is_none());
    }

    // Test 6: Threshold of 0 → everything is non-unique.
    #[test]
    fn threshold_zero_everything_non_unique() {
        let k = 3;
        let start = ck(b"ACG");
        let end = ck(b"GTA");
        // Even count=0 (absent from index) is >= 0, so non-unique.
        let index = make_index(&[(start, 0), (end, 0)]);
        let result = validate_anchors(b"ACGGTA", k, &index, 0).unwrap();
        assert!(!result.start_unique);
        assert!(!result.end_unique);
        assert!(result.warning.is_some());
    }

    // Test 7: Anchor at exactly threshold → unique (threshold is exclusive upper bound).
    #[test]
    fn anchor_exactly_at_threshold_is_unique() {
        let k = 3;
        let start = ck(b"ACG");
        let end = ck(b"GTA");
        let index = make_index(&[(start, 100), (end, 100)]);
        // count=100 with threshold=101 → unique (100 < 101)
        let result = validate_anchors(b"ACGGTA", k, &index, 101).unwrap();
        assert!(
            result.start_unique,
            "count=100 < threshold=101 should be unique"
        );
        assert!(
            result.end_unique,
            "count=100 < threshold=101 should be unique"
        );
        assert!(result.warning.is_none());

        // count=100 with threshold=100 → NOT unique (100 is not < 100)
        let result2 = validate_anchors(b"ACGGTA", k, &index, 100).unwrap();
        assert!(
            !result2.start_unique,
            "count=100 >= threshold=100 should be non-unique"
        );
        assert!(
            !result2.end_unique,
            "count=100 >= threshold=100 should be non-unique"
        );
    }
}
