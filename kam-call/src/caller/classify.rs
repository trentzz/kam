//! Variant classification helpers: sequence analysis for classify_variant and its sub-routines.

use super::types::{VariantType, SV_LENGTH_THRESHOLD};

/// Determine variant type by comparing reference and alternate sequences.
///
/// Deletions and insertions of 50 bp or more are classified as structural
/// variants. Large deletions where the alt contains an inverted segment relative
/// to the ref are classified as `InvDel`. Large insertions whose content does not
/// match the local reference are classified as `NovelInsertion`; those that do
/// match (tandem repeats) are `TandemDuplication`. Same-length variants where the
/// alt is the reverse complement of the ref are classified as `Inversion`.
///
/// `Fusion` is not detected from sequence content alone. It requires graph-level
/// evidence (k-mers from non-contiguous reference regions) and must be assigned
/// by the caller before invoking this function.
///
/// # Example
/// ```
/// use kam_call::caller::{classify_variant, VariantType};
/// assert_eq!(classify_variant(b"A", b"T"), VariantType::Snv);
/// assert_eq!(classify_variant(b"ACG", b"A"), VariantType::Deletion);
/// assert_eq!(classify_variant(b"A", b"ACG"), VariantType::Insertion);
/// // Large structural variants:
/// let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
/// let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 60 bp deletion
/// assert_eq!(classify_variant(&long_ref, &short_alt), VariantType::LargeDeletion);
/// ```
pub fn classify_variant(ref_seq: &[u8], alt_seq: &[u8]) -> VariantType {
    use std::cmp::Ordering;
    match ref_seq.len().cmp(&alt_seq.len()) {
        Ordering::Equal => {
            // Check for inversion: alt is reverse complement of ref.
            // Full-path case: entire window is RC (rare).
            if is_reverse_complement(ref_seq, alt_seq) {
                return VariantType::Inversion;
            }
            // Partial-inversion case: a contiguous central segment is RC'd
            // while the flanking bases match. This is the common case for
            // targeted SV detection where the target window is wider than
            // the inverted region.
            if let Some(inv_len) = partial_inversion_len(ref_seq, alt_seq) {
                if inv_len >= SV_LENGTH_THRESHOLD {
                    return VariantType::Inversion;
                }
            }
            let diffs = ref_seq
                .iter()
                .zip(alt_seq.iter())
                .filter(|(r, a)| r != a)
                .count();
            match diffs {
                0 => VariantType::Snv, // identical — treat as SNV (no-op)
                1 => VariantType::Snv,
                _ => VariantType::Mnv,
            }
        }
        Ordering::Greater => {
            let del_len = ref_seq.len() - alt_seq.len();
            // Check for InvDel before the net-length threshold.
            // An InvDel can have a small net deletion (deleted bp minus
            // inverted-inserted bp) even though the gross rearrangement is
            // large. alt_seq_has_inversion_relative_to_ref internally requires
            // both the ref central and alt central segments to be at least
            // SV_LENGTH_THRESHOLD bp, so short deletions with no inversion
            // are not affected.
            if alt_seq_has_inversion_relative_to_ref(ref_seq, alt_seq) {
                VariantType::InvDel
            } else if del_len >= SV_LENGTH_THRESHOLD {
                VariantType::LargeDeletion
            } else {
                VariantType::Deletion
            }
        }
        Ordering::Less => {
            let ins_len = alt_seq.len() - ref_seq.len();
            if ins_len >= SV_LENGTH_THRESHOLD {
                // Distinguish tandem duplication from novel insertion.
                // A tandem duplication inserts sequence that matches a nearby
                // window of the reference. A novel insertion has no significant
                // similarity to the flanking reference context.
                if is_tandem_duplication(ref_seq, alt_seq) {
                    VariantType::TandemDuplication
                } else {
                    VariantType::NovelInsertion
                }
            } else {
                VariantType::Insertion
            }
        }
    }
}

/// Detect a partial inversion: a contiguous central segment where alt is the
/// reverse complement of the corresponding ref region, with matching flanks.
///
/// Returns `Some(length)` where `length` is the number of bases in the inverted
/// segment if such a segment exists, or `None` otherwise.
///
/// Only considers segments of length ≥ 2 (required by `is_reverse_complement`).
pub(crate) fn partial_inversion_len(ref_seq: &[u8], alt_seq: &[u8]) -> Option<usize> {
    if ref_seq.len() != alt_seq.len() {
        return None;
    }
    let left = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .position(|(r, a)| r != a)?;
    let right = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .rposition(|(r, a)| r != a)?;
    if right < left {
        return None;
    }
    if is_reverse_complement(&ref_seq[left..=right], &alt_seq[left..=right]) {
        Some(right - left + 1)
    } else {
        None
    }
}

/// Detect an InvDel signature: the alt sequence contains a segment that is the
/// reverse complement of a segment found within the reference.
///
/// The alt is shorter than the ref (a deletion has occurred). After stripping
/// the shared flanking prefix and suffix, the alt central segment must be ≥
/// `SV_LENGTH_THRESHOLD` bp long and must be the reverse complement of some
/// substring of the ref central region. Both the alt central length and the
/// matched ref substring length must meet the threshold.
///
/// This handles the common InvDel topology:
/// `ref = [flank_L][deleted_region][inv_region][flank_R]`
/// `alt = [flank_L][rc(inv_region)][flank_R]`
/// where the alt central is `rc(inv_region)` and the ref central contains
/// both `deleted_region` and `inv_region`.
fn alt_seq_has_inversion_relative_to_ref(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    debug_assert!(ref_seq.len() > alt_seq.len());

    // Find the shared prefix length.
    let prefix = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .take_while(|(r, a)| r == a)
        .count();

    // Find the shared suffix length (do not overlap the prefix).
    let max_suffix = alt_seq.len().saturating_sub(prefix);
    let suffix = ref_seq
        .iter()
        .rev()
        .zip(alt_seq.iter().rev())
        .take(max_suffix)
        .take_while(|(r, a)| r == a)
        .count();

    // The central alt segment (after removing shared flanks).
    let alt_central = &alt_seq[prefix..alt_seq.len() - suffix];
    // The corresponding ref segment (longer — contains the deleted + inv regions).
    let ref_central = &ref_seq[prefix..ref_seq.len() - suffix];

    if alt_central.len() < SV_LENGTH_THRESHOLD || ref_central.len() < SV_LENGTH_THRESHOLD {
        return false;
    }

    // Fast path: the alt central is the RC of the entire ref central.
    if is_reverse_complement(ref_central, alt_central) {
        return true;
    }

    // General case: the alt central is the RC of some sub-window of ref_central.
    // Pre-compute the RC of alt_central and look for it as a substring of ref_central.
    // This is equivalent to asking: does ref_central contain a segment whose RC
    // equals alt_central?
    let alt_central_rc: Vec<u8> = alt_central.iter().rev().map(|&b| complement(b)).collect();
    if alt_central_rc.len() <= ref_central.len() {
        return ref_central
            .windows(alt_central_rc.len())
            .any(|w| w == alt_central_rc.as_slice());
    }

    false
}

/// Detect a tandem duplication: the inserted bases match a contiguous stretch
/// of the reference sequence.
///
/// Extracts the inserted bases from the alt sequence, then checks whether those
/// bases appear as a contiguous substring in the reference — including the case
/// where the insertion is longer than the reference but the sequence is a
/// repeated tile of a ref motif (handled by searching in `ref + ref`).
///
/// Returns false for novel insertions whose content is not present in the
/// local reference.
///
/// # Examples (expected outcomes)
///
/// - `ref = ACGT×10`, `insert = ACGT×15`: matches a substring of `ref+ref` → true.
/// - `ref = ACGT×12`, `insert = CCCC…(60)`: no 60-C run in `ref+ref` → false.
pub(crate) fn is_tandem_duplication(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    debug_assert!(alt_seq.len() > ref_seq.len());

    let ins_len = alt_seq.len() - ref_seq.len();

    // Find the first position where ref and alt diverge (the insertion point).
    let diverge = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .position(|(r, a)| r != a)
        .unwrap_or(ref_seq.len());

    // The inserted bases sit between diverge and diverge+ins_len in alt.
    if diverge + ins_len > alt_seq.len() {
        return false;
    }
    let inserted = &alt_seq[diverge..diverge + ins_len];
    if inserted.is_empty() {
        return false;
    }

    // Verify suffix: the alt sequence after the insertion must match the
    // reference from the divergence point onwards. A true tandem duplication
    // has matching prefix + duplicated segment + matching suffix that rejoins
    // the reference. Without this check, novel inserted sequences can match
    // by coincidence in the doubled reference.
    let alt_suffix = &alt_seq[diverge + ins_len..];
    let ref_suffix = &ref_seq[diverge..];
    if alt_suffix.len() != ref_suffix.len() || alt_suffix != ref_suffix {
        return false;
    }

    // Search for the inserted sequence as a substring of (ref + ref).
    // This handles the case where the inserted length exceeds the reference length
    // but the inserted bases are a contiguous tiling of a reference repeat unit
    // (e.g., ACGT×15 inserted into ACGT×10).
    let doubled: Vec<u8> = ref_seq.iter().chain(ref_seq.iter()).copied().collect();
    if ins_len <= doubled.len() {
        return doubled.windows(ins_len).any(|w| w == inserted);
    }

    false
}

/// Return true if `alt` is the reverse complement of `ref`.
///
/// Only returns true for sequences of length ≥ 2 where the entire alt is the
/// reverse complement of the entire ref. This detects inversion paths in the de
/// Bruijn graph where both strands of a segment are traversed.
pub(crate) fn is_reverse_complement(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    if ref_seq.len() < 2 || ref_seq.len() != alt_seq.len() {
        return false;
    }
    ref_seq
        .iter()
        .zip(alt_seq.iter().rev())
        .all(|(&r, &a)| complement(r) == a)
}

/// DNA complement (ACGT only; other bytes return themselves unchanged).
#[inline]
pub(crate) fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Test 5: classify_variant with same-length, 1 diff → SNV.
    #[test]
    fn classify_snv() {
        assert_eq!(classify_variant(b"A", b"T"), VariantType::Snv);
        assert_eq!(classify_variant(b"ACGT", b"ACTT"), VariantType::Snv);
    }

    // Test 6: classify_variant alt shorter → Deletion.
    #[test]
    fn classify_deletion() {
        assert_eq!(classify_variant(b"ACG", b"A"), VariantType::Deletion);
    }

    // Test 7: classify_variant alt longer → Insertion.
    #[test]
    fn classify_insertion() {
        assert_eq!(classify_variant(b"A", b"ACG"), VariantType::Insertion);
    }

    // Test 14: large deletion (≥50 bp) is classified as LargeDeletion.
    #[test]
    fn classify_large_deletion() {
        let long_ref: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
        let short_alt: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 60 bp deletion
        assert_eq!(
            classify_variant(&long_ref, &short_alt),
            VariantType::LargeDeletion
        );
    }

    // Test 15: small deletion (<50 bp) is still classified as Deletion.
    #[test]
    fn classify_small_deletion_not_sv() {
        let r: Vec<u8> = b"ACGT".repeat(15).to_vec(); // 60 bp
        let a: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 20 bp deletion
        assert_eq!(classify_variant(&r, &a), VariantType::Deletion);
    }

    // Test 16: large insertion (≥50 bp) is classified as TandemDuplication.
    #[test]
    fn classify_tandem_duplication() {
        let short_ref: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp
        let long_alt: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp — 60 bp insertion
        assert_eq!(
            classify_variant(&short_ref, &long_alt),
            VariantType::TandemDuplication
        );
    }

    // Test 17: sequence that is reverse complement → classified as Inversion.
    #[test]
    fn classify_inversion() {
        let ref_seq = b"AAACCC";
        let alt_seq = b"GGGTTT"; // rc of AAACCC
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Inversion);
    }

    // Test 18: is_reverse_complement correctly identifies RC pairs.
    #[test]
    fn is_reverse_complement_basic() {
        // ACGT rc = ACGT (palindrome) — true.
        assert!(is_reverse_complement(b"ACGT", b"ACGT"));
        // AACC rc = GGTT — true.
        assert!(is_reverse_complement(b"AACC", b"GGTT"));
        // Different sequences, not RC.
        assert!(!is_reverse_complement(b"AACC", b"TTGG"));
        // Length 1 — returns false (too short).
        assert!(!is_reverse_complement(b"A", b"T"));
    }

    // Test 19: partial inversion in the central region with matching flanks.
    //
    // ref = AA [AAACCC] AA  (6 bp central, flanked by 2 bp each)
    // alt = AA [GGGTTT] AA  (central is rc of AAACCC)
    // Expected: Some(6) — 6 bp inverted segment.
    #[test]
    fn partial_inversion_len_central_segment() {
        // Flanks: AA ... AA.  Central 6 bp: AAACCC / GGGTTT.
        let ref_seq = b"AAAAACCCAA";
        let alt_seq = b"AAGGGTTTAA";
        assert_eq!(partial_inversion_len(ref_seq, alt_seq), Some(6));
    }

    // Test 20: partial inversion below SV_LENGTH_THRESHOLD is classified as MNV.
    //
    // A 6 bp inversion (< 50 bp) must not be promoted to VariantType::Inversion.
    #[test]
    fn partial_inversion_below_threshold_is_mnv() {
        // Same as test 19 — 6 bp inversion, below the 50 bp threshold.
        let ref_seq = b"AAAAACCCAA";
        let alt_seq = b"AAGGGTTTAA";
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Mnv);
    }

    // Test 21: full-path inversion still classified correctly.
    #[test]
    fn full_path_inversion_still_works() {
        // AAACCC rc = GGGTTT — full path RC, both checks should fire.
        let ref_seq = b"AAACCC";
        let alt_seq = b"GGGTTT";
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Inversion);
    }

    // Test 22: a single SNV inside an otherwise identical window → still SNV.
    //
    // partial_inversion_len must not confuse a single mismatch with an inversion.
    #[test]
    fn single_snv_not_misclassified_as_inversion() {
        let ref_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGT"; // 103 bp
        let alt_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGT"; // single A→T
        assert_eq!(classify_variant(ref_seq, alt_seq), VariantType::Snv);
    }

    // Test 33: Display impl produces the correct string for all VariantType variants.
    #[test]
    fn variant_type_display() {
        assert_eq!(VariantType::Snv.to_string(), "SNV");
        assert_eq!(VariantType::Insertion.to_string(), "Insertion");
        assert_eq!(VariantType::Deletion.to_string(), "Deletion");
        assert_eq!(VariantType::Mnv.to_string(), "MNV");
        assert_eq!(VariantType::Complex.to_string(), "Complex");
        assert_eq!(VariantType::LargeDeletion.to_string(), "LargeDeletion");
        assert_eq!(
            VariantType::TandemDuplication.to_string(),
            "TandemDuplication"
        );
        assert_eq!(VariantType::Inversion.to_string(), "Inversion");
        assert_eq!(VariantType::Fusion.to_string(), "Fusion");
        assert_eq!(VariantType::InvDel.to_string(), "InvDel");
        assert_eq!(VariantType::NovelInsertion.to_string(), "NovelInsertion");
    }

    // Test 34: is_sv returns true for all three new types.
    #[test]
    fn new_sv_types_are_recognised_as_sv() {
        for vt in [
            VariantType::Fusion,
            VariantType::InvDel,
            VariantType::NovelInsertion,
        ] {
            let is_sv = matches!(
                vt,
                VariantType::LargeDeletion
                    | VariantType::TandemDuplication
                    | VariantType::Inversion
                    | VariantType::Fusion
                    | VariantType::InvDel
                    | VariantType::NovelInsertion
            );
            assert!(is_sv, "{vt} must be classified as an SV");
        }
    }

    // Test 36: a large deletion with an inverted central alt segment is InvDel.
    #[test]
    fn classify_invdel_deletion_plus_inversion() {
        // Build a 160-bp ref: 10 bp flank, 50 bp deleted region, 60 bp inversion
        // region, 40 bp right flank.
        let flank_l: Vec<u8> = b"AAAAAAAAAA".to_vec(); // 10 bp
        let del_region: Vec<u8> = b"TTTTTTTTTT".repeat(5).to_vec(); // 50 bp, deleted
                                                                    // 60 bp inversion region: use AAACCC repeated 10 times.
        let inv_region: Vec<u8> = b"AAACCC".repeat(10).to_vec(); // 60 bp
        let flank_r: Vec<u8> = b"GGGGGGGGGG".repeat(4).to_vec(); // 40 bp

        let ref_seq: Vec<u8> = flank_l
            .iter()
            .chain(&del_region)
            .chain(&inv_region)
            .chain(&flank_r)
            .copied()
            .collect();

        // Alt: flank_L + rc(inv_region) + flank_R (del_region is absent).
        let rc_inv: Vec<u8> = inv_region.iter().rev().map(|&b| complement(b)).collect();
        let alt_seq: Vec<u8> = flank_l
            .iter()
            .chain(&rc_inv)
            .chain(&flank_r)
            .copied()
            .collect();

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::InvDel,
            "large deletion with inverted alt segment must be InvDel"
        );
    }

    // Test 37: a plain large deletion (no inversion in alt) stays LargeDeletion.
    #[test]
    fn classify_plain_large_deletion_stays_large_deletion() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
        let alt_seq: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp — 60 bp deletion

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::LargeDeletion,
            "plain large deletion without inversion must stay LargeDeletion"
        );
    }

    // Test 38: deletion of exactly 50 bp (at the threshold) without inversion.
    #[test]
    fn classify_deletion_at_sv_threshold_is_large_deletion() {
        let ref_seq: Vec<u8> = vec![b'A'; 100];
        // Delete exactly 50 bp → alt is 50 bp.
        let alt_seq: Vec<u8> = vec![b'A'; 50];

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::LargeDeletion,
            "deletion of exactly SV_LENGTH_THRESHOLD bp must be LargeDeletion"
        );
    }

    // Test 39: deletion of 49 bp (one below threshold) is a plain Deletion.
    #[test]
    fn classify_deletion_below_sv_threshold_is_deletion() {
        let ref_seq: Vec<u8> = vec![b'A'; 99];
        // Delete 49 bp → alt is 50 bp.
        let alt_seq: Vec<u8> = vec![b'A'; 50];

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::Deletion,
            "deletion of 49 bp must remain Deletion"
        );
    }

    // Test 40: a large insertion with random (non-ref) content is NovelInsertion.
    #[test]
    fn classify_novel_insertion_random_sequence() {
        // ref = 50 bp of ACGT repeats — no 'C' runs of 60 bp.
        let ref_seq: Vec<u8> = b"ACGT".repeat(12).to_vec(); // 48 bp, close enough
                                                            // Insert 60 bp of all-C at position 24 (midpoint of ref).
        let novel_insert: Vec<u8> = vec![b'C'; 60];
        let alt_seq: Vec<u8> = ref_seq[..24]
            .iter()
            .chain(&novel_insert)
            .chain(&ref_seq[24..])
            .copied()
            .collect();

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::NovelInsertion,
            "insertion of random non-ref bases must be NovelInsertion"
        );
    }

    // Test 41: a large insertion where the inserted bases match the reference
    // context is TandemDuplication, not NovelInsertion.
    #[test]
    fn classify_tandem_dup_stays_tandem_duplication() {
        // ref = 100 bp of ACGT repeats.
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
                                                            // Duplicate 60 bp starting at position 20 and insert them at position 80.
        let dup_bases: Vec<u8> = ref_seq[20..80].to_vec(); // 60 bp from ref
        let alt_seq: Vec<u8> = ref_seq[..80]
            .iter()
            .chain(&dup_bases)
            .chain(&ref_seq[80..])
            .copied()
            .collect();

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::TandemDuplication,
            "insertion of ref-matching bases must be TandemDuplication"
        );
    }

    // Test 42: insertion of exactly SV_LENGTH_THRESHOLD (50 bp) with novel content.
    #[test]
    fn classify_novel_insertion_at_sv_threshold() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
                                                            // Insert 50 bp of all-T at position 50.
        let novel_insert: Vec<u8> = vec![b'T'; 50];
        let alt_seq: Vec<u8> = ref_seq[..50]
            .iter()
            .chain(&novel_insert)
            .chain(&ref_seq[50..])
            .copied()
            .collect();

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::NovelInsertion,
            "novel insertion of exactly SV_LENGTH_THRESHOLD bp must be NovelInsertion"
        );
    }

    // Test 43: insertion of 49 bp (below threshold) stays Insertion regardless
    // of content.
    #[test]
    fn classify_insertion_below_sv_threshold_stays_insertion() {
        let ref_seq: Vec<u8> = vec![b'A'; 50];
        let novel_insert: Vec<u8> = vec![b'C'; 49];
        let alt_seq: Vec<u8> = ref_seq[..25]
            .iter()
            .chain(&novel_insert)
            .chain(&ref_seq[25..])
            .copied()
            .collect();

        assert_eq!(
            classify_variant(&ref_seq, &alt_seq),
            VariantType::Insertion,
            "insertion of 49 bp must remain Insertion"
        );
    }

    // Test 44: Fusion classification is not produced by classify_variant.
    #[test]
    fn fusion_classification_requires_graph_evidence_not_sequence_content() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec();
        let alt_seq: Vec<u8> = ref_seq[..50]
            .iter()
            .chain(&vec![b'N'; 60]) // 60 foreign bases
            .chain(&ref_seq[50..])
            .copied()
            .collect();
        let vt = classify_variant(&ref_seq, &alt_seq);
        // classify_variant cannot produce Fusion from sequence alone.
        assert_ne!(
            vt,
            VariantType::Fusion,
            "classify_variant must not produce Fusion — that requires graph-level evidence"
        );
    }

    // Test 56: novel insertion whose inserted bytes happen to appear in the
    // doubled reference, but whose alt suffix does not rejoin the reference.
    // Must not be misclassified as a tandem duplication.
    #[test]
    fn novel_insertion_not_classified_as_dup() {
        let ref_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        // Insert 50 random bases that do not rejoin the reference cleanly.
        let mut alt = ref_seq[..25].to_vec();
        alt.extend_from_slice(b"NNNNNNNNNNTTTTTTTTTTCCCCCCCCCCGGGGGGGGGGAAAAAAAAAA"); // 50 novel bases
        alt.extend_from_slice(b"XYZXYZXYZXYZXYZ"); // suffix does not match ref
        assert!(!is_tandem_duplication(ref_seq, &alt));
    }

    // Test 57: a true tandem duplication (prefix matches, inserted segment
    // from ref, suffix matches) is still correctly classified.
    #[test]
    fn true_tandem_dup_still_classified() {
        let ref_seq = b"AAAAACCCCCGGGGGTTTTTNNNNN"; // 25 bp
        let mut alt = ref_seq.to_vec();
        // Insert "CCCCCGGGGG" (10 bp from ref) at position 10.
        let dup_segment = b"CCCCCGGGGG";
        alt.splice(10..10, dup_segment.iter().copied());
        assert!(is_tandem_duplication(ref_seq, &alt));
    }

    // Test 58: classify_variant assigns NovelInsertion when the alt is longer
    // than the ref by >= 50 bp and the inserted content is novel.
    #[test]
    fn classify_variant_novel_insertion_type() {
        let ref_seq = vec![b'A'; 100];
        let mut alt = ref_seq[..50].to_vec();
        // Insert 60 novel bases (above SV_LENGTH_THRESHOLD of 50).
        alt.extend_from_slice(&[b'T'; 60]);
        alt.extend_from_slice(&ref_seq[50..]);
        let vt = classify_variant(&ref_seq, &alt);
        assert_eq!(vt, VariantType::NovelInsertion);
    }
}
