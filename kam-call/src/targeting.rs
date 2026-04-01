//! Tumour-informed filtering for monitoring mode.
//!
//! When a set of expected somatic variants is known in advance (e.g. from a
//! matched tumour biopsy or a reference standard truth set), `--target-variants`
//! restricts output to calls that match one of those expected alleles exactly.
//! All other PASS calls are relabelled [`VariantFilter::NotTargeted`].
//!
//! This replicates the tumour-informed filtering used in alignment-based
//! personalised ctDNA monitoring pipelines and reduces false positives to near
//! zero by ignoring background biological cfDNA variants.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::caller::{VariantCall, VariantFilter, VariantType};

/// A set of expected somatic variant keys: (chrom, pos, ref_allele, alt_allele).
///
/// The `pos` field uses whatever coordinate convention the VCF source uses
/// (parsed directly from the integer in column 2).  The variant extraction
/// from called paths uses the same convention, so matching is internally
/// consistent.
pub type TargetVariantSet = HashSet<(String, i64, String, String)>;

/// A set of (chrom, pos) pairs derived from a target VCF.
///
/// Used for position-based tumour-informed matching when `--ti-position-tolerance`
/// is non-zero.  A call passes if its (chrom, pos) is within the tolerance of any
/// entry in this set, regardless of REF/ALT.
pub type TargetPositionSet = HashSet<(String, i64)>;

/// Return true if any target entry is a BND record (ALT contains `]` or `[`).
fn targets_has_bnd(targets: &TargetVariantSet) -> bool {
    targets
        .iter()
        .any(|(_, _, _, alt)| alt.contains(']') || alt.contains('['))
}

/// Load expected variants from a VCF file.
///
/// Parses CHROM, POS, REF, and ALT from each non-header line.
/// Multi-allelic ALT columns are split on `','` and each allele is added
/// separately.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or any line is malformed.
///
/// # Example
///
/// ```no_run
/// use kam_call::targeting::load_target_variants;
/// let set = load_target_variants("truth.vcf").unwrap();
/// assert!(!set.is_empty());
/// ```
pub fn load_target_variants<P: AsRef<Path>>(
    path: P,
) -> Result<TargetVariantSet, Box<dyn std::error::Error>> {
    let file = File::open(path.as_ref()).map_err(|e| {
        format!(
            "cannot open --target-variants VCF '{}': {e}",
            path.as_ref().display()
        )
    })?;
    let reader = BufReader::new(file);
    let mut set = TargetVariantSet::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.splitn(6, '\t').collect();
        if parts.len() < 5 {
            return Err(format!(
                "--target-variants VCF line {}: expected ≥5 tab-separated fields, got {}",
                line_no + 1,
                parts.len()
            )
            .into());
        }
        let chrom = parts[0].to_string();
        let pos: i64 = parts[1].parse().map_err(|e| {
            format!(
                "--target-variants VCF line {}: cannot parse POS '{}': {e}",
                line_no + 1,
                parts[1]
            )
        })?;
        let reference = parts[3].to_string();
        // Split multi-allelic ALT and add each separately.
        for alt in parts[4].split(',') {
            set.insert((chrom.clone(), pos, reference.clone(), alt.to_string()));
        }
    }

    Ok(set)
}

/// Extract the variant key from a called path pair.
///
/// Given the `target_id` (e.g. `"chr17:7572800-7572900"`), the full reference
/// path sequence, and the full alternate path sequence, this function computes
/// the (chrom, pos, ref_allele, alt_allele) key using the same algorithm as
/// the benchmarking script's `extract_called_variants`.
///
/// The `pos` field is 1-based (standard VCF convention) so that keys produced
/// here can be compared directly against entries loaded by `load_target_variants`.
/// Target FASTA headers use 0-based coordinates, so `pos = target_start + offset + 1`.
///
/// Returns `None` if the sequences are identical (no variant) or if the
/// target_id cannot be parsed.
///
/// # Example
///
/// ```
/// use kam_call::targeting::extract_variant_key;
/// let ref_seq = b"ACGTACGT";
/// let alt_seq = b"ACTTACGT";
/// let key = extract_variant_key("chr1:100-108", ref_seq, alt_seq);
/// // SNV at 0-based window offset 2, target_start=100 (0-based) → VCF pos 103.
/// assert_eq!(key, Some(("chr1".to_string(), 103, "G".to_string(), "T".to_string())));
/// ```
pub fn extract_variant_key(
    target_id: &str,
    ref_seq: &[u8],
    alt_seq: &[u8],
) -> Option<(String, i64, String, String)> {
    // Parse "chrN:start-end" — start is whatever coordinate the target FASTA uses.
    let (chrom, target_start) = parse_target_id(target_id)?;

    // Find leftmost differing position.
    let diff_pos = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .position(|(r, a)| r != a)?;

    // Trim common suffix to get the minimal differing region.
    let ref_rem = &ref_seq[diff_pos..];
    let alt_rem = &alt_seq[diff_pos..];
    let common_suffix = ref_rem
        .iter()
        .rev()
        .zip(alt_rem.iter().rev())
        .take_while(|(r, a)| r == a)
        // Do not consume the entire remaining region — leave at least one base.
        .take(ref_rem.len().min(alt_rem.len()).saturating_sub(1))
        .count();

    let ref_trimmed = if common_suffix > 0 {
        &ref_seq[diff_pos..ref_seq.len() - common_suffix]
    } else {
        &ref_seq[diff_pos..]
    };
    let alt_trimmed = if common_suffix > 0 {
        &alt_seq[diff_pos..alt_seq.len() - common_suffix]
    } else {
        &alt_seq[diff_pos..]
    };

    if ref_seq.len() == alt_seq.len() {
        // SNV / MNV.  Convert 0-based target offset to 1-based VCF position.
        let ref_allele = bytes_to_str(ref_trimmed);
        let alt_allele = bytes_to_str(alt_trimmed);
        let genomic_pos = target_start + diff_pos as i64 + 1;
        Some((chrom, genomic_pos, ref_allele, alt_allele))
    } else {
        indel_key(
            chrom,
            target_start,
            diff_pos,
            ref_seq,
            alt_seq,
            ref_trimmed,
            alt_trimmed,
        )
    }
}

/// Apply tumour-informed filter to a list of calls.
///
/// Calls with `filter == Pass` that do not match any entry in `targets` are
/// relabelled [`VariantFilter::NotTargeted`].  All other calls are unchanged.
///
/// # Example
///
/// ```
/// use std::collections::HashSet;
/// use kam_call::caller::{VariantCall, VariantFilter, VariantType};
/// use kam_call::targeting::apply_target_filter;
///
/// let mut call = VariantCall {
///     target_id: "chr1:0-100".to_string(),
///     variant_type: VariantType::Snv,
///     ref_sequence: b"ACGTACGT".to_vec(),
///     alt_sequence: b"ACTTACGT".to_vec(),
///     vaf: 0.01,
///     vaf_ci_low: 0.005,
///     vaf_ci_high: 0.02,
///     n_molecules_ref: 99,
///     n_molecules_alt: 1,
///     n_duplex_alt: 0,
///     n_simplex_alt: 1,
///     n_simplex_fwd_alt: 1,
///     n_simplex_rev_alt: 0,
///     n_duplex_ref: 0,
///     n_simplex_ref: 99,
///     mean_alt_error_prob: 0.001,
///     min_variant_specific_duplex: 0,
///     mean_variant_specific_molecules: 1.0,
///     confidence: 0.99,
///     strand_bias_p: 0.5,
///     filter: VariantFilter::Pass,
/// };
///
/// let targets = HashSet::new(); // empty — nothing targeted
/// apply_target_filter(&mut [call], &targets);
/// ```
pub fn apply_target_filter(calls: &mut [VariantCall], targets: &TargetVariantSet) {
    let has_bnd_targets = targets_has_bnd(targets);
    for call in calls.iter_mut() {
        if call.filter != VariantFilter::Pass {
            continue;
        }
        // Fusion calls use junction FASTA headers as target_id; these cannot be
        // parsed as chrN:start-end coordinates. When the target VCF contains BND
        // records, treat all Fusion calls as targeted (pass through).
        if call.variant_type == VariantType::Fusion && has_bnd_targets {
            continue;
        }
        let key = extract_variant_key(&call.target_id, &call.ref_sequence, &call.alt_sequence);
        let is_targeted = key.is_some_and(|k| targets.contains(&k));
        if !is_targeted {
            call.filter = VariantFilter::NotTargeted;
        }
    }
}

/// Apply tumour-informed filter with optional position-based fallback.
///
/// When `tolerance > 0`, a PASS call is retained if:
/// 1. Its (CHROM, POS, REF, ALT) matches an entry in `targets` exactly, OR
/// 2. Its (CHROM, POS) is within `tolerance` bp of any target position.
///
/// When `tolerance == 0`, this is equivalent to `apply_target_filter`.
///
/// The position fallback is needed for large SVs where kam reports partial
/// alleles (e.g., a 1bp junction insertion) that cannot match a full SV truth
/// allele by exact REF/ALT comparison.
///
/// # Example
///
/// ```
/// use std::collections::HashSet;
/// use kam_call::caller::{VariantCall, VariantFilter, VariantType};
/// use kam_call::targeting::{apply_target_filter_with_tolerance, TargetVariantSet};
///
/// let call = VariantCall {
///     target_id: "chr1:490-600".to_string(),
///     variant_type: VariantType::Snv,
///     ref_sequence: b"ACGTACGT".to_vec(),
///     alt_sequence: b"ACTCACGT".to_vec(),
///     vaf: 0.01,
///     vaf_ci_low: 0.005,
///     vaf_ci_high: 0.02,
///     n_molecules_ref: 99,
///     n_molecules_alt: 1,
///     n_duplex_alt: 0,
///     n_simplex_alt: 1,
///     n_simplex_fwd_alt: 1,
///     n_simplex_rev_alt: 0,
///     n_duplex_ref: 0,
///     n_simplex_ref: 99,
///     mean_alt_error_prob: 0.001,
///     min_variant_specific_duplex: 0,
///     mean_variant_specific_molecules: 1.0,
///     confidence: 0.99,
///     strand_bias_p: 0.5,
///     filter: VariantFilter::Pass,
/// };
///
/// // Truth VCF has a DUP at pos 500; call is at pos 493 (partial allele).
/// // With tolerance=10, the call is within 10bp of pos 500 → PASS.
/// let mut targets = TargetVariantSet::new();
/// targets.insert(("chr1".to_string(), 500, "C".to_string(), "CDUP100".to_string()));
/// let mut calls = vec![call];
/// apply_target_filter_with_tolerance(&mut calls, &targets, 10);
/// assert_eq!(calls[0].filter, VariantFilter::Pass);
/// ```
pub fn apply_target_filter_with_tolerance(
    calls: &mut [VariantCall],
    targets: &TargetVariantSet,
    tolerance: i64,
) {
    let has_bnd_targets = targets_has_bnd(targets);
    // Build a (chrom, pos) position set from the target VCF for fast proximity checks.
    let positions: TargetPositionSet = targets
        .iter()
        .map(|(chrom, pos, _, _)| (chrom.clone(), *pos))
        .collect();

    for call in calls.iter_mut() {
        if call.filter != VariantFilter::Pass {
            continue;
        }
        // Fusion calls: pass through when BND targets exist.
        if call.variant_type == VariantType::Fusion && has_bnd_targets {
            continue;
        }

        let key = extract_variant_key(&call.target_id, &call.ref_sequence, &call.alt_sequence);
        if key.is_none() {
            call.filter = VariantFilter::NotTargeted;
            continue;
        }
        let (call_chrom, call_pos, ref_allele, alt_allele) = key.unwrap();

        // 1. Exact match.
        if targets.contains(&(call_chrom.clone(), call_pos, ref_allele, alt_allele)) {
            continue; // already PASS
        }

        // 2. Position-based fallback (only when tolerance > 0).
        if tolerance > 0 {
            let near_target = positions.iter().any(|(t_chrom, t_pos)| {
                t_chrom == &call_chrom && (call_pos - t_pos).abs() <= tolerance
            });
            if near_target {
                continue; // PASS via position proximity
            }
        }

        call.filter = VariantFilter::NotTargeted;
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Parse "chrN:start-end" into (chrom, start).
pub(crate) fn parse_target_id(target_id: &str) -> Option<(String, i64)> {
    let (chrom, rest) = target_id.split_once(':')?;
    let (start_str, _end_str) = rest.split_once('-')?;
    let start: i64 = start_str.parse().ok()?;
    Some((chrom.to_string(), start))
}

/// Convert a byte slice to a UTF-8 string (lossily).
fn bytes_to_str(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).into_owned()
}

/// Compute the variant key for an indel call, with left-normalisation.
///
/// Mirrors the indel branch of the Python `extract_called_variants` function.
#[allow(clippy::too_many_arguments)]
fn indel_key(
    chrom: String,
    target_start: i64,
    diff_pos: usize,
    ref_seq: &[u8],
    alt_seq: &[u8],
    ref_trimmed: &[u8],
    alt_trimmed: &[u8],
) -> Option<(String, i64, String, String)> {
    // Remove inner common suffix between ref_trimmed and alt_trimmed.
    let inner_cs = ref_trimmed
        .iter()
        .rev()
        .zip(alt_trimmed.iter().rev())
        .take_while(|(r, a)| r == a)
        .count();

    let ref_min = if inner_cs > 0 {
        &ref_trimmed[..ref_trimmed.len() - inner_cs]
    } else {
        ref_trimmed
    };
    let alt_min = if inner_cs > 0 {
        &alt_trimmed[..alt_trimmed.len() - inner_cs]
    } else {
        alt_trimmed
    };

    // Remove inner common prefix.
    let inner_cp = ref_min
        .iter()
        .zip(alt_min.iter())
        .take_while(|(r, a)| r == a)
        .count();
    let ref_min = &ref_min[inner_cp..];
    let alt_min = &alt_min[inner_cp..];
    let indel_start = diff_pos + inner_cp;

    if ref_seq.len() > alt_seq.len() {
        // Deletion: ref_min is the deleted sequence.
        deletion_key(chrom, target_start, indel_start, ref_seq, ref_min)
    } else {
        // Insertion: alt_min is the inserted sequence.
        insertion_key(chrom, target_start, indel_start, ref_seq, alt_min)
    }
}

fn deletion_key(
    chrom: String,
    target_start: i64,
    indel_start: usize,
    ref_seq: &[u8],
    del_seq_orig: &[u8],
) -> Option<(String, i64, String, String)> {
    let mut del_seq = del_seq_orig.to_vec();
    let mut anchor_pos = indel_start as i64 - 1;

    // Left-normalise: shift anchor left while last base of del_seq matches anchor.
    while anchor_pos > 0
        && !del_seq.is_empty()
        && ref_seq[anchor_pos as usize]
            == *del_seq.last().expect("del_seq non-empty in while guard")
    {
        // Rotate right: move last element to front.
        let last = *del_seq.last().expect("del_seq non-empty in while guard");
        del_seq.rotate_right(1);
        del_seq[0] = last;
        anchor_pos -= 1;
    }

    if anchor_pos >= 0 {
        let anchor = ref_seq[anchor_pos as usize];
        let mut ref_allele = vec![anchor];
        ref_allele.extend_from_slice(&del_seq);
        let alt_allele = vec![anchor];
        // Convert 0-based anchor offset to 1-based VCF position.
        let genomic_pos = target_start + anchor_pos + 1;
        Some((
            chrom,
            genomic_pos,
            bytes_to_str(&ref_allele),
            bytes_to_str(&alt_allele),
        ))
    } else {
        let genomic_pos = target_start + indel_start as i64 + 1;
        Some((chrom, genomic_pos, bytes_to_str(&del_seq), String::new()))
    }
}

fn insertion_key(
    chrom: String,
    target_start: i64,
    indel_start: usize,
    ref_seq: &[u8],
    ins_seq_orig: &[u8],
) -> Option<(String, i64, String, String)> {
    let mut ins_seq = ins_seq_orig.to_vec();
    let mut anchor_pos = indel_start as i64 - 1;

    // Left-normalise: shift anchor left while last base of ins_seq matches anchor.
    while anchor_pos > 0
        && !ins_seq.is_empty()
        && ref_seq[anchor_pos as usize]
            == *ins_seq.last().expect("ins_seq non-empty in while guard")
    {
        let last = *ins_seq.last().expect("ins_seq non-empty in while guard");
        ins_seq.rotate_right(1);
        ins_seq[0] = last;
        anchor_pos -= 1;
    }

    if anchor_pos >= 0 {
        let anchor = ref_seq[anchor_pos as usize];
        let ref_allele = vec![anchor];
        let mut alt_allele = vec![anchor];
        alt_allele.extend_from_slice(&ins_seq);
        // Convert 0-based anchor offset to 1-based VCF position.
        let genomic_pos = target_start + anchor_pos + 1;
        Some((
            chrom,
            genomic_pos,
            bytes_to_str(&ref_allele),
            bytes_to_str(&alt_allele),
        ))
    } else {
        let genomic_pos = target_start + indel_start as i64 + 1;
        Some((chrom, genomic_pos, String::new(), bytes_to_str(&ins_seq)))
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Test 1: SNV at 0-based window offset 2, target starts at 100 (0-based).
    // 1-based VCF pos = 100 + 2 + 1 = 103.
    #[test]
    fn extract_snv_key() {
        let ref_seq = b"ACGTACGT";
        let alt_seq = b"ACTTACGT";
        let key = extract_variant_key("chr1:100-108", ref_seq, alt_seq);
        assert_eq!(
            key,
            Some(("chr1".to_string(), 103, "G".to_string(), "T".to_string()))
        );
    }

    // Test 2: identical sequences → no variant.
    #[test]
    fn extract_no_variant_returns_none() {
        let seq = b"ACGTACGT";
        assert_eq!(extract_variant_key("chr1:100-108", seq, seq), None);
    }

    // Test 3: parse_target_id round-trip.
    #[test]
    fn parse_target_id_ok() {
        let r = parse_target_id("chr17:7572800-7572900");
        assert_eq!(r, Some(("chr17".to_string(), 7_572_800)));
    }

    // Test 4: apply_target_filter marks un-targeted calls.
    #[test]
    fn apply_target_filter_marks_not_targeted() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        let call = VariantCall {
            target_id: "chr1:100-108".to_string(),
            variant_type: VariantType::Snv,
            ref_sequence: b"ACGTACGT".to_vec(),
            alt_sequence: b"ACTTACGT".to_vec(),
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules_ref: 99,
            n_molecules_alt: 1,
            n_duplex_alt: 0,
            n_simplex_alt: 1,
            n_simplex_fwd_alt: 1,
            n_simplex_rev_alt: 0,
            n_duplex_ref: 0,
            n_simplex_ref: 99,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 1.0,
            confidence: 0.99,
            strand_bias_p: 0.5,
            filter: VariantFilter::Pass,
        };

        let empty_targets = TargetVariantSet::new();
        apply_target_filter(&mut [call.clone()], &empty_targets);
        // Call was Pass, not in empty set → NotTargeted.
        apply_target_filter(&mut [call.clone()], &empty_targets);

        let mut calls = vec![call.clone()];
        apply_target_filter(&mut calls, &empty_targets);
        assert_eq!(calls[0].filter, VariantFilter::NotTargeted);
    }

    // Test 5: apply_target_filter passes targeted calls.
    #[test]
    fn apply_target_filter_passes_match() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        let call = VariantCall {
            target_id: "chr1:100-108".to_string(),
            variant_type: VariantType::Snv,
            ref_sequence: b"ACGTACGT".to_vec(),
            alt_sequence: b"ACTTACGT".to_vec(),
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules_ref: 99,
            n_molecules_alt: 1,
            n_duplex_alt: 0,
            n_simplex_alt: 1,
            n_simplex_fwd_alt: 1,
            n_simplex_rev_alt: 0,
            n_duplex_ref: 0,
            n_simplex_ref: 99,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 1.0,
            confidence: 0.99,
            strand_bias_p: 0.5,
            filter: VariantFilter::Pass,
        };

        let mut targets = TargetVariantSet::new();
        targets.insert(("chr1".to_string(), 103, "G".to_string(), "T".to_string()));

        let mut calls = vec![call];
        apply_target_filter(&mut calls, &targets);
        assert_eq!(calls[0].filter, VariantFilter::Pass);
    }

    // Test 6: non-Pass calls are not touched by apply_target_filter.
    #[test]
    fn apply_target_filter_ignores_non_pass() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        let call = VariantCall {
            target_id: "chr1:100-108".to_string(),
            variant_type: VariantType::Snv,
            ref_sequence: b"ACGTACGT".to_vec(),
            alt_sequence: b"ACTTACGT".to_vec(),
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules_ref: 99,
            n_molecules_alt: 1,
            n_duplex_alt: 0,
            n_simplex_alt: 1,
            n_simplex_fwd_alt: 1,
            n_simplex_rev_alt: 0,
            n_duplex_ref: 0,
            n_simplex_ref: 99,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 1.0,
            confidence: 0.5,
            strand_bias_p: 0.5,
            filter: VariantFilter::LowConfidence,
        };

        let empty = TargetVariantSet::new();
        let mut calls = vec![call];
        apply_target_filter(&mut calls, &empty);
        assert_eq!(calls[0].filter, VariantFilter::LowConfidence);
    }

    // Test 8: position tolerance passes nearby calls regardless of REF/ALT.
    #[test]
    fn apply_target_filter_with_tolerance_passes_nearby_call() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        // Call is at pos 493 (partial 1bp insertion); truth is at pos 500 (full DUP).
        // Target: (chr1, 500, "C", "CDUP").
        // Called key for chr1:490-600, ref=b"ACGT", alt=b"ACCT":
        //   diff at offset 2 → SNV pos = 490 + 2 + 1 = 493.
        let call = VariantCall {
            target_id: "chr1:490-600".to_string(),
            variant_type: VariantType::Snv,
            ref_sequence: b"ACGT".to_vec(),
            alt_sequence: b"ACCT".to_vec(),
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules_ref: 99,
            n_molecules_alt: 1,
            n_duplex_alt: 0,
            n_simplex_alt: 1,
            n_simplex_fwd_alt: 1,
            n_simplex_rev_alt: 0,
            n_duplex_ref: 0,
            n_simplex_ref: 99,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 1.0,
            confidence: 0.99,
            strand_bias_p: 0.5,
            filter: VariantFilter::Pass,
        };

        let mut targets = TargetVariantSet::new();
        // Truth is at pos 500 with a different REF/ALT — exact match will fail.
        targets.insert(("chr1".to_string(), 500, "C".to_string(), "CDUP".to_string()));

        // With tolerance=10: call at pos 493, target at pos 500 → |493-500|=7 ≤ 10 → Pass.
        let mut calls = vec![call.clone()];
        apply_target_filter_with_tolerance(&mut calls, &targets, 10);
        assert_eq!(calls[0].filter, VariantFilter::Pass);

        // With tolerance=5: |493-500|=7 > 5 → NotTargeted.
        let mut calls2 = vec![call];
        apply_target_filter_with_tolerance(&mut calls2, &targets, 5);
        assert_eq!(calls2[0].filter, VariantFilter::NotTargeted);
    }

    // Test 9: exact match still works when tolerance > 0.
    #[test]
    fn apply_target_filter_with_tolerance_exact_match_still_works() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        let call = VariantCall {
            target_id: "chr1:100-108".to_string(),
            variant_type: VariantType::Snv,
            ref_sequence: b"ACGTACGT".to_vec(),
            alt_sequence: b"ACTTACGT".to_vec(),
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules_ref: 99,
            n_molecules_alt: 1,
            n_duplex_alt: 0,
            n_simplex_alt: 1,
            n_simplex_fwd_alt: 1,
            n_simplex_rev_alt: 0,
            n_duplex_ref: 0,
            n_simplex_ref: 99,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 1.0,
            confidence: 0.99,
            strand_bias_p: 0.5,
            filter: VariantFilter::Pass,
        };

        let mut targets = TargetVariantSet::new();
        // Exact match: chr1, pos 103, G→T.
        targets.insert(("chr1".to_string(), 103, "G".to_string(), "T".to_string()));

        let mut calls = vec![call];
        apply_target_filter_with_tolerance(&mut calls, &targets, 0);
        assert_eq!(calls[0].filter, VariantFilter::Pass);
    }

    // Test 10: Fusion calls pass through when BND targets exist.
    #[test]
    fn apply_target_filter_passes_fusion_with_bnd_targets() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        let call = VariantCall {
            target_id: "BCR_ABL1__chr1:150-200__chr1:900-950__fusion".to_string(),
            variant_type: VariantType::Fusion,
            ref_sequence: b"ACGT".to_vec(),
            alt_sequence: b"ACGT".to_vec(),
            vaf: 0.40,
            vaf_ci_low: 0.38,
            vaf_ci_high: 0.42,
            n_molecules_ref: 600,
            n_molecules_alt: 400,
            n_duplex_alt: 0,
            n_simplex_alt: 400,
            n_simplex_fwd_alt: 200,
            n_simplex_rev_alt: 200,
            n_duplex_ref: 0,
            n_simplex_ref: 600,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 400.0,
            confidence: 0.99,
            strand_bias_p: 0.5,
            filter: VariantFilter::Pass,
        };

        // Target VCF has BND record (ALT contains ']').
        let mut targets = TargetVariantSet::new();
        targets.insert((
            "chr1".to_string(),
            200,
            "A".to_string(),
            "A]chr1:900]".to_string(),
        ));

        let mut calls = vec![call];
        apply_target_filter(&mut calls, &targets);
        // Fusion call should remain Pass when BND targets exist.
        assert_eq!(calls[0].filter, VariantFilter::Pass);
    }

    // Test 11: Fusion calls still get NotTargeted when no BND targets exist.
    #[test]
    fn apply_target_filter_not_targeted_fusion_without_bnd_targets() {
        use crate::caller::{VariantCall, VariantFilter, VariantType};

        let call = VariantCall {
            target_id: "BCR_ABL1__chr1:150-200__chr1:900-950__fusion".to_string(),
            variant_type: VariantType::Fusion,
            ref_sequence: b"ACGT".to_vec(),
            alt_sequence: b"ACGT".to_vec(),
            vaf: 0.40,
            vaf_ci_low: 0.38,
            vaf_ci_high: 0.42,
            n_molecules_ref: 600,
            n_molecules_alt: 400,
            n_duplex_alt: 0,
            n_simplex_alt: 400,
            n_simplex_fwd_alt: 200,
            n_simplex_rev_alt: 200,
            n_duplex_ref: 0,
            n_simplex_ref: 600,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 0,
            mean_variant_specific_molecules: 400.0,
            confidence: 0.99,
            strand_bias_p: 0.5,
            filter: VariantFilter::Pass,
        };

        // Target VCF has only SNV records (no BND).
        let mut targets = TargetVariantSet::new();
        targets.insert(("chr1".to_string(), 103, "G".to_string(), "T".to_string()));

        let mut calls = vec![call];
        apply_target_filter(&mut calls, &targets);
        // No BND targets → fusion call should be NotTargeted.
        assert_eq!(calls[0].filter, VariantFilter::NotTargeted);
    }

    // Test 7: load_target_variants parses a minimal VCF string.
    #[test]
    fn load_target_variants_parses_vcf() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let vcf_path = dir.path().join("truth.vcf");
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.1").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT").unwrap();
        writeln!(f, "chr17\t7572812\t.\tG\tT").unwrap();
        writeln!(f, "chr7\t55242468\t.\tAGTT\tA").unwrap();

        let set = load_target_variants(&vcf_path).unwrap();
        assert!(set.contains(&(
            "chr17".to_string(),
            7_572_812,
            "G".to_string(),
            "T".to_string()
        )));
        assert!(set.contains(&(
            "chr7".to_string(),
            55_242_468,
            "AGTT".to_string(),
            "A".to_string()
        )));
        assert_eq!(set.len(), 2);
    }
}
