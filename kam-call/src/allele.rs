//! Minimal VCF allele extraction and left-normalisation.
//!
//! Converts full-sequence path representations (100+ bp ref and alt) to the
//! minimal VCF allele form required for standard variant comparison. For SNVs
//! this trims the flanking shared context. For indels it additionally applies
//! VCF left-normalisation so that deletions or insertions in repetitive regions
//! are always anchored as far left as possible — matching the convention used by
//! bcftools, GATK, and standard truth sets.
//!
//! # Coordinate convention
//!
//! Target IDs in the form `chrN:START-END` use a 1-based closed interval where
//! START is the genomic position of the first base in the stored sequence.
//! A base at sequence index `i` (0-based) is therefore at VCF position
//! `START + i` (1-based).

use std::fmt;

// ─── Public types ─────────────────────────────────────────────────────────────

/// A minimal, left-normalised variant allele ready for VCF output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MinimalAllele {
    /// Chromosome or sequence name.
    pub chrom: String,
    /// 1-based genomic position of the first base of the REF allele.
    pub pos: u64,
    /// Reference allele string.
    pub ref_allele: String,
    /// Alternate allele string.
    pub alt_allele: String,
}

impl fmt::Display for MinimalAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{} {}->{}",
            self.chrom, self.pos, self.ref_allele, self.alt_allele
        )
    }
}

// ─── Public functions ─────────────────────────────────────────────────────────

/// Extract a minimal left-normalised VCF allele from full target sequences.
///
/// `target_id` must match the format `chrN:START-END` where START is the
/// 1-based genomic position of the first base in `ref_seq`. Returns `None`
/// if `target_id` cannot be parsed or if `ref_seq` and `alt_seq` are identical.
///
/// # Algorithm
///
/// 1. Find the first position where `ref_seq` and `alt_seq` differ.
/// 2. Trim the common suffix of both full sequences to isolate the minimal
///    differing region.
/// 3. For indels: strip any residual inner common prefix and suffix of the
///    trimmed region to isolate the inserted or deleted bases.
/// 4. Left-normalise: shift the anchor left while the preceding base in the
///    reference matches the last base of the inserted/deleted sequence.
/// 5. Prepend the anchor base as required by the VCF specification.
///
/// # Example
///
/// ```
/// use kam_call::allele::extract_minimal_allele;
///
/// // Deletion of one G from a "GG" run.
/// // ref: ...T(50) G(51) G(52) T(53)...
/// // alt: ...T(50) G(51) T(52)...     (G at index 51 deleted)
/// // First diff is at index 52; left-normalisation shifts to anchor at index 50.
/// let ref_seq: Vec<u8> = (0u8..101).map(|i| if i == 51 || i == 52 { b'G' }
///     else if i == 50 || i == 53 { b'T' } else { b'A' }).collect();
/// let alt_seq: Vec<u8> = {
///     let mut v = ref_seq[..51].to_vec();
///     v.extend_from_slice(&ref_seq[52..]);
///     v
/// };
/// let a = extract_minimal_allele("chr2:1000-1100", &ref_seq, &alt_seq).unwrap();
/// assert_eq!(a.pos, 1050);
/// assert_eq!(a.ref_allele, "TG");
/// assert_eq!(a.alt_allele, "T");
/// ```
pub fn extract_minimal_allele(
    target_id: &str,
    ref_seq: &[u8],
    alt_seq: &[u8],
) -> Option<MinimalAllele> {
    let (chrom, target_start) = parse_target_id(target_id)?;

    // Step 1: find the first differing position.
    let diff_pos = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .position(|(r, a)| r != a)?;

    // Step 2: trim the outer common suffix (must not consume diff_pos).
    let mut outer_cs = 0usize;
    while outer_cs < ref_seq.len().saturating_sub(diff_pos + 1)
        && outer_cs < alt_seq.len().saturating_sub(diff_pos + 1)
        && ref_seq[ref_seq.len() - 1 - outer_cs] == alt_seq[alt_seq.len() - 1 - outer_cs]
    {
        outer_cs += 1;
    }
    let ref_trimmed = &ref_seq[diff_pos..ref_seq.len() - outer_cs];
    let alt_trimmed = &alt_seq[diff_pos..alt_seq.len() - outer_cs];

    if ref_seq.len() == alt_seq.len() {
        // SNV / MNV: no indel normalisation needed.
        Some(MinimalAllele {
            chrom,
            pos: target_start + diff_pos as u64,
            ref_allele: bytes_to_string(ref_trimmed),
            alt_allele: bytes_to_string(alt_trimmed),
        })
    } else {
        let is_deletion = ref_seq.len() > alt_seq.len();
        indel_allele(
            chrom,
            target_start,
            ref_seq,
            ref_trimmed,
            alt_trimmed,
            diff_pos,
            is_deletion,
        )
    }
}

// ─── Private helpers ──────────────────────────────────────────────────────────

/// Parse `chrN:START-END` → `(chrom, start)`.
///
/// `start` is the integer before the dash, treated as the 1-based position of
/// the first base in the target sequence (see module-level coordinate note).
fn parse_target_id(id: &str) -> Option<(String, u64)> {
    let colon = id.find(':')?;
    let rest = &id[colon + 1..];
    let dash = rest.find('-')?;
    let chrom = id[..colon].to_string();
    let start: u64 = rest[..dash].parse().ok()?;
    Some((chrom, start))
}

/// Compute the left-normalised indel allele from trimmed region slices.
fn indel_allele(
    chrom: String,
    target_start: u64,
    ref_seq: &[u8],
    ref_trimmed: &[u8],
    alt_trimmed: &[u8],
    diff_pos: usize,
    is_deletion: bool,
) -> Option<MinimalAllele> {
    // Step 3: strip residual inner common suffix (allows alt_min to be empty).
    let mut inner_cs = 0usize;
    while inner_cs < ref_trimmed.len()
        && inner_cs < alt_trimmed.len()
        && ref_trimmed[ref_trimmed.len() - 1 - inner_cs]
            == alt_trimmed[alt_trimmed.len() - 1 - inner_cs]
    {
        inner_cs += 1;
    }
    let ref_min = &ref_trimmed[..ref_trimmed.len() - inner_cs];
    let alt_min = &alt_trimmed[..alt_trimmed.len() - inner_cs];

    // Strip inner common prefix.
    let inner_cp = ref_min
        .iter()
        .zip(alt_min.iter())
        .take_while(|(r, a)| r == a)
        .count();
    let event_start = diff_pos + inner_cp;

    if is_deletion {
        // Deleted bases: portion of ref_min beyond the shared prefix.
        let mut del_seq: Vec<u8> = ref_min[inner_cp..].to_vec();
        normalise_and_anchor(
            chrom,
            target_start,
            ref_seq,
            &mut del_seq,
            event_start,
            true,
        )
    } else {
        // Inserted bases: portion of alt_min beyond the shared prefix.
        let mut ins_seq: Vec<u8> = alt_min[inner_cp..].to_vec();
        normalise_and_anchor(
            chrom,
            target_start,
            ref_seq,
            &mut ins_seq,
            event_start,
            false,
        )
    }
}

/// Shift `event_seq` left and build the anchor-prefixed alleles.
///
/// For a deletion: REF = anchor + del_seq, ALT = anchor.
/// For an insertion: REF = anchor, ALT = anchor + ins_seq.
///
/// Left-normalisation: while `anchor_pos > 0` and the preceding reference base
/// equals the last base of `event_seq`, rotate `event_seq` right by one and
/// decrement `anchor_pos`.
fn normalise_and_anchor(
    chrom: String,
    target_start: u64,
    ref_seq: &[u8],
    event_seq: &mut Vec<u8>,
    event_start: usize,
    is_deletion: bool,
) -> Option<MinimalAllele> {
    // The anchor position is the base immediately before the event.
    let mut anchor_pos = event_start.saturating_sub(1);

    // Step 4: left-normalise.
    while anchor_pos > 0 && !event_seq.is_empty() && ref_seq[anchor_pos] == *event_seq.last()? {
        let last = event_seq.remove(event_seq.len() - 1);
        event_seq.insert(0, last);
        anchor_pos -= 1;
    }

    // Step 5: prepend the anchor base.
    let anchor = ref_seq[anchor_pos] as char;
    let (ref_allele, alt_allele) = if is_deletion {
        (
            format!("{anchor}{}", bytes_to_string(event_seq)),
            anchor.to_string(),
        )
    } else {
        (
            anchor.to_string(),
            format!("{anchor}{}", bytes_to_string(event_seq)),
        )
    };

    Some(MinimalAllele {
        chrom,
        pos: target_start + anchor_pos as u64,
        ref_allele,
        alt_allele,
    })
}

fn bytes_to_string(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).into_owned()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn check(target_id: &str, ref_seq: &[u8], alt_seq: &[u8], pos: u64, r: &str, a: &str) {
        let result = extract_minimal_allele(target_id, ref_seq, alt_seq)
            .unwrap_or_else(|| panic!("expected Some for {target_id}"));
        assert_eq!(
            result.pos, pos,
            "pos mismatch: got {}, want {}",
            result.pos, pos
        );
        assert_eq!(
            result.ref_allele, r,
            "REF mismatch: got {}, want {r}",
            result.ref_allele
        );
        assert_eq!(
            result.alt_allele, a,
            "ALT mismatch: got {}, want {a}",
            result.alt_allele
        );
    }

    // Test 1: SNV at index 50 of a 101-bp target.
    #[test]
    fn snv_correct_position() {
        let mut r = vec![b'A'; 101];
        let mut a = r.clone();
        r[50] = b'C';
        a[50] = b'T';
        check("chr1:1000-1100", &r, &a, 1050, "C", "T");
    }

    // Test 2: Identical sequences → None.
    #[test]
    fn identical_returns_none() {
        let seq = vec![b'A'; 101];
        assert!(extract_minimal_allele("chr1:1000-1100", &seq, &seq).is_none());
    }

    // Test 3: Invalid target_id → None.
    #[test]
    fn invalid_id_returns_none() {
        let r = vec![b'A'; 5];
        let mut a = r.clone();
        a[2] = b'T';
        assert!(extract_minimal_allele("notvalid", &r, &a).is_none());
    }

    // Test 4: Deletion in a GG run is left-normalised to the T anchor.
    //
    // ref: ...T(50) G(51) G(52) T(53)...
    // alt: ...T(50) G(51) T(52)...       (G at index 51 deleted)
    // First diff at index 52 (because ref[51]==alt[51] after the shift);
    // left-normalise → anchor at index 50 (T).
    // Expected: pos=target_start+50, REF="TG", ALT="T".
    #[test]
    fn deletion_left_normalised_in_gg_run() {
        let mut r = vec![b'A'; 101];
        r[50] = b'T';
        r[51] = b'G';
        r[52] = b'G';
        r[53] = b'T';
        // Delete r[51] (first G).
        let mut a: Vec<u8> = r[..51].to_vec();
        a.extend_from_slice(&r[52..]);
        check("chr2:1000-1100", &r, &a, 1050, "TG", "T");
    }

    // Test 5: Simple 1-bp insertion with no repeat — no left shift.
    //
    // ref: ...A(50) A(51)... (all A)
    // alt: ...A(50) G A(51)... (G inserted)
    // anchor = A at index 50; A ≠ G, so no shift.
    // Expected: pos=target_start+50, REF="A", ALT="AG".
    #[test]
    fn insertion_no_repeat_anchor_correct() {
        let r: Vec<u8> = vec![b'A'; 100];
        let mut a: Vec<u8> = r[..51].to_vec();
        a.push(b'G');
        a.extend_from_slice(&r[51..]);
        check("chr3:500-599", &r, &a, 550, "A", "AG");
    }

    // Test 6: MNV (two adjacent SNVs) — same length, no normalisation.
    #[test]
    fn mnv_same_length() {
        let mut r = vec![b'A'; 10];
        let mut a = r.clone();
        r[3] = b'C';
        r[4] = b'G';
        a[3] = b'T';
        a[4] = b'T';
        check("chrX:100-109", &r, &a, 103, "CG", "TT");
    }

    // Test 7: Multi-base deletion with distinct flanking bases (no repeat).
    // Delete "GCT" at indices 40-42.
    #[test]
    fn multi_base_deletion_no_repeat() {
        let mut r: Vec<u8> = vec![b'A'; 103];
        r[39] = b'C'; // anchor
        r[40] = b'G';
        r[41] = b'C';
        r[42] = b'T';
        r[43] = b'T'; // base after deletion (distinct to prevent shift)
        let mut a: Vec<u8> = r[..40].to_vec();
        a.extend_from_slice(&r[43..]);
        // anchor=r[39]='C', del_seq="GCT", del_seq[-1]='T'; r[39]='C' ≠ 'T' → no shift
        check("chr5:500-602", &r, &a, 539, "CGCT", "C");
    }

    // Test 8: parse_target_id handles standard coordinates.
    #[test]
    fn parse_target_id_standard() {
        assert_eq!(
            parse_target_id("chr2:29318295-29318395"),
            Some(("chr2".to_string(), 29318295))
        );
    }

    // Test 9: parse_target_id rejects malformed inputs.
    #[test]
    fn parse_target_id_malformed() {
        assert!(parse_target_id("no_colon").is_none());
        assert!(parse_target_id("chr1:abc-100").is_none());
        assert!(parse_target_id("chr1:100").is_none()); // no dash
    }

    // Test 10: Insertion with repeat context is left-normalised.
    //
    // Insert T into a "TT" run.
    // ref: ...A(49) T(50) T(51) A(52)...
    // alt: ...A(49) T(50) T(51) T(52) A(53)... (T inserted after index 51)
    // First diff at index 52 (ref[52]=A, alt[52]=T after insert shifts everything).
    // ins_seq = "T"; ref[52]='A'... actually wait, the first diff is where the
    // sequence shifts, which is at index 52. ins_seq[-1]='T', ref[52-1]=ref[51]='T'
    // → shift: anchor_pos=50. ref[50]='T', ins_seq[-1]='T' → shift: anchor_pos=49.
    // ref[49]='A' ≠ 'T' → stop.
    // Expected: pos=target_start+49, REF="A", ALT="AT".
    #[test]
    fn insertion_left_normalised_in_tt_run() {
        let mut r: Vec<u8> = vec![b'A'; 100];
        r[50] = b'T';
        r[51] = b'T';
        // Insert T after index 51.
        let mut a: Vec<u8> = r[..52].to_vec();
        a.push(b'T');
        a.extend_from_slice(&r[52..]);
        check("chr7:1000-1099", &r, &a, 1049, "A", "AT");
    }
}
