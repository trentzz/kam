//! Endpoint fingerprinting for UMI collision detection.
//!
//! Computes a 64-bit fingerprint from the template endpoints of a read pair.
//! Used as a secondary grouping signal to detect UMI collisions — two different
//! molecules that happen to share the same 5bp random UMI.
//!
//! # Encoding
//!
//! Each base is 2-bit encoded: A=00, C=01, G=10, T=11. N is treated as A (00).
//! The fingerprint packs 8 bases from each of 4 endpoints:
//!
//! ```text
//! bits [63:48]  R1 first 8bp
//! bits [47:32]  R1 last  8bp
//! bits [31:16]  R2 first 8bp
//! bits [15:0]   R2 last  8bp
//! ```
//!
//! 4 endpoints × 8 bases × 2 bits = 64 bits, exactly fitting a `u64`.

/// Number of bases used from each endpoint.
pub const FINGERPRINT_BASES: usize = 8;

/// Maximum bit differences (popcount of XOR) allowed for two fingerprints to be
/// considered compatible (same molecule). Tolerates ~2 base sequencing errors
/// across all 4 endpoints at Q20 quality.
pub const FINGERPRINT_MAX_DIFF: u32 = 4;

/// Encode a single base to its 2-bit representation.
///
/// - A → 0b00
/// - C → 0b01
/// - G → 0b10
/// - T → 0b11
/// - N (or anything else) → 0b00 (treated as A)
#[inline]
fn encode_base(base: u8) -> u64 {
    match base.to_ascii_uppercase() {
        b'C' => 0b01,
        b'G' => 0b10,
        b'T' => 0b11,
        _ => 0b00, // A, N, or any unrecognised base → 00
    }
}

/// Encode up to `FINGERPRINT_BASES` bases from `seq` into a 16-bit field packed
/// into the low bits of a `u64`. If `seq` has fewer than `FINGERPRINT_BASES`
/// bases, the remaining high bits are left as zero (padding).
///
/// Bases are packed MSB-first so that the first base occupies the highest bits
/// of the field, preserving lexicographic ordering.
#[inline]
fn encode_endpoint(seq: &[u8]) -> u64 {
    let len = seq.len().min(FINGERPRINT_BASES);
    let mut field: u64 = 0;
    for &base in seq.iter().take(len) {
        field = (field << 2) | encode_base(base);
    }
    // Shift left to fill the full 16-bit field regardless of how many bases were
    // encoded, so that short sequences are zero-padded on the right.
    let remaining = FINGERPRINT_BASES - len;
    field <<= remaining * 2;
    field
}

/// Compute a 64-bit fingerprint from template read-pair endpoints.
///
/// Uses the first and last [`FINGERPRINT_BASES`] (8) bases of each template
/// read, 2-bit encoded and packed into a single `u64`.
///
/// Two read pairs from the same molecule will have nearly identical fingerprints
/// (within sequencing-error tolerance). Two read pairs from different molecules
/// at the same locus that share a UMI collision will have different fingerprints.
///
/// # Short templates
///
/// If a template is shorter than `2 * FINGERPRINT_BASES` (16 bp), the first and
/// last windows may overlap or be the same slice. No panic occurs — the full
/// available sequence is used for both endpoints, and any unused bits are zero.
///
/// # N bases
///
/// N bases (and any unrecognised character) are encoded as A (00).
///
/// # Examples
///
/// ```
/// use kam_assemble::fingerprint::compute_endpoint_fingerprint;
///
/// let r1 = b"ACGTACGTTTTTTTTTACGT";
/// let r2 = b"TTTTTTTTACGTACGTTTTT";
/// let fp = compute_endpoint_fingerprint(r1, r2);
/// assert_eq!(fp, compute_endpoint_fingerprint(r1, r2)); // deterministic
/// ```
pub fn compute_endpoint_fingerprint(template_r1: &[u8], template_r2: &[u8]) -> u64 {
    // R1 endpoints
    let r1_first = encode_endpoint(template_r1.get(..FINGERPRINT_BASES).unwrap_or(template_r1));
    let r1_last_start = template_r1.len().saturating_sub(FINGERPRINT_BASES);
    let r1_last = encode_endpoint(&template_r1[r1_last_start..]);

    // R2 endpoints
    let r2_first = encode_endpoint(template_r2.get(..FINGERPRINT_BASES).unwrap_or(template_r2));
    let r2_last_start = template_r2.len().saturating_sub(FINGERPRINT_BASES);
    let r2_last = encode_endpoint(&template_r2[r2_last_start..]);

    // Pack into u64: [r1_first 16b][r1_last 16b][r2_first 16b][r2_last 16b]
    (r1_first << 48) | (r1_last << 32) | (r2_first << 16) | r2_last
}

/// Check whether two fingerprints are compatible — i.e., could plausibly come
/// from the same molecule.
///
/// Compatibility is determined by the popcount (number of set bits) of the XOR
/// of the two fingerprints. A threshold of [`FINGERPRINT_MAX_DIFF`] (4) bit
/// differences tolerates approximately 2 base sequencing errors distributed
/// across all 4 endpoints.
///
/// # Examples
///
/// ```
/// use kam_assemble::fingerprint::{compute_endpoint_fingerprint, fingerprints_compatible};
///
/// let r1 = b"ACGTACGTTTTTTTTTACGT";
/// let r2 = b"TTTTTTTTACGTACGTTTTT";
/// let fp = compute_endpoint_fingerprint(r1, r2);
/// assert!(fingerprints_compatible(fp, fp));
/// ```
pub fn fingerprints_compatible(fp1: u64, fp2: u64) -> bool {
    (fp1 ^ fp2).count_ones() <= FINGERPRINT_MAX_DIFF
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------------------------------------------------------------------------
    // Test 1: Identical templates produce identical fingerprints
    // ---------------------------------------------------------------------------
    #[test]
    fn test_identical_templates_identical_fingerprints() {
        let r1 = b"ACGTACGTNNNNNNNNTTTTACGT";
        let r2 = b"TTTTTTTTACGTACGTCCCCCCCC";
        let fp1 = compute_endpoint_fingerprint(r1, r2);
        let fp2 = compute_endpoint_fingerprint(r1, r2);
        assert_eq!(fp1, fp2);
    }

    // ---------------------------------------------------------------------------
    // Test 2: Templates differing by 1 base produce compatible fingerprints
    // ---------------------------------------------------------------------------
    #[test]
    fn test_single_base_difference_is_compatible() {
        let r1a: &[u8] = b"ACGTACGTTTTTTTTTACGT";
        let r2a: &[u8] = b"TTTTTTTTACGTACGTTTTT";

        // Introduce one base change in the first position of r1 (A→C)
        let r1b: &[u8] = b"CCGTACGTTTTTTTTTACGT";
        let r2b: &[u8] = b"TTTTTTTTACGTACGTTTTT";

        let fp_a = compute_endpoint_fingerprint(r1a, r2a);
        let fp_b = compute_endpoint_fingerprint(r1b, r2b);

        // A→C is 00→01: 1 bit difference — well within the threshold of 4
        assert!(fingerprints_compatible(fp_a, fp_b));
    }

    // ---------------------------------------------------------------------------
    // Test 3: Completely different templates produce incompatible fingerprints
    // ---------------------------------------------------------------------------
    #[test]
    fn test_completely_different_templates_incompatible() {
        let r1a: &[u8] = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let r2a: &[u8] = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        let r1b: &[u8] = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let r2b: &[u8] = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

        let fp_a = compute_endpoint_fingerprint(r1a, r2a);
        let fp_b = compute_endpoint_fingerprint(r1b, r2b);

        // A=00, T=11 — every 2-bit pair differs, so XOR is all 1s (64 bits set)
        assert!(!fingerprints_compatible(fp_a, fp_b));
    }

    // ---------------------------------------------------------------------------
    // Test 4: Short templates (< 16 bp) don't panic
    // ---------------------------------------------------------------------------
    #[test]
    fn test_short_templates_no_panic() {
        // Template shorter than 16bp
        let short_r1: &[u8] = b"ACGT";
        let short_r2: &[u8] = b"TTTT";
        let _fp = compute_endpoint_fingerprint(short_r1, short_r2);

        // Empty template
        let _fp2 = compute_endpoint_fingerprint(b"", b"");

        // Single-base template
        let _fp3 = compute_endpoint_fingerprint(b"A", b"C");

        // Exactly 8 bases (boundary — first == last endpoint)
        let _fp4 = compute_endpoint_fingerprint(b"ACGTACGT", b"TTTTTTTT");
    }

    // ---------------------------------------------------------------------------
    // Test 5: N bases are handled (treated as A)
    // ---------------------------------------------------------------------------
    #[test]
    fn test_n_bases_treated_as_a() {
        // A sequence of all A should produce the same fingerprint as all N
        let r1_a: &[u8] = b"AAAAAAAAAAAAAAAAAAAAAAAA";
        let r2_a: &[u8] = b"AAAAAAAAAAAAAAAAAAAAAAAA";

        let r1_n: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNN";
        let r2_n: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNN";

        let fp_a = compute_endpoint_fingerprint(r1_a, r2_a);
        let fp_n = compute_endpoint_fingerprint(r1_n, r2_n);

        assert_eq!(fp_a, fp_n, "N should be encoded identically to A");
    }

    // ---------------------------------------------------------------------------
    // Test 6: Fingerprint is deterministic (same input always same output)
    // ---------------------------------------------------------------------------
    #[test]
    fn test_fingerprint_deterministic() {
        let r1: &[u8] = b"GCTAGCTAGCTAGCTAGCTA";
        let r2: &[u8] = b"TAGCTAGCTAGCTAGCTAGC";

        let results: Vec<u64> = (0..10)
            .map(|_| compute_endpoint_fingerprint(r1, r2))
            .collect();

        assert!(
            results.windows(2).all(|w| w[0] == w[1]),
            "fingerprint must be deterministic across all calls"
        );
    }
}
