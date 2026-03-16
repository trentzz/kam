//! 2-bit k-mer encoding for DNA sequences.
//!
//! Each base is represented as 2 bits: A=00, C=01, G=10, T=11.
//! Bases are packed left-to-right: the first base occupies the highest bit
//! pair in the u64. This allows direct numeric comparison for canonical form
//! selection.
//!
//! Only k≤31 is supported; k=31 uses 62 of the 64 bits available in a u64.

/// Encode a single base as 2 bits: A=0b00, C=0b01, G=0b10, T=0b11.
///
/// Returns `None` for `N` or any unrecognised byte.
///
/// # Examples
///
/// ```
/// use kam_index::encode::encode_base;
/// assert_eq!(encode_base(b'A'), Some(0b00));
/// assert_eq!(encode_base(b'C'), Some(0b01));
/// assert_eq!(encode_base(b'G'), Some(0b10));
/// assert_eq!(encode_base(b'T'), Some(0b11));
/// assert_eq!(encode_base(b'N'), None);
/// ```
pub fn encode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

/// Decode a 2-bit value back to its uppercase base character.
///
/// Only the two least-significant bits of `encoded` are examined.
///
/// # Examples
///
/// ```
/// use kam_index::encode::decode_base;
/// assert_eq!(decode_base(0b00), b'A');
/// assert_eq!(decode_base(0b01), b'C');
/// assert_eq!(decode_base(0b10), b'G');
/// assert_eq!(decode_base(0b11), b'T');
/// ```
pub fn decode_base(encoded: u8) -> u8 {
    match encoded & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => unreachable!(),
    }
}

/// Encode a k-mer sequence into a `u64`.
///
/// Bases are packed left-to-right: the first base occupies the two
/// most-significant bits used (i.e., bits `2*(k-1)` and `2*(k-1)+1` from the
/// right, at position `(k-1)` in 2-bit units from the bottom).
///
/// Returns `None` if:
/// - the sequence contains `N` or any unrecognised base, or
/// - `k > 31` (would require more than 62 bits).
///
/// # Examples
///
/// ```
/// use kam_index::encode::{encode_kmer, decode_kmer};
/// let enc = encode_kmer(b"ACGT").unwrap();
/// assert_eq!(decode_kmer(enc, 4), b"ACGT");
/// ```
pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    let k = seq.len();
    if k > 31 {
        return None;
    }
    let mut val: u64 = 0;
    for &base in seq {
        let bits = encode_base(base)? as u64;
        val = (val << 2) | bits;
    }
    Some(val)
}

/// Decode a `u64` back into a k-mer sequence of length `k`.
///
/// # Examples
///
/// ```
/// use kam_index::encode::{encode_kmer, decode_kmer};
/// let enc = encode_kmer(b"GATTACA").unwrap();
/// assert_eq!(decode_kmer(enc, 7), b"GATTACA");
/// ```
pub fn decode_kmer(encoded: u64, k: usize) -> Vec<u8> {
    (0..k)
        .map(|i| {
            // Base i is at bit position 2*(k-1-i) from the right.
            let shift = 2 * (k - 1 - i);
            let bits = ((encoded >> shift) & 0b11) as u8;
            decode_base(bits)
        })
        .collect()
}

/// Compute the reverse complement of an encoded k-mer.
///
/// The complement of each base is obtained by flipping both bits (XOR with
/// `0b11`). The order is then reversed so the last base becomes the first.
///
/// # Examples
///
/// ```
/// use kam_index::encode::{encode_kmer, reverse_complement};
/// // ACGT is a palindrome: its reverse complement is itself.
/// let enc = encode_kmer(b"ACGT").unwrap();
/// assert_eq!(reverse_complement(enc, 4), enc);
///
/// // AAAA → TTTT
/// let aaaa = encode_kmer(b"AAAA").unwrap();
/// let tttt = encode_kmer(b"TTTT").unwrap();
/// assert_eq!(reverse_complement(aaaa, 4), tttt);
/// ```
pub fn reverse_complement(encoded: u64, k: usize) -> u64 {
    // Complement all bases: XOR the 2k used bits with all-ones mask.
    let mask = if k == 0 { 0 } else { (1u64 << (2 * k)) - 1 };
    let complemented = (encoded ^ mask) & mask;

    // Reverse the order of the k bases.
    let mut reversed: u64 = 0;
    let mut remaining = complemented;
    for _ in 0..k {
        reversed = (reversed << 2) | (remaining & 0b11);
        remaining >>= 2;
    }
    reversed
}

/// Compute the canonical form of a k-mer: `min(kmer, reverse_complement(kmer))`.
///
/// # Examples
///
/// ```
/// use kam_index::encode::{encode_kmer, canonical};
/// let fwd = encode_kmer(b"AAAC").unwrap();
/// let rev = encode_kmer(b"GTTT").unwrap(); // rev-comp of AAAC
/// assert_eq!(canonical(fwd, 4), fwd.min(rev));
/// ```
pub fn canonical(encoded: u64, k: usize) -> u64 {
    let rc = reverse_complement(encoded, k);
    encoded.min(rc)
}

/// Indicates which strand a k-mer was observed on when canonicalised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerStrand {
    /// The original k-mer was already the canonical (smaller) form.
    Forward,
    /// The reverse complement was smaller; the canonical form is the rev-comp.
    Reverse,
}

/// Canonicalise a k-mer and return both the canonical value and the strand.
///
/// Returns `(canonical_kmer, KmerStrand::Forward)` when the original k-mer is
/// already the smaller value, and `(canonical_kmer, KmerStrand::Reverse)` when
/// the reverse complement is smaller.
///
/// # Examples
///
/// ```
/// use kam_index::encode::{encode_kmer, canonical_with_strand, KmerStrand};
/// // ACGT is self-complementary; the original is returned as Forward.
/// let enc = encode_kmer(b"ACGT").unwrap();
/// let (can, strand) = canonical_with_strand(enc, 4);
/// assert_eq!(can, enc);
/// assert_eq!(strand, KmerStrand::Forward);
/// ```
pub fn canonical_with_strand(encoded: u64, k: usize) -> (u64, KmerStrand) {
    let rc = reverse_complement(encoded, k);
    if encoded <= rc {
        (encoded, KmerStrand::Forward)
    } else {
        (rc, KmerStrand::Reverse)
    }
}

/// Iterator over all k-mers in a DNA sequence.
///
/// Uses a sliding window for O(1) per k-mer after the first: the window is
/// shifted left by 2 bits and the new base's 2-bit encoding is OR-ed in.
/// When an `N` (or unrecognised base) is encountered the window is reset and
/// the iterator skips forward until a full window of valid bases is available
/// again.
///
/// Each item is `(position, encoded_kmer)` where `position` is the 0-based
/// index of the *first* base of the k-mer in the original sequence.
///
/// # Examples
///
/// ```
/// use kam_index::encode::KmerIterator;
/// let seq = b"ACGTA";
/// let kmers: Vec<_> = KmerIterator::new(seq, 3).collect();
/// assert_eq!(kmers.len(), 3);
/// assert_eq!(kmers[0].0, 0); // first k-mer starts at position 0
/// ```
pub struct KmerIterator<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
    /// Current sliding window value.
    window: u64,
    /// Number of valid bases currently loaded in the window (0..=k).
    loaded: usize,
    /// Bitmask for the lowest 2*k bits.
    mask: u64,
}

impl<'a> KmerIterator<'a> {
    /// Create a new `KmerIterator` over `seq` with k-mer length `k`.
    ///
    /// # Panics
    ///
    /// Panics if `k == 0` or `k > 31`.
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!(k > 0 && k <= 31, "k must be in 1..=31, got {k}");
        let mask = (1u64 << (2 * k)) - 1;
        KmerIterator {
            seq,
            k,
            pos: 0,
            window: 0,
            loaded: 0,
            mask,
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.pos >= self.seq.len() {
                return None;
            }

            let base = self.seq[self.pos];
            self.pos += 1;

            match encode_base(base) {
                Some(bits) => {
                    self.window = ((self.window << 2) | bits as u64) & self.mask;
                    self.loaded += 1;

                    if self.loaded >= self.k {
                        // The start position of this k-mer in seq.
                        let kmer_start = self.pos - self.k;
                        return Some((kmer_start, self.window));
                    }
                }
                None => {
                    // Reset on N or unrecognised base.
                    self.window = 0;
                    self.loaded = 0;
                }
            }
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Test 1: encode_base for each canonical base and N.
    #[test]
    fn encode_base_values() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b'N'), None);
    }

    // Test 2: decode_base roundtrips correctly for all 4 values.
    #[test]
    fn decode_base_roundtrip() {
        for bits in 0u8..4 {
            let base = decode_base(bits);
            assert_eq!(encode_base(base), Some(bits));
        }
    }

    // Test 3: encode_kmer + decode_kmer roundtrip for k=5 and k=31.
    #[test]
    fn encode_decode_kmer_roundtrip_k5() {
        let seq = b"ACGTA";
        let enc = encode_kmer(seq).unwrap();
        assert_eq!(decode_kmer(enc, 5), seq.to_vec());
    }

    #[test]
    fn encode_decode_kmer_roundtrip_k31() {
        // 31-base sequence using all four bases.
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        assert_eq!(seq.len(), 31);
        let enc = encode_kmer(seq).unwrap();
        assert_eq!(decode_kmer(enc, 31), seq.to_vec());
    }

    // Test 4: encode_kmer returns None for sequence with N.
    #[test]
    fn encode_kmer_none_on_n() {
        assert_eq!(encode_kmer(b"ACNGT"), None);
        assert_eq!(encode_kmer(b"NAAAA"), None);
        assert_eq!(encode_kmer(b"AAAAN"), None);
    }

    // Test 5: encode_kmer returns None for k > 31.
    #[test]
    fn encode_kmer_none_on_k_gt_31() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bases
        assert_eq!(seq.len(), 32);
        assert_eq!(encode_kmer(seq), None);
    }

    // Test 6: reverse_complement("ACGT") is self-complementary.
    #[test]
    fn reverse_complement_acgt_palindrome() {
        let enc = encode_kmer(b"ACGT").unwrap();
        assert_eq!(reverse_complement(enc, 4), enc);
    }

    // Test 7: reverse_complement("AAAA") = encode("TTTT").
    #[test]
    fn reverse_complement_aaaa_is_tttt() {
        let aaaa = encode_kmer(b"AAAA").unwrap();
        let tttt = encode_kmer(b"TTTT").unwrap();
        assert_eq!(reverse_complement(aaaa, 4), tttt);
    }

    // Test 8: canonical always returns the smaller of kmer and revcomp.
    #[test]
    fn canonical_returns_smaller() {
        let seqs: &[&[u8]] = &[b"ACGT", b"AAAC", b"GTTT", b"CCCC", b"TTTT"];
        for &seq in seqs {
            let enc = encode_kmer(seq).unwrap();
            let k = seq.len();
            let rc = reverse_complement(enc, k);
            let can = canonical(enc, k);
            assert_eq!(can, enc.min(rc), "failed for {:?}", seq);
        }
    }

    // Test 9: canonical_with_strand returns Forward when kmer is already canonical.
    #[test]
    fn canonical_with_strand_forward_when_already_canonical() {
        // AAAC < GTTT (its rev-comp), so AAAC is already canonical.
        let enc = encode_kmer(b"AAAC").unwrap();
        let rc = reverse_complement(enc, 4);
        // Ensure our assumption holds.
        assert!(enc < rc, "test assumption: AAAC should be < its rev-comp");
        let (can, strand) = canonical_with_strand(enc, 4);
        assert_eq!(can, enc);
        assert_eq!(strand, KmerStrand::Forward);
    }

    // Test 10: canonical_with_strand returns Reverse when revcomp is smaller.
    #[test]
    fn canonical_with_strand_reverse_when_revcomp_smaller() {
        // TTTT > AAAA (its rev-comp), so the rev-comp is canonical.
        let enc = encode_kmer(b"TTTT").unwrap();
        let rc = reverse_complement(enc, 4);
        assert!(enc > rc, "test assumption: TTTT should be > its rev-comp");
        let (can, strand) = canonical_with_strand(enc, 4);
        assert_eq!(can, rc);
        assert_eq!(strand, KmerStrand::Reverse);
    }

    // Test 11: KmerIterator yields L-k+1 k-mers for a clean sequence of length L.
    #[test]
    fn kmer_iterator_correct_count() {
        let seq = b"ACGTACGTACGT"; // length 12
        let k = 5;
        let kmers: Vec<_> = KmerIterator::new(seq, k).collect();
        assert_eq!(kmers.len(), seq.len() - k + 1);
    }

    // Test 12: KmerIterator skips positions with N, resetting the window.
    #[test]
    fn kmer_iterator_skips_n() {
        // "ACGTNACGT" — the N splits this into two windows of 4 bases each,
        // but we need k=3 to get k-mers from both halves.
        let seq = b"ACGTNACGT";
        let k = 3;
        let kmers: Vec<_> = KmerIterator::new(seq, k).collect();

        // Before N: ACG, CGT (positions 0, 1)
        // After N: ACG, CGT (positions 5, 6)
        assert_eq!(kmers.len(), 4);
        assert_eq!(kmers[0].0, 0);
        assert_eq!(kmers[1].0, 1);
        assert_eq!(kmers[2].0, 5);
        assert_eq!(kmers[3].0, 6);

        // Encoded values should match naively-encoded k-mers.
        assert_eq!(kmers[0].1, encode_kmer(b"ACG").unwrap());
        assert_eq!(kmers[1].1, encode_kmer(b"CGT").unwrap());
        assert_eq!(kmers[2].1, encode_kmer(b"ACG").unwrap());
        assert_eq!(kmers[3].1, encode_kmer(b"CGT").unwrap());
    }

    // Test 13: KmerIterator sliding window gives same results as naive re-encoding.
    #[test]
    fn kmer_iterator_matches_naive_encoding() {
        let seq = b"GATTACAGATTACA";
        let k = 5;
        let iter_results: Vec<_> = KmerIterator::new(seq, k).collect();

        // Naive: encode each window from scratch.
        let naive: Vec<(usize, u64)> = (0..=seq.len() - k)
            .map(|i| (i, encode_kmer(&seq[i..i + k]).unwrap()))
            .collect();

        assert_eq!(iter_results, naive);
    }
}
