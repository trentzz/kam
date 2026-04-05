# Task 010: 2-Bit K-mer Encoding

## Location
`kam-index/src/encode.rs`

## What to implement

2-bit encoding for DNA k-mers: pack bases into u64 for k≤31, compute canonical form (min of kmer and reverse complement), and extract k-mers from sequences.

## Interface

```rust
/// Encode a single base as 2 bits: A=0b00, C=0b01, G=0b10, T=0b11
/// Returns None for N or invalid bases.
pub fn encode_base(base: u8) -> Option<u8>;

/// Decode a 2-bit value back to a base
pub fn decode_base(encoded: u8) -> u8;

/// Encode a k-mer sequence into a u64. Returns None if sequence contains N or is longer than 31.
pub fn encode_kmer(seq: &[u8]) -> Option<u64>;

/// Decode a u64 back into a k-mer sequence of length k
pub fn decode_kmer(encoded: u64, k: usize) -> Vec<u8>;

/// Compute the reverse complement of an encoded k-mer
pub fn reverse_complement(encoded: u64, k: usize) -> u64;

/// Compute canonical form: min(kmer, reverse_complement(kmer))
pub fn canonical(encoded: u64, k: usize) -> u64;

/// Track which strand a k-mer was observed on
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerStrand {
    Forward,   // kmer == canonical (the original was already the smaller)
    Reverse,   // kmer != canonical (the reverse complement was smaller)
}

/// Canonicalize and return strand information
pub fn canonical_with_strand(encoded: u64, k: usize) -> (u64, KmerStrand);

/// Iterator over k-mers in a sequence, skipping positions with N
pub struct KmerIterator<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
}

impl<'a> KmerIterator<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Self;
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = (usize, u64);  // (position, encoded kmer)
}
```

## Implementation notes

- 2-bit encoding: A=00, C=01, G=10, T=11
- Pack left-to-right: first base in highest bits
- Reverse complement: reverse the order AND complement each base (A↔T, C↔G → flip both bits)
- For k=31, uses 62 bits of the u64 — mask off upper 2 bits
- KmerIterator should use a sliding window (shift + add new base) for efficiency, not re-encode each k-mer from scratch. Reset window when encountering N.

## Tests required

1. encode_base: A→0, C→1, G→2, T→3, N→None
2. decode_base roundtrips correctly
3. encode_kmer + decode_kmer roundtrip for k=5, k=31
4. encode_kmer returns None for sequences with N
5. encode_kmer returns None for k>31
6. reverse_complement("ACGT") = encode("ACGT") (self-complementary palindrome)
7. reverse_complement("AAAA", k=4) = encode("TTTT")
8. canonical always returns the smaller of kmer and revcomp
9. canonical_with_strand returns Forward when kmer is already canonical
10. canonical_with_strand returns Reverse when revcomp is smaller
11. KmerIterator yields correct count: seq of length L → L-k+1 k-mers
12. KmerIterator skips positions with N (resets window)
13. KmerIterator sliding window gives same results as naive re-encoding

## Definition of done

- `cargo test -p kam-index` passes
- `cargo clippy -p kam-index -- -D warnings` passes
- Doc comments with examples
