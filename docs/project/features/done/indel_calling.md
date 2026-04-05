# Indel Calling

## Status

Done.

## Summary

Indel variants (insertions and deletions) are now called and scored correctly.
This required two changes: proper VCF allele minimisation in the Rust output
layer, and left-normalisation in both the Rust code and the benchmarking
scoring script.

## Problem

Before this feature, indel sensitivity was zero across all tested depths and
VAF levels. The root causes were:

1. **No VCF coordinate extraction.** The `write_vcf` function in
   `kam-call/src/output.rs` wrote the raw target_id string as CHROM and
   hardcoded POS=1. Truth VCFs use proper genomic coordinates, so no indel
   (or SNV) could ever match.

2. **Broken minimisation in scoring script.** The inner common-suffix loop in
   `run_titration_batch.py` used `< len(alt_trimmed) - 1` as its bound. This
   prevented stripping the final character of `alt_trimmed`, so a pure
   deletion (where `alt_min` should reduce to empty) was never minimised
   correctly, and the extracted deleted sequence was wrong.

3. **No left-normalisation.** Standard truth VCFs are left-aligned. Indels in
   repeat regions are ambiguous (e.g., deleting one G from "GGG" can be
   anchored at any of the three positions). Without left-normalisation, called
   variants and truth variants differed in position even when they described
   the same biological event.

## Solution

### `kam-call/src/allele.rs` (new module)

A self-contained module that converts full-sequence path representations
(ref_seq / alt_seq of 100+ bp) to minimal left-normalised VCF alleles.

Public API:

```rust
pub fn extract_minimal_allele(
    target_id: &str,   // "chrN:START-END" (1-based closed interval)
    ref_seq: &[u8],
    alt_seq: &[u8],
) -> Option<MinimalAllele>

pub struct MinimalAllele {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}
```

Algorithm:
1. Find the first differing position (`diff_pos`).
2. Trim the outer common suffix from both sequences.
3. For SNVs/MNVs (equal length): return the trimmed region directly.
4. For indels: strip any residual inner common suffix and prefix to isolate
   the inserted or deleted bases.
5. Left-normalise: while `anchor_pos > 0` and the preceding reference base
   equals the last base of the event sequence, rotate the event sequence
   right by one and decrement `anchor_pos`.
6. Prepend the anchor base as required by the VCF specification.

Coordinate note: `chrN:START-END` uses a 1-based closed interval where START
is the genomic position of the first base in the stored sequence. A base at
sequence index `i` (0-based) maps to VCF position `START + i` (1-based).

### `kam-call/src/output.rs`

`write_vcf` now calls `extract_minimal_allele` to convert each call to
proper CHROM/POS/REF/ALT before writing the VCF record. Falls back to the
alignment-free format (target_id as CHROM, POS=1) when the target_id cannot
be parsed.

### `benchmarking/scripts/run_titration_batch.py`

Fixed the `extract_called_variants` function to:
- Use `< len(alt_trimmed)` (not `- 1`) in the inner common-suffix loop.
- Apply the same left-normalisation algorithm as the Rust code.

## Tests

`kam-call/src/allele.rs` contains 10 unit tests and 1 doc test:

- SNV at correct genomic position
- Identical sequences return `None`
- Invalid target_id returns `None`
- Deletion in GG run left-normalised to T anchor
- Simple insertion with no repeat (no shift)
- MNV (equal-length, no normalisation)
- Multi-base deletion with distinct flanks (no shift)
- `parse_target_id` with standard coordinates
- `parse_target_id` rejects malformed inputs
- Insertion in TT run left-normalised

`kam-call/src/output.rs` test 8 (`vcf_uses_genomic_coordinates_for_snv`)
verifies that the VCF writer emits correct CHROM, POS, REF, and ALT for a
`chrN:START-END` target_id.

## Benchmark Results

Indel sensitivity at 250k reads, 5 ng input, 1% VAF: ~0.27.
Indel sensitivity at 500k reads, 5 ng input, 1% VAF: ~0.31.
Full per-depth results in `benchmarking/results/tables/`.
