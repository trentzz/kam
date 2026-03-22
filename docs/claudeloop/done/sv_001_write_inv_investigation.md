# Task SV-001: Document Inversion Detection Root Cause

## What to do

Write `docs/research/inversion_detection.md` documenting the investigation into
why inversions are classified as MNV instead of Inversion.

## Findings (already established)

**Symptom**: SV benchmark produces MNV calls on the INV target instead of
Inversion calls, even at 5% VAF where 23 supporting molecules are found.

**Root cause**: `classify_variant` in `kam-call/src/caller.rs` (line 362) checks
`is_reverse_complement(ref_seq, alt_seq)` against the **entire** target window
path. A real inversion only inverts a central segment — the flanking sequences
remain unchanged. So the full-path RC check always fails.

**Evidence**: The MNV call at 5% VAF on `chr1:849-1049_INV_100bp` has ref and
alt sequences that share identical prefixes and suffixes (~51 bp each), with
only the central 100 bp differing. That central alt segment IS the RC of the
corresponding ref segment — confirming the graph found the correct path.

**Fix required**: Add partial-inversion detection to `classify_variant`. Find the
leftmost and rightmost differing positions, extract both central segments, and
check `is_reverse_complement(ref_mid, alt_mid)`. If true and segment length ≥ 50
bp, return `VariantType::Inversion`.

**No synthetic path needed**: Unlike DUP, the graph traversal already succeeds
using the junction k-mers. The problem is purely in post-hoc classification.

## Definition of done

- `docs/research/inversion_detection.md` exists with symptom, diagnostic,
  root cause, fix, and expected outcome sections.
