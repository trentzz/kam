# SV-EXP-002: InvDel Classification

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_001_new_variant_types.md
**Status**: todo

## Goal

Classify a path as `InvDel` instead of `LargeDeletion` when the differing
region between the ref and alt path sequences contains a reverse-complement
segment. Done looks like: a synthetic path with a deletion flanking an
inverted region is classified as `InvDel`; a plain deletion is still
classified as `LargeDeletion`.

## Steps

1. Read `kam-call/src/classify.rs` (or the file containing `classify_variant`)
   in full.

2. Locate the `LargeDeletion` classification branch. It currently triggers when
   the alt path is shorter than the ref path by ≥50 bp (or similar threshold).

3. Before emitting `LargeDeletion`, check for an inverted sub-segment:
   - Extract the differing region from the ref path sequence.
   - Compute the reverse complement of substrings of that region.
   - If any substring of length ≥20 bp is present in the alt path sequence
     as its reverse complement, classify as `InvDel`.

4. Add a helper function `contains_rc_segment(ref_region: &[u8], alt_seq: &[u8],
   min_len: usize) -> bool` in `classify.rs`.

5. Write unit tests:
   - Ref path: `AAAA[CCCCGTGT]AAAA` (deletion of 8 bp). Alt: shorter.
     → `LargeDeletion`.
   - Ref path with inverted flanking: alt contains RC of a ref sub-segment.
     → `InvDel`.
   - Edge case: inverted segment shorter than `min_len` → still `LargeDeletion`.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- The RC check is an additional heuristic on top of the existing length check.
  It does not change the threshold for calling `LargeDeletion`; it only
  re-classifies after the deletion is detected.
- `min_len = 20` is a reasonable starting value. Make it a named constant.
- The exact sequences available in the classifier depend on the path walking
  output. Read the classify function carefully to understand what ref/alt
  sequences are accessible.
