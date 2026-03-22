# Inversion Detection: Investigation and Fix

## Symptom

SV benchmark results show no `Inversion`-classified calls on the INV target
(`chr1:849-1049_INV_100bp`) at any VAF level (0.5%–5%). Instead, the target
produces a single high-evidence `MNV` call alongside many low-evidence `SNV`
calls. At 5% VAF, the MNV call has 23 supporting molecules, a confidence of
1.0, and a PASS filter — all indicating a real variant is detected. The variant
type classification is wrong.

## Diagnostic Data

Inspecting the 5% VAF MNV call:

- **ref_seq**: 200 bp reference window (`chr1:849-1049`)
- **alt_seq**: 200 bp sequence sharing an identical 51 bp prefix and 49 bp
  suffix with the reference; only the central 100 bp differs
- The central alt segment is the reverse complement of the corresponding
  central ref segment — confirmed by manual comparison

The graph traversal via junction k-mers therefore **succeeded**: it found the
inverted alt path and scored it with 23 molecules. Classification failed.

## Root Cause

`classify_variant` in `kam-call/src/caller.rs` checks inversion with:

```rust
if is_reverse_complement(ref_seq, alt_seq) {
    return VariantType::Inversion;
}
```

`is_reverse_complement` compares the **entire** ref path against the **entire**
alt path. For a partial inversion — where only a central segment is inverted
and the flanks remain unchanged — this check always fails. The function then
falls through to count differing positions, finding many (the full 100 bp
inverted region), and returns `MNV`.

No synthetic path construction is needed: the graph already finds the correct
inversion path using the junction k-mers from `sv_junctions.fa`. The only
failure is post-hoc classification.

## Fix

Add `partial_inversion_len(ref_seq, alt_seq) -> Option<usize>` to find the
contiguous differing region between two equal-length sequences and test whether
that region is a reverse complement:

1. Find `left` = first index where `ref[i] != alt[i]`.
2. Find `right` = last index where `ref[i] != alt[i]`.
3. Extract `ref_mid = ref[left..=right]`, `alt_mid = alt[left..=right]`.
4. Return `Some(right - left + 1)` if `is_reverse_complement(ref_mid, alt_mid)`.

In `classify_variant`, insert this check after the existing full-path RC check:

```rust
if let Some(inv_len) = partial_inversion_len(ref_seq, alt_seq) {
    if inv_len >= SV_LENGTH_THRESHOLD {
        return VariantType::Inversion;
    }
}
```

This correctly classifies the existing 23-molecule MNV call as an Inversion
without changing any indexing or path-finding logic.

## Expected Outcome

After the fix:

| VAF  | INV detected | Expected filter |
|------|-------------|----------------|
| 5%   | Yes         | PASS            |
| 2%   | Yes         | PASS or LowConfidence |
| 1%   | Uncertain   | LowConfidence   |
| 0.5% | Uncertain   | LowConfidence or not detected |

## What Was Implemented

See task SV-002 for the implementation details and test results.

## Open Questions

- At 0.5% and 1% VAF, the graph may not find the inversion path if junction
  k-mer coverage is too sparse. If so, a synthetic path construction analogous
  to `synthesize_dup_alt_path` could be added as a fallback.
- The threshold for partial inversion length is currently `SV_LENGTH_THRESHOLD`
  (50 bp). This should be reviewed for smaller inversions.
