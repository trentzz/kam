# Task SV-002: Fix Inversion Classifier in classify_variant

## Location

`kam-call/src/caller.rs` — `classify_variant` function (line 357) and helpers.

## What to implement

Add `partial_inversion_len(ref_seq, alt_seq) -> Option<usize>` that detects
whether the differing region between two equal-length sequences is a reverse
complement of the corresponding reference segment.

Algorithm:
1. Find `left` = first index where `ref_seq[i] != alt_seq[i]`.
2. Find `right` = last index where `ref_seq[i] != alt_seq[i]`.
3. If no differing positions, return None.
4. Check `is_reverse_complement(&ref_seq[left..=right], &alt_seq[left..=right])`.
5. If true, return `Some(right - left + 1)`.

In `classify_variant`, after the existing full-path RC check, add:
```rust
if let Some(inv_len) = partial_inversion_len(ref_seq, alt_seq) {
    if inv_len >= SV_LENGTH_THRESHOLD {
        return VariantType::Inversion;
    }
}
```

## Tests to add

In `#[cfg(test)]` block:
- Partial inversion detected: ref=`AAAA[ACGT]TTTT`, alt=`AAAA[ACGT_RC]TTTT`
  where the central 4+ bases are RC'd and flanks match.
- Partial inversion below threshold (< 50 bp) → MNV, not Inversion.
- Full-path inversion still classified as Inversion (existing test must pass).
- SNV inside an otherwise identical window → still SNV.

## Definition of done

- `cargo test -p kam-call` passes with new tests.
- `cargo clippy -p kam-call -- -D warnings` passes.
- `cargo fmt -- --check` passes.
