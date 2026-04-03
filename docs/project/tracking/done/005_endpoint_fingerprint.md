# Task 005: Endpoint Fingerprint for UMI Collision Detection

## Location
`kam-assemble/src/fingerprint.rs`

## What to implement

A function that computes a 64-bit fingerprint from the template endpoints of a read pair. This is used as a secondary grouping signal to detect UMI collisions — two different molecules that happen to share the same 5bp random UMI.

## Context

With only 1,024 possible UMIs (Twist 5bp random), collisions are non-trivial at high depth. Two reads with identical canonical UMI pairs but clearly different template sequences are likely from different molecules. The endpoint fingerprint provides a fast check.

## Interface

```rust
/// Compute a 64-bit fingerprint from template read pair endpoints.
/// Uses the first and last FINGERPRINT_BASES (8) of each template read.
///
/// Two read pairs from the same molecule will have nearly identical
/// fingerprints (within sequencing error tolerance).
/// Two read pairs from different molecules at the same locus with
/// a UMI collision will have different fingerprints.
pub fn compute_endpoint_fingerprint(
    template_r1: &[u8],
    template_r2: &[u8],
) -> u64;

/// Check if two fingerprints are compatible (could be from the same molecule).
/// Allows for some bit differences due to sequencing errors.
pub fn fingerprints_compatible(fp1: u64, fp2: u64) -> bool;
```

## Implementation notes

- Use first 8bp and last 8bp of template_r1, and first 8bp and last 8bp of template_r2
- If template is shorter than 16bp, use entire template
- 2-bit encode bases and pack into u64
- For compatibility check, count bit differences (popcount of XOR) — threshold at ~4 bit differences (accounts for 1-2 base sequencing errors across 4 endpoints)

## Tests required

1. Identical templates produce identical fingerprints
2. Templates differing by 1 base produce compatible fingerprints
3. Completely different templates produce incompatible fingerprints
4. Short templates (< 16bp) don't panic
5. Templates with N bases handled (treat N as A or use 5th encoding)
6. Fingerprint is deterministic (same input always same output)

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
- Doc comments with examples
