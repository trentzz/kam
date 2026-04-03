# Task 004: CanonicalUmiPair Implementation

## Location
`kam-core/src/molecule.rs`

## What to implement

The `CanonicalUmiPair` constructor and strand determination. Two read pairs from the same molecule (one from each strand) must produce the same `CanonicalUmiPair`, and we need to know which strand a given R1 came from.

## Interface (implement exactly this)

```rust
impl CanonicalUmiPair {
    /// Create a canonical UMI pair from R1 and R2 UMIs.
    /// The lexicographically smaller UMI is always stored as umi_a.
    /// This ensures forward and reverse strand reads from the same molecule
    /// hash to the same key.
    ///
    /// # Example
    /// ```
    /// let pair = CanonicalUmiPair::new(*b"ACGTA", *b"TGCAT");
    /// assert_eq!(pair.umi_a, *b"ACGTA");
    /// assert_eq!(pair.umi_b, *b"TGCAT");
    /// ```
    pub fn new(umi_r1: [u8; 5], umi_r2: [u8; 5]) -> Self { ... }

    /// Determine which strand R1 came from based on whether its UMI
    /// is the canonical "a" (Forward) or "b" (Reverse).
    pub fn strand_of_r1(&self, umi_r1: &[u8; 5]) -> Strand { ... }
}
```

## Tests required

1. `new()` with umi_r1 < umi_r2 lexicographically — umi_a should be umi_r1
2. `new()` with umi_r1 > umi_r2 — umi_a should be umi_r2 (swapped)
3. `new()` with umi_r1 == umi_r2 — should work (same UMI, possible self-complementary)
4. `strand_of_r1()` returns Forward when umi_r1 matches umi_a
5. `strand_of_r1()` returns Reverse when umi_r1 matches umi_b
6. Two read pairs from same molecule (UMIs swapped) produce identical `CanonicalUmiPair`
7. `CanonicalUmiPair` implements Hash consistently (equal pairs have equal hashes)
8. All-A UMI vs all-T UMI orders correctly

## Definition of done

- `cargo test -p kam-core` passes
- `cargo clippy -p kam-core -- -D warnings` passes
- Doc comments with examples
