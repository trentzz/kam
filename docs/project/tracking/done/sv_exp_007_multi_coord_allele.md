# SV-EXP-007: Multi-Coordinate Allele Extraction for Fusions

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_006_fusion_classify.md
**Status**: todo

## Goal

Extend allele extraction to handle fusion events where breakpoint coordinates
exist on two different target regions. Done looks like: a `Fusion` variant
record contains breakpoint positions for both partner A and partner B, suitable
for BND VCF notation.

## Steps

1. Read the current `MinimalAllele` or equivalent allele struct in
   `kam-call/src/` and understand the coordinate fields.

2. Determine whether `MinimalAllele` can hold two positions or whether a new
   `FusionAllele` struct is needed:
   - If `MinimalAllele` holds a single `(chrom, pos, ref, alt)`, extend it
     with optional `partner: Option<(chrom, pos)>`.
   - Or create a `FusionAllele { a_chrom, a_pos, b_chrom, b_pos, junction_seq }`.

3. In the allele extraction step for `Fusion` calls:
   - Extract the junction sequence from the alt path.
   - Determine the breakpoint position within each partner from the k-mer
     annotations (last k-mer fully in A, first k-mer fully in B).
   - Populate both coordinate fields.

4. Pass the fusion allele coordinates to the VCF writer (used in SV-EXP-009).

5. Write unit tests:
   - Fusion allele extraction produces correct breakpoint coordinates for a
     known junction sequence.
   - Coordinates are consistent with the junction k-mers (no off-by-one).

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- Breakpoint precision is limited by k-mer size. The reported position is the
  last reference position fully covered by the path before the junction, which
  may be `k - 1` bp before the true breakpoint.
- Document this limitation in a comment on the struct.
