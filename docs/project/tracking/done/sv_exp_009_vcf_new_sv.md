# SV-EXP-009: VCF Output for New SV Types

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_007_multi_coord_allele.md, sv_exp_008_sv_thresholds.md
**Status**: todo

## Goal

Update VCF output to correctly represent `InvDel`, `NovelInsertion`, and
`Fusion`. Done looks like: VCF records for InvDel include `SVTYPE=INVDEL`;
NovelInsertion records include `SVTYPE=INS;NOVEL=1`; Fusion records use BND
notation with two records per event.

## Steps

1. Read the VCF output code in `kam-call/src/output.rs` or `kam/src/output.rs`
   in full.

2. For `InvDel`:
   - Add `SVTYPE=INVDEL` to the INFO field.
   - REF and ALT follow standard SV notation (symbolic allele or sequence
     allele, whichever is current convention).

3. For `NovelInsertion`:
   - Set `SVTYPE=INS` and add `NOVEL=1` to INFO.
   - Add `NOVEL` to the VCF header `INFO` definitions.

4. For `Fusion`:
   - Emit two BND records per event (standard VCF 4.3 BND notation).
   - Record 1: breakpoint on partner A, ALT contains partner B coordinate.
   - Record 2: breakpoint on partner B, ALT contains partner A coordinate.
   - Add `SVTYPE=BND` to INFO.
   - Both records share the same EVENT INFO tag (e.g. the fusion name or a
     generated ID).

5. Update the VCF header writer to include new INFO and FORMAT definitions for
   the new fields (`NOVEL`, `EVENT`).

6. Write tests:
   - VCF output for an `InvDel` call contains `SVTYPE=INVDEL`.
   - VCF output for a `NovelInsertion` call contains `SVTYPE=INS;NOVEL=1`.
   - VCF output for a `Fusion` call produces exactly two BND records with
     correct ALT notation.

7. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- BND notation: `ALT = T]chr2:321682]` for a translocation to chr2:321682
  on the forward strand. Follow VCF 4.3 spec section 5.4.
- If the current VCF writer does not support symbolic ALT alleles or multi-
  record events, extend it minimally. Do not refactor the whole writer.
- The TSV output format (if supported) can use `InvDel`, `NovelInsertion`,
  and `Fusion` as plain strings in the type column.
