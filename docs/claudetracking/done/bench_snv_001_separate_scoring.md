# BENCH-SNV-001: Separate SNV and Indel Scoring

**Epic**: BENCH-SNV (`overallplans/BENCH-SNV.md`)
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

Update `score_variants.py` to classify each variant as SNV or indel and report
separate TP/FP/FN counts for each type. The overall metrics remain, and
per-type metrics are added alongside them.

Done when: `score_variants.py --truth ... --called ...` outputs a single TSV
row with both `snv_tp`, `snv_sensitivity`, `indel_tp`, `indel_sensitivity`
columns alongside the existing overall columns.

## Steps

1. Add a helper to classify (REF, ALT) as snv, indel, or other.
   - SNV: len(REF)==1 and len(ALT)==1
   - Indel: otherwise
2. Update `load_vcf_variants` to return a dict keyed by type (snv/indel/all).
   The truth VCF has `TYPE=SNP` or `TYPE=INDEL` in INFO; use that when present.
   Fall back to length-based classification.
3. Compute `compute_metrics` separately for each subset.
4. Add per-type columns to COLUMNS and the output row.
5. Update the stderr summary line to print per-type sensitivity.

## Notes

- The truth VCF at `scripts/truth_variants.vcf` uses `TYPE=SNP` / `TYPE=INDEL`
  in INFO. Use this rather than length inference where possible.
- Keep the existing overall columns intact so downstream scripts don't break.
- Called VCF may not have TYPE tags; use length-based classification for those.
