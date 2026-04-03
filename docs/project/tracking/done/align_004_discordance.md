# ALIGN-004: Discordance Analysis

**Epic**: ALIGN-COMPARE (docs/claudetracking/overallplans/ALIGN-COMPARE.md)
**Priority**: high
**Depends on**: align_003_concordance.md
**Status**: todo

## Goal

Categorise every discordant (variant × sample) pair and identify whether
discordance is random (stochastic at low VAF) or systematic (a class of
variants one method consistently misses). Done looks like: a discordance
analysis report in `docs/benchmarking/align_compare/discordance_analysis.md`
with a categorised table of discordant cases.

## Steps

1. Extract discordant rows from `concordance_table.tsv`:
   - `kam_only`: detected by kam but not alignment.
   - `alignment_only`: detected by alignment but not kam.

2. For each discordant case, assign a category:

   **VAF-below-threshold**: true VAF ≤ 0.25% and the missing method has no
   call at all (not just filtered). Likely stochastic low-VAF miss.

   **Filter-mismatch**: variant was called by both methods but filtered
   differently (e.g. alignment PASS, kam LowConfidence). Record both filters.

   **Systematic-type**: >80% of samples for this variant are discordant in
   the same direction. Indicates a structural issue (e.g. a repeat region, a
   GC-rich context, an SV not covered by the k-mer target).

   **Isolated**: discordant in only 1–2 samples with no pattern. Flag for
   manual inspection.

3. Write `docs/benchmarking/align_compare/scripts/categorise_discordance.py`
   to automate steps 1–2 and output
   `docs/benchmarking/align_compare/discordance_table.tsv`.

4. Write `discordance_analysis.md`:
   - Counts per category.
   - Top 10 systematic discordant variants with explanation.
   - Root cause hypothesis for each systematic case.

## Notes

- Systematic cases are the most interesting for the paper. A variant that is
  always missed by alignment but always detected by kam (or vice versa) reveals
  a fundamental difference in the two approaches.
- Manual inspection of systematic cases may require reading alignment BAMs or
  kam QC logs. Document tools used.
