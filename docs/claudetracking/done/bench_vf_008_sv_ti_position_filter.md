# BENCH-VF-008: Position-based tumour-informed filter for large SVs

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/bench_varforge.md)
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

Implement a position-based tumour-informed filter mode for large structural variants.
Currently, tumour-informed mode uses exact REF/ALT matching against target variants.
For large SVs (≥20bp), kam detects only partial alleles (e.g., a 1bp insertion at the
DUP junction instead of the full 100bp duplication). No match is possible against the
full truth allele, giving 0% tumour-informed sensitivity for all SVs.

Done looks like: a `--ti-position-tolerance N` flag (or similar) that enables position-based
tumour-informed matching; SV tumour-informed sensitivity matches or exceeds discovery sensitivity.

## Steps

1. Read `docs/research/sv_benchmark_investigation.md` for context.

2. In `kam-call/src/` (or `kam/src/`), add a position-based matching mode to the
   tumour-informed filter:
   - If `--ti-position-tolerance N` is set, a call PASSES tumour-informed if its position
     is within N bp of any target variant position, regardless of REF/ALT.
   - Default: 10bp (matching the SV scoring tolerance).

3. Update `docs/benchmarking/sv/scripts/run_sv_suite.sh` to use the new flag when running
   in tumour-informed mode for SV datasets.

4. Re-run SV tumour-informed mode; re-score; update figures.

## Notes

- This is a new CLI flag; check CLAUDE.md architecture decisions before adding.
- The position-based filter is less strict than REF/ALT matching. Document the trade-off.
- For SNVs and small indels, the REF/ALT exact match is preferred (keeps FPs low).
  Consider making the position tolerance opt-in rather than the default.
