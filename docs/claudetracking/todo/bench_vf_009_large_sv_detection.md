# BENCH-VF-009: Investigate large SV detection reliability

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/bench_varforge.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Understand why large SVs (100bp inversions, 99bp insertions) are detected unreliably
in discovery mode. The DUP at pos 500 is consistently detected (50% sensitivity = 1/2
variants), but the INV at pos 900 is detected only at some VAF levels. INS and INVDEL
are near-zero until ≥1% VAF.

Done looks like: a documented root cause with evidence, and either a fix or a clear
statement that large SV detection requires architectural changes beyond the current scope.

## Steps

1. Read `docs/research/sv_benchmark_investigation.md`.

2. Compare the k-mer graph structure for DUP vs INV targets:
   - Run kam with debug logging enabled on sim_sv_vaf0050_a and sim_sv_vaf0100_a
   - Compare the path walk output for the DUP target (chr1:399-699) vs INV target (chr1:799-1099)
   - Identify why the INV path is inconsistently found

3. Check if the INV detection correlates with specific simulation seeds (replicate a vs b):
   - Do both replicates at the same VAF always give the same pattern, or is it random?
   - If random: the junction k-mer for 100bp INV survives path walking stochastically

4. Test with smaller SV sizes: generate a 20bp, 40bp, 60bp, 80bp inversion test case
   and check where detection becomes unreliable.

5. Write up findings in `docs/research/` and create targeted follow-up tasks if a fix exists.

## Notes

- Do not attempt to fix junctions: `--sv-junctions` causes DFS explosion (see investigation).
- The current pathfind approach (de Bruijn + DFS) may be fundamentally unsuited for
  large SVs. The investigation should honestly assess this.
- If detection requires architectural work (e.g., junction-indexed graph without DFS),
  this is a new epic, not a bug fix.
