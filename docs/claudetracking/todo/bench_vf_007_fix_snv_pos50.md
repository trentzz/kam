# BENCH-VF-007: Fix SNV benchmark edge position

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/bench_varforge.md)
**Priority**: medium
**Depends on**: BENCH-VF-006
**Status**: todo

## Goal

Move the SNV at position 50 to a position further from the chromosome edge.
Currently, the chr1:0-150 target region starts at reference position 0. The SNV at pos 50
is only 50bp from the chromosome start. Varforge may generate fewer reads covering this
position due to edge effects (fragments cannot extend before position 1). This causes
systematic 80% rather than 100% SNV sensitivity at all VAF levels.

Done looks like: SNV sensitivity reaches 100% at moderate VAF levels (≥0.25%); benchmark
targets updated; all 50 SNV datasets re-simulated and re-scored.

## Steps

1. Update `SNV_ROWS` in `make_snv_indel_suite.py`:
   - Change pos 50 to pos 300 (currently used for the first target, but that target will
     be adjusted too). Use a position that is ≥200bp from either chromosome end.
   - Ensure the new positions do not conflict with indel positions (100, 300, 550, 750, 1000).
   - Suggested new positions: 150, 350, 600, 850, 1100 (well-separated and away from edges).
   - Read the actual reference bases at these positions and use correct REF alleles.

2. Update `snvindel_targets.fa`:
   - Adjust or add target regions to cover the new SNV positions.
   - Keep 100bp flanking on each side; use 0-based start in FASTA headers.

3. Re-run `python3 docs/benchmarking/snvindel/scripts/make_snv_indel_suite.py` to regenerate
   truth VCFs.

4. Delete and re-simulate all 50 SNV datasets:
   ```
   rm -rf docs/benchmarking/snvindel/results/sim_snv_vaf*_[ab]
   for cfg in docs/benchmarking/snvindel/configs/snv_vaf*.yaml; do
       varforge simulate --config "$cfg"
   done
   ```

5. Delete and re-run kam on SNV datasets:
   ```
   rm -rf docs/benchmarking/snvindel/results/kam_snv_vaf*_[ab]
   bash docs/benchmarking/snvindel/scripts/run_snv_indel_suite.sh
   ```
   (This will re-run SNVs only if indel results already exist and are complete.)

   Actually, the run script runs both SNV and indel sequentially. Consider skipping
   already-done indel datasets using the `--force false` (default) logic.

6. Re-score: `python3 docs/benchmarking/snvindel/scripts/score_snv_indel_suite.py`

7. Regenerate figures: `python3 docs/benchmarking/snvindel/scripts/plot_varforge_sensitivity.py`

8. Update `docs/research/snvindel_benchmark_scoring_investigation.md` with the fix.

## Notes

- The indel datasets do NOT need to be re-run; only SNVs are affected.
- The run script has a `--force` flag; use default (no --force) to skip existing indel results.
- Reference sequence at new candidate positions:
  - pos 150: check with `python3 -c "ref=open(...).read()..."`
  - Avoid positions that are repetitive or low-complexity (check GC content).
