# SV-STRAND-003: Re-run SV benchmark suite after strand bias fix

**Epic**: SV-STRAND (docs/claudetracking/overallplans/SV-STRAND.md)
**Priority**: medium
**Depends on**: SV-STRAND-001, SV-STRAND-002
**Status**: todo

## Goal

Re-run the full SV benchmark suite with the updated binary and update the
summary figures and tables. Done looks like: `sv_summary.tsv` shows INV
sensitivity at moderate and high VAF (2–10%) increased substantially from
the StrandBias-filtered baseline; figures are regenerated.

## Steps

1. Build the release binary:
   ```bash
   cargo build --release
   ```

2. Re-run the SV suite with `--force`:
   ```bash
   bash docs/benchmarking/sv/scripts/run_sv_suite.sh --force
   ```

3. Re-score and update the summary:
   ```bash
   python3 docs/benchmarking/sv/scripts/score_sv_suite.py
   ```

4. Regenerate figures:
   ```bash
   python3 docs/benchmarking/sv/scripts/plot_sv_sensitivity.py
   ```

5. Document findings in `docs/research/sv_large_detection_investigation.md`
   under a new "SV-STRAND results" section. Note how INV sensitivity changes
   across the VAF sweep.

## Notes

- Expected outcome: INV calls at 2–10% VAF that were StrandBias-filtered
  should now PASS, increasing sensitivity from near 0% to >50%.
- At low VAF (0.05–1.25%), the min_molecules bottleneck is the dominant
  failure mode. Improvement there requires SV-SCORE.
- If SV-SCORE-003 was run before this task, re-run with `--force` again so
  results reflect both fixes.
