# SV-SCORE-003: Re-run SV benchmark suite after scoring fix

**Epic**: SV-SCORE (docs/claudetracking/overallplans/SV-SCORE.md)
**Priority**: medium
**Depends on**: SV-SCORE-001, SV-SCORE-002
**Status**: todo

## Goal

Re-run the full SV benchmark suite with the updated binary and update the
summary figures and tables. Done looks like: `sv_summary.tsv` and
`sv_sensitivity.png` reflect results from the new scoring model; the INV
NALT values are substantially higher than 1–3 at moderate VAF.

## Steps

1. Build the release binary:
   ```bash
   cargo build --release
   ```

2. Re-run the SV suite with `--force` to overwrite existing results:
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

5. Check `sv_summary.tsv` for improvement in INV sensitivity at low VAF.
   Document findings in `docs/research/sv_large_detection_investigation.md`
   under a new "SV-SCORE results" section.

## Notes

- Expected outcome: INV NALT at 1% VAF should increase from 1–3 to a value
  proportional to the true molecule count (~10 at 1% VAF, 5000x coverage).
- If INV sensitivity at low VAF does not improve, check that
  `mean_variant_specific_molecules` is being populated correctly and that the
  variant is classified as `Inversion`.
- StrandBias will still filter INV calls at ≥2% VAF. That is addressed in the
  SV-STRAND epic and is expected in these results.
