# Task SV-005: Write SV Sensitivity Report

## What to do

After tasks SV-003 and SV-004 are complete, write
`docs/research/sv_sensitivity_report.md` summarising:

1. **Results table**: TP/FP/FN, filter status, and molecule count for each
   SV type at each VAF level in both discovery and tumour-informed mode.

2. **Analysis**:
   - At what VAF does each SV type become detectable?
   - Does tumour-informed mode eliminate all FPs?
   - How do molecule counts compare to SNV/indel performance at equivalent VAF?
   - What are the remaining limitations?

3. **Open questions** for future work:
   - Inversion sensitivity at 0.5% and 1% VAF.
   - Whether synthetic path construction would improve DUP at very low VAF.
   - Translocation and foreign insertion detection (Phase 2).

## Definition of done

- `docs/research/sv_sensitivity_report.md` exists with the results table
  and analysis sections filled in.
- No placeholder sections remain.
