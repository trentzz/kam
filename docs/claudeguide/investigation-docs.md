# Investigation Documentation

Whenever a diagnostic investigation uncovers a non-obvious finding, write it
up in `docs/research/` before moving on. This applies to: root causes of bugs,
unexpected benchmark results, failed hypotheses, surprising measurements.

## Required sections

1. **Symptom** — what went wrong or looked unexpected, and the diagnostic data
   that revealed it.
2. **Hypothesis** — what was tested and whether it was confirmed or disproved.
3. **Root cause** — the actual cause, with supporting evidence.
4. **Fix** — what would fix it and the expected impact.
5. **Result** — what was implemented and what the measured outcome was.

## Existing write-ups (for style reference)

- `docs/research/sensitivity_investigation.md`
- `docs/research/anchor_missing_investigation.md`
- `docs/research/sv_threshold_investigation.md`
- `docs/research/inversion_detection.md`

A fix without documentation is half the work. The investigation log is the
primary context for future sessions and paper writing.
