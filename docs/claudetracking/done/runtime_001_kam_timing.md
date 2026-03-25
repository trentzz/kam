# RUNTIME-001: Time kam on Titration Samples

**Epic**: RUNTIME (docs/claudetracking/overallplans/RUNTIME.md)
**Priority**: medium
**Depends on**: align_002_run_titration.md
**Status**: todo

## Goal

Record per-stage wall-clock times for `kam run` on a representative subset of
the titration samples. Done looks like: a TSV at
`docs/benchmarking/runtime/kam_times.tsv` with one row per (sample × run),
columns for each pipeline stage and total time.

## Steps

1. Select a representative subset of 4–8 samples at different tumour fractions
   (e.g. 1%, 5%, 10%, 50%) to characterise how runtime scales with coverage.

2. For each selected sample, run:
   ```bash
   /usr/bin/time -v kam run --r1 ... --r2 ... --targets ... \
       --output /dev/null 2> timing_raw/${sample}_time.txt
   ```
   Also save the kam QC JSON for per-stage breakdowns.

3. Parse the QC JSON for per-stage times (fields such as `assemble_wall_s`,
   `index_wall_s`, `pathfind_wall_s`, `call_wall_s` — confirm field names by
   reading an actual QC JSON output).

4. Parse `timing_raw/*.txt` for total wall-clock and peak RSS from
   `/usr/bin/time -v`.

5. Write results to `docs/benchmarking/runtime/kam_times.tsv`:
   - Columns: sample, read_pairs, total_wall_s, assemble_wall_s,
     index_wall_s, pathfind_wall_s, call_wall_s, peak_rss_mb.

6. Run each sample 3 times to check for variance; record the median.

## Notes

- Reuse the runs from ALIGN-002 if timing was recorded there. If not, re-run
  with `/usr/bin/time -v` wrapping.
- Per-stage times from QC JSON are more accurate than timing the binary
  externally. Use QC JSON values as primary; `/usr/bin/time -v` as secondary.
- Record the machine spec (CPU, RAM, disk type) in
  `docs/benchmarking/runtime/MACHINE_INFO.md`.
