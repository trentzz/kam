# RUNTIME-002: Time Alignment Pipeline on Same Samples

**Epic**: RUNTIME (docs/claudetracking/overallplans/RUNTIME.md)
**Priority**: medium
**Depends on**: runtime_001_kam_timing.md
**Status**: todo

## Goal

Record per-stage wall-clock times for the alignment-based pipeline (HUMID +
Jellyfish + km + kmtools) on the same samples used in RUNTIME-001. Done looks
like: a TSV at `docs/benchmarking/runtime/alignment_times.tsv` with the same
schema as `kam_times.tsv`.

## Steps

1. Use the same 4–8 samples from RUNTIME-001.

2. Write `docs/benchmarking/runtime/scripts/run_alignment_timed.sh`:
   - Stage 1 (HUMID): time HUMID execution on the FASTQ pair.
   - Stage 2 (Jellyfish count + dump): time Jellyfish separately.
   - Stage 3 (km + kmtools): time combined.
   - Wrap each stage with `/usr/bin/time -v`.

3. Run the script on each sample. Run 3 times; record the median per stage.

4. Write results to `docs/benchmarking/runtime/alignment_times.tsv`:
   - Columns: sample, read_pairs, total_wall_s, humid_wall_s,
     jellyfish_wall_s, km_wall_s, peak_rss_mb.

5. If the full alignment pipeline is not available (e.g. HUMID/km not installed
   on the same machine), document the limitation and use a single BWA + GATK
   HaplotypeCaller run as an alternative alignment pipeline baseline.

## Notes

- Run the alignment pipeline on the same machine as RUNTIME-001 to ensure
  fair comparison (same CPU, same disk).
- Both pipelines should process the same FASTQ files. Do not use pre-sorted
  or pre-deduplicated inputs for the alignment pipeline.
- If HUMID is not available, document which tools were used as the alignment
  pipeline equivalent and what their roles are.
