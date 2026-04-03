# RUNTIME: Speed Comparison vs Alignment-Based Pipeline

**Status**: todo
**Priority**: medium
**Branch**: epic/RUNTIME

## Goal

Measure and compare wall-clock runtime for kam vs the alignment-based pipeline
(HUMID + Jellyfish + km + kmtools, or BWA + GATK) on the same titration
samples. Produce a per-stage bar chart showing where time is spent in each
pipeline.

## Motivation

Speed is a secondary selling point for kam: alignment-free processing should be
faster, especially on small targeted panels. Quantifying this difference gives
a concrete runtime claim for the methods paper and helps identify bottlenecks
for future optimisation.

## Design

### Samples

Run on the same 24 titration samples used in ALIGN-COMPARE, or a representative
subset (e.g. 4 samples at 1%, 5%, 10%, 50% tumour fraction).

### Measurement

Record wall-clock time using `time` or `/usr/bin/time -v` per stage. For kam,
use per-stage times from the QC JSON output (`stage_wall_seconds` fields).
For the alignment pipeline, wrap each tool invocation.

### Stages

| kam stage | Alignment stage equivalent |
|-----------|---------------------------|
| assemble | HUMID |
| index | Jellyfish count |
| pathfind | km + kmtools |
| call | (no equivalent, included in km) |

Report total pipeline time and per-stage breakdown.

### Figure

Bar chart: x = pipeline stage (grouped by pipeline), y = wall-clock seconds.
One panel per sample size (reads count). Include error bars if multiple samples
are measured.

## Child tasks

| ID | File | Status |
|----|------|--------|
| RUNTIME-001 | todo/runtime_001_kam_timing.md | todo |
| RUNTIME-002 | todo/runtime_002_alignment_timing.md | todo |
| RUNTIME-003 | todo/runtime_003_runtime_figure.md | todo |

## Dependencies

- ALIGN-002 (kam run on titration samples) — reuse the same runs with timing.
- Access to the alignment pipeline tools on the same HPC node.

## Scope

- `docs/benchmarking/runtime/` — timing scripts, raw timing TSVs, figures
- `docs/benchmarking/runtime/scripts/` — shell wrappers for timing

## Out of scope

- Profiling internal bottlenecks within kam (a separate perf epic)
- Parallelism scaling experiments
- Memory usage comparison
