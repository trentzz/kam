# BENCH-SV: SV Benchmarking Expansion

**Status**: active
**Priority**: medium

## Goal

Expand the SV benchmarking suite to cover more structural variant types. The
current benchmark tests deletion (DEL), tandem duplication (DUP), and inversion
(INV). Additional types to add:

- Insertion (INS): de novo sequence inserted, no reference anchor at end
- Large deletion (>500 bp): tests whether coverage spans the junction
- Breakend / translocation (TRA): inter-chromosomal or long-range rearrangement
- Complex inversion with flanking deletions (INVDEL)

Each type gets a truth VCF and a varforge config committed to the repo.

## Motivation

The current benchmark only tests three SV types. Insertions are structurally
different from deletions (no reference sequence at the alt allele) and have not
been tested. Large deletions and translocations probe different failure modes in
the path-walking stage.

## Child tasks

| ID | File | Status |
|----|------|--------|
| BENCH-SV-001 | done/bench_sv_001_new_sv_types.md | done |

## Scope

- `docs/benchmarking/sv/data/` — truth VCFs for new SV types
- `docs/benchmarking/sv/configs/` — varforge configs for new SV types

## Out of scope

- Changes to the kam Rust code
- Re-running the existing SV benchmark
