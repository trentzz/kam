# BENCH-SV-NEW: Varforge Benchmark Suite for New SV Types

**Status**: done
**Priority**: high
**Branch**: epic/BENCH-SV-NEW

## Goal

Extend the varforge benchmark suite to cover the three new SV types added in
SV-EXPAND: InvDel, NovelInsertion, and Fusion. Each type gets 50 configs
(25 VAF levels × 2 replicates), matching the existing BENCH-VARFORGE design.
Results are scored and visualised as sensitivity curves per type.

## Motivation

SV-EXPAND adds new variant types to the classifier. Without benchmarks, there
is no evidence that classification is correct or that sensitivity is acceptable.
The existing `sv_vaf*` suite covers DEL, DUP, and INV only. This epic closes
the gap.

## Design

### VAF levels

Same 25 levels as BENCH-VARFORGE:

```
0.05%, 0.075%, 0.10%, 0.125%, 0.15%, 0.20%, 0.25%, 0.30%, 0.40%, 0.50%,
0.60%, 0.75%, 1.00%, 1.25%, 1.50%, 2.00%, 2.50%, 3.00%, 4.00%, 5.00%,
6.00%, 7.00%, 8.00%, 9.00%, 10.00%
```

Two replicates (seeds A and B) at each level.

### Variant classes

| Class | Prefix | Varforge support |
|-------|--------|-----------------|
| InvDel | `invdel_vaf*` | Add to varforge engine |
| NovelInsertion | `novins_vaf*` | Add to varforge engine |
| Fusion | `fusion_vaf*` | Requires SV-EXP-005 (junction FASTA) |

### Output per dataset

Same as BENCH-VARFORGE: `calls_discovery.vcf` and `calls_tumour_informed.vcf`.

### Result aggregation

Per-type sensitivity curves and a combined summary table with sensitivity at
key VAFs (0.25%, 0.5%, 1%, 2%, 5%) for each mode.

## Child tasks

| ID | File | Status |
|----|------|--------|
| BENCH-SVN-001 | done/bench_svn_001_invdel_configs.md | done |
| BENCH-SVN-002 | done/bench_svn_002_novel_ins_configs.md | done |
| BENCH-SVN-003 | done/bench_svn_003_fusion_configs.md | done |
| BENCH-SVN-004 | done/bench_svn_004_run_score.md | done |
| BENCH-SVN-005 | done/bench_svn_005_sv_figures.md | done |

## Dependencies

- SV-EXP-001 through SV-EXP-003 (new variant types and classification)
- SV-EXP-005 through SV-EXP-006 (fusion junctions and classification, for BENCH-SVN-003)

## Scope

- `docs/benchmarking/sv_new/configs/` — InvDel, NovelInsertion, Fusion configs
- `docs/benchmarking/sv_new/data/` — truth VCFs
- `docs/benchmarking/sv_new/scripts/` — generation and scoring scripts
- `docs/benchmarking/sv_new/results/` — sensitivity tables and figures

## Out of scope

- Changes to the kam Rust code
- Re-running existing DEL/DUP/INV benchmarks
