# BENCH-VARFORGE: Extensive Varforge Benchmarking Suite

**Status**: active
**Priority**: high
**Branch**: epic/BENCH-VARFORGE

## Goal

Build a comprehensive synthetic benchmarking suite with ~50 configs per
variant class. The suite covers SNVs, indels, and SVs across a dense VAF
sweep and multiple replicates. Every dataset is run in both discovery and
tumour-informed modes. Results are recorded per-dataset and aggregated into
sensitivity curves and summary tables.

## Motivation

The current synthetic benchmark has 4 VAF levels per type (0.5%, 1%, 2%,
5%). This is enough to confirm detection works but not enough to characterise
the sensitivity curve, find detection thresholds, or compare discovery against
tumour-informed mode across a meaningful range. Fifty configs per type gives
dense enough coverage to draw smooth curves and identify the VAF at which
sensitivity drops below clinically relevant thresholds (e.g. 50%, 90%).

## Design

### VAF levels

50 configs per type = 25 VAF levels × 2 replicates (different seeds).

VAF levels (purity = VAF × 2 for diploid):

```
0.05%, 0.075%, 0.10%, 0.125%, 0.15%, 0.20%, 0.25%, 0.30%, 0.40%, 0.50%,
0.60%, 0.75%, 1.00%, 1.25%, 1.50%, 2.00%, 2.50%, 3.00%, 4.00%, 5.00%,
6.00%, 7.00%, 8.00%, 9.00%, 10.00%
```

Two replicates at each level (seeds A and B) to quantify stochastic
variance at low VAF.

### Variant classes

| Class | Dataset | Truth VCF | Config prefix |
|-------|---------|-----------|---------------|
| SNV | 5 SNVs on 2000 bp ref | truth_snvs_vafXXX_{a,b}.vcf | snv_vafXXX_{a,b}.yaml |
| Indel | 5 indels on 2000 bp ref | truth_indels_vafXXX_{a,b}.vcf | indel_vafXXX_{a,b}.yaml |
| SV (DEL+DUP+INV) | 3 SVs on 2000 bp ref | truth_svs_vafXXX_{a,b}.vcf | sv_vafXXX_{a,b}.yaml |
| INS | 1 insertion on 2000 bp ref | truth_ins_vafXXX_{a,b}.vcf | ins_vafXXX_{a,b}.yaml |
| INVDEL | 1 INVDEL on 2000 bp ref | truth_invdel_vafXXX_{a,b}.vcf | invdel_vafXXX_{a,b}.yaml |

SNV, indel, and SV (DEL+DUP+INV) get the full 50-config sweep. INS and
INVDEL get 50 configs each as well. Large DEL is excluded until the varforge
engine.rs panic is resolved.

### Outputs per dataset

Every dataset must produce both:
- `calls_discovery.vcf` — kam run without `--target-variants`
- `calls_tumour_informed.vcf` — kam run with `--target-variants`

See `docs/claudeguide/benchmarking.md` for the full output naming convention.

### Result aggregation

For each variant class, produce:
1. A per-dataset sensitivity table (one row per dataset: VAF, replicate,
   discovery sensitivity, tumour-informed sensitivity).
2. A sensitivity curve plot: x = VAF, y = sensitivity, two lines (discovery
   and tumour-informed), with replicate scatter.
3. A summary table: sensitivity at key VAFs (0.1%, 0.25%, 0.5%, 1%, 2%, 5%)
   for each mode.

## Child tasks

| ID | File | Status |
|----|------|--------|
| BENCH-VF-001 | todo/bench_vf_001_snv_suite.md | todo |
| BENCH-VF-002 | todo/bench_vf_002_indel_suite.md | todo |
| BENCH-VF-003 | todo/bench_vf_003_sv_suite.md | todo |
| BENCH-VF-004 | todo/bench_vf_004_run_all.md | todo |
| BENCH-VF-005 | todo/bench_vf_005_score_aggregate.md | todo |
| BENCH-VF-006 | todo/bench_vf_006_figures.md | todo |

## Scope

- `docs/benchmarking/snvindel/configs/` — SNV and indel configs
- `docs/benchmarking/snvindel/data/` — SNV and indel truth VCFs
- `docs/benchmarking/sv/configs/` — SV configs
- `docs/benchmarking/sv/data/` — SV truth VCFs
- New scoring/aggregation scripts in each benchmark's `scripts/` directory
- New figure output in each benchmark's `results/` directory

## Out of scope

- Changes to the kam Rust code
- Large DEL benchmark (blocked on varforge bug)
- Combined SNV+indel configs (existing suite covers this)
