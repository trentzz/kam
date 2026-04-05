# BENCH-SVN-004: Run All New SV Benchmarks and Score

**Epic**: BENCH-SV-NEW (docs/claudetracking/overallplans/BENCH-SV-NEW.md)
**Priority**: high
**Depends on**: bench_svn_001_invdel_configs.md, bench_svn_002_novel_ins_configs.md, bench_svn_003_fusion_configs.md
**Status**: todo

## Goal

Run kam on all new SV benchmark FASTQ pairs (InvDel, NovelInsertion, Fusion)
in both discovery and tumour-informed modes, then score each run against its
truth VCF. Done looks like: per-dataset sensitivity TSVs exist for all 150
datasets (50 per type), and a summary TSV aggregates sensitivity at key VAF
levels.

## Steps

1. Write `docs/benchmarking/sv_new/scripts/run_sv_new_suite.sh`:
   - Loop over all 50 InvDel, 50 NovelInsertion, and 50 Fusion configs.
   - For each config, run:
     ```bash
     kam run --r1 sim_*/R1.fastq.gz --r2 sim_*/R2.fastq.gz \
         --targets ... --config ... \
         --output results/{type}_vaf{tag}_{rep}/calls_discovery.vcf
     kam run --r1 ... --r2 ... --targets ... --target-variants truth.vcf \
         --output results/{type}_vaf{tag}_{rep}/calls_tumour_informed.vcf
     ```
   - Parallelise with `xargs -P 8` or GNU parallel.

2. Write `docs/benchmarking/sv_new/scripts/score_sv_new_suite.py`:
   - For each dataset, compare calls VCF to truth VCF.
   - A variant is a true positive if chrom, pos, ref, alt match (or within
     ±5 bp for SV breakpoints).
   - Compute sensitivity = TP / (TP + FN) per dataset.
   - Output: `results/sensitivity_invdel.tsv`, `sensitivity_novins.tsv`,
     `sensitivity_fusion.tsv`.
   - Columns: vaf_tag, rep, discovery_sensitivity, ti_sensitivity.

3. Run the suite and scoring scripts. Verify all 150 datasets produce non-empty
   result files.

4. Produce a summary table at key VAF levels (0.25%, 0.5%, 1%, 2%, 5%).

## Notes

- The `±5 bp` breakpoint tolerance for SV matching is standard for SV
  benchmarking (see SURVIVOR and TRUVARI conventions).
- If any SV type has zero detections across all VAF levels, flag as a
  blocking issue before proceeding to figures.
- Both discovery and tumour-informed modes must be run and scored. Results
  from tumour-informed mode are expected to be the primary headline metric.
