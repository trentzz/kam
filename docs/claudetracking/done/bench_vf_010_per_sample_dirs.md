# BENCH-VF-010: Per-sample benchmark result directories

**Epic**: BENCH-VARFORGE (docs/claudetracking/overallplans/BENCH-VARFORGE.md)
**Priority**: high
**Depends on**: BENCH-VF-004
**Status**: todo

## Goal

Reorganise benchmark results into per-sample directories, each containing the
varforge config, truth TSV, kam discovery TSV, and kam tumour-informed TSV.
Push these result files to GitHub (remove them from .gitignore). Done looks
like: every benchmarked sample has a `samples/{type}_vaf{tag}_{rep}/`
directory committed to the repo, with four files per sample.

## Steps

1. Update `.gitignore`:
   - Remove the two blanket ignores for `docs/benchmarking/*/results/`.
   - Add `docs/benchmarking/*/results/sim_*/` (FASTQ simulation outputs, large).
   - Add `docs/benchmarking/*/results/kam_*/` (old flat VCF result dirs).
   - The new `samples/` directories are not ignored.

2. Write `docs/benchmarking/build_sample_dirs.py`:
   - For each suite (`sv`, `snvindel`):
     - For each sample (all type/tag/rep combinations from the suite scripts):
       - Create `docs/benchmarking/{suite}/samples/{type}_vaf{tag}_{rep}/`
       - Copy `configs/{type}_vaf{tag}_{rep}.yaml` → `config.yaml`
       - Write `varforge_cmd.txt`:
         `varforge run docs/benchmarking/{suite}/configs/{type}_vaf{tag}_{rep}.yaml`
       - Convert truth VCF from `data/` → `truth.tsv`
         (columns: chrom, pos, ref, alt, vaf, type)
       - Convert `results/kam_{type}_vaf{tag}_{rep}/calls_discovery.vcf` →
         `discovery.tsv`
         (columns: chrom, pos, ref, alt, filter, vaf, vaf_lo, vaf_hi,
          nref, nalt, ndupalt, nsimalt, sbp, conf)
       - Convert `results/kam_{type}_vaf{tag}_{rep}/calls_tumour_informed.vcf`
         → `tumour_informed.tsv`

3. Run the script:
   ```bash
   python3 docs/benchmarking/build_sample_dirs.py
   ```
   Confirm that sample directories are created for all expected samples.

4. Update `run_sv_suite.sh` and `run_snv_indel_suite.sh` to also write to the
   `samples/` directories going forward (so future re-runs stay in sync).

5. Commit the new `samples/` directories and updated `.gitignore`.

## Notes

- Truth TSV source: `docs/benchmarking/{suite}/data/truth_{type}_{tag}_{rep}.vcf`
  for snvindel; `docs/benchmarking/sv/data/truth_{svtype}_{tag}_{rep}.vcf`
  for sv suite.
- For samples where only VCF exists (all pre-TSV runs), convert VCF to TSV in
  the script using VCF INFO field parsing.
- FASTQ files (`sim_*` directories) stay gitignored — they are large (~500MB
  each) and reproducible from the config + varforge.
- The `kam_*` flat directories can stay for now; they are ignored by git.
