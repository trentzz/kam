# Benchmarking Handoff

Instructions for running the titration benchmarks on a machine with sufficient RAM.

## Memory requirements

The titration batch runs 24 samples (Twist cfDNA Standard v2, 15ng and 30ng at 8 VAF levels each).
At 1M read pairs per sample:

- Most samples: 1.2–1.6 GB peak RSS.
- Low-VAF and low-concentration samples (especially 15ng): up to 5+ GB peak RSS due to
  de Bruijn graph path explosion at high read depth.

**Minimum RAM recommended: 16 GB.** 32 GB preferred to run all 24 samples without skips.

The OOM kill threshold is set at `PEAK_RSS_LIMIT_MB = 5000` (5 GB) in
`benchmarking/scripts/run_titration_batch.py`. Raise this if the machine has more RAM.

## Setup

1. **Build the binary:**

   ```bash
   cargo build --release
   ```

   The binary lands at `target/release/kam`.

2. **Set the FASTQ directory:**

   ```bash
   export KAM_FASTQ_DIR=/path/to/titration-nondedup/fastqs
   ```

   The directory must contain files matching:
   `TWIST_STDV2_{ng}ng_VAF_{vaf}pc_*_R{1,2}.fastq.gz`

   where `{ng}` is `15ng` or `30ng`, and `{vaf}` uses `p` as the decimal separator
   (e.g. `0p1` for 0.1%, `2` for 2%).

3. **Verify truth variants file exists:**

   `benchmarking/scripts/truth_variants.vcf` should already exist (checked in to the repo).
   It contains 375 variants (205 SNV, 170 indel). If it is missing or needs regeneration:

   ```bash
   export KAM_PROBES_TSV=/path/to/titration.probes.QC.pass.tsv
   python3 benchmarking/scripts/tsv_to_vcf.py
   ```

4. **Install Python dependencies:**

   ```bash
   pip install psutil
   ```

## Running the batch

```bash
python3 benchmarking/scripts/run_titration_batch.py \
  --fastq-dir /path/to/titration-nondedup/fastqs \
  2>&1 | tee benchmarking/results/titration_batch.log
```

All paths are overridable via CLI flags:

```
--fastq-dir PATH       Directory with TWIST_STDV2_* FASTQ pairs (required on new machine)
--truth-vcf PATH       Truth VCF (default: benchmarking/scripts/truth_variants.vcf)
--targets PATH         Target FASTA (default: benchmarking/scripts/targets_100bp.fa)
--kam-binary PATH      kam binary (default: target/release/kam)
--results-dir PATH     Output directory (default: benchmarking/results/tables)
--reads N              Read pairs per sample (default: 1000000)
--rss-limit-mb N       OOM kill threshold in MB (default: 5000)
```

Example with raised limits for a 32 GB machine:

```bash
python3 benchmarking/scripts/run_titration_batch.py \
  --fastq-dir /data/titration/fastqs \
  --reads 2000000 \
  --rss-limit-mb 24000 \
  2>&1 | tee benchmarking/results/titration_batch.log
```

Results are written to `benchmarking/results/tables/titration_results.tsv` after each sample
completes, so partial runs are safe to resume (just re-run; completed entries are not
re-processed in the current design, though the file is overwritten on restart).

## Output columns

The TSV has these column groups:

| Group | Columns |
|---|---|
| Identity | `sample`, `ng`, `vaf` |
| Assembly | `molecules`, `duplex`, `duplex_pct` |
| Called | `variants_called` |
| Overall metrics | `tp fp fn tn sensitivity precision specificity fpr fdr fnr f1` |
| SNV metrics | `snv_tp snv_fp snv_fn snv_sensitivity snv_precision snv_specificity snv_fpr` |
| Indel metrics | `indel_tp indel_fp indel_fn indel_sensitivity indel_precision indel_specificity indel_fpr` |
| Stage timing | `t_assemble_ms t_index_ms t_pathfind_ms t_call_ms t_output_ms` |
| Stage peak RSS | `rss_assemble_mb rss_index_mb rss_pathfind_mb rss_call_mb rss_output_mb` |
| Stage peak CPU | `cpu_assemble_pct cpu_index_pct cpu_pathfind_pct cpu_call_pct cpu_output_pct` |
| Runtime | `peak_rss_mb wall_time_s exit_code` |

Skipped samples (OOM kill) have empty metric columns and `exit_code = -1`.

## Tuning for RAM

If you still see OOM kills, reduce `READS_PER_SAMPLE` in the script:

```python
READS_PER_SAMPLE = 500_000   # try 500K if 1M causes OOM
PEAK_RSS_LIMIT_MB = 16_000   # raise limit if machine has >16 GB free
```

At 200K reads per sample, no sample exceeds 760 MB. At 1M reads, most samples use
~1.6 GB but some 15ng samples use >5 GB.

## After the batch

Once `titration_results.tsv` is complete:

1. Copy it back to the repo at `benchmarking/results/tables/titration_results.tsv`.
2. Run the plot script to regenerate figures:

   ```bash
   python3 benchmarking/scripts/plot_results.py
   ```

   Figures land in `docs/paper/figures/`. Rebuild the paper after updating figures:

   ```bash
   cd docs/paper
   for fmt in main main_normal main_normal_2col main_2p; do
     pdflatex -interaction=nonstopmode "$fmt.tex"
     biber "$fmt"
     pdflatex -interaction=nonstopmode "$fmt.tex"
     pdflatex -interaction=nonstopmode "$fmt.tex"
   done
   cp main.pdf pdfs/kam_v0.1_big.pdf
   cp main_normal.pdf pdfs/kam_v0.1_normal.pdf
   cp main_normal_2col.pdf pdfs/kam_v0.1_normal_2col.pdf
   cp main_2p.pdf pdfs/kam_v0.1_2p.pdf
   ```

3. Commit and push.
