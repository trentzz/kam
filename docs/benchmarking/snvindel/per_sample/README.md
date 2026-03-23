# Per-sample variant TSVs: SNV and indel titration

This directory contains raw kam output TSVs for all 24 titration samples, in both
discovery and tumour-informed modes.

## File naming

```
{sample}_discovery.tsv    — discovery mode (no --target-variants)
{sample}_monitoring.tsv   — tumour-informed mode (--target-variants truth_variants.vcf)
```

## Generating the files

Run from the repository root with FASTQ data available:

```bash
KAM_FASTQ_DIR=/data/titration-nondedup/fastqs \
python docs/benchmarking/scripts/generate_per_sample_tsvs.py
```

Options:
- `--reads N`: read pairs per sample (default 2,000,000 — matches the canonical v10 run)
- `--fastq-dir PATH`: override FASTQ directory
- `--output-dir PATH`: override output directory (default: this directory)

The script skips samples whose TSVs already exist, so it is safe to rerun after partial
completion.

## Contents of each TSV

All calls (PASS and filtered), one row per candidate variant path per target. Columns:

| Column | Description |
|--------|-------------|
| `target_id` | Target region (e.g. `chr17:7674220-7674320`) |
| `variant_type` | SNV / MNV / Insertion / Deletion / Complex |
| `ref_seq` | Full reference path sequence |
| `alt_seq` | Full alt path sequence |
| `vaf` | VAF point estimate |
| `vaf_ci_low` | 95% CI lower bound |
| `vaf_ci_high` | 95% CI upper bound |
| `n_molecules_ref` | Reference-supporting molecules |
| `n_molecules_alt` | Alt-supporting molecules |
| `n_duplex_alt` | Alt-specific duplex molecules |
| `n_simplex_alt` | Alt simplex molecules |
| `strand_bias_p` | Fisher strand bias p-value |
| `confidence` | Posterior confidence |
| `filter` | PASS / LowConfidence / StrandBias / LowDuplex / HighVaf / NotTargeted |

## Discovery vs tumour-informed

**Discovery**: `filter=PASS` means the call passed all quality filters. At 2M reads,
these samples produce PASS calls only for variants that genuinely exist in the sample.
Background biology (germline, clonal haematopoiesis) does not accumulate enough alt
molecules at 2M reads to pass the default filters.

**Monitoring**: `filter=PASS` means the call passed quality filters AND matches an entry
in `truth_variants.vcf`. All other quality-passing calls are relabelled `NotTargeted`.
At full read depth (not 2M reads), the 15ng 0% VAF sample produces 62 background biology
PASS calls in discovery mode; all 62 become NotTargeted in tumour-informed mode.

## Sample list (24 samples)

| Sample | DNA | VAF |
|--------|-----|-----|
| Sample_5ng_VAF_0pc | 5ng | 0% (negative control) |
| Sample_5ng_VAF_0p001pc | 5ng | 0.001% |
| Sample_5ng_VAF_0p01pc | 5ng | 0.01% |
| Sample_5ng_VAF_0p1pc | 5ng | 0.1% |
| Sample_5ng_VAF_0p25pc | 5ng | 0.25% |
| Sample_5ng_VAF_0p5pc | 5ng | 0.5% |
| Sample_5ng_VAF_1pc | 5ng | 1% |
| Sample_5ng_VAF_2pc | 5ng | 2% |
| Sample_15ng_VAF_0pc | 15ng | 0% (negative control) |
| Sample_15ng_VAF_0p001pc | 15ng | 0.001% |
| Sample_15ng_VAF_0p01pc | 15ng | 0.01% |
| Sample_15ng_VAF_0p1pc | 15ng | 0.1% |
| Sample_15ng_VAF_0p25pc | 15ng | 0.25% |
| Sample_15ng_VAF_0p5pc | 15ng | 0.5% |
| Sample_15ng_VAF_1pc | 15ng | 1% |
| Sample_15ng_VAF_2pc | 15ng | 2% |
| Sample_30ng_VAF_0pc | 30ng | 0% (negative control) |
| Sample_30ng_VAF_0p001pc | 30ng | 0.001% |
| Sample_30ng_VAF_0p01pc | 30ng | 0.01% |
| Sample_30ng_VAF_0p1pc | 30ng | 0.1% |
| Sample_30ng_VAF_0p25pc | 30ng | 0.25% |
| Sample_30ng_VAF_0p5pc | 30ng | 0.5% |
| Sample_30ng_VAF_1pc | 30ng | 1% |
| Sample_30ng_VAF_2pc | 30ng | 2% |
