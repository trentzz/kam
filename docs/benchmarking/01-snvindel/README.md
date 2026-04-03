# SNV/Indel Benchmark

Evaluation of kam SNV and indel detection on the Twist cfDNA Pan-Cancer
Reference Standard v2 titration series.

## Dataset

- **Sample set**: 24 samples — 3 input masses (5 ng, 15 ng, 30 ng) × 8 VAF
  levels (0%, 0.001%, 0.01%, 0.1%, 0.25%, 0.5%, 1%, 2%).
- **Truth set**: 375 validated somatic variants (205 SNVs, 170 indels) across
  84 cancer genes on GRCh38.
- **Chemistry**: Twist duplex UMI (5 bp UMI, 2 bp skip, inline).
- **Read depths tested**: 250K, 500K, 1M, 2M read pairs per sample.

## Key results (2M read pairs, tumour-informed monitoring mode)

| Metric | Value |
|--------|-------|
| Sensitivity at 2% VAF | 52–61% (varies by input mass) |
| SNV sensitivity at 2% VAF | 69–80% |
| Indel sensitivity at 2% VAF | 31–39% |
| Precision (monitoring mode) | 1.0 at all VAF levels |
| False positives | 0 (monitoring mode) |

See [`docs/project/investigations/snvindel_results_v10_2m_reads_k31.md`](../../project/investigations/snvindel_results_v10_2m_reads_k31.md)
for full results tables and figures.

## Folder structure

```
01-snvindel/
├── README.md              — this file
├── HANDOFF.md             — how to re-run the benchmark on a new machine
├── design/                — benchmark design documents
│   ├── 00_overview.md
│   ├── 01_synthetic_twist.md
│   └── 02_titration.md
├── scripts/               — all scripts for running and analysing the benchmark
│   ├── run_titration_batch.py   — main batch runner (24 samples)
│   ├── run_titration_benchmarks.sh
│   ├── score_variants.py        — TP/FP/FN scoring against truth VCF
│   ├── plot_results.py          — generate sensitivity/precision figures
│   ├── generate_per_sample_tsvs.py
│   └── ...
├── per_sample/            — per-sample TSV outputs and full pipeline outputs
│   ├── README.md
│   └── full_outputs/      — assembly/index/pathfind/call QC JSON + variants TSV
├── sweeps/                — read depth and k-mer size sweep results
├── summary/               — aggregated figures and summary tables
├── results_aggregate/     — aggregate result graphs across all conditions
├── summary_aggregate/     — aggregate CSV result tables
└── logs/                  — batch run logs (gitignored)
```

Large generated files (varforge FASTQs, per-sample BAMs) live in `bigdata/benchmarking/01-snvindel/`.

## Re-running the benchmark

See [`HANDOFF.md`](HANDOFF.md) for full instructions including memory
requirements and FASTQ directory setup.

Quick start:

```bash
cargo build --release
export FASTQ_DIR=/path/to/twist_titration_fastqs
python3 docs/benchmarking/01-snvindel/scripts/run_titration_batch.py \
    --fastq-dir "$FASTQ_DIR" \
    --output-dir docs/benchmarking/01-snvindel/per_sample \
    --depths 2000000
```

## Known limitations

- Indel sensitivity (31–39% at 2% VAF) is substantially lower than SNV
  sensitivity (69–80%). The primary cause is end-anchor coverage gaps: the DFS
  walk requires an exact k-mer match at the reference end, and for indels the
  end anchor is often absent from low-coverage molecules. This is documented in
  `docs/project/investigations/anchor_missing_investigation.md`.
- Sensitivity at 5 ng input is lower than at 15–30 ng (51.7% vs 59–61% at
  2% VAF) due to fewer starting molecules and higher stochastic dropout.
