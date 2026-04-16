# Benchmark 07 — SNV/Indel ML Boost v1

Evaluates ML models on the real Twist UMI titration dataset.
Compares four conditions:

| Condition | ML model | TI |
|-----------|----------|----|
| Discovery (no ML, no TI) | — | — |
| ML v1 | twist-duplex-v1 | — |
| ML v2 | twist-duplex-v2 | — |
| TI only | — | ✓ |

`twist-duplex-v2` is the current bundled model (AUPRC 0.973).
`twist-duplex-v1` is kept for comparison (AUPRC 0.606).

See `discovery-precision-analysis.md` for the full ML model failure analysis
and `report.md` for the TI vs discovery comparison.

## Dataset
- 24 samples: 3 DNA inputs (5 ng, 15 ng, 30 ng) × 8 VAF levels (0 %–2 %)
- FASTQs: `/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs/`
- Targets: `../01-snvindel/scripts/targets_100bp.fa` (375 loci, 100 bp)
- Truth VCF: `../01-snvindel/scripts/truth_variants.vcf` (375 variants: 205 SNV, 170 indel)

## Running

Install psutil if needed: `pip install psutil`

### Discovery (no ML, no TI)
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --output titration_2mreads_disc.tsv
```

### Discovery + ML v1
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --ml-model twist-duplex-v1 \
  --output titration_2mreads_disc_ml.tsv
```

### Discovery + ML v2
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --ml-model twist-duplex-v2 \
  --output titration_2mreads_disc_ml_v2.tsv
```

### TI only (no ML)
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --target-variants docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
  --output titration_2mreads_ti.tsv
```

### TI + duplex confirmation filter
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --target-variants docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
  --min-alt-duplex 1 \
  --output titration_2mreads_ti_minduplex1.tsv
```

## Outputs
- `results/titration_2mreads_disc.tsv` — aggregate metrics, discovery run
- `results/titration_2mreads_disc_ml.tsv` — aggregate metrics, discovery + ML v1
- `results/titration_2mreads_disc_ml_v2.tsv` — aggregate metrics, discovery + ML v2
- `results/titration_2mreads_ti.tsv` — aggregate metrics, TI-only run
- `results/titration_2mreads_ti_minduplex1.tsv` — aggregate metrics, TI + `--min-alt-duplex 1`
- `results/per_sample/{sample}.targets.tsv` — per-sample, per-variant detail
- `results/vcfs/{sample}.monitoring.vcf` — per-sample called variants
- `report.md` — analysis and conclusions (TI vs discovery)
- `duplex-filter-analysis.md` — sensitivity cost of `--min-alt-duplex 1` vs TI-only
- `discovery-precision-analysis.md` — discovery mode precision study: 5 filter conditions, ML v1/v2 failure analysis
