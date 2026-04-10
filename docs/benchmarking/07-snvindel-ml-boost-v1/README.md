# Benchmark 07 — SNV/Indel ML Boost v1

Evaluates the `twist-duplex-v1` ML model on the real Twist UMI titration dataset.
Runs in tumour-informed (TI) mode and compares three conditions:

| Condition | ML | TI |
|-----------|----|----|
| Baseline (no ML, no TI) | — | — |
| TI only | — | ✓ |
| ML + TI | twist-duplex-v1 | ✓ |

## Dataset
- 24 samples: 3 DNA inputs (5 ng, 15 ng, 30 ng) × 8 VAF levels (0 %–2 %)
- FASTQs: `/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs/`
- Targets: `../01-snvindel/scripts/targets_100bp.fa` (375 loci, 100 bp)
- Truth VCF: `../01-snvindel/scripts/truth_variants.vcf` (375 variants: 205 SNV, 170 indel)

## Running

Install psutil if needed: `pip install psutil`

### TI only (no ML)
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --target-variants docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
  --output titration_2mreads_ti.tsv
```

### ML + TI
```bash
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --target-variants docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
  --ml-model twist-duplex-v1 \
  --save-vcfs docs/benchmarking/07-snvindel-ml-boost-v1/results/vcfs \
  --per-sample-dir docs/benchmarking/07-snvindel-ml-boost-v1/results/per_sample \
  --output titration_2mreads_ml_twist_duplex_ti.tsv
```

## Outputs
- `results/titration_2mreads_ti.tsv` — aggregate metrics, TI-only run
- `results/titration_2mreads_ml_twist_duplex_ti.tsv` — aggregate metrics, ML+TI run
- `results/per_sample/{sample}.targets.tsv` — per-sample, per-variant detail
- `results/vcfs/{sample}.monitoring.vcf` — per-sample called variants
- `report.md` — analysis and conclusions
