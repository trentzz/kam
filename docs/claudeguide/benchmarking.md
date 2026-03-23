# Benchmarking Guide

## Overview

kam has two synthetic benchmark suites. Both use varforge to simulate Twist
duplex UMI FASTQ data from a truth VCF, then run kam and score the output.

| Suite | Location | What it tests |
|-------|----------|---------------|
| SNV/indel | `docs/benchmarking/snvindel/` | SNV and indel sensitivity across VAF levels |
| SV | `docs/benchmarking/sv/` | Structural variant detection across SV types and VAF levels |

---

## Varforge datasets

### SNV/indel (`docs/benchmarking/snvindel/`)

Configs in `configs/`, truth VCFs in `data/`, generated FASTQs in `results/`
(gitignored).

| Dataset | VAF levels | Variants |
|---------|-----------|----------|
| SNV only | 0.5%, 1%, 2%, 5% | 5 SNVs |
| Indel only | 0.5%, 1%, 2%, 5% | 5 indels (2–5 bp) |
| Combined | 0.5%, 1%, 2%, 5% | 10 (5 SNV + 5 indel) |

Config naming: `snv_vaf010.yaml`, `indel_vaf005.yaml`, `combined_vaf020.yaml`
(three-digit suffix: 005 = 0.5%, 010 = 1%, 020 = 2%, 050 = 5%).

To regenerate a dataset:
```
varforge simulate --config docs/benchmarking/snvindel/configs/snv_vaf010.yaml
```

### SV (`docs/benchmarking/sv/`)

Configs in `configs/`, truth VCFs in `data/`, generated FASTQs in `results/`
(gitignored).

| Dataset | VAF levels | Notes |
|---------|-----------|-------|
| DEL + DUP + INV | 0.5%, 1%, 2%, 5% | Standard SV types |
| INS (99 bp) | 0.5%, 1%, 2%, 5% | No reference anchor at alt end |
| Large DEL (600 bp) | 0.5%, 1%, 2%, 5% | Configs committed; varforge panics (engine.rs bug) |
| INVDEL | 0.5%, 1%, 2%, 5% | Complex rearrangement |
| DEL + DUP + INV | 0.25% | Reserved; varforge panics |

Config naming: `sim_vaf010.yaml`, `sim_ins_vaf005.yaml`,
`sim_largedel_vaf020.yaml`.

Purity formula (diploid): `purity = VAF × 2`. Seeds follow the pattern:
standard SV 1xxx, INS 30xx, large DEL 31xx, INVDEL 32xx.

---

## Running kam on a dataset

Every benchmarking run must produce **both** outputs:

1. **Discovery mode** — no prior knowledge of variants.
2. **Monitoring (tumour-informed) mode** — filters to known variants via
   `--target-variants`.

Name the outputs consistently:
```
kam_<dataset>/calls_discovery.vcf
kam_<dataset>/calls_monitoring.vcf
```

Example for an SV dataset at 1% VAF:
```bash
# Discovery
kam run \
  --r1 docs/benchmarking/sv/results/sim_vaf010/SV_VAF010_R1.fastq.gz \
  --r2 docs/benchmarking/sv/results/sim_vaf010/SV_VAF010_R2.fastq.gz \
  --targets docs/benchmarking/sv/data/sv_targets.fa \
  --sv-junctions docs/benchmarking/sv/data/sv_junctions.fa \
  --output-dir docs/benchmarking/sv/results/kam_vaf010 \
  --output-format vcf

# Monitoring
kam run \
  --r1 docs/benchmarking/sv/results/sim_vaf010/SV_VAF010_R1.fastq.gz \
  --r2 docs/benchmarking/sv/results/sim_vaf010/SV_VAF010_R2.fastq.gz \
  --targets docs/benchmarking/sv/data/sv_targets.fa \
  --sv-junctions docs/benchmarking/sv/data/sv_junctions.fa \
  --output-dir docs/benchmarking/sv/results/kam_vaf010_monitoring \
  --target-variants docs/benchmarking/sv/data/truth_svs_vaf010.vcf \
  --output-format vcf
```

Never report results from only one mode. The discovery vs tumour-informed comparison
is a primary result.

---

## Scoring

Use `docs/benchmarking/snvindel/scripts/score_variants.py` to compute
TP/FP/FN, sensitivity, precision, and F1. It classifies variants as SNV or
INDEL and reports both separately.

```
python docs/benchmarking/snvindel/scripts/score_variants.py \
  --truth <truth_vcf> \
  --calls <kam_output_vcf>
```

Matching is by (REF, ALT) sequence only — no position matching.

---

## Adding new datasets

1. Add a truth VCF to `data/` with the correct VAF in the INFO field.
2. Add a varforge config to `configs/` pointing to the truth VCF. Follow the
   naming convention and seed pattern above.
3. Run varforge to generate the FASTQs.
4. Add the results directory to `.gitignore` if it is not already covered by
   the existing wildcard entries.
5. Commit the config and truth VCF. Do not commit the FASTQs.
