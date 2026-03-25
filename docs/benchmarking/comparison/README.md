# Thesis Alignment Baseline Comparison

This directory holds the standardised baseline derived from the thesis
alignment-based variant detection results, used to compare against kam output.

## Files

| File | Description |
|------|-------------|
| `alignment_baseline.csv` | Per-(sample, variant) detection results from the alignment pipeline |
| `truth_variants.csv` | Truth variant set (375 variants) |
| `scripts/parse_thesis_results.py` | Script that produced the two CSV files above |

## Data Provenance

**Source:** `/mnt/tzeng-local/tzeng-thesis/` (local thesis data, not in repo)

**Truth set:** `titration.probes.QC.pass.tsv`
375 variants from the TWIST Standard v2 panel: 205 SNVs and 170 indels.
These passed a QC filter applied in the thesis.

**Alignment pipeline:** HUMID v2 deduplication + RASCALL v3 calling.
This is the primary result reported in the thesis for the titration experiment.

SNV results come from:
```
benchmark-dataframes/dedup-v2-generate-v3/
```

Indel results come from:
```
benchmark-dataframes-indel/dedup-v2/
```

Each directory contains 24 per-sample TSV files: 3 input amounts (5 ng, 15 ng,
30 ng) x 8 VAF levels (0%, 0.001%, 0.01%, 0.1%, 0.25%, 0.5%, 1%, 2%).

## alignment_baseline.csv columns

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | string | Source filename stem |
| `input_ng` | int | Input DNA amount in nanograms |
| `vaf_level` | string | VAF label from filename (e.g. `2pc`, `0p1pc`) |
| `chromosome` | string | Chromosome (e.g. `chr1`) |
| `position` | int | 1-based genomic position |
| `ref` | string | Reference allele |
| `alt` | string | Alternate allele |
| `variant_type` | string | `SNP` or `INDEL` |
| `alignment_detected` | bool | Whether RASCALL reported detection |
| `alignment_vaf` | float or empty | RASCALL-reported VAF (`RASCALL rVAF`); empty if not detected |
| `alignment_depth` | int or empty | Minimum coverage at locus (`RASCALL Min_coverage`); empty if not detected |

VAF label to numeric fraction mapping:

| Label | Fraction |
|-------|----------|
| `0pc` | 0.0 |
| `0p001pc` | 0.00001 |
| `0p01pc` | 0.0001 |
| `0p1pc` | 0.001 |
| `0p25pc` | 0.0025 |
| `0p5pc` | 0.005 |
| `1pc` | 0.01 |
| `2pc` | 0.02 |

## truth_variants.csv columns

| Column | Type | Description |
|--------|------|-------------|
| `chromosome` | string | Chromosome |
| `position` | int | 1-based genomic position |
| `ref` | string | Reference allele |
| `alt` | string | Alternate allele |
| `variant_type` | string | `SNP` or `INDEL` |

## Regenerating the files

```bash
python3 docs/benchmarking/comparison/scripts/parse_thesis_results.py
```

The script accepts `--data-root` and `--out-dir` to override the default paths.
Run `python3 ... --help` for details.

## Notes

- The alignment pipeline did not report per-allele depth separately, only minimum
  coverage across the 70 bp target sequence. Use `alignment_depth` as an
  approximation of locus depth, not a precise AD/DP ratio.
- The 0% VAF samples (negative controls) are included. Any detections there are
  false positives in the alignment pipeline.
- Ref and Alt for each variant are joined from the truth set by (chromosome,
  position). The per-sample TSV files record only chromosome and position.
