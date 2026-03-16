# Output Format Specifications

## Annotated FASTQ Headers (After Molecule Assembly)

kam-assemble outputs consensus reads as FASTQ with enriched headers encoding molecule provenance. This format is compatible with any tool that reads FASTQ (extra tags are ignored) while providing full metadata for downstream kam stages.

### Tag Definitions

Following SAM tag naming conventions for compatibility if later converted to BAM:

| Tag | Type | Description | Example |
|-----|------|-------------|---------|
| `MI` | String | Molecule ID (hash of canonical UMI pair) | `MI:Z:a3f8b2c1` |
| `RX` | String | Raw UMI pair, dash-separated | `RX:Z:ACGTA-TGCAT` |
| `FS` | Integer | Total family size (all reads contributing) | `FS:i:12` |
| `DS` | Integer | Duplex support count (reads on both strands) | `DS:i:8` |
| `SF` | Integer | Simplex forward count | `SF:i:3` |
| `SR` | Integer | Simplex reverse count | `SR:i:1` |
| `FT` | String | Family type: `DUPLEX`, `SIMPLEX_FWD`, `SIMPLEX_REV`, `SINGLETON` | `FT:Z:DUPLEX` |
| `CQ` | Float | Mean consensus quality (probability of error) | `CQ:f:0.0003` |
| `CP` | Float | Estimated UMI collision probability for this family | `CP:f:0.012` |
| `SK` | String | Skip bases observed (for Twist QC) | `SK:Z:AG` |

### Header Format
```
@READ_001 MI:Z:a3f8b2c1 RX:Z:ACGTA-TGCAT FS:i:12 DS:i:8 SF:i:3 SR:i:1 FT:Z:DUPLEX CQ:f:0.0003 CP:f:0.012 SK:Z:AG
ACGTACGTACGTACGTACGT...
+
IIIIIIIIIIIIIIIIIII...
```

### Quality Line
The quality line encodes the **consensus quality** per base, not the raw read quality. Phred-encoded as standard FASTQ. For duplex consensus positions, quality reflects the duplex error suppression (~Q70 theoretical max, capped at Q60 in practice since that exceeds Illumina's measurement range).

### Compatibility Notes
- Tags use SAM-style `XX:T:value` format for easy BAM conversion
- `MI` and `RX` are standard SAM tags already used by fgbio
- `FS`, `DS`, `SF`, `SR` are custom but follow SAM naming conventions (uppercase two-letter)
- Any tool that ignores header annotations treats this as normal FASTQ
- `RX` format matches fgbio convention (dash-separated UMI pair)

## QC JSON (Per-Stage Output)

Each kam stage outputs a structured JSON file for Nextflow pipeline validation.

### Molecule Assembly QC
```json
{
  "stage": "molecule_assembly",
  "version": "0.1.0",
  "sample_id": "PATIENT_001",
  "timestamp": "2026-03-14T02:17:43Z",
  "git_sha": "a3f8b2c1",
  "chemistry": "twist-umi-duplex",
  "metrics": {
    "n_input_read_pairs": 4821044,
    "n_molecules": 612847,
    "n_duplex_families": 489478,
    "n_simplex_fwd_families": 73927,
    "n_simplex_rev_families": 49442,
    "n_singletons": 0,
    "duplex_fraction": 0.798,
    "mean_family_size": 7.2,
    "median_family_size": 6,
    "max_family_size": 47,
    "low_qual_umi_discarded": 14463,
    "low_qual_umi_fraction": 0.003,
    "estimated_umi_collision_rate": 0.012,
    "n_umi_groups_before_clustering": 98234,
    "n_umi_groups_after_clustering": 89012,
    "skip_base_consistency": 0.997,
    "mean_consensus_quality": 0.0004
  },
  "sanity_checks": {
    "duplex_fraction_above_50pct": true,
    "low_qual_umi_below_5pct": true,
    "collision_rate_below_5pct": true,
    "sufficient_molecules": true,
    "skip_bases_consistent": true
  },
  "warnings": [],
  "passed": true
}
```

### K-mer Indexing QC
```json
{
  "stage": "kmer_indexing",
  "metrics": {
    "n_target_kmers_in_allowlist": 31248,
    "n_kmers_indexed": 28901,
    "n_kmers_with_duplex_support": 27234,
    "mean_molecule_depth_per_kmer": 891.2,
    "n_zero_coverage_target_kmers": 47,
    "memory_usage_mb": 412
  }
}
```

### Graph Walking QC
```json
{
  "stage": "graph_walking",
  "metrics": {
    "n_targets_queried": 1247,
    "n_targets_with_variant_paths": 23,
    "n_anchors_flagged_non_unique": 3,
    "mean_paths_per_target": 1.8,
    "n_targets_with_boundary_issues": 0
  }
}
```

### Variant Calling QC
```json
{
  "stage": "variant_calling",
  "metrics": {
    "n_variants_called": 7,
    "n_variants_duplex_supported": 6,
    "n_variants_simplex_only": 1,
    "n_variants_with_strand_bias": 0,
    "mean_vaf": 0.023,
    "min_vaf": 0.001,
    "max_vaf": 0.142
  }
}
```

## Variant Call Output

### TSV (Primary — km-compatible)

Tab-separated, one variant per line. Headers:

```
target_id    variant_type    ref_path    alt_path    vaf    vaf_ci_low    vaf_ci_high    n_molecules_ref    n_molecules_alt    n_duplex_alt    n_simplex_alt    strand_bias_p    confidence    filter
```

| Column | Description |
|--------|-------------|
| `target_id` | Target sequence identifier |
| `variant_type` | SNV, INS, DEL, MNV, COMPLEX |
| `ref_path` | Reference k-mer path sequence |
| `alt_path` | Variant k-mer path sequence |
| `vaf` | Point estimate of variant allele frequency |
| `vaf_ci_low` | Lower bound of 95% credible interval |
| `vaf_ci_high` | Upper bound of 95% credible interval |
| `n_molecules_ref` | Molecules supporting reference path |
| `n_molecules_alt` | Molecules supporting variant path |
| `n_duplex_alt` | Of alt molecules, how many are duplex-confirmed |
| `n_simplex_alt` | Of alt molecules, simplex-only count |
| `strand_bias_p` | Fisher's exact test p-value for strand bias |
| `confidence` | Posterior probability that variant is real (0-1) |
| `filter` | PASS, STRAND_BIAS, LOW_CONFIDENCE, LOW_DUPLEX, COLLISION_RISK |

### VCF (Secondary — Standard Format)

Standard VCF 4.3 with custom INFO and FORMAT fields:

```vcf
##INFO=<ID=NMOL,Number=1,Type=Integer,Description="Number of molecules supporting variant">
##INFO=<ID=NDUP,Number=1,Type=Integer,Description="Number of duplex molecules supporting variant">
##INFO=<ID=SBP,Number=1,Type=Float,Description="Strand bias p-value (Fisher exact)">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
##FORMAT=<ID=CI,Number=2,Type=Float,Description="95% credible interval for VAF">
##FORMAT=<ID=CONF,Number=1,Type=Float,Description="Posterior probability variant is real">
```

Note: VCF requires chromosome and position coordinates, which kam doesn't have (alignment-free). Options:
1. Use target_id as CHROM and relative position within target as POS
2. If targets have genomic coordinates in metadata, use those
3. Output a "pseudo-VCF" that uses target coordinates

### JSON (For Programmatic Consumption)

Full structured output including all evidence details, per-k-mer molecule counts, path details. Used by `kam generate-report` for the morning report.

## UMI Collision Probability Reporting

### Per-Family Estimate
For each molecule family, estimate the probability that at least one read in the family is a collision (different original molecule with same UMI):

```
P(collision) = 1 - (1 - 1/N_umi)^(n_molecules_at_position - 1)
```

Where:
- `N_umi` = UMI diversity (1,024 for Twist 5bp)
- `n_molecules_at_position` = estimated number of molecules at this genomic position

Since we don't have genomic position (alignment-free), estimate from:
- Total molecules / estimated panel coverage area
- Or from the local k-mer depth as a proxy

### Per-Sample Summary
Report in QC JSON:
- Estimated overall collision rate
- Number of families where P(collision) > 0.05
- Recommendation on whether depth exceeds UMI diversity limits

## Filter Expression Syntax (Future)

For configurable per-family filtering, support expressions like:
```
--filter "duplex_size >= 3 AND simplex_size >= 1"
--filter "family_size >= 5 OR (duplex_size >= 2 AND consensus_quality < 0.001)"
--downgrade-rule "duplex_size < 2 AND simplex_size >= 5 => SIMPLEX"
```

Initial implementation: fixed threshold parameters (min-family-size, min-duplex-size). Filter DSL as a later enhancement.

## References

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of the `RX` (raw UMI) and `MI` (molecule ID) SAM tag conventions adopted here, and the dash-separated UMI pair format.]

The SAM/BAM Format Specification Working Group (2024) Sequence Alignment/Map Format Specification. Available at: https://samtools.github.io/hts-specs/SAMv1.pdf (Accessed: March 2026). [Normative reference for SAM tag naming conventions (`XX:T:value` format) and the VCF 4.3 standard cited for the secondary variant output format.]

Twist Biosciences (2024) Twist UMI Adapter System. Product documentation. Available at: https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system (Accessed: March 2026). [Source of the 1,024 UMI diversity figure used in the per-family collision probability formula.]

Bourgey, M. et al. (2019) 'km: find-mutation', GitHub repository. Available at: https://github.com/iric-soft/km (Accessed: March 2026). [Source of the km-compatible TSV output format adopted as the primary variant call output.]
