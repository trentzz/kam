# Configuration Examples

This page describes each example configuration file in the `examples/` directory.
Run any example with:

```bash
kam run --config examples/<name>.toml
```

CLI flags always override values set in the config file, so you can point an
example at your own files without editing the TOML:

```bash
kam run --config examples/minimal.toml \
        --r1 my_R1.fq.gz \
        --r2 my_R2.fq.gz \
        --targets my_panel.fa \
        --output-dir results/
```

---

## Basic usage

### minimal.toml

The smallest valid configuration. Provides only the four required fields
(`r1`, `r2`, `targets`, `output_dir`) and lets every other parameter fall back
to the built-in default. Use this as a starting point or for quick exploratory
runs where the defaults are acceptable.

```bash
kam run --config examples/minimal.toml
```

### discovery-mode.toml

Standard discovery mode for ctDNA panel sequencing. Detects SNVs, indels, and
SVs across all target regions without prior knowledge of expected variants.
Outputs both TSV and VCF with a QC JSON. The thresholds (confidence 0.99,
min 2 alt molecules) are appropriate for most Twist UMI duplex panel runs.

```bash
kam run --config examples/discovery-mode.toml
```

### tumour-informed.toml

Tumour-informed monitoring mode. Supply a VCF of known somatic variants from
the matched tumour biopsy (`target_variants`). Only calls matching that truth
set within the position tolerance are marked PASS. All other calls are emitted
with a `TumourInformed` filter tag so they remain visible but do not pollute
the signal. Near-zero false-positive rate for serial ctDNA monitoring.

```bash
kam run --config examples/tumour-informed.toml
```

### high-sensitivity.toml

Relaxed thresholds for maximum variant detection. Lowers the confidence
threshold to 0.90, accepts single-molecule calls, relaxes the strand bias
filter, and accepts UMI bases at Q15. Expect a higher false-positive rate.
Use for initial discovery, wet-lab validation, or when downstream filtering
will handle FP removal.

```bash
kam run --config examples/high-sensitivity.toml
```

### high-specificity.toml

Strict thresholds for zero false positives. Requires a posterior probability
of 0.999, at least three alt-supporting molecules, one duplex molecule on the
alt allele, and a strict strand bias filter. Also blocks germline-compatible
VAF (max VAF 0.35). Appropriate for clinical research reporting.

```bash
kam run --config examples/high-specificity.toml
```

---

## Chemistry

### twist-umi-duplex.toml

Full reference configuration for the Twist Biosciences duplex UMI panel.
Read structure is 5bp UMI + 2bp skip + template on both R1 and R2. All fields
are documented with their purpose and recommended values. Use this as a
reference when adapting to a different chemistry.

```bash
kam run --config examples/twist-umi-duplex.toml
```

### simplex-umi-12bp.toml

Configuration for protocols using a 12bp UMI with no skip/spacer and simplex
(single-strand) sequencing. Sets `duplex = false` and relaxes the strand bias
threshold to account for the structural asymmetry of simplex reads.

```bash
kam run --config examples/simplex-umi-12bp.toml
```

### simplex-umi-9bp.toml

Configuration for protocols using a 9bp UMI (262,144 possible values) with no
skip and simplex sequencing. Collision probability is low at typical depths;
the assembler corrects single-base UMI errors via Hamming-1 clustering.

```bash
kam run --config examples/simplex-umi-9bp.toml
```

### simplex-umi-8bp.toml

Configuration for protocols using an 8bp UMI (65,536 possible values) with
simplex sequencing. Notes the collision risk at high depth and recommends
switching to a longer UMI if molecule count inflation is observed.

```bash
kam run --config examples/simplex-umi-8bp.toml
```

### no-umi.toml

Configuration for standard (non-UMI) sequencing libraries. Sets
`umi_length = 0` and `skip_length = 0`. Without UMIs, duplicate collapse is
not available and sensitivity at low VAF is limited by sequencing error.
Suitable for high-depth amplicon panels where k-mer based calling is the goal.

```bash
kam run --config examples/no-umi.toml
```

---

## SV detection

### sv-detection.toml

SV mode for detecting large deletions, tandem duplications, novel insertions,
inversions, and InvDel events alongside SNVs and indels. Large deletions,
tandem duplications, and novel insertions are detected automatically without
any extra input. Inversions and InvDel events require an SV junction FASTA
(`sv_junctions`) whose k-mers are added to the allowlist so that
breakpoint-spanning reads are captured. SV thresholds are relaxed relative to
SNV thresholds because breakpoint k-mers are naturally rare.

```bash
kam run --config examples/sv-detection.toml
```

### fusion-detection.toml

Fusion gene detection using synthetic fusion target sequences. Each entry in
`fusion_targets` encodes both partner coordinates in its ID and spans the
breakpoint. Uses `kmer_size = 35` to reduce off-target k-mer hits across
homologous regions near fusion partners.

```bash
kam run --config examples/fusion-detection.toml
```

### sv-high-sensitivity.toml

Maximally permissive SV thresholds for low-VAF SV detection (0.5–2%). Lowers
`sv_min_confidence` to 0.70 and accepts Q15 UMI bases. Intended for
benchmarking varforge SV sweeps or research discovery. Expect elevated false
positives; all output should be reviewed before reporting.

```bash
kam run --config examples/sv-high-sensitivity.toml
```

---

## Panel and application

### ctdna-monitoring.toml

Production-grade ctDNA monitoring configuration. Combines tumour-informed
filtering, strict confidence (0.999), duplex confirmation on the alt allele,
and a germline VAF cutoff. Designed for serial plasma time points using a
matched tumour biopsy VCF. Use the same config and same `target_variants` file
across all time points for a given patient.

```bash
kam run --config examples/ctdna-monitoring.toml \
        --r1 plasma_T1_R1.fq.gz \
        --r2 plasma_T1_R2.fq.gz \
        --output-dir results/plasma_T1/
```

### research-discovery.toml

Full-verbosity research mode. Enables all four output formats (TSV, VCF, CSV,
JSON), a comprehensive set of detailed log channels, and a slightly relaxed
confidence threshold (0.95) to surface low-support variants for investigation.
Detailed logs include every dropped read, UMI clustering decision, family
summary, and filtered variant. Runtime and disk usage are significantly higher
than production configs.

```bash
kam run --config examples/research-discovery.toml
```

### large-panel.toml

Optimised for panels covering hundreds of genes (> 500 kb target space). Uses
`kmer_size = 35` to maintain specificity across the wide target space, and
raises the thread count to 16 to handle the larger indexing and calling
workload. All other thresholds are standard.

```bash
kam run --config examples/large-panel.toml
```

### low-input.toml

Adjusted for low-input samples (5–10 ng DNA). Accepts UMI bases at Q15 to
recover more reads from low-yield libraries, drops very short templates
aggressively (min 40 bp) to remove adapter dimers, and accepts singleton
families because UMI family sizes are small at low input. Enables the
`reads_dropped` and `families` log channels to monitor library quality.

```bash
kam run --config examples/low-input.toml
```

---

## Performance

### fast-mode.toml

Speed-optimised configuration. Uses `kmer_size = 25` for faster indexing,
outputs TSV only, suppresses the QC JSON, and disables all logging. Suitable
for CI integration tests, rapid QC screening, or preliminary passes before a
full-accuracy run.

```bash
kam run --config examples/fast-mode.toml
```

### thorough-mode.toml

Maximum-accuracy configuration for runs where runtime is not a constraint.
Uses `kmer_size = 39`, strict quality filters, duplex confirmation, all output
formats, and the full logging suite. Every pipeline decision is recorded.
Use for validation cohorts, reference samples, or paper result runs.

```bash
kam run --config examples/thorough-mode.toml
```
