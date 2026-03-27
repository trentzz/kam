# Tumour-Informed Pipeline Example

This example shows all three kam pipeline modes for ctDNA liquid biopsy analysis
using Twist UMI duplex chemistry:

1. **Discovery** — finds all variants across a panel, no filtering.
2. **Tumour-informed** — restricts PASS calls to a known mutation set from a matched biopsy.
3. **Alt-seq** — detects large SVs using junction sequences generated from alt allele contexts.

All three modes use [multiseqex](https://github.com/trentzz/multiseqex) to convert
genomic coordinates into the FASTA target files that kam requires.

---

## Prerequisites

- **multiseqex >= 0.2.1** — generates target FASTA files from a reference genome.

  ```bash
  # From crates.io
  cargo install multiseqex

  # From source (to get the latest version)
  git clone https://github.com/trentzz/multiseqex
  cargo install --path multiseqex
  ```

- **kam** — the variant detection tool.

  ```bash
  cargo install kam-bio
  ```

- **Reference FASTA** (e.g. `hg38.fa`). multiseqex builds a `.fai` index automatically
  if one is not present alongside the FASTA file.

---

## Input files

| File | Description | Format |
|------|-------------|--------|
| `inputs/mutations.vcf` | Somatic mutations from the matched tumour biopsy | VCF 4.2 |
| `inputs/sv_mutations.vcf` | SV mutations (deletions, InvDels) from the matched biopsy | VCF 4.2 with full REF sequences |
| `inputs/panel.bed` | Panel target regions for discovery mode | BED (0-based half-open) |

### mutations.vcf

Standard VCF 4.2. The ID column must be `.` so that multiseqex produces clean
`chr:start-end` FASTA headers. The `AF` INFO field records allele frequency from
the matched biopsy and is optional but informative.

```
#CHROM  POS        ID  REF      ALT  QUAL  FILTER  INFO
chr7    55181378   .   T        A    .     PASS    AF=0.42
chr17   7674220    .   C        T    .     PASS    AF=0.38
chr12   25398284   .   C        A    .     PASS    AF=0.61
chr3    178952085  .   A        G    .     PASS    AF=0.29
chr1    115256530  .   GGAATTT  G    .     PASS    AF=0.33
```

### sv_mutations.vcf

Same format as `mutations.vcf`, but REF alleles must be real nucleotide sequences
(not `N` or symbolic alleles). multiseqex `--alt-seq` works by substituting the ALT
sequence into the extracted reference window; it needs the real REF bases to do this.
These sequences come from the variant caller that produced the SV calls.

```
#CHROM  POS       ID  REF                           ALT                       FILTER  INFO
chr13   32337671  .   TATGCCTGACATGTTAAGCAGTGG...   T                         PASS    SVTYPE=DEL;SVLEN=-50;AF=0.27
chr17   41276045  .   ATGCTAGCATGCATGCTAGCATGC...   ATGCTAGCATGCATGCTCGATACG  PASS    SVTYPE=INVDEL;SVLEN=-60;AF=0.18
```

### panel.bed

0-based half-open coordinates, tab-separated. Regions are wider windows (~500-1000 bp)
covering the target hotspots. Used for discovery mode to define where to look across
the panel.

```
# CHROM  START      END        NAME
chr7     55181100   55181700   EGFR_exon20
chr17    7673800    7674700    TP53_exon7_8
chr12    25397800   25398800   KRAS_exon2
chr3     178951700  178952500  PIK3CA_exon20
chr13    32337300   32338100   BRCA2_region
```

---

## How it works

### Discovery mode

multiseqex reads `panel.bed` and extracts one FASTA sequence per region from the
reference. The `--name-template "{chr}:{start}-{end}"` flag ensures headers are in
the exact `chr:start-end` format that kam requires.

```bash
multiseqex hg38.fa \
    --bed inputs/panel.bed \
    --name-template "{chr}:{start}-{end}" \
    -o targets.fa
```

The resulting `targets.fa` has headers like:

```
>chr7:55181100-55181700
ATGCTAGCATGC...
```

kam indexes k-mers from the target sequences and the plasma reads, walks the de
Bruijn graph to find alt paths, then calls variants where alt k-mers are supported
by enough molecules. Without `--target-variants`, every variant passing quality
thresholds is marked PASS.

### Tumour-informed mode

multiseqex reads `mutations.vcf` and extracts one 401 bp window per variant (±200 bp
around the mutation site).

```bash
multiseqex hg38.fa \
    --vcf inputs/mutations.vcf \
    --flank 200 \
    --name-template "{chr}:{start}-{end}" \
    -o targets.fa
```

The resulting `targets.fa` has headers like:

```
>chr7:55181178-55181578
ATGCTAGCATGC...
```

kam runs with `--target-variants inputs/mutations.vcf`. Only variants whose position
matches an entry in the VCF (within `ti_position_tolerance` bp) are marked PASS. All
other variants are emitted with a `TumourInformed` filter tag and remain in the output
for review.

**Coordinate offset.** multiseqex uses the 1-based POS coordinate from the VCF to
label the FASTA header start. kam treats the header start as 0-based, producing a
systematic 1 bp offset: a variant at VCF POS 55181378 is called at position 55181379
in kam's output. Setting `ti_position_tolerance = 1` absorbs this offset. Do not use
0 when targets are generated from a VCF.

### Alt-seq mode

Alt-seq mode detects large structural variants (large deletions, InvDels) using two
multiseqex runs.

**Step 1** — extract reference windows (used for depth measurement and REF k-mer
building):

```bash
multiseqex hg38.fa \
    --vcf inputs/sv_mutations.vcf \
    --flank 200 \
    --name-template "{chr}:{start}-{end}" \
    -o sv_targets.fa
```

**Step 2** — generate junction sequences (alt allele contexts for the k-mer allowlist):

```bash
multiseqex hg38.fa \
    --vcf inputs/sv_mutations.vcf \
    --alt-seq \
    --flank 100 \
    -o sv_junctions.fa
```

`--alt-seq` replaces the REF bases in each extracted window with the ALT sequence,
producing a sequence that spans the breakpoint as it appears in a tumour read. kam
adds k-mers from `sv_junctions.fa` to the allowlist so that reads crossing the
breakpoint are captured during indexing. The headers in `sv_junctions.fa` are ignored
by kam; only the sequences matter.

SVs use `ti_position_tolerance = 10` because partial allele representations and
breakpoint imprecision can produce offsets larger than 1 bp.

---

## Running manually

### Discovery

```bash
# 1. Extract panel targets
multiseqex hg38.fa \
    --bed inputs/panel.bed \
    --name-template "{chr}:{start}-{end}" \
    -o results/targets.fa

# 2. Run kam
kam run \
    --config configs/discovery.toml \
    --r1 plasma_R1.fq.gz \
    --r2 plasma_R2.fq.gz \
    --targets results/targets.fa \
    --output-dir results/discovery/
```

### Tumour-informed

```bash
# 1. Extract per-mutation target windows
multiseqex hg38.fa \
    --vcf inputs/mutations.vcf \
    --flank 200 \
    --name-template "{chr}:{start}-{end}" \
    -o results/targets.fa

# 2. Run kam
kam run \
    --config configs/tumour-informed.toml \
    --r1 plasma_R1.fq.gz \
    --r2 plasma_R2.fq.gz \
    --targets results/targets.fa \
    --target-variants inputs/mutations.vcf \
    --output-dir results/tumour-informed/
```

### Alt-seq

```bash
# 1. Extract SV reference windows
multiseqex hg38.fa \
    --vcf inputs/sv_mutations.vcf \
    --flank 200 \
    --name-template "{chr}:{start}-{end}" \
    -o results/sv_targets.fa

# 2. Generate SV junction sequences
multiseqex hg38.fa \
    --vcf inputs/sv_mutations.vcf \
    --alt-seq \
    --flank 100 \
    -o results/sv_junctions.fa

# 3. Run kam
kam run \
    --config configs/alt-seq.toml \
    --r1 plasma_R1.fq.gz \
    --r2 plasma_R2.fq.gz \
    --targets results/sv_targets.fa \
    --sv-junctions results/sv_junctions.fa \
    --target-variants inputs/sv_mutations.vcf \
    --output-dir results/alt-seq/
```

Expected output structure:

```
results/
├── targets.fa                  # intermediate: target sequences for discovery/TI
├── discovery/
│   ├── variants.tsv
│   ├── variants.vcf
│   ├── qc.json
│   └── logs/
├── tumour-informed/
│   ├── variants.tsv
│   ├── variants.vcf
│   ├── qc.json
│   └── logs/
└── alt-seq/
    ├── sv_targets.fa           # intermediate: SV reference windows
    ├── sv_junctions.fa         # intermediate: SV junction sequences
    ├── variants.tsv
    ├── variants.vcf
    ├── qc.json
    └── logs/
```

---

## Running with the script

```bash
# Discovery
./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ discovery

# Tumour-informed
./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ tumour-informed

# Alt-seq (SV detection)
./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ alt-seq
```

Override tool paths if they are not in `PATH`:

```bash
KAM=/opt/tools/kam MULTISEQEX=~/.cargo/bin/multiseqex \
    ./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ tumour-informed
```

---

## Output

Each kam run produces:

| File | Description |
|------|-------------|
| `variants.tsv` | All candidate variants with molecule counts, VAF, confidence, and filter status |
| `variants.vcf` | Same variants in VCF 4.2 format |
| `qc.json` | Per-stage QC metrics (molecule counts, UMI collision rate, family size distribution) |
| `logs/` | Per-target assembly and drop logs |

**PASS vs filtered records.** In discovery mode, any variant passing quality thresholds
(`min_confidence`, `min_alt_molecules`, `strand_bias_threshold`, `max_vaf`) is marked
PASS. In tumour-informed mode, only variants matching a position in the `target_variants`
VCF (within `ti_position_tolerance`) are marked PASS. All other detected variants are
labelled `TumourInformed` and remain visible in the output but do not appear as PASS
calls. This is the mechanism that gives tumour-informed mode its near-zero false positive
rate for serial monitoring.
