# Fusion and Translocation Detection

## What fusions are

A fusion (or translocation) joins sequences from two non-adjacent genomic loci into a single
chimeric molecule. In cancer, gene fusions often arise from chromosomal rearrangements. Common
examples: BCR-ABL1 (CML), EML4-ALK (NSCLC), TMPRSS2-ERG (prostate cancer).

In kam, a fusion is detected when molecules span a synthetic sequence that joins two distinct
reference loci. Wild-type molecules do not span both loci; only tumour-derived rearranged
molecules do.

---

## Simpler input: junction sequences from BAM

If you have observed a fusion in a BAM file (e.g. via IGV), the simplest approach is
`--junction-sequences`. Copy the chimeric read sequence spanning the fusion junction and
save it in a FASTA file with any header:

```
>BCR_ABL1_from_igv
CAGAGGAAGAGCTGCAGACCATGCGTCTGTGGCCGCTGATCCTGTCGGAGCCATAATTAGCAGACCGGTACCTGAGGACTG
```

```bash
kam run \
  --r1 plasma_R1.fastq.gz \
  --r2 plasma_R2.fastq.gz \
  --targets panel_targets.fa \
  --junction-sequences junction_from_bam.fa \
  --output-dir results/ \
  --output-format tsv,vcf
```

Key benefits over `--fusion-targets`:

- **No coordinates needed.** Any FASTA header works.
- **Strand orientation is automatic.** The observed sequence already encodes the correct
  orientation.
- **Inserted nucleotides are included.** Random bases at the junction from DNA repair are
  part of the copied sequence.

The trade-off: `--junction-sequences` uses total library depth as the VAF denominator, not
partner-locus depth. For accurate partner-depth VAF, use `--fusion-targets` with coordinate
headers. For monitoring (is the fusion present or absent), `--junction-sequences` is
sufficient and far simpler.

See `guides/patient-sv-monitoring.md` for a complete walkthrough.

---

## How fusion detection works

Fusion detection reuses the standard walk/score/call pipeline. The difference is in the VAF
denominator and the input FASTA format.

1. **Fusion targets are synthetic sequences**: each entry in the `--fusion-targets` FASTA is a
   concatenation of a segment from partner locus A and a segment from partner locus B. The
   header encodes both loci.

2. **K-mers are added to the allowlist**: fusion target k-mers are added to the allowlist
   alongside normal target k-mers. This ensures that reads spanning the fusion junction are
   captured during indexing.

3. **The reference path is the fusion sequence itself**: for a fusion target, the "reference"
   in the graph walk is the synthetic junction sequence. Wild-type molecules from either locus
   appear as background. The path that covers the full fusion junction is the alt path.

4. **Partner depth provides the VAF denominator**: wild-type molecules at locus A and locus B
   are counted separately. Their average depth is the denominator for fusion VAF estimation.
   This is correct because a fused molecule replaces one copy from each locus.

5. **BND VCF notation**: fusions are written as paired BND records in VCF output, one record
   per breakpoint partner, following VCF 4.3 BND conventions.

---

## Fusion target FASTA format

The header for each entry must follow this exact format:

```
>{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion
```

Where:
- `{name}` is the fusion name (e.g. `BCR_ABL1`). Underscores are allowed; double underscores
  are not (they are used as field separators).
- `{chromA}:{startA}-{endA}` is the partner A breakpoint region (0-based, half-open interval).
- `{chromB}:{startB}-{endB}` is the partner B breakpoint region (0-based, half-open interval).
- `__fusion` is a required literal suffix that tells kam this is a fusion target.

The sequence is the partner A segment followed directly by the partner B segment. The `endA -
startA` bases are the partner A contribution; the `endB - startB` bases are partner B.

Example header:

```
>BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion
```

Here:
- `name = BCR_ABL1`
- Partner A: `chr22:23632500-23632550` (50 bp, the BCR breakpoint region)
- Partner B: `chr9:130854000-130854050` (50 bp, the ABL1 breakpoint region)
- Breakpoint position in the sequence: 50 (the first base of partner B)

A complete example entry:

```
>BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion
CAGAGGAAGAGCTGCAGACCATGCGTCTGTGGCCGCTGAAGGTGAGGAAAGGATTTTTGTGGAGCCATAATTAGCAGACC
```

The sequence has 80 bp: 40 bp from chr22 (BCR) and 40 bp from chr9 (ABL1). At least k=31 bases
on each side of the breakpoint are needed to produce junction-spanning k-mers.

---

## How to create fusion targets for your panel

1. Identify the breakpoint coordinates from a known fusion annotation or prior sequencing.

2. Extract 75 bp ending at the left partner breakpoint and 75 bp starting at the right partner
   breakpoint:

```bash
PARTNER_A=$(samtools faidx hg38.fa chr22:23632425-23632500 | grep -v '^>' | tr -d '\n')
PARTNER_B=$(samtools faidx hg38.fa chr9:130854000-130854075 | grep -v '^>' | tr -d '\n')
echo ">BCR_ABL1__chr22:23632425-23632500__chr9:130854000-130854075__fusion"
echo "${PARTNER_A}${PARTNER_B}"
```

3. The header coordinates must match the extracted sequence lengths. In this example:
   - Partner A: `chr22:23632425-23632500` = 75 bp
   - Partner B: `chr9:130854000-130854075` = 75 bp
   - `breakpoint_pos = endA - startA = 23632500 - 23632425 = 75`

4. Verify by counting characters: the sequence should be exactly `(endA - startA) +
   (endB - startB)` bases.

---

## Recommended settings

Default settings work for fusion detection. No fusion-specific thresholds are currently exposed.
The fusion caller uses `sv_min_confidence` (default 0.95) and `sv_min_alt_molecules` (default 1)
from the main caller config.

Typical command:

```
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --fusion-targets fusion_targets.fa \
  --output-dir results/ \
  --output-format tsv,vcf
```

Fusion targets and normal panel targets are processed on the same run. Pass normal targets to
`--targets` and fusion targets to `--fusion-targets`. They are kept separate internally: fusion
targets are called after normal targets so that partner depths are available for VAF estimation.

---

## What kam is good at

**Targeted fusion monitoring**: if your panel includes the breakpoint regions for known fusions,
kam detects them with the same sensitivity as SNVs once VAF exceeds 0.25%.

**Paired BND records**: VCF output includes both breakpoints as BND records with `MATEID`
linking them, compatible with standard SV tools.

**Simultaneous detection with SNVs and SVs**: run normal variant detection and fusion detection
in a single pass.

---

## What kam is not good at

**De novo fusion discovery**: kam is a targeted caller. It detects fusions at breakpoints
specified in `--fusion-targets`. It does not scan for novel fusions outside the panel.

**Fusions with inexact or variable breakpoints**: the partner depth denominator assumes that the
breakpoint is at the position encoded in the header. If the actual breakpoint differs from the
expected position by more than k bases, sensitivity drops.

**Fusions involving highly repetitive sequences**: breakpoint regions in repetitive elements
produce ambiguous k-mers that pollute the allowlist. Avoid designing fusion targets within
repetitive regions.

---

## BND VCF output format

Fusions produce two VCF records. For BCR-ABL1:

```
##ALT=<ID=BND,Description="Breakend (fusion junction)">
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of the mate breakend">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr22	23632500	BCR_ABL1_bnd1	N	N[chr9:130854000[	.	PASS	VAF=0.010000;...;SVTYPE=BND;MATEID=BCR_ABL1_bnd2
chr9	130854000	BCR_ABL1_bnd2	N	]chr22:23632500]N	.	PASS	VAF=0.010000;...;SVTYPE=BND;MATEID=BCR_ABL1_bnd1
```

The INFO fields carry the same statistical evidence as normal variant calls: `VAF`, `VAF_LO`,
`VAF_HI`, `NALT`, `NDUPALT`, `CONF`, etc.

**Note on fusion VAF**: the `VAF` field in a fusion VCF record is computed using the partner
depth as the denominator, not the fusion target depth. The fusion target depth reflects only
the rearranged molecules, which is not a valid denominator for VAF. Partner depth (the average
molecule count at the two partner loci) is the appropriate reference.
