# Patient SV Monitoring: From BAM Observation to ctDNA Tracking

## Overview

This guide covers the workflow for monitoring a known structural variant or fusion in serial
cfDNA (cell-free DNA) samples. The typical scenario: a tissue biopsy reveals a fusion or
rearrangement, you see the chimeric reads in IGV, and you want to track that event in the
patient's blood over time.

This is the primary use case for `--junction-sequences`. You provide the observed junction
sequence directly. kam adds its k-mers to the allowlist, walks the junction as a target, and
reports whether supporting molecules exist in each cfDNA sample.

**Who this is for**: clinical researchers, bioinformaticians, and molecular pathologists who
have identified a structural event in a tumour and want to monitor it in liquid biopsy samples.

**What you need**:

- A BAM file (or IGV session) showing the chimeric reads at the structural variant
- Paired-end cfDNA FASTQ files to analyse (one or more time points)
- A panel target FASTA for the standard variant targets

---

## Prerequisites

| Requirement | Description |
|---|---|
| BAM with the SV | Tumour or tissue biopsy BAM showing chimeric reads spanning the junction |
| cfDNA FASTQs | Paired R1/R2 FASTQ files from plasma samples (Twist UMI duplex chemistry) |
| Panel targets | FASTA of your standard panel target windows (`--targets`) |
| kam binary | Version 0.3.0 or later (junction sequence support) |

You do not need a reference genome. kam is alignment-free.

You do not need to know the exact genomic coordinates of the junction. The sequence itself is
sufficient.

---

## Understanding junction sequence input

### Why the observed sequence handles orientation automatically

When you see a chimeric read in IGV, the sequence you observe already encodes the correct
strand orientation. If partner A is on the forward strand and partner B is reverse-complemented,
the read sequence shows exactly that arrangement. By copying the sequence as displayed, you
capture the orientation without needing to specify it separately.

This is a key advantage over coordinate-based input (`--fusion-targets`), which assumes
forward-forward orientation by default and requires explicit handling for other orientations.

### Why random nucleotides at the junction are not a problem

Many structural rearrangements have 1-5 random nucleotides inserted at the breakpoint during
DNA repair (non-templated insertion). These bases are part of the observed sequence. When you
copy the chimeric read, the inserted bases are included. kam will index and walk the sequence
as provided, including any inserted bases. No special handling is needed.

### When to use each input mode

| Flag | Input format | When to use |
|---|---|---|
| `--junction-sequences` | Any FASTA, any header | You have the observed sequence from BAM/IGV. No coordinates needed. |
| `--fusion-targets` | Coordinate-encoded header | You have exact breakpoint coordinates and want partner-depth VAF. |
| `--sv-junctions` | Any FASTA, any header | You have inversion or InvDel junction sequences for intra-chromosomal SVs within a panel target. |

`--junction-sequences` is the simplest option. Use it when you have the sequence and want to
monitor it. Use `--fusion-targets` only when you need the partner-depth VAF denominator (which
requires coordinates to identify the two partner loci).

`--sv-junctions` is for augmenting the k-mer allowlist so that inversion and InvDel
breakpoint-spanning reads are captured during indexing. It works alongside `--targets`, not as
a standalone target. `--junction-sequences` entries are walked as standalone targets.

---

## Step-by-step walkthrough

### Step 1: Find the chimeric read in IGV

Open the tumour BAM in IGV and navigate to the region of interest.

Look for:
- Reads with soft-clipped segments (coloured clips in IGV)
- Split reads spanning two chromosomes or distant loci
- Discordant read pairs pointing to different chromosomes

Select a chimeric read that spans the junction. In IGV, right-click the read and choose
"Copy read sequence". Alternatively, use "Copy details" and extract the SEQ field.

You need at least 60 bases of sequence spanning the junction: 30 bases on each side of the
breakpoint. Longer is better. 80-120 bases gives more robust k-mer coverage.

If multiple chimeric reads span the same junction, pick the longest one. If they disagree on
the junction sequence (different inserted bases, different breakpoint positions), use the
consensus or the most common variant.

### Step 2: Create the junction FASTA

Create a FASTA file with one entry per junction you want to monitor. The header can be any
descriptive string. It does not need coordinates.

```
>BCR_ABL1_patient_junction
CAGAGGAAGAGCTGCAGACCATGCGTCTGTGGCCGCTGATCCTGTCGGAGCCATAATTAGCAGACCGGTACCTGAGGACTG
>EML4_ALK_patient_junction
TGGAGCCTTGAGATCTTCAACTGCTGAAGACTTCTAGTCCTACAGATCAACAATTTACCCAG
```

Header format does not matter. These are all valid:

```
>my_fusion_1
>patient_007_BCR-ABL1_v1
>junction_from_igv_chr22_chr9
>any_text_you_want
```

Save this file as `junction_sequences.fa`.

### Step 3: Run kam

```bash
kam run \
  --r1 plasma_R1.fastq.gz \
  --r2 plasma_R2.fastq.gz \
  --targets panel_targets.fa \
  --junction-sequences junction_sequences.fa \
  --output-dir results/ \
  --output-format tsv,vcf
```

This runs the full pipeline. Junction sequence k-mers are added to the allowlist alongside
your panel targets. Each junction entry is walked as a standalone target.

### Step 4: Interpret the output

Open `results/variants.tsv`. Look for rows where `target_id` matches your junction FASTA
headers.

Key columns:

| Column | What to look for |
|---|---|
| `target_id` | Your junction header (e.g. `BCR_ABL1_patient_junction`) |
| `variant_class` | `Fusion` or the detected SV type |
| `vaf` | Variant allele frequency in the cfDNA sample |
| `n_molecules_alt` | Number of molecules supporting the junction |
| `n_duplex_alt` | Duplex-confirmed molecules (strongest evidence) |
| `confidence` | Posterior probability the call is real |
| `filter` | `PASS` or the reason for filtering |

A PASS call with `n_molecules_alt >= 2` and `n_duplex_alt >= 1` is strong evidence that the
junction is present in the cfDNA. A call with `n_molecules_alt = 1` and `n_duplex_alt = 0`
is weaker but may still be real at very low VAF.

If no row appears for your junction target, it means zero junction-spanning k-mers were found
in the cfDNA sample. The variant is either absent or below the detection limit.

---

## Worked example: monitoring BCR-ABL1 in CML

### Scenario

A patient with chronic myeloid leukaemia (CML) has a confirmed BCR-ABL1 fusion from bone
marrow biopsy. The chimeric read sequence from the tumour BAM is:

```
CAGAGGAAGAGCTGCAGACCATGCGTCTGTGGCCGCTGATCCTGTCGGAGCCATAATTAGCAGACCGGTACCTGAGGACTG
```

This 80 bp sequence spans the BCR-ABL1 junction. The first 40 bases are from BCR (chr22), the
last 40 bases are from ABL1 (chr9), with 2 inserted bases (`TC`) at position 39-40 from
non-templated repair.

### Junction FASTA

Create `bcr_abl1_junction.fa`:

```
>BCR_ABL1_CML_patient
CAGAGGAAGAGCTGCAGACCATGCGTCTGTGGCCGCTGATCCTGTCGGAGCCATAATTAGCAGACCGGTACCTGAGGACTG
```

### Run command

```bash
kam run \
  --r1 plasma_T1_R1.fastq.gz \
  --r2 plasma_T1_R2.fastq.gz \
  --targets cml_panel_targets.fa \
  --junction-sequences bcr_abl1_junction.fa \
  --output-dir results/plasma_T1/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

### Expected output

In `results/plasma_T1/variants.tsv`, a positive result looks like:

```
target_id                vaf      n_molecules_ref  n_molecules_alt  n_duplex_alt  confidence  filter
BCR_ABL1_CML_patient     0.0034   112000           385              72            0.9999      PASS
```

This indicates 385 molecules spanning the BCR-ABL1 junction at 0.34% VAF, with 72 duplex
confirmations. The call is high confidence.

### Serial monitoring

Run the same command on each time point, changing only the FASTQ paths and output directory:

```bash
for tp in T1 T2 T3 T4; do
  kam run \
    --r1 "plasma_${tp}_R1.fastq.gz" \
    --r2 "plasma_${tp}_R2.fastq.gz" \
    --targets cml_panel_targets.fa \
    --junction-sequences bcr_abl1_junction.fa \
    --output-dir "results/plasma_${tp}/" \
    --max-vaf 0.35 \
    --output-format tsv,vcf
done
```

Extract the fusion VAF across time points:

```bash
for tp in T1 T2 T3 T4; do
  echo -n "${tp}: "
  awk -F'\t' '$1 == "BCR_ABL1_CML_patient" {print $3}' "results/plasma_${tp}/variants.tsv"
done
```

---

## Using tumour-informed mode for maximum precision

If you want zero false positives across the entire call set (not just the junction), combine
`--junction-sequences` with `--target-variants`. The junction sequence provides the k-mers for
detection. The target variants VCF restricts PASS calls to only the expected variants.

```bash
kam run \
  --r1 plasma_R1.fastq.gz \
  --r2 plasma_R2.fastq.gz \
  --targets panel_targets.fa \
  --junction-sequences bcr_abl1_junction.fa \
  --target-variants patient_somatic.vcf \
  --output-dir results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

The `patient_somatic.vcf` should contain the BND record for the fusion alongside any SNV and
indel targets from the tumour biopsy. Calls not matching a target variant entry are labelled
`NotTargeted`.

This combination gives both: the sequence-level detection from `--junction-sequences` and the
specificity of tumour-informed filtering from `--target-variants`.

---

## Troubleshooting

### No junction k-mers found

**Symptom**: no row appears for your junction target in the output.

**Possible causes**:

1. **Sequence too short.** The junction sequence must be at least `2 * k` bases (default 62 bp)
   to produce k-mers spanning the breakpoint. Extend to at least 80 bp by including more
   flanking sequence from the chimeric read.

2. **Sequence from the wrong strand.** If you copied the sequence from a reverse-strand read,
   it may not match the forward-strand k-mers in the library. Try reverse-complementing the
   sequence. In practice, kam canonicalises k-mers, so this should not be an issue for standard
   k-mer sizes.

3. **The variant is absent from this sample.** If the patient is in remission or the cfDNA
   fraction is very low, there may be zero junction-spanning molecules.

### Low molecule support

**Symptom**: the call has `n_molecules_alt = 1` or `n_molecules_alt = 2` with no duplex.

**What it means**: the VAF may be very low (< 0.1%). At low VAF, stochastic sampling means
few molecules carry the variant allele.

**What to do**:

- Check the `n_molecules_ref` column. If total depth is low (< 50,000 molecules), increased
  sequencing depth would help.
- If monitoring over time, a consistent low-level signal across multiple time points is more
  credible than a single observation.
- Consider lowering `--sv-min-confidence` to 0.90 for exploratory analysis (increases
  sensitivity but also increases false positives).

### Call filtered NotTargeted

**Symptom**: the junction call has `filter = NotTargeted`.

**Cause**: you passed `--target-variants` but the junction call does not match any entry in
the VCF. The junction call coordinates (derived from the FASTA header) do not match the VCF
CHROM/POS/REF/ALT.

**Fix**: either add the junction as a BND entry in the target variants VCF, or increase
`--ti-position-tolerance` to allow fuzzy matching, or remove `--target-variants` if you do
not need tumour-informed filtering for this run.

### Many false positives from discovery mode

**Symptom**: many PASS calls that are not real variants.

**What to do**: add `--target-variants` pointing to a VCF of known somatic variants from the
tumour biopsy. This restricts PASS calls to expected variants only. Discovery mode typically
produces 35-72 background calls per sample. Monitoring mode produces zero false positives
for a correct truth set.

---

## Input mode comparison

| Feature | `--junction-sequences` | `--fusion-targets` | `--sv-junctions` |
|---|---|---|---|
| Header format | Any | Must follow `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion` | Any |
| Coordinates needed | No | Yes | No |
| Walked as standalone target | Yes | Yes | No (augments allowlist only) |
| VAF denominator | Total library depth | Partner locus depth | N/A (paired with `--targets`) |
| Handles all orientations | Yes (sequence encodes it) | FF only (currently) | Yes (sequence encodes it) |
| Handles inserted bases at junction | Yes (in the sequence) | No (concatenates reference segments) | Yes (in the sequence) |
| Best for | BAM/IGV observations, monitoring | Known coordinate-based fusions | Inversion/InvDel junction k-mers |
| VCF output | Uses header as target_id | Paired BND records with coordinates | SV calls on parent target |

---

## Tips

- **Copy more sequence than you think you need.** 80-120 bp gives robust k-mer coverage.
  31 bp on each side of the junction is the minimum.

- **Use the same junction FASTA across all time points.** Consistency ensures that VAF
  trends are comparable.

- **Combine with standard panel targets.** `--junction-sequences` adds to, not replaces,
  the panel targets. You get SNV/indel calls from `--targets` and junction monitoring from
  `--junction-sequences` in the same run.

- **One FASTA for all junctions.** Put all junctions (multiple fusions, multiple SVs) in a
  single FASTA file. Each entry is walked independently.

- **Check the sequence length.** If a junction call never appears despite high depth, the
  sequence may be too short. Aim for at least 60 bases, ideally 80-120.
