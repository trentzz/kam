# Structural Variant Detection

## What SVs kam can detect

kam detects five structural variant types:

| Type | VCF symbolic allele | Description |
|------|--------------------|------------------------------------|
| `LargeDeletion` | `<DEL>` | Deletion of ≥ 50 bp |
| `TandemDuplication` | `<DUP>` | Insertion where the alt is a tandem repeat of nearby reference sequence, ≥ 50 bp longer than ref |
| `Inversion` | `<INV>` | Same-length rearrangement where the alt is the reverse complement of the corresponding ref segment |
| `NovelInsertion` | `<INS>` | Insertion of ≥ 50 bp that is not a tandem repeat of nearby reference sequence |
| `InvDel` | `<INVDEL>` | Inversion combined with a flanking deletion (compound SV) |

The threshold between short indel and SV is 50 bp (`SV_LENGTH_THRESHOLD`).

---

## How SV detection works

SV detection reuses the same de Bruijn graph walk as SNV/indel detection, with one addition:
**junction k-mers**.

For most SVs, the breakpoint sequence is not present in the target reference. A deletion removes
sequence, so the k-mers spanning the deletion breakpoint do not appear in the unmodified
reference window. An inversion introduces a reverse-complement junction that is absent from the
forward-strand reference.

The `--sv-junctions` FASTA provides these breakpoint sequences. When loaded:

1. Junction k-mers are added to the allowlist. This ensures reads spanning the breakpoint are
   captured during indexing.
2. The path walker's alt-path budget is extended (+250 k-mers) to accommodate long SV paths.
3. The caller classifies paths as SV types when the length difference or reverse-complement test
   exceeds the 50 bp threshold.

For SV types, the alt molecule count uses the mean molecule count over variant-specific k-mers
(k-mers present in the alt path but absent from the reference path), not the minimum. This
avoids undercounting: long SV paths (70–150 k-mers) have very few variant-specific k-mers per
molecule, and the minimum bottlenecks at 1–3 even at moderate VAF.

---

## Benchmark sensitivity

On synthetic titration data (varforge VAF sweeps, 2 replicates each):

| VAF | SV sensitivity (discovery mode) |
|-----|--------------------------------|
| 0.05% | 25% |
| 0.10% | 75% |
| 0.15% | 75% |
| 0.20% | 75% |
| 0.25% | 100% |
| 0.30% | 100% |
| ≥0.40% | 100% |

SVs reach 100% sensitivity at 0.25% VAF and maintain it above that threshold.

---

## Required inputs: --sv-junctions FASTA

The `--sv-junctions` flag accepts a FASTA file containing the breakpoint-spanning sequences for
each SV target. These sequences are synthetic: they represent the sequence that appears at the
junction in the tumour genome.

Each junction entry needs two things:
- A sequence that spans the breakpoint. Provide at least k (default 31) bases on each side of
  the breakpoint, so at least 62 bp total. Longer is better: 100–200 bp gives more k-mers for
  robustly capturing junction-spanning reads.
- A FASTA header identifying the target. The header can be any string; it does not need to
  follow a special format.

---

## Junction FASTA format

```
>TP53_del_exon5-6_junction
CTTGGGCTTGCAGCCCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCT
>EGFR_del_exon19_junction
AAGGTGAGGAAAGGATTTTTGTGGAGCCATAATTAGCAGACCGGTACCTGAGG
```

Each sequence must be long enough that at least one k-mer (31 bp by default) falls entirely
within the junction, spanning both sides of the breakpoint.

---

## Creating junction sequences for your panel

To create a junction sequence for a known deletion:

1. Identify the breakpoint coordinates (chr, left breakpoint, right breakpoint).
2. Extract 75 bp of reference sequence ending at the left breakpoint.
3. Extract 75 bp of reference sequence starting at the right breakpoint.
4. Concatenate the two segments. The result is the 150 bp junction sequence.

Example using samtools:

```bash
# Deletion breakpoints: chr17:7674100 (left), chr17:7675200 (right)
LEFT=$(samtools faidx hg38.fa chr17:7673925-7674100 | grep -v '^>' | tr -d '\n')
RIGHT=$(samtools faidx hg38.fa chr17:7675200-7675375 | grep -v '^>' | tr -d '\n')
echo ">TP53_large_del_junction"
echo "${LEFT}${RIGHT}"
```

For inversions, the junction sequence on the alt haplotype contains a reverse-complemented
segment. Construct this synthetically:

```bash
# Inversion: chr17:7674100–7674500 is inverted
LEFT=$(samtools faidx hg38.fa chr17:7673925-7674100 | grep -v '^>' | tr -d '\n')
INV=$(samtools faidx hg38.fa chr17:7674100-7674500 | grep -v '^>' | tr -d '\n' | rev | tr ACGTacgt TGCAtgca)
RIGHT=$(samtools faidx hg38.fa chr17:7674500-7674675 | grep -v '^>' | tr -d '\n')
echo ">TP53_inversion_junction"
echo "${LEFT}${INV:0:75}${RIGHT:0:75}"  # Only need the breakpoint ends
```

For tandem duplications, the junction is the right end of the duplicated segment joined to the
left end:

```bash
# TandemDup: chr17:7674100–7674300 is duplicated
RIGHT_END=$(samtools faidx hg38.fa chr17:7674225-7674300 | grep -v '^>' | tr -d '\n')
LEFT_START=$(samtools faidx hg38.fa chr17:7674100-7674175 | grep -v '^>' | tr -d '\n')
echo ">TP53_tandem_dup_junction"
echo "${RIGHT_END}${LEFT_START}"
```

---

## Recommended settings

| Setting | Default | Recommended for SV |
|---------|---------|-------------------|
| `sv_min_confidence` | 0.95 | Keep at 0.95 (lower than SNV threshold is appropriate for SVs) |
| `sv_min_alt_molecules` | 1 | Keep at 1 in monitoring mode; raise to 2 for discovery |
| `sv_strand_bias_threshold` | 1.0 (disabled) | Keep disabled for inversions; inversions are structurally strand-biased |
| `--min-confidence` | 0.99 | SNV/indel threshold; does not apply to SV types |

The SV-specific thresholds are intentionally lower than the SNV thresholds. A 150 bp deletion
supported by two molecules is qualitatively stronger evidence than two molecules supporting a
single-base change: the probability that a random sequencing error mimics a 150 bp deletion is
negligible. The posterior model accounts for this, so the threshold can safely be lower.

---

## What kam is good at

**SV detection from 0.25% VAF**: 100% sensitivity at 0.25% VAF in benchmarks, with 75%
sensitivity at 0.1% VAF.

**Zero false positives in tumour-informed mode**: monitoring mode with `--target-variants`
suppresses all calls not matching a pre-specified allele. Combined with SV junction k-mers,
this provides clean monitoring of known SVs.

**All five SV types on the same run**: pass all junction sequences in a single
`--sv-junctions` file alongside normal targets. No separate run is needed.

---

## What kam is not good at

**InvDel (compound SVs)**: InvDel requires both a length change and a reverse-complement
segment above the 50 bp threshold. These are rare in panel sequencing and the evidence is
fragmented across multiple variant-specific k-mers. Treat InvDel calls as provisional and
validate with orthogonal methods.

**SVs without junction k-mers**: if you do not provide `--sv-junctions`, kam will only detect
SVs that happen to be fully contained within a target window (uncommon for deletions >50 bp).
Junction k-mers are required for reliable SV detection.

**De novo SV discovery**: kam is a targeted caller. It detects SVs at specified target windows.
It does not scan the genome for novel SVs outside the panel.

**Very long SVs (>500 bp)**: the path walker budget (default +250 k-mers for SV paths) may not
accommodate very long tandem duplications or insertions. Use the `_maxpathN` suffix in the
target ID to increase the budget for specific targets.

---

## Example command

```
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --sv-junctions sv_junctions.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

For monitoring mode with a known SV:

```
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --sv-junctions sv_junctions.fa \
  --target-variants biopsy_somatic.vcf \
  --output-dir results/ \
  --output-format tsv,vcf
```

The biopsy VCF should contain the SV call using the same `target_id` as CHROM and the symbolic
alt notation (`<DEL>`, `<DUP>`, `<INV>`, `<INS>`, `<INVDEL>`).

---

## SV-specific output fields

SV calls in the VCF carry additional INFO tags:

| Tag | Meaning |
|-----|---------|
| `SVTYPE` | `DEL`, `DUP`, `INV`, `INS`, or `INVDEL` |
| `SVLEN` | Length difference: negative for deletions, positive for insertions, 0 for inversions |

The REF field in the VCF is the anchor base (single base). The ALT field is the symbolic
allele (`<DEL>`, etc.).

In TSV output, `ref_seq` and `alt_seq` contain the full reconstructed path sequences. For a
150 bp deletion, `ref_seq` is ~150 bp longer than `alt_seq`.

The `n_molecules_alt` field for SV calls reflects the mean molecule count over variant-specific
k-mers, not the minimum. This is noted in `call_qc.json` but is transparent in the output.
