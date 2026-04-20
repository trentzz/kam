# Alt-as-Ref Mode

## What alt-as-ref mode is

In standard discovery mode, kam identifies variants by finding paths through the de Bruijn graph that diverge from the reference. The graph is built from k-mers supported by at least 2 molecules. At low VAF, some variant-specific k-mers fall below this threshold by chance, breaking the alt path and causing a miss.

Alt-as-ref mode adds ALT allele sequences for known variants — SNVs and indels — directly to the k-mer allowlist. These alt k-mers are indexed regardless of their raw molecule support, giving the graph more signal to reconstruct the alt path at low VAF.

It is most useful when you have a list of known variants from a prior tissue biopsy and want to maximise sensitivity for those exact positions in a liquid biopsy sample.

---

## How it works

When `--alt-as-ref alt_seqs.fa` is passed, kam:

1. Reads each sequence from the FASTA and extracts its canonical k-mers.
2. Adds these k-mers to the allowlist alongside the reference target k-mers.
3. During two-pass indexing, any molecule whose reads contain an alt k-mer is captured for full indexing — even if it contains no reference k-mers.
4. During graph construction, alt k-mers with molecule support contribute to the de Bruijn graph, strengthening the alt path signal.

The effect is that variants whose alt k-mers would otherwise be filtered out by the molecule-support threshold are instead retained, allowing the path walker to find and score them.

This is distinct from `--sv-junctions`, which provides junction sequences for structural variant breakpoints (inversions, InvDel events). Use `--alt-as-ref` for SNVs and short indels; use `--sv-junctions` for SVs.

---

## Alt-allele sequence file format

The `--alt-as-ref` flag takes a standard FASTA file. **Any header format is accepted** — headers are only used for logging and are never parsed for coordinates. Only the sequences matter.

The one requirement is that each sequence is at least `k` bp long (default k=31). For short alt alleles (SNVs, short indels), include flanking reference context on each side so the junction-spanning k-mers can be extracted.

---

## Workflow A — from a VCF (have coordinates)

Use this when you have a VCF of known variants from a tissue biopsy.

**Prerequisites:** a GRCh38 reference FASTA indexed with `samtools faidx`, and [`multiseqex`](https://github.com/trentzz/multiseqex) (`cargo install multiseqex`).

### 1. Generate alt-allele sequences

```bash
multiseqex hg38.fa \
  --vcf tumour_biopsy_mutations.vcf \
  --alt-seq --flank 100 \
  -o alt_seqs.fa
```

The `--flank 100` argument adds 100 bp of reference sequence on each side of the ALT allele, ensuring junction-spanning k-mers are generated even for 1 bp SNVs and short indels. Without flanking, most alt sequences are shorter than k=31 and produce zero k-mers.

Verify the output:

```bash
awk '/^>/{next} {print length($0)}' alt_seqs.fa | sort -n | head -5
# All values should be well above 31
```

### 2. Run kam

```bash
kam run \
  --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa \
  --output-dir results/ \
  --alt-as-ref alt_seqs.fa \
  --target-variants tumour_biopsy_mutations.vcf \
  --output-format-override tsv,vcf
```

The `--target-variants` flag is optional but recommended. It applies tumour-informed filtering so only variants matching the known list are marked PASS.

### 3. Inspect results

```bash
grep -v "^#" results/variants.vcf | awk '$7=="PASS"' | wc -l
awk -F'\t' 'NR==1 || $7=="PASS"' results/variants.tsv | head -20
```

---

## Workflow B — from raw sequences (no coordinates)

Use this when you have an observed alt sequence from IGV, a BAM pileup, or km/kmtools output, but do not have precise genomic coordinates.

### 1. Create a FASTA with your sequences

Any header is accepted. The sequence must be long enough to contain k-mers spanning the variant junction — a ~70 bp window around the alt allele is sufficient at k=31:

```
>observed_patient_042_tp53
GCAATGCCCTCAGAATCTGTTCAGTTTGTACTTCTGATGACAGAAAGAGCCTCAGAACATCCCCAAACT
>igv_copy_egfr_exon19
TAGCTGGATTGTCAGTGCGCTTTTCCCAACACCACTTGCTCCAACCACCACCAGTTTGTACTCAGTCATTT
```

If your sequence is shorter than 31 bp, manually add flanking reference bases from the genome or use Workflow A instead.

### 2. Run kam

```bash
kam run \
  --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa \
  --output-dir results/ \
  --alt-as-ref raw_alt_seqs.fa \
  --output-format-override tsv,vcf
```

Without `--target-variants`, all quality-passing variants are marked PASS — standard discovery mode with alt k-mer augmentation.

---

## Combined with TI rescue

For maximum sensitivity at very low VAF, combine `--alt-as-ref` with `--ti-rescue`. Rescue probing queries the k-mer index directly for each undetected TI variant, recovering sub-threshold evidence that the path walker missed:

```bash
kam run \
  --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa \
  --output-dir results/ \
  --alt-as-ref alt_seqs.fa \
  --target-variants tumour_biopsy_mutations.vcf \
  --ti-rescue \
  --min-alt-molecules 1 \
  --output-format-override tsv,vcf
```

---

## When to use alt-as-ref vs other modes

| Scenario | Recommended mode |
|---|---|
| No prior variant list, panel-wide discovery | Standard discovery (no extra flags) |
| Known variants from tissue biopsy, standard VAF (≥0.5%) | `--target-variants` alone |
| Known variants, low VAF (<0.5%), short indels | `--alt-as-ref` + `--target-variants` |
| Raw observed sequences without coordinates | `--alt-as-ref` with Workflow B |
| Structural variant breakpoints (inversions, InvDel) | `--sv-junctions` |
| No known variants but want raw junction sequences | `--junction-sequences` |

---

## Limitations

- **Sequence length:** each sequence must be ≥ k bp (default 31). Short raw alt alleles without flanking produce zero k-mers and have no effect.
- **Large indels:** indels larger than ~30 bp create a long chain of alt-specific k-mers. If molecule coverage is very low, many k-mers in the chain may still fall below the graph construction threshold. This is an inherent graph-based detection limit, not a bug.
- **SNV impact:** for SNVs, the benefit of `--alt-as-ref` is smaller than for indels. SNVs produce k-mers in the target reference window already, so molecule capture is rarely the bottleneck.
- **Not a substitute for depth:** alt-as-ref lowers the effective detection threshold but cannot recover variants with zero supporting molecules. At very low VAF, sequencing depth remains the primary lever.
- **Symbolic alleles (Workflow A only):** `multiseqex --alt-seq` requires real nucleotide REF and ALT sequences. `<DEL>`, `<INS>`, and `<CNV>` are not supported.
