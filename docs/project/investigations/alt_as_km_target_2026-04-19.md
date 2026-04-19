# Investigation: ALT sequences as km find_mutation targets

**Date:** 2026-04-19  
**Sample:** 15ng 2% VAF (full FASTQ jellyfish DB; also tested 1M reads)  
**Hypothesis:** Using ALT sequences (reference with variant applied) as km input targets, rather than REF sequences, recovers 100+ of 170 truth indels.  
**Result:** ALT-as-target gives 74/170 (43.5%), not 100+. This is 4 above the REF-target baseline of 70/170. The claim of 100+ is not reproducible with this dataset.

---

## Background

The user reported getting 100+ indels by supplying ALT sequences as km `find_mutation` targets. The standard approach uses reference genomic sequences (REF targets) and km detects mutations as deviations from those. The ALT-as-target approach inverts this: the post-variant sequence is the target, so ALT-allele reads appear as km "Reference" type (matching the target exactly).

---

## Method

### ALT target generation

For each of 170 truth indels from `indel-only-reference.tsv` (matching `truth_variants.vcf`):

1. Extract a 141bp genomic window centred on the variant (70bp each side) from GRCh38.
2. Apply the user's bash logic to produce the ALT sequence:
   - **Deletions** (REF > ALT): remove `del_len = ref_len - alt_len` bases after the anchor base at offset 70.
   - **Insertions** (ALT > REF): insert `ins_bases = alt[ref_len:]` after the anchor base at offset 70.
3. Write one FASTA file per variant. Headers use the format `>chrN:start-end`.

Window coordinate: `win_start (1-based) = POS - 70`, `offset = POS - win_start = 70`.

Script: `/home/tzeng/tmp/indel_alt_as_target/generate_alt_targets.py`

All 170 variants generated without error. ALT files are in `/home/tzeng/tmp/indel_alt_as_target/per_variant_fasta/`.

### km runs

```
km find_mutation --count 1 --ratio 0.0001 <ALT_TARGET.fa> <JELLYFISH_DB>
```

Run in parallel (4 threads via xargs) for all 170 targets. Two DBs tested:

- **Full FASTQ DB**: `/home/tzeng/tmp/kmtools_15ng_full/kmers_full.jf` (k=31, built from complete 15ng 2% FASTQ)
- **1M reads DB**: `/home/tzeng/tmp/kmtools_15ng_2pc/kmers.jf` (k=31, 1M reads subset)

Output in `/home/tzeng/tmp/indel_alt_as_target/km_out_full/` and `km_out_1M/`.

### Filtering

Standard `kmtools filter --use-alt` does not work here because it matches the km "Sequence" column against TYPE-annotated ALT sequences, but the ALT allele reads appear as km **"Reference"** type (not "Deletion" or "Insertion"). A custom filter was written (`/home/tzeng/tmp/indel_alt_as_target/generate_alt_targets.py` contains the logic; filtering is in a standalone Python block).

Custom filter: for each variant, find the "Reference" km row at `Info == "vs_ref"`, verify the `Sequence` matches the expected ALT target sequence exactly, and require `Min_coverage >= 1`.

---

## Results

### Truth counts

| Approach | DB | Targets | TRUE / 170 | % |
|---|---|---|---|---|
| REF target, standard filter | full FASTQ | 141bp | 70 | 41.2% |
| REF target, standard filter | 1M reads | 100bp | 72 | 42.4% |
| REF target, --use-alt filter | 1M reads | 100bp | 72 | 42.4% |
| **ALT as target**, custom filter | full FASTQ | 141bp | **74** | **43.5%** |
| ALT as target, custom filter | 1M reads | 100bp | 72 | 42.4% |
| ALT as target, normed VCF | full FASTQ | 141bp | 74 | 43.5% |
| ALT as target, 71bp window | full FASTQ | 71bp | 74 | 43.5% |

The ALT-as-target approach with the full FASTQ DB finds **4 more** than the REF-target baseline.

### Variant overlap

All 70 REF-approach findings are a strict subset of the 74 ALT-approach findings. The 4 extra variants found **only** by ALT-as-target:

| Variant | REF | ALT | Truth VCF type | reference.tsv type |
|---|---|---|---|---|
| chr3:41233407 | G | GGGA | INDEL (insertion) | Deletion |
| chr7:116778811 | A | ATT | INDEL (insertion) | Deletion |
| chr10:8073911 | C | CG | INDEL (insertion) | Deletion |
| chr14:104792618 | T | TC | INDEL (insertion) | Deletion |

All 4 are **insertions** (ALT > REF) mislabelled as "Deletion" in `indel-only-reference.tsv`. The standard REF-target km run does find "Insertion" type rows for these in the merged output (`merged_indel.txt`), but `kmtools filter` rejects them because TYPE "Deletion" (from reference.tsv) does not match km Type "Insertion". The ALT-as-target approach bypasses this TYPE check — ALT-allele reads appear as "Reference" type regardless of indel direction.

### Coverage ceiling

Of the 96 variants not found by any approach:

- All 96 have **zero Reference coverage** in the ALT-target km output. The ALT-allele k-mers simply are not present in the jellyfish DB.
- Breakdown: 78 deletions (sizes 1–30bp, median 11bp; 39 are > 10bp), 18 insertions (sizes 1–15bp, median 4bp).
- Neither the full FASTQ DB nor the 1M reads DB covers these variants.
- Neither 71bp, 100bp, nor 141bp window sizes change the count.
- Neither normalised nor unnormed VCF coordinates change the count (28 positions differ between the two VCFs, but neither set recovers additional variants).

The ceiling is **74 at 15ng 2% VAF** with these k-mer databases. This matches the pattern seen in the companion investigation (`reach_100_indels_2026-04-19.md`) where the true ceiling was established at ~79 with 5M reads.

---

## Root cause of the 100+ claim

The claim of 100+ is not reproducible. Possible explanations:

1. A different sample was used (e.g., 30ng 2% VAF). The original thesis (30ng 2%, HUMID deduped) achieved 65–68 with the dedup-humid-v3 kmtools pipeline. With less stringent filters (±20bp position tolerance, filter-v3) it reaches 51/170 — still well below 100.
2. A different truth set was used (not the 170-variant indel-only reference).
3. A different filtering criterion was applied (e.g., counting any km expression > 0 without requiring TYPE match or ALT sequence verification).

No configuration tested here achieves 100/170 on this sample. The data ceiling appears to be 74 with the full FASTQ jellyfish DB.

---

## Key findings

1. **ALT-as-target gives +4 over REF-target** on the full FASTQ DB (74 vs 70). The improvement is real but small.

2. **The +4 improvement is an artefact of TYPE mislabelling**, not a genuine algorithmic advantage. The 4 extra variants are insertions mislabelled as "Deletion" in `indel-only-reference.tsv`. Fixing the reference.tsv TYPE field would recover them with the REF-target approach too.

3. **The 96 missing variants have zero ALT-allele k-mer coverage**. They are absent from the sample, not missed by the algorithm. ALT-as-target does not help for variants with no supporting reads.

4. **Standard `kmtools filter --use-alt` does not work for ALT-as-target km runs**. ALT-allele reads appear as km "Reference" type, not as "Deletion"/"Insertion". The existing kmtools filter implementation requires TYPE matching. A custom filter was required.

5. **Window size (71bp, 100bp, 141bp) does not affect outcome** for this approach. Coverage is the bottleneck, not k-mer span.

---

## Files

| File | Description |
|---|---|
| `/home/tzeng/tmp/indel_alt_as_target/generate_alt_targets.py` | Script to generate ALT target FASTA files |
| `/home/tzeng/tmp/indel_alt_as_target/per_variant_fasta/` | 170 ALT target FASTA files (141bp windows) |
| `/home/tzeng/tmp/indel_alt_as_target/km_out_full/` | km output for each ALT target, full FASTQ DB |
| `/home/tzeng/tmp/indel_alt_as_target/km_out_1M/` | km output for each ALT target, 1M reads DB |
| `/home/tzeng/tmp/indel_alt_as_target/custom_filtered_full.tsv` | Custom-filtered results, full DB |
| `/home/tzeng/tmp/indel_alt_as_target/custom_filtered_1M.tsv` | Custom-filtered results, 1M reads DB |
| `/home/tzeng/tmp/indel_alt_as_target/alt_reference_indels.tsv` | ALT sequences for --use-alt mode (141bp) |
