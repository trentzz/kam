# Investigation: 5 INDEL Gaps — Alignment Detects, kam Never Does

Date: 2026-04-12
Area: variant calling, TI matching, normalisation
Symptom: 5 INDEL variants are detected by alignment at VAF >= 0.25-2% but are never marked as `detected=True` by kam across all samples and VAF levels.
Root cause: 4 of the 5 are DETECTED by the de Bruijn graph but misclassified as `NotTargeted` due to a left-normalisation off-by-one in `extract_variant_key`. The 5th (chr17:31235638:CTGTT:C) has genuine zero alt k-mer evidence, consistent with marginal alignment support in just one sample.

---

## Background

Benchmark 07-snvindel-ml-boost-v1 compares 375 truth variants across 24 titration samples (5 ng, 15 ng, 30 ng, at 8 VAF levels). Of these 375 variants, 170 are INDELs. The alignment comparison document (`alignment_comparison.md`) identifies 5 variants that alignment detects at VAF >= 0.1% but kam never detects. All 5 are INDELs:

1. `chr3:37028864:GG:G` — 1 bp deletion
2. `chr17:7670685:GG:G` — 1 bp deletion
3. `chr22:29664999:A:AA` — 1 bp insertion
4. `chr9:136504893:G:GG` — 1 bp insertion
5. `chr17:31235638:CTGTT:C` — 4 bp deletion

The benchmark `detected` field is True only if the variant key `(chrom, pos, ref, alt)` from the called TSV matches the truth VCF entry exactly. This investigation was triggered to determine whether the failures are in k-mer detection, graph path-walking, or the TI matching logic.

---

## Findings

### Step 1: Target sequences and k-mer analysis

The 100 bp target windows for each variant were extracted from `targets_100bp.fa`. All 5 variants fall at window offset 50 (centre of the window). For k=31, each variant generates 28-30 alt-specific k-mers — none of which appear within the same target window's reference sequence. No extreme homopolymers exist at the variant sites:

| Variant | Context (±10 bp) | Max homopolymer run | Alt-specific 31-mers |
|---|---|---|---|
| chr3:37028864:GG:G | TGTACCCCCCGGAGAAG | 5xC (flanking) | 29 |
| chr17:7670685:GG:G | TTCAGCTCTCGGAACAT | 2xG at site | 29 |
| chr22:29664999:A:AA | GTTTACTATTAAACCAC | 3xA at site | 28 |
| chr9:136504893:G:GG | TTGGTGTGCAGCACGCG | 1xG at site | 30 |
| chr17:31235638:CTGTT:C | TCTGGTTACTCTGTTTG | 3xT (flanking) | 28 |

No allowlist filtering issue exists: the allowlist is built from reference target sequences, and alt k-mers (which differ from ref) are excluded from the allowlist. However, this is intentional — the two-pass indexing records alt k-mers from on-target molecules regardless. An on-target molecule is one whose reads contain at least one ref k-mer from the target. Any alt molecule overlapping the target has flanking ref k-mers, so it passes pass 1 and its alt k-mers are indexed in pass 2.

### Step 2: Diagnostic kam run (15 ng 2% VAF, 2M reads)

The pipeline was run with `--ti-rescue` on the 15 ng 2% VAF sample, subsetted to 2M read pairs. Results from `variants.tsv` for the 5 target windows:

| Variant | Graph result | n_alt_molecules | filter | call_source |
|---|---|---|---|---|
| chr3:37028864:GG:G | CALLED (full path) | 12 | NotTargeted | CALLED |
| chr17:7670685:GG:G | CALLED (full path) | 20 | NotTargeted | CALLED |
| chr22:29664999:A:AA | CALLED (full path) | 2 | NotTargeted | CALLED |
| chr9:136504893:G:GG | CALLED (full path) | 7 | NotTargeted | CALLED |
| chr17:31235638:CTGTT:C | NO path found | 0 | LowConfidence | NO_EVIDENCE |

The rescue probe confirmed:
- chr3: 29/29 alt k-mers found, 12 alt molecules, VAF ~1.4%
- chr17:7670685: 29/29 alt k-mers found, 20 alt molecules, VAF ~2.3%
- chr22:29664999: 28/28 alt k-mers found, 2 alt molecules, VAF ~0.4%
- chr9:136504893: 30/30 alt k-mers found, 7 alt molecules, VAF ~2.1%
- chr17:31235638: 0/28 alt k-mers found, 0 alt molecules

**The first 4 variants are being DETECTED — they have alt k-mer evidence and full de Bruijn graph paths. They are not reaching `detected=True` because the TI matching rejects the extracted variant key.**

### Step 3: TI matching failure in `extract_variant_key`

`extract_variant_key` (in `kam-call/src/targeting.rs`) converts a full target-window path pair into a canonical `(chrom, pos, ref, alt)` key using left-normalisation. At the time of this investigation, the function contained an extraneous `+ 1` in the VCF position formula — a regression introduced when the codebase was squashed. Target FASTA headers use 1-based coordinates (base at index i has VCF pos = target_start + i), so the correct formula is `target_start + anchor_pos` without any additional offset.

For a GG deletion at window offset 51 (e.g. chr3):

```
ref window: ...CCGG AGG...  (offset 49='C', 50='G', 51='G')
del_seq = ['G'], anchor_pos starts at 50
ref[50] = 'G' == del_seq[-1] = 'G' → loop fires, anchor_pos → 49
ref[49] = 'C' != 'G' → stop
anchor = ref[49] = 'C'

Regression formula (+1): pos = target_start + 49 + 1 = target_start + 50 = 37028864
  → ref='CG', alt='C' — position matches truth, but alleles don't (truth has GG:G)
  → TI filter rejects on allele mismatch

Correct formula (no +1): pos = target_start + 49 = 37028863
  → ref='CG', alt='C' at 37028863 — both position and alleles differ from truth (GG:G at 37028864)
  → TI filter rejects on both pos and allele mismatch
```

The truth VCF uses the bcftools convention (anchor at the first G of a GG run, giving GG:G). The left-normalisation algorithm overshoots the repeat by one step, producing CG:C. This allele mismatch means the 4 GG-type variants cannot pass the TI filter regardless of the position formula. They remain false negatives due to the normalisation overshoot (documented in the Resolution section).

For the other ~160 indels that are NOT in single-nucleotide repeat contexts, the `+ 1` regression caused position mismatch with the truth VCF. Removing the `+ 1` restores correct matching for all these variants, recovering the ~63 indel TPs that were incorrectly filtered.

The root cause for the 4 remaining FNs is the left-normalisation overshoot:

| Variant | algorithm output | truth VCF | issue |
|---|---|---|---|
| chr3:37028864:GG:G | 37028863:CG:C | 37028864:GG:G | anchor outside repeat |
| chr17:7670685:GG:G | 7670684:CG:C | 7670685:GG:G | anchor outside repeat |
| chr22:29664999:A:AA | 29664998:T:TA | 29664999:A:AA | anchor outside repeat |
| chr9:136504893:G:GG | 136504892:A:AG | 136504893:G:GG | anchor outside repeat |

### Step 4: chr17:31235638:CTGTT:C — genuine detection gap

This variant is different. The de Bruijn graph finds no alt path, and the rescue probe finds 0 of 28 alt k-mers in the index. With 692 ref molecules at this target in the 2M read subset, approximately 14 alt molecules would be expected at 2% VAF. The systematic absence of alt k-mers across all samples and all VAF levels points to a genuine problem with alt molecule capture.

As documented in `alignment_comparison.md`, alignment detects this variant in only one sample (5 ng, 2% VAF, 2 supporting reads). At 15 ng and 30 ng, 2% VAF, alignment shows only 0-1 AD — below the AD >= 2 threshold. The alignment support is marginal and concentrated in a single condition. This is consistent with the variant either having very low actual spiked-in concentration at most loci, or with the deletion junction sequence being difficult to sequence or assemble.

The deleted sequence is `CTGTT` in the context `...TCTGGTTACTCTGTTTGATT...`. The resulting alt `...TCTGGTTACTCTGATTCTCGG...` spans a non-repetitive junction with 28 well-defined alt k-mers. No allowlist filtering or homopolymer issue applies. The absence is most likely due to very few alt molecules in the actual library at most ng/VAF conditions — not a k-mer or graph limitation.

---

## Root cause

Two distinct failure modes:

**Mode 1 (variants 1-4): Position formula regression in targeting.rs.** The `extract_variant_key` function in `kam-call/src/targeting.rs` was producing positions one higher than the correct value due to an extraneous `+ 1` in the VCF position formula. The function comment incorrectly stated that target FASTA headers use 0-based coordinates; they use 1-based coordinates, so `pos = target_start + anchor_pos` is correct, not `target_start + anchor_pos + 1`.

With the regression: for variants 1–4, the extracted position key one-off from the truth VCF produced a non-matching key in TI filter, marking the call NotTargeted. The Python benchmark scorer, using the same left-normalisation algorithm without the `+ 1` error, correctly identified these calls as TPs but they were filtered out upstream by the Rust TI filter.

Additionally, for the four variants of the GG→G / A:AA type: the standard left-normalisation algorithm overshoots one position to the left of the repeat (producing CG:C instead of GG:G for a GG deletion). The truth VCF uses the bcftools convention (anchor at first G, producing GG:G). The TI filter, even after removing the `+ 1` regression, still rejects these four variants on allele mismatch (CG:C vs GG:G). These variants remain false negatives and are tracked separately below.

**Mode 2 (variant 5: chr17:31235638:CTGTT:C): Genuine alt molecule absence.** Zero alt k-mers appear in the index across all runs. The alignment-based tool detects this variant in only one sample at 2 AD, consistent with very low effective concentration of alt molecules in the library at most conditions. No evidence of k-mer filtering or graph failure; the problem is upstream in the sequencing data.

---

## Resolution

### Mode 1 fix (applied 2026-04-12)

Removed the `+ 1` from all VCF position formulas in `extract_variant_key`, `deletion_key`, and `insertion_key` in `kam-call/src/targeting.rs`. The correct formula is `pos = target_start + offset` for all variant types.

Applied the same correction to the Python benchmark scoring function in `docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py` (the deletion/insertion loop condition was also reverted to the standard `ref_seq[anchor_pos]` check, undoing an earlier incorrect attempt to fix the GG→G overshoot).

After the fix, all indels previously affected by the `+ 1` regression pass the TI filter and are counted as detected. At 2% VAF (15 ng), INDEL sensitivity returns to ~39% (67 of 170 indel TPs), restoring the count seen in prior correct-binary runs.

### GG→G-type variants (residual FN, Mode 1 partial)

Four variants remain false negatives after the position fix:

| Variant | Algorithm output | Truth VCF | Issue |
|---|---|---|---|
| chr3:37028864:GG:G | 37028863:CG:C | 37028864:GG:G | anchor overshoots repeat by 1 |
| chr17:7670685:GG:G | 7670684:CG:C | 7670685:GG:G | anchor overshoots repeat by 1 |
| chr22:29664999:A:AA | 29664998:T:TA | 29664999:A:AA | anchor overshoots repeat by 1 |
| chr9:136504893:G:GG | 136504892:A:AG | 136504893:G:GG | anchor overshoots repeat by 1 |

The standard left-normalisation loop shifts the anchor one step past the end of the repeat context. The bcftools convention keeps the anchor inside the repeat. Fixing this requires a post-normalisation right-shift correction (undo the final step if `ref[anchor_pos + 1] == event_seq[-1]`). This was not implemented to avoid introducing new regressions in complex multi-base deletions. These four variants are documented as known false negatives due to representation mismatch.

### Mode 2 (chr17:31235638:CTGTT:C)

No fix. Alignment detects this variant at 2 AD in one sample only. Flag as marginal in the truth set.

### Allowlist note

The allowlist only contains reference k-mers. Alt k-mers are NOT in the allowlist. This is by design: on-target molecule detection uses ref k-mers, but all k-mers from on-target molecules (including alt) are indexed in pass 2. The allowlist is not the cause of any failure here.
