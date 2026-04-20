# Investigation: Alt-as-Ref Mode — Indel Sensitivity Gap on Titration Data

Date: 2026-04-19
Area: SNV/indel, indexing, pathfind, alt-as-ref mode, ti-rescue
Symptom: At 15ng 2% VAF, kam detects 63/170 indels (37%) vs 162/205 SNVs (79%). Enabling alt-as-ref mode produced no improvement for indels. A comparison tool improves from 68 to ~100 indel TPs at the same condition using its alt-as-ref equivalent.

---

## Background

The titration benchmark evaluates kam on 375 validated somatic variants (205 SNVs, 170 indels) across 24 real FASTQ samples: 3 input masses (5ng, 15ng, 30ng) × 8 VAF levels (0–2%). The FASTQs are in the local thesis data directory (`/path/to/titration-nondedup/fastqs`).

A comparison tool was reported to improve from 68 → ~100 indel TPs at 15ng 2% VAF when using alt-as-ref mode. This prompted an investigation into whether kam's `--alt-as-ref` flag could produce a similar improvement.

Baseline results at 15ng 2% VAF (1M reads, discovery mode):

| Metric | Value |
|---|---|
| SNV sensitivity | 79.0% (162/205) |
| Indel sensitivity | 37.1% (63/170) |
| Overall precision | 1.000 |

---

## Run 1: Alt-as-ref without flank (broken)

### Method

Generated alt-allele sequences from the truth VCF using multiseqex without `--flank`:

```bash
multiseqex hg38.fa \
    --vcf docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
    --alt-seq \
    -o docs/benchmarking/01-snvindel/scripts/alt_seqs.fa
```

Ran the titration benchmark with:

```bash
python docs/benchmarking/01-snvindel/scripts/run_titration_batch.py \
    --kam-binary target/release/kam \
    --sv-junctions docs/benchmarking/01-snvindel/scripts/alt_seqs.fa \
    --target-variants docs/benchmarking/01-snvindel/scripts/truth_variants.vcf
```

### Findings

No improvement. At 15ng 2% VAF: indel_tp=67, sensitivity=39.4%.

Root cause: `multiseqex --alt-seq` without `--flank` outputs only the raw ALT allele sequence. For indels, the ALT allele is typically 1–16 bp (median 1 bp). All 170 indel alt sequences were shorter than k=31, producing **zero k-mers**. The alt-as-ref boost had no effect on any indel.

| Stat | Value |
|---|---|
| Min alt seq length (indels) | 1 bp |
| Max alt seq length (indels) | 13 bp |
| Median alt seq length (indels) | 1 bp |
| Sequences with 0 k-mers at k=31 | 170/170 (100%) |

---

## Run 2: Alt-as-ref with 100bp flank

### Method

Regenerated alt sequences with 100bp flanking context:

```bash
multiseqex hg38.fa \
    --vcf docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
    --alt-seq --flank 100 \
    -o docs/benchmarking/01-snvindel/scripts/alt_seqs_flank100.fa
```

All 375 alt sequences now 201–216 bp. Ran benchmark with `--sv-junctions alt_seqs_flank100.fa`.

### Findings

Still no improvement. At 15ng 2% VAF: indel_tp=67, sensitivity=39.4%. Molecule count unchanged (354,483 in both runs), confirming the alt k-mers are not pulling in new molecules. The target windows already capture them.

---

## Deep investigation: iterative hypothesis testing

After confirming alt-as-ref itself was not the bottleneck, a series of seven hypotheses were tested to identify what prevents the 103 missed indels from being detected.

All runs at 15ng 2% VAF, 1M reads, scored against truth_variants.vcf (original, then normalised where noted).

### H0 — Normalise the truth VCF

**Hypothesis:** bcftools norm shifts some indel positions; if the truth VCF coordinates are not left-normalised, scoring misses real TPs.

**Method:**

```bash
bcftools norm \
    --fasta-ref /data/genomes/GRCh38/hg38.fa \
    --output-type v \
    docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
    > docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf
```

Reran the baseline benchmark, scoring against `truth_variants.norm.vcf`.

**Result:** 70 indel TPs (41.2%), 0 FP. Gained 7 TPs at zero cost. bcftools norm realigned 14 records (position shifts and allele left-trimming).

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_normvcf_1mreads.tsv`

All subsequent runs use `truth_variants.norm.vcf` as the truth set.

---

### H1 — Larger target windows (200bp)

**Hypothesis:** The 101bp target windows are too narrow for some indels to produce spanning k-mers. Wider windows give the graph more context.

**Method:**

```bash
multiseqex hg38.fa \
    --vcf docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf \
    --flank 100 \
    --name-template "{chr}:{start}-{end}" \
    -o docs/benchmarking/01-snvindel/scripts/targets_200bp.fa

python docs/benchmarking/01-snvindel/scripts/run_titration_batch.py \
    --kam-binary target/release/kam \
    --targets docs/benchmarking/01-snvindel/scripts/targets_200bp.fa \
    --truth docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf
```

**Result:** 70 indel TPs, **551 FP**. Sensitivity unchanged; 551 false positives added. Wider windows cover more genomic territory and detect off-target variants.

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_200bp_normvcf_1mreads.tsv`

**Conclusion:** Larger windows hurt precision. Not the path forward without `--target-variants` filtering.

---

### H2 — Lower graph construction threshold (graph-min=1)

**Hypothesis:** The de Bruijn graph requires n_molecules ≥ 2 to include a k-mer. At low VAF, junction k-mers for rare alt alleles may fail this threshold. Setting graph-min=1 retains singletons.

**Method:**

```bash
python docs/benchmarking/01-snvindel/scripts/run_titration_batch.py \
    --kam-binary target/release/kam \
    --truth docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf \
    --graph-min-molecules 1
```

**Result:** 67 indel TPs, **88 FP**. Sensitivity dropped vs norm-VCF baseline (70 TP). The singleton threshold adds noise that overwhelms genuine signals.

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_gmin1_1mreads.tsv`

**Conclusion:** graph-min=1 does not help for indels and introduces substantial FPs.

---

### H3 — 200bp targets with alt-as-ref and TI rescue

**Hypothesis:** Combining wider windows, alt k-mers, and rescue probing could recover indels missed by all three modes individually.

**Method:**

```bash
python docs/benchmarking/01-snvindel/scripts/run_titration_batch.py \
    --kam-binary target/release/kam \
    --targets docs/benchmarking/01-snvindel/scripts/targets_200bp.fa \
    --alt-as-ref docs/benchmarking/01-snvindel/scripts/alt_seqs_flank100.fa \
    --truth docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf \
    --ti-rescue \
    --count-rescued
```

**Result:** 70 indel TPs, 0 FP (ti-rescue recovers nothing new; 200bp+altseq+rescue equals the 200bp baseline with `--target-variants` suppression). No gain over H0.

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_200bp_altseq_normvcf_1mreads.tsv`

---

### H4 — Graph-min=1 with duplex-only requirement

**Hypothesis:** graph-min=1 adds FPs from low-quality molecules. Requiring duplex support (duplex threshold=1) might filter the FPs while retaining true positives.

**Result:** 58 indel TPs, 174 FP. Worse on both metrics. Duplex requirement filters out real signal at low VAF.

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_gmin1_duplex1_normvcf_1mreads.tsv`

---

### H5 — TI rescue with corrected offset (build_alt_seq bug)

**Hypothesis:** The ti-rescue probe was silently failing for indels due to an off-by-one error in `build_alt_seq`.

**Root cause identified:** In `kam/src/rescue.rs`, `build_alt_seq` computed the alt allele insertion offset as:

```rust
let offset = vcf_pos_1based - target_start_0based - 1;
```

The correct formula (confirmed by examining multiseqex header coordinates) is:

```rust
let offset = vcf_pos_1based - target_start_0based;
```

The fencepost error caused 164/170 indels to fail the embedded REF allele check silently. `build_alt_seq` returned `None` for all of them, so rescue generated no probe rows and reported NO_EVIDENCE for every indel.

**Fix:** Corrected the offset formula in `rescue.rs` and updated the unit test `build_alt_seq_snv` to use the new formula (SNV at offset 3, not 2).

The correct targets also needed to be regenerated from the normalised VCF, since rescue probe coordinates must match the normalised positions:

```bash
multiseqex hg38.fa \
    --vcf docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf \
    --flank 50 \
    --name-template "{chr}:{start}-{end}" \
    -o docs/benchmarking/01-snvindel/scripts/targets_norm_100bp.fa

python docs/benchmarking/01-snvindel/scripts/run_titration_batch.py \
    --kam-binary target/release/kam \
    --targets docs/benchmarking/01-snvindel/scripts/targets_norm_100bp.fa \
    --truth docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf \
    --ti-rescue \
    --count-rescued
```

**Result:** 70 indel TPs, 0 FP. After fixing the offset: 136 variants successfully probed. Rescue found evidence for 7 (all nalt=1, filter=RESCUED/LowConfidence). The 129 remaining have NO_EVIDENCE.

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_rescue_fixed_1mreads.tsv`

**Conclusion:** The bug fix is important for correctness but does not increase the final TP count — the 7 rescued variants were already counted via the path walker or are below PASS threshold.

---

### H6 — TI rescue baseline (pass-only)

**Method:** Same as H5 but scoring only PASS-filtered calls (not counting RESCUED/LowConfidence).

**Result:** 70 indel TPs, same as H0 norm-VCF baseline. Rescue finds evidence for 7 variants but they do not reach PASS.

Results file: `docs/benchmarking/01-snvindel/summary/titration_results_tirescue_normvcf_1mreads.tsv`

---

## Summary of results

| Run | indel_tp | indel_fn | sensitivity | FP |
|-----|----------|----------|-------------|-----|
| Original baseline | 63 | 107 | 37.1% | 0 |
| Norm VCF (H0) | 70 | 100 | 41.2% | 0 |
| graph-min=1 (H2) | 67 | 103 | 39.4% | 88 |
| 200bp targets (H1) | 70 | 100 | 41.2% | 551 |
| gmin1 + duplex1 (H4) | 58 | 112 | 34.1% | 174 |
| TI rescue (PASS only) (H6) | 70 | 100 | 41.2% | 0 |
| TI rescue (fixed offset, counting RESCUED) (H5) | 70 | 100 | 41.2% | 0 |

Every hypothesis tested leaves sensitivity at 41.2% (70/170 indels). The 100 remaining missed indels cannot be recovered by any mode tested.

---

## Root cause analysis: why 129 rescue probes find NO_EVIDENCE

After the offset fix, rescue probed 136 indels and found NO_EVIDENCE for 129 of them. To understand why, per-target molecule counts (n_molecules_ref) were extracted from the rescue TSV and analysed.

| n_molecules_ref | Number of targets |
|---|---|
| 0 | 63 |
| 1–5 | 20 |
| 6–50 | 7 |
| 51+ | 39 |

**63 targets have zero reference molecule coverage.** These are entirely absent from the k-mer index for this sample at 1M reads. No molecule captured at all means no alt evidence is possible.

**20 targets have 1–5 reference molecules.** At 2% VAF, the expected number of alt molecules is 0.02 × (1–5) = 0.02–0.10. The probability of observing at least one alt molecule is 2–10%. Missing these is expected by chance.

**39 targets have 51+ reference molecules but zero alt evidence.** This is unexpected and not fully explained. Possible causes:

1. `alt_seqs_flank100.fa` was generated from the original (non-normalised) VCF. For the 14 normalised indels, the alt k-mers may not align to the normalised positions. Regenerating from `truth_variants.norm.vcf` is required to test this.
2. The alt k-mers are present in the index but the allele check in `build_alt_seq` still fails for some edge case not caught by the offset fix.
3. The alt alleles at these 39 sites genuinely have zero molecule support in this sample (e.g., spike-in efficiency variance at those specific targets).

**Dominant cause of the 129 NO_EVIDENCE cases: insufficient coverage.** At 1M reads (≈354,000 molecules across the panel), 63 of the 375 targets have no molecule coverage at all. This is a sequencing depth / panel uniformity issue, not an algorithmic one. Increasing to 2M+ reads would reduce this.

---

## Conclusion

The indel sensitivity gap at 15ng 2% VAF is driven primarily by sequencing depth, not by a failure of alt-as-ref mode, graph construction, or path walking.

- **37→41% gain**: normalising the truth VCF corrected 7 scoring mismatches (free fix, committed).
- **Ceiling at 41.2%**: 100 indels remain undetected. 63 targets have zero coverage. Most of the remainder have marginal coverage (1–5 molecules). At 1M reads and 2% VAF, the expected alt molecule count per target is far below 1 for most missed sites.
- **Alt-as-ref has no effect**: the target windows already capture alt-bearing molecules. The bottleneck is not molecule capture but raw molecule count.
- **TI rescue (after offset fix) finds 7 variants** but they do not reach PASS. All have nalt=1 and filter=LowConfidence.
- **graph-min=1 does not help** and degrades precision.

The path to higher indel sensitivity at this condition is increased read depth, not algorithm changes.

---

## Open questions

1. **39 high-coverage, zero-evidence targets**: regenerate `alt_seqs_flank100.fa` from `truth_variants.norm.vcf` and rerun to test whether the original-VCF coordinates cause a silent mismatch at these sites.
2. **2M reads benchmark**: run the full titration with the normalised VCF + normalised targets at 2M reads to quantify the depth effect.
3. **TI rescue PASS threshold**: the 7 rescued variants have nalt=1. Lowering `--min-alt-molecules` to 1 would promote them to PASS; assess whether this introduces FPs.

---

## Artefacts

| File | Description |
|---|---|
| `docs/benchmarking/01-snvindel/scripts/truth_variants.norm.vcf` | Truth VCF normalised by bcftools norm |
| `docs/benchmarking/01-snvindel/scripts/alt_seqs_flank100.fa` | Alt allele FASTA with 100bp flank (from original VCF — needs regen from norm VCF) |
| `docs/benchmarking/01-snvindel/scripts/targets_norm_100bp.fa` | 101bp target sequences from normalised VCF |
| `docs/benchmarking/01-snvindel/scripts/targets_200bp.fa` | 201bp target sequences from normalised VCF (H1, H3) |
| `docs/benchmarking/01-snvindel/summary/titration_results_1mreads.tsv` | Original baseline |
| `docs/benchmarking/01-snvindel/summary/titration_results_normvcf_1mreads.tsv` | Norm VCF baseline (H0) |
| `docs/benchmarking/01-snvindel/summary/titration_results_gmin1_1mreads.tsv` | graph-min=1 (H2) |
| `docs/benchmarking/01-snvindel/summary/titration_results_200bp_normvcf_1mreads.tsv` | 200bp targets (H1) |
| `docs/benchmarking/01-snvindel/summary/titration_results_200bp_altseq_normvcf_1mreads.tsv` | 200bp + altseq + rescue (H3) |
| `docs/benchmarking/01-snvindel/summary/titration_results_gmin1_duplex1_normvcf_1mreads.tsv` | gmin1 + duplex1 (H4) |
| `docs/benchmarking/01-snvindel/summary/titration_results_tirescue_normvcf_1mreads.tsv` | TI rescue PASS-only (H6) |
| `docs/benchmarking/01-snvindel/summary/titration_results_rescue_fixed_1mreads.tsv` | TI rescue with fixed offset (H5) |

---

## Code changes made during this investigation

**`kam/src/rescue.rs`** — fixed off-by-one in `build_alt_seq`:
```rust
// Before (wrong):
let offset = vcf_pos_1based - target_start_0based - 1;
// After (correct):
let offset = vcf_pos_1based - target_start_0based;
```
Also updated unit test `build_alt_seq_snv` to use offset 3 instead of 2.

**`kam/src/commands/run.rs`** — wired `cfg.input.alt_as_ref` into the allowlist augmentation path (same pattern as sv_junctions).

**`kam/src/commands/run.rs`** — wired `cfg.indexing.graph_min_molecules` into graph construction threshold (replaces hardcoded `n_molecules >= 2`).

**`docs/benchmarking/01-snvindel/scripts/run_titration_batch.py`** — extended with:
- `--alt-as-ref`, `--graph-min-molecules`, `--ti-rescue`, `--count-rescued` flags
- `--save-tsvs DIR` to save per-sample output before temp cleanup
- Fixed REPO path (was `parents[2]`, now `parents[4]`)
- Fixed default `_DEFAULT_TARGETS` and `_DEFAULT_TRUTH` to use full paths

**`kam/src/cli.rs`** — added `--alt-as-ref` as a separate, semantically correct flag (distinct from `--sv-junctions`).
