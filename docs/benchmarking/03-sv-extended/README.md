# sv_new Benchmark Suite

Benchmarks for new SV types: inversions with deletion (InvDel), novel insertions
(NovIns), and fusion/translocation detection.

---

## Fusion Benchmark

### Design

The fusion benchmark tests kam's ability to detect a synthetic BCR-ABL1-like
breakpoint joining two regions on the 2000 bp chr1 reference:

- Gene A: chr1:100-300 (0-based), breakpoint at position 200 (anchor base: A).
- Gene B: chr1:800-1000 (0-based), breakpoint at position 900 (anchor base: T).

The fusion target for kam (`data/fusion_targets.fa`) is a 100 bp synthetic
sequence: the 50 bp ending at gene A's breakpoint followed by the 50 bp
starting at gene B's breakpoint.

```
>BCR_ABL1__chr1:150-200__chr1:900-950__fusion
[50bp from chr1:150-200][50bp from chr1:900-950]
```

### Why varforge cannot simulate fusions directly

varforge parses variant VCFs to spike mutations into reads generated from a
single reference sequence. Its VCF parser explicitly skips BND records (the
VCF notation for breakends). The `Translocation` variant type exists in
varforge's internal model but is not reachable from a VCF input. Passing a
BND VCF to varforge produces no fusion reads.

### Fused reference approach

Because varforge cannot generate fusion-spanning reads from a BND VCF, the
benchmark uses a fused reference strategy:

1. A synthetic reference (`data/fusion_ref.fa`) is constructed by
   concatenating the gene A side of the breakpoint (chr1:0-200) with the
   gene B side (chr1:900-2000). The breakpoint falls at position 200 in this
   fused reference.

2. varforge simulates reads from `fusion_ref.fa` at reduced coverage. Every
   read that spans position 200 in `fusion_ref.fa` is a fusion-spanning read
   containing the exact junction k-mers that kam's fusion detection uses.

3. A second varforge run simulates wild-type background reads from the normal
   `ref.fa` at the remaining coverage.

4. The two FASTQ sets are concatenated to produce the final mixed sample.

Coverage split at VAF `v`:

```
fusion_coverage = 5000 × v
wt_coverage     = 5000 × (1 - v)
```

### Files

| File | Description |
|------|-------------|
| `data/fusion_ref.fa` | Fused reference: chr1:0-200 + chr1:900-2000 (1300 bp) |
| `data/fusion_targets.fa` | kam fusion target FASTA (100 bp, one entry) |
| `data/truth_fusion_vaf<T>_<rep>.vcf` | BND truth VCF for each replicate |
| `configs/fusion_vaf<T>_<rep>.yaml` | varforge config for fusion reads |
| `configs/fusion_vaf<T>_<rep>_wt.yaml` | varforge config for WT background reads |
| `scripts/make_fusion_suite.py` | Generator script for all configs and VCFs |

### Generating the suite

Run from the repo root:

```bash
python3 docs/benchmarking/03-sv-extended/scripts/make_fusion_suite.py
```

This produces 50 fusion configs, 50 WT background configs, and 50 truth VCFs
(25 VAF levels × 2 replicates each). Seeds start at 130000+ to avoid
collisions with other suites (InvDel uses 90000+, NovIns uses 110000+).

### Running a simulation

For each VAF level and replicate, run two varforge jobs and then concatenate:

```bash
# Step 1: Simulate fusion reads.
varforge simulate -c configs/fusion_vaf0100_a.yaml

# Step 2: Simulate wild-type background reads.
varforge simulate -c configs/fusion_vaf0100_a_wt.yaml

# Step 3: Concatenate R1 and R2 FASTQ files.
cat results/sim_fusion_vaf0100_a/sample_R1.fastq.gz \
    results/sim_fusion_vaf0100_a_wt/sample_R1.fastq.gz \
    > results/sim_fusion_vaf0100_a_mixed/sample_R1.fastq.gz

cat results/sim_fusion_vaf0100_a/sample_R2.fastq.gz \
    results/sim_fusion_vaf0100_a_wt/sample_R2.fastq.gz \
    > results/sim_fusion_vaf0100_a_mixed/sample_R2.fastq.gz
```

### Running kam

Pass the fusion targets file via `--fusion-targets`:

Use `fusion_partner_targets.fa` (two 200 bp windows around each breakpoint)
as the primary `--targets` FASTA. Do NOT use `ref.fa` (2000 bp) — processing
the full reference causes each run to take 30+ minutes rather than ~20 seconds.

`fusion_partner_targets.fa` is pre-generated in `data/` and contains:
- `fusion_partner_A__chr1:100-300` (200 bp around gene A breakpoint)
- `fusion_partner_B__chr1:800-1000` (200 bp around gene B breakpoint)

```bash
kam run \
  --r1 results/sim_fusion_vaf0100_a_mixed/sample_R1.fastq.gz \
  --r2 results/sim_fusion_vaf0100_a_mixed/sample_R2.fastq.gz \
  --targets docs/benchmarking/03-sv-extended/data/fusion_partner_targets.fa \
  --fusion-targets docs/benchmarking/03-sv-extended/data/fusion_targets.fa \
  --output-dir results/kam_fusion_vaf0100_a/ \
  --output-format-override vcf,tsv
```

### Truth evaluation

A call is a true positive if kam emits a BND record matching `BCR_ABL1` with
breakpoints at chr1:200 and chr1:900 (within ±5 bp). Use the scoring script:

```bash
python3 docs/benchmarking/03-sv-extended/scripts/score_sv_new_suite.py \
  --variant-type fusion
```

---

## InvDel Benchmark

### Design

The InvDel truth variant is a 120 bp InvDel at chr1:300 on the 2000 bp
synthetic reference: 79 bp deleted, 60 bp inverted-inserted, net length
change 19 bp. Configs are generated by `scripts/make_invdel_suite.py`
and use seed range 90000+. 25 VAF levels × 2 replicates = 50 samples.

### Misclassification bug and fix

InvDel variants were originally misclassified as `Deletion` with a
`StrandBias` filter at all VAF levels. Root cause: `classify_variant` in
`kam-call/src/caller.rs` gated InvDel detection behind a net-length-change
check (`if del_len >= SV_LENGTH_THRESHOLD`). The benchmark InvDel has a net
length change of only 19 bp (79 deleted − 60 inverted-inserted), which is
below the 50 bp threshold, so the inversion test was never reached and the
variant was returned as `Deletion`. `assign_filter` then applied the
SNV/indel strand bias threshold (0.01) rather than the disabled SV
threshold (1.0), filtering all calls as `StrandBias`.

The fix (commit `0a9406f`) moves the `alt_seq_has_inversion_relative_to_ref`
call before the net-length gate. If the inversion pattern is detected, the
function returns `InvDel` immediately regardless of net length change.

### Results (current binary, post-fix)

| VAF  | Result | Notes |
|------|--------|-------|
| 0.5% | PASS | Pre-fix: Deletion/StrandBias |
| 1%   | PASS | Pre-fix: Deletion/StrandBias |
| 2%   | PASS | Pre-fix: Deletion/StrandBias |
| 5%   | PASS | Pre-fix: Deletion/StrandBias |

Any committed benchmark results predating commit `0a9406f` show the wrong
classification and should be treated as stale.

See `docs/project/investigations/invdel_misclassification_investigation.md`
for the full root-cause analysis.

---

## NovIns Benchmark

### Design

The NovelInsertion truth variant is an 80 bp novel sequence inserted at
chr1:400 on the 2000 bp synthetic reference. The inserted sequence is
verified to share no 20-mer with the reference or its reverse complement.
Configs are generated by `scripts/make_novins_suite.py` and use seed range
110000+. 25 VAF levels × 2 replicates = 50 samples.

### Stale committed results

**The committed benchmark results for NovIns are stale.** Results in
`results/sim_*/discovery/variants.tsv` show `TandemDuplication` as the
variant type because they were generated with a binary built before the
`NovelInsertion` classification existed (pre-commit `dad04a4`). In that
earlier code, all large insertions were unconditionally classified as
`TandemDuplication` regardless of sequence content.

The current code is correct: `is_tandem_duplication()` checks whether the
inserted sequence is a substring of the doubled reference, and the random
novel sequence fails that check, producing `NovelInsertion`. The detection
itself (alt path found, molecule count, confidence) is unaffected by the
classification fix — only the `variant_type` column changes.

To get correct results, regenerate by running the benchmark with the
current binary.

See `docs/project/investigations/novins_misclassification_investigation.md`
for the full analysis including commit timeline and code traces.
