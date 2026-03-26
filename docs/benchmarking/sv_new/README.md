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
python3 docs/benchmarking/sv_new/scripts/make_fusion_suite.py
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
  --targets docs/benchmarking/sv_new/data/fusion_partner_targets.fa \
  --fusion-targets docs/benchmarking/sv_new/data/fusion_targets.fa \
  --output-dir results/kam_fusion_vaf0100_a/ \
  --output-format-override vcf,tsv
```

### Truth evaluation

A call is a true positive if kam emits a BND record matching `BCR_ABL1` with
breakpoints at chr1:200 and chr1:900 (within ±5 bp). Use the standard
`scripts/score_benchmark.py` with `--variant-type fusion`.

---

## InvDel Benchmark

See `scripts/make_invdel_suite.py` for the generator.

The InvDel truth variant is a 120 bp InvDel at chr1:300 on the 2000 bp
synthetic reference. Configs use the seed range 90000+.

---

## NovIns Benchmark

See `scripts/make_novins_suite.py` for the generator.

Novel insertion benchmarks use the seed range 110000+.
