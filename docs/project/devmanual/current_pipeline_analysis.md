# Current Pipeline Analysis: Information Loss at Each Stage

## Pipeline Flow

```
FASTQ (duplex UMI sequenced)
  ↓
HUMID — reference-free deduplication/clustering
  ↓
Jellyfish — k-mer counting → .jf database
  ↓
km (find_mutation) — k-mer graph walking against target sequences → TSV output
  ↓
kmtools (chunk/merge/filter/plot) — parallelisation wrapper, filtering, visualisation
```

`multiseqex` handles batching reference sequence extraction for target inputs to km.
The pipeline is containerised via Docker.

## Information Flow Compression

```
Duplex reads with UMIs
[molecule identity, strand identity, base quality per position, read pair linkage]
    ↓ HUMID (loses duplex strand pairing, no consensus calling)
Deduplicated FASTQ
[sequence only, no molecule provenance]
    ↓ Jellyfish -C (loses strand info, quality info, read-of-origin)
k-mer count table {kmer_sequence → integer}
    ↓ km graph walking
Variant ratio [float]
    ↓ kmtools filter
Pass/fail against threshold
```

## Per-Stage Loss Analysis

### Stage 1: HUMID
- **Lost:** Duplex strand pairing. No concept of complementary UMI pairs.
- **Consequence:** Duplex read families handled incorrectly — either collapsed to single simplex read (losing strand info for duplex consensus) or reads from complementary strand families misclustered.
- **Impact:** By the time data hits Jellyfish, the duplex error-correction benefit paid for in library prep is partially or wholly lost.

### Stage 2: Jellyfish
- **Canonical k-mer collapsing (`-C` flag):** Erases strand of observation. A+B and B+A count merged. For duplex data, asymmetry between forward and reverse k-mer counts IS signal (true variant = both strands; error = one strand). Lost.
- **Quality score discarding:** No quality-weighted counting. K-mer seen 10× at Q30 looks identical to 10× at Q10. For 0.1% VAF discrimination, this matters enormously.
- **Read-of-origin discarding:** Cannot trace k-mer back to which read/pair/UMI family. 5 appearances from 5 families (strong) vs 5 PCR duplicates of one molecule (weak) — indistinguishable.
- **What Jellyfish does well:** Extremely fast, memory-efficient. Hard to beat for pure counting throughput.

### Stage 3: km (find_mutation)
- **`--ratio` parameter:** Single global threshold (e.g., 0.0001) without accounting for coverage depth. A ratio of 0.0001 at 10,000× is very different evidence from same ratio at 100×.
- **No statistical model:** Hard cutoff, not Bayesian or significance-based.
- **Calibrated for RNA-seq:** Original context, not ultra-deep ctDNA panel sequencing.

### Stage 4: kmtools
- Exists primarily to work around km's lack of multithreading (chunks targets, runs parallel km processes, merges).
- `filter` matches output against known variant reference.
- `plot` provides VAF histograms, type distributions, per-sample summaries.

## What's Lost By the End

By the time km runs, you have a single integer count per k-mer. Thrown away:
1. Which UMI family contributed this k-mer (critical for duplex consensus)
2. Which strand it came from (critical for duplex validation)
3. Quality scores of bases in the k-mer (critical for error discrimination at 0.1% VAF)
4. Whether k-mers co-occur on the same read (adjacent k-mer support = stronger evidence)

The pipeline bets that after dedup + counting, true mutant k-mers will have counts reflecting true VAF. But without strand-separated counting or quality weighting, the noise floor is set by Jellyfish's flat integer count rather than by the duplex chemistry's actual error rate.
