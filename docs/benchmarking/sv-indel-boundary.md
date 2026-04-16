# The SV/INDEL Boundary in kam: Design, Benchmarking, and Classification

This document answers a question that comes up naturally when evaluating kam's
large-variant detection: what do we actually mean by "SV" in these benchmarks,
and where does a long INDEL end and a structural variant begin?

---

## The short answer

In kam, the boundary is 50 bp. Variants with a length difference between the
reference and alternate allele of less than 50 bp are classified as INDELs;
variants with a length difference of 50 bp or more are classified as large
variants (LargeDeletion, TandemDuplication, NovelInsertion, Inversion, InvDel,
or Fusion). The 50 bp threshold is a convention borrowed from the VCF
specification and the broader SV field. It is not a mechanistic boundary.

The critical point, which the early benchmark documentation obscured, is that
two of the six "SV" types — LargeDeletion and NovelInsertion, and to a large
extent TandemDuplication — are detected by exactly the same de Bruijn graph
path walk as a 10 bp or 30 bp INDEL. They are long INDELs that happen to cross
the 50 bp threshold. The remaining three types — Inversion, InvDel, and Fusion
— are genuinely structurally distinct and require junction k-mer injection to
be detectable at all.

---

## What the benchmarks actually test

### Benchmark 02-sv-core

Three variant types embedded in a 2000 bp synthetic reference (chr1):

| Variant | Position | Size | Detected via |
|---------|----------|------|--------------|
| LargeDeletion | chr1:200-299 | 100 bp removed | Standard path walk |
| TandemDuplication | chr1:500-599 | 100 bp inserted in tandem | Standard path walk |
| Inversion | chr1:900-999 | 100 bp RC segment | Junction k-mer injection |

VAF levels: 0.5%, 1%, 2%, 5%. Simulated with varforge at 5000× raw read
depth (~1000 molecule families per target), Twist duplex UMI chemistry. All
simulated molecules were simplex (no duplex pairs produced by the simulator at
any VAF level).

### Benchmark 03-sv-extended

Three additional types:

| Variant | Size | Detected via |
|---------|------|--------------|
| InvDel | 120 bp inversion + net deletion | Junction k-mer injection |
| NovelInsertion | 99 bp novel sequence | Standard path walk |
| Fusion (BCR-ABL1-like) | Cross-locus breakpoint at chr1:200 + chr1:900 | Fusion junction injection |

25 VAF levels (0.05%–10%) × 2 replicates each per type.

---

## Why long deletions and insertions really are just long INDELs

### The detection mechanism

kam builds a de Bruijn graph from all k-mers in the molecule index. For a
target region, the reference path is the sequence of overlapping k-mers tiling
the reference. An alternate allele creates a diverging path that rejoins the
reference after the variant site.

For a 10 bp deletion, the alt path skips 10 reference bases and rejoins. For a
100 bp deletion, the alt path skips 100 reference bases and rejoins. The graph
walk (depth-first search from the left anchor k-mer to the right anchor k-mer)
finds both. The key requirement is that the target window is wide enough to
contain anchor k-mers on both sides of the variant. For a 100 bp deletion in a
100 bp target, the anchors are tight; for a 100 bp deletion in a 200 bp target,
there is comfortable flanking sequence.

### What changes at large sizes

Two things change as variant size increases, but neither is mechanistic:

1. **Evidence scoring.** For short variants, `min_molecules` (the minimum
   molecule count across all k-mers on the alt path) is a reliable signal.
   For a 100-k-mer alt path, `min_molecules` bottlenecks at 1–3 regardless
   of true VAF, because every position in a long path must be independently
   covered. For large variants (LargeDeletion, TandemDuplication, Inversion),
   kam uses `mean_variant_specific_molecules` — the mean molecule count across
   only the k-mers unique to the alt allele — which provides a more stable
   estimate.

2. **Calling thresholds.** The default confidence threshold is 0.99 for
   SNVs and INDELs and 0.95 for large variants. The minimum alt molecule count
   is 2 for SNVs/INDELs and 1 for large variants. The rationale: a 100 bp
   structural signature supported by a single molecule is qualitatively
   stronger evidence than a single molecule supporting a single base change,
   because background sequencing error cannot independently produce a 100 bp
   structural event.

### The 50 bp threshold is a convention

Nothing in the biology or the detection algorithm changes abruptly at 50 bp.
A 49 bp deletion goes through the INDEL code path; a 51 bp deletion goes through
the large-variant code path with relaxed thresholds. In practice this means
a 51 bp deletion is *more* likely to pass filters than a 49 bp deletion, because
it benefits from lower confidence requirement and single-molecule support.

The 50 bp boundary is where the VCF specification recommends switching from
inline REF/ALT representation to symbolic alleles (`<DEL>`, `<INS>`, etc.).
kam follows this convention in its VCF output but applies it to the calling
logic as well, which is a design decision that could reasonably have been made
differently.

---

## Classification at the boundary: what actually happens

When kam recovers an alt path, it computes `ref_len` and `alt_len` from the
path sequences. The classification logic in `kam-call/src/caller.rs`:

```
if |alt_len - ref_len| < 50:
    → SNV, MNV, Insertion, Deletion, or Complex (INDEL class)
else:
    → check for inversion signal
    → if net deletion: LargeDeletion or InvDel
    → if net insertion:
        → if inserted sequence matches nearby ref window: TandemDuplication
        → else: NovelInsertion
    → if same length and central segment is RC: Inversion
```

A variant exactly at the boundary will be classified as INDEL (the check is
`< 50`, not `<= 50`). A 50 bp deletion is classified as `LargeDeletion`.

In the 02-sv-core benchmark, the embedded deletions are 100 bp, well above the
threshold. But in a real panel sample, a deletion could fall anywhere from 1 bp
to several kilobases. Any deletion of 50 bp or more will be classified and
reported as `LargeDeletion`, triggering the large-variant thresholds regardless
of whether the user intended to call "SVs".

---

## What the benchmarks do and do not show

### What they show

- A 100 bp deletion is reliably detected at 0.5% VAF (PASS, 14 alt molecules,
  confidence 1.000) using only the standard path walk, no junction sequences.
- A 100 bp tandem duplication is detected at 0.5% VAF (PASS, 2 alt molecules,
  confidence 0.982) with the SV-specific confidence threshold of 0.95.
- A 100 bp inversion is detected at 0.5% VAF (PASS, 3 alt molecules,
  confidence 0.999) using junction k-mer injection.
- At 1% VAF, tandem duplication coverage is 1 molecule (LowConf), reflecting
  the limited number of k-mers unique to the head-to-tail junction in a
  simulated dataset.
- Monitoring mode (truth VCF provided) achieves precision 1.0 at all VAF
  levels with zero false positives.

### What they do not show

- **Real deletion/insertion spectra.** The benchmark embeds exactly one 100 bp
  deletion, one 100 bp tandem duplication, and one 100 bp inversion. Real
  panels may have deletions of 50–500 bp, some of which are clinically
  characterised and some of which are not. The benchmark does not stress-test
  the classifier across the size range.
- **Variants near the 50 bp threshold.** A 48 bp deletion goes through the
  INDEL path; a 52 bp deletion goes through the large-variant path with
  relaxed thresholds. This boundary is untested.
- **Duplex evidence.** varforge produced only simplex molecules at all VAF
  levels. All detections in the benchmark are simplex-based. The duplex
  evidence scoring path (`min_variant_specific_duplex`) is not exercised.
- **Real somatic large indels.** The titration dataset (benchmarks 01 and 07)
  has a truth set of 375 SNVs and indels, of which the indels are up to ~20 bp.
  No large deletions or tandem duplications are in the titration truth set.
  The SV benchmarks are entirely synthetic.

---

## Practical implication: long INDELs in discovery mode

On a real panel sample in discovery mode, any variant with a 50 bp or larger
length difference will be called as a large variant type and will use the
relaxed thresholds. This means:

- A 60 bp germline deletion in the target region will appear as `LargeDeletion`
  at high confidence and high VAF. It will not be filtered by tumour-informed
  mode unless it is in the target VCF.
- A 55 bp tandem repeat expansion will appear as `TandemDuplication`.

If the intent is to call only SNVs and short indels, the `--max-sv-size` flag
(if implemented) or post-call filtering on `SVTYPE` in the VCF can suppress
large-variant calls.

In tumour-informed monitoring mode, this is not a concern: large variants not
in the target VCF receive `NotTargeted` and do not inflate false positive counts.

---

## Summary

The "SV" label in these benchmarks covers two mechanistically different things:

1. **Long INDELs** (LargeDeletion, NovelInsertion, TandemDuplication): detected
   by the standard de Bruijn graph path walk, same as short INDELs. The 50 bp
   threshold triggers adjusted evidence scoring and relaxed confidence thresholds,
   not a different detection mode. These could equally be called "large INDELs"
   or "long-range INDELs" without loss of precision.

2. **True structural rearrangements** (Inversion, InvDel, Fusion): require
   junction k-mer injection because their breakpoints create k-mers absent from
   the panel target sequences. Detection would fail without the `--sv-junctions`
   or `--fusion-targets` input.

The synthetic benchmarks deliberately include both classes. The 100 bp deletion
and 100 bp tandem duplication in benchmark 02 are testing whether kam's standard
path walk scales to longer events. The 100 bp inversion is testing whether
junction injection works. These are separate questions bundled under a shared
"SV benchmark" label, which is worth being explicit about when interpreting the
results.
