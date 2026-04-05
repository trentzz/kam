# Endpoint Fingerprinting for UMI Collision Detection

## The Problem

Twist's 5bp random UMIs have only 1,024 possible values. At high sequencing depths, different original DNA molecules can receive the same UMI pair by chance — a UMI collision. Without detection, these collisions merge unrelated molecules into the same family, corrupting the consensus and potentially creating false variants or masking real ones.

Alignment-based tools (fgbio) solve this trivially: reads with the same UMI but different mapping coordinates are from different molecules. In an alignment-free pipeline, we need a different signal.

## The Solution: Template Endpoint Fingerprinting

Different DNA molecules are different fragments of the genome. Even if they overlap, they almost certainly start and end at different positions (cfDNA fragments have a distribution of lengths and breakpoints). By comparing the beginning and end of the template sequences from R1 and R2, we can detect when two read pairs with the same UMI are from different fragments.

### What Goes Into the Fingerprint

For a paired-end read pair, there are four "endpoints" — sequences near the physical ends of the original DNA fragment:

```
Original cfDNA fragment:
5'─────────────────────────────────3'
   ←── R2 template ───|skip|UMI]     (R2 reads the reverse strand)
   [UMI|skip|── R1 template ──→      (R1 reads the forward strand)

Endpoints used:
  R1 first N bp  = fragment start (forward)
  R1 last N bp   = where R1 reaches toward fragment middle
  R2 first N bp  = fragment end (reverse complement)
  R2 last N bp   = where R2 reaches toward fragment middle
```

The first bases of R1 and R2 are the most informative — they correspond to the physical breakpoints of the cfDNA fragment, which are essentially random. Two different fragments are overwhelmingly likely to differ here.

### 2-Bit Encoding

Each base is encoded as 2 bits: A=00, C=01, G=10, T=11. N is treated as A (00). This packs bases tightly into an integer for fast comparison.

## Tradeoff Analysis: How Many Bases Per Endpoint

The fingerprint is packed into a `u64` (64 bits). With 4 endpoints and 2 bits per base, the budget is:

```
4 endpoints × N bases × 2 bits/base ≤ 64 bits
N ≤ 8 bases per endpoint
```

### Option A: 4 bases per endpoint (32 bits total)

```
Fingerprint: [R1_first_4bp][R1_last_4bp][R2_first_4bp][R2_last_4bp]
             8 bits        8 bits       8 bits        8 bits = 32 bits
```

**Discriminating power:** 4^16 = ~4.3 billion possible fingerprints. Sounds like a lot, but at any given UMI group the comparison is between fragments at the same locus. Fragments starting within 4bp of each other (common at highly covered positions) will share the first-4bp of R1, reducing effective discrimination to the other 3 endpoints.

**Pros:**
- Only uses 32 of 64 bits — remaining 32 bits could store other metadata
- Very tolerant of sequencing errors (fewer bits to go wrong)

**Cons:**
- At a hotspot with 1,000 molecules, some unrelated fragments will share 4bp endpoints by chance
- Estimated collision rate within a UMI group: non-negligible at depths >5,000×

### Option B: 8 bases per endpoint (64 bits total) — RECOMMENDED

```
Fingerprint: [R1_first_8bp][R1_last_8bp][R2_first_8bp][R2_last_8bp]
             16 bits       16 bits      16 bits       16 bits = 64 bits
```

**Discriminating power:** 4^32 = ~1.8 × 10^19 possible fingerprints. Two random fragments sharing all 32 bases is astronomically unlikely.

**Pros:**
- Exact encoding, no hashing needed — comparison is a single `u64 == u64` or XOR + popcount
- Extremely low false-merge rate even at very high depth
- Fast: integer comparison, no string operations

**Cons:**
- Uses the full 64 bits — no room for extra metadata
- More sensitive to sequencing errors (more bits that could differ between reads from the same molecule)

### Option C: 12 bases per endpoint (96 bits — needs hashing)

```
Raw: 4 endpoints × 12bp × 2 bits = 96 bits → doesn't fit in u64
Must hash down to 64 bits
```

**Pros:**
- Maximum discrimination

**Cons:**
- Hashing is lossy — two different fingerprints can hash to the same value (hash collision on top of UMI collision — confusing to debug)
- Hashing is slower than direct encoding
- Marginal benefit over 8bp — the discriminating power of 8bp is already far beyond what's needed

### Comparison Table

| Bases/endpoint | Total bits | Fits u64 | Discrimination | Error tolerance | Recommended |
|---------------|------------|----------|----------------|-----------------|-------------|
| 4 | 32 | Yes (spare room) | Moderate | High | No — too coarse at high depth |
| 6 | 48 | Yes (spare room) | Good | Good | Acceptable alternative |
| 8 | 64 | Exactly | Excellent | Moderate | **Yes — best balance** |
| 12 | 96 | No (needs hash) | Overkill | Lower | No — unnecessary complexity |

## Tradeoff Analysis: Bit Difference Threshold

Reads from the same molecule won't have perfectly identical fingerprints because of sequencing errors in the endpoint bases. The bit difference threshold controls how many bits can differ before two fingerprints are considered "incompatible" (different molecules).

### How Sequencing Errors Affect Fingerprints

A single base sequencing error flips 1-2 bits in the fingerprint (e.g., A→C is 00→01, one bit; A→T is 00→11, two bits). With 32 bases in the fingerprint and a ~1% per-base error rate (Q20):

- Expected errors in 32 bases: ~0.32 (usually 0, sometimes 1)
- Each error flips 1-2 bits
- Expected bit differences between two reads of the same molecule: 0-4 bits (from 0-2 errors in each read)

### Threshold Options

**Strict: 2 bit differences**

- Tolerates ~1 base error across all endpoints
- Very few false merges — almost any two different fragments will differ by more than 2 bits
- Risk: splits real families when 2+ endpoint bases have errors
- Estimated family splitting rate: ~5% at Q20 quality

**Moderate: 4 bit differences — RECOMMENDED**

- Tolerates ~2 base errors across all endpoints
- Low false merge rate — two random fragments sharing 30/32 bases is extremely unlikely
- Handles the common case of 0-1 sequencing errors per read comfortably
- Estimated family splitting rate: <1% at Q20 quality

**Lenient: 6 bit differences**

- Tolerates ~3 base errors
- Keeps families together even with poor quality endpoint bases
- Risk: fragments that share a common start position (same restriction site, same cfDNA breakpoint hotspot) may be falsely merged
- Relevant at nucleosome positioning hotspots where cfDNA fragments preferentially break at specific positions

**Very lenient: 8 bit differences**

- Tolerates ~4 base errors
- Essentially only catches completely unrelated fragments
- At this threshold, the fingerprint is barely discriminating

### Threshold Comparison

| Max bit diffs | Base errors tolerated | Family split rate (Q20) | False merge risk | Recommended |
|--------------|----------------------|------------------------|-----------------|-------------|
| 2 | ~1 | ~5% | Very low | No — too aggressive |
| 4 | ~2 | <1% | Low | **Yes — best balance** |
| 6 | ~3 | <0.1% | Moderate at hotspots | Acceptable for low-quality data |
| 8 | ~4 | <0.01% | Higher | No — too permissive |

## Edge Cases

### Short Templates

If a template is shorter than 8bp, we can't extract 8bp from the start and 8bp from the end without overlapping. Strategy:

- Template < 16bp: use the full template for that read's portion (first half + second half), pad remaining bits with zeros
- Template < 8bp: use all bases, pad heavily — the fingerprint is less discriminating but that's a short fragment that's probably low-quality anyway
- Template = 0bp: fingerprint is all zeros for that read — rely on the other read's contribution

### N Bases

Treat N as A (00). This means an N→A error is invisible in the fingerprint, but N bases are already low quality and typically caught by the UMI quality filter upstream.

### Identical Fingerprints from Different Molecules

Even with 8bp/4-bit-diff, two different molecules can have matching fingerprints if they happen to start and end at the same positions (within 2bp). This happens when:

- cfDNA fragments break at the same nucleosome boundary (real biological phenomenon — cfDNA has preferred breakpoints at ~10bp periodicity around nucleosomes)
- Restriction enzyme digestion creates fixed endpoints

This is inherently undetectable from sequence alone — it's the same limitation that alignment-based tools face when two fragments have identical start/end coordinates. The UMI is the only signal that distinguishes them, and if the UMI also collides, they're indistinguishable. The collision probability estimation (reported in QC) accounts for this.

## Recommendation

**Default: 8 bases per endpoint, 4 bit difference threshold.** Configurable via:

```
--fingerprint-bases 8       # bases per endpoint (4, 6, or 8)
--fingerprint-max-diff 4    # max bit differences for "compatible" (0-8)
```

This is the right default because:
1. 8bp exactly fills a u64 — no hashing, no wasted space, single-instruction comparison
2. 4 bit diffs tolerates realistic sequencing error rates while catching genuine collisions
3. At Twist's 1,024 UMI diversity, collision detection genuinely matters at depths >2,000×
4. The false-merge and false-split rates are both <1% at Q20 quality

Users who know they have very high quality data can tighten to `--fingerprint-max-diff 2`. Users with poor quality or FFPE data can loosen to `--fingerprint-max-diff 6`.

## References

Twist Biosciences (2024) Twist UMI Adapter System. Product documentation. Twist Biosciences. Available at: https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system (Accessed: March 2026). [Source of the 5bp random UMI design with 1,024 possible values and the 2bp skip/spacer structure `5M2S+T`.]

Snyder, M.W., Kircher, M., Hill, A.J., Daza, R.M. and Shendure, J. (2016) 'Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin', Cell, 164(1–2), pp. 57–68. doi: 10.1016/j.cell.2015.11.050. [Background on cfDNA fragment length distributions and preferred breakpoints at nucleosome boundaries, relevant to the discussion of biological fragmentation and endpoint discrimination.]

Illumina (2021) Quality Scores for Next-Generation Sequencing. Technical Note. Available at: https://www.illumina.com/content/dam/illumina/gcs/inherited-assets/working-with-us/support/documents/illumina_sequencing_introduction.pdf. [General source for Phred quality score definitions (Q20 = 1% error rate, Q60 as practical measurement ceiling) used in bit-difference threshold analysis.]
