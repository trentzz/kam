# Twist Biosciences UMI Duplex Chemistry

## Read Structure

Twist UMI adapters use read structure `5M2S+T` for both R1 and R2 (fgbio notation):

- `5M` — 5 bases of molecular barcode (the UMI)
- `2S` — 2 bases to skip (monotemplate/spacer sequence, not biological)
- `+T` — all remaining bases are template (genomic sequence)

```
R1: [UMI_A 5bp][SKIP 2bp][genomic sequence...]
R2: [UMI_B 5bp][SKIP 2bp][genomic sequence...]
```

## UMI Tags

- `ZA` — 5' UMI (from R1 first 5bp)
- `ZB` — 3' UMI (from R2 first 5bp)
- `RX` — string containing both UMIs separated by a dash

## Duplex Pairing

R1 and R2 carry UMIs from opposite strands of the same original cfDNA molecule:
- Read pair 1: ZA=ACGTA, ZB=TGCAT → molecule X, forward strand
- Read pair 2: ZA=TGCAT, ZB=ACGTA → molecule X, reverse strand (same molecule!)

Canonical form: `min(ZA+ZB, ZB+ZA)` lexicographically → groups both strand orientations.

## Twist vs IDT: Critical Differences

| Property | Twist | IDT |
|----------|-------|-----|
| UMI length | 5bp random | 8bp fixed |
| UMI diversity | 4^5 = 1,024 | From known set |
| Error correction | Hamming distance only (no whitelist) | Whitelist-based correction possible |
| Collision risk | Significant at high depth | Much lower |

## The 5bp Random UMI Collision Problem

With only 1,024 possible UMIs, collision probability at high depth is non-trivial. At 10,000× coverage, two different molecules at the same genomic position can share the same UMI by chance. This is a source of false negatives (genuinely different molecules collapsed into one family).

**Mitigation:** Use endpoint fingerprinting as a secondary grouping signal. Two reads with matching canonical UMI pairs but very different template endpoint sequences are likely UMI collisions at different genomic positions — split them into separate families.

## Skip Bases as QC Signal

The 2bp skip sequence is NOT random — it's a fixed monotemplate from the ligation chemistry, consistent within an adapter lot. Uses:
- Validate skip base sequence matches expected value
- Flag reads where skip quality is below threshold
- Use skip sequence consistency across a read family as additional confidence signal

## Concrete Parsing Logic

```
Position 0-4:  UMI (5bp)     → extract, validate quality
Position 5-6:  Skip (2bp)    → extract, validate against expected sequence
Position 7+:   Template      → genomic sequence for downstream analysis
```

Minimum read length: 7 + MIN_TEMPLATE_LEN (need UMI + skip + enough template for k-mers)

## Grouping Strategy

1. **Primary key:** canonical_umi_pair with Hamming distance ≤ 1 tolerance
2. **Secondary collision filter:** within each primary group, sub-group by endpoint_fingerprint (hash of first+last 8bp of template R1 and R2)
3. **Family classification:**
   - n_fwd > 0 && n_rev > 0 → duplex family (highest confidence)
   - n_fwd > 0 && n_rev == 0 → simplex forward only
   - n_fwd == 0 && n_rev > 0 → simplex reverse only
   - Single read total → singleton (typically discard or flag)

## Configurable Parameters

- Minimum UMI base quality (Phred threshold) — tradeoff between sensitivity and specificity
- Hamming distance tolerance for UMI grouping (default: 1)
- Endpoint fingerprint edit distance threshold for collision detection

## References

Twist Biosciences (2024) Twist UMI Adapter System. Product documentation. Available at: https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system (Accessed: March 2026). [Primary source for all read structure details: `5M2S+T` notation, 5bp random UMI, 2bp skip/spacer, `ZA`/`ZB`/`RX` tag names, duplex strand-pairing logic, and the 4^5 = 1,024 UMI diversity figure.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of the fgbio read structure notation (`5M2S+T`), the `RX` dash-separated UMI format convention, and the canonical-UMI lexicographic min() pairing convention used in `GroupReadsByUmi`.]

Smith, T., Heger, A. and Sudbery, I. (2017) 'UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy', Genome Research, 27(3), pp. 491–499. doi: 10.1101/gr.209601.116. [Background reference for Hamming distance tolerance in UMI grouping and the general context of sequencing-error-tolerant UMI matching.]
