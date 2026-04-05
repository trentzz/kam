# Seven Differentiated Value Areas for kam

These are the gaps where a Rust tool bridging HUMID's reference-free FASTQ approach with fgbio's duplex semantics would have no direct open-source competitor.

## 1. Native Duplex UMI Strand-Pairing

No open-source tool does this cleanly at the FASTQ level. The core algorithm: when you see a UMI pair A+B, canonicalise it (so A+B and B+A map to the same molecule), track which strand each read came from, and group them into a duplex family. Requires a different data model than simplex — grouping complementary pairs, not just identical/similar UMIs.

## 2. Configurable Per-Family Filtering with Transparency

fgbio's `min-reads` parameter is a blunt instrument. What's needed: per-family reporting of (total reads, forward strand reads, reverse strand reads, mean base quality, UMI distance from canonical), plus flexible filter expressions. Users should be able to say "keep duplex families with ≥3 reads on each strand" or "downgrade to simplex consensus if duplex family size < 2 but simplex ≥ 5." Currently one global threshold forces a blind tradeoff.

## 3. Quality-Weighted Consensus Calling at the FASTQ Level

Most tools either pick a representative read (HUMID, UMI-tools dedup mode) or require alignment first (fgbio). Quality-weighted consensus calling directly on FASTQ — using Phred scores per base position across the family — eliminates the alignment requirement for error correction.

## 4. Mixed Simplex/Duplex Output with Per-Read Provenance Tags

Encode rich per-consensus metadata: family size, strand support, consensus confidence score per base. Downstream variant callers can make better use of this information than existing tools provide.

## 5. UMI Collision Detection and Reporting

At high sequencing depths (≥10,000×), UMI space saturation becomes significant — two different molecules can get the same UMI by chance. None of the open tools model this or report collision probability estimates. For ctDNA at 0.1% VAF this genuinely matters, especially with Twist's 5bp random UMIs (only 1,024 possible values).

## 6. Reference-Free Operation with BAM Output

HUMID's key insight (reference-free, FASTQ-in) combined with fgbio's output model (BAM with RX/MI/cD tags) — nobody does both. A tool that takes FASTQ and emits tagged FASTQ or pre-alignment uBAM with full UMI metadata.

## 7. Streaming/Memory-Efficient Processing

HUMID's trie has the stack overflow issue. fgbio loads everything into memory grouped by genomic coordinate. For 10,000× deep ctDNA panels with millions of reads from a small target region, this is a real problem. A streaming, memory-bounded approach using sorted UMI buckets handles this better.

## References

Langedijk, R. et al. (no formal publication date) HUMID: fast reference-free duplicate removal for FASTQ files with UMIs. GitHub repository. Available at: https://github.com/fengsong77/HUMID (Accessed: March 2026). [Source of the reference-free FASTQ deduplication trie algorithm and the stack overflow limitation on large datasets.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of the `GroupReadsByUmi` / `CallDuplexConsensusReads` pipeline, `min-reads` parameter behaviour, and RX/MI/cD SAM tag conventions.]

Smith, T., Heger, A. and Sudbery, I. (2017) 'UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy', Genome Research, 27(3), pp. 491–499. doi: 10.1101/gr.209601.116. [Source of network-based UMI deduplication methods (directional, adjacency, cluster) referenced as prior art with no duplex support.]

Twist Biosciences (2024) Twist UMI Adapter System. Product documentation. Available at: https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system (Accessed: March 2026). [Source of the 5bp random UMI design, 1,024 UMI diversity figure, and the collision risk at high depth.]
