# HUMID: Analysis and Gaps

## What It Is

HUMID is a C++ tool for reference-free FASTQ deduplication using a trie (prefix tree) data structure rather than the network/graph methods used by UMI-tools. It's a small project (5 GitHub stars, 2 forks) from Leiden University Medical Center.

### Core Algorithm: The Trie

When a read arrives, its UMI (and optionally part of the read sequence) is traversed character-by-character down the trie. Deduplication uses Hamming distance lookups within the trie rather than pairwise O(n²) comparison.

- Trie initialized with ACGTN alphabet (A=0, C=1, etc.)
- Nodes sized to contain the character with highest index
- Fast-fail Hamming distance — gives up quickly when distance threshold exceeded
- Clustering methods (from fastqdedup, inspired by HUMID's trie):
  - `highest_count`: selects one read from a cluster
  - `adjacency`: from highest-count read, selects all within specified distance
  - `directional`: uses counts to determine PCR/sequencing artefact vs genuine difference

### What HUMID Can Do

- Reference-free deduplication from FASTQ (no BAM needed)
- UMI-aware dedup where UMIs are in the read header (via fastp or BCL Convert)
- Hamming distance-based fuzzy UMI matching for sequencing error correction
- Paired-end support
- Multiple clustering strategies

### What HUMID Cannot Do (Gaps for kam)

1. **No duplex UMI support** — treats UMIs as opaque strings, no concept of complementary strand pairing. Would wrongly split duplex pairs into two families or miss strand-pairing entirely.

2. **No consensus sequence generation** — picks one representative read, discards others. Does not collapse a family into a quality-weighted consensus sequence (the core operation for ctDNA error correction).

3. **No family size tracking/filtering** — doesn't output metadata about how many raw reads supported each consensus. fgbio encodes this in cD/cE SAM tags; HUMID has nothing equivalent.

4. **No quality-weighted decisions** — doesn't use base quality scores when selecting which read to keep from a duplicate cluster.

5. **No inline UMI extraction** — requires fastp to first move UMI from read sequence to header. No native inline UMI trimming + extraction.

6. **Stack overflow on large datasets** — recursive trie traversal can overflow the call stack with maximum clustering method on large datasets.

7. **No BAM output with RX/MI tags** — FASTQ-only output. Duplex workflows need BAM with RX (UMI) and MI (molecule ID) tags for downstream tools.

## References

Langedijk, R. et al. (no formal publication date) HUMID: fast reference-free duplicate removal for FASTQ files with UMIs. GitHub repository. Available at: https://github.com/fengsong77/HUMID (Accessed: March 2026). [Primary subject of this analysis. Source of the trie-based deduplication algorithm, ACGTN alphabet encoding, Hamming distance fast-fail traversal, and stack overflow limitation.]

fastqdedup (no formal publication date) fastqdedup: FASTQ deduplication tool inspired by HUMID. GitHub repository. Available at: https://github.com/LUMC/fastqdedup (Accessed: March 2026). [Source of the clustering method implementations — `highest_count`, `adjacency`, and `directional` — cited as inspired by HUMID's trie.]

Smith, T., Heger, A. and Sudbery, I. (2017) 'UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy', Genome Research, 27(3), pp. 491–499. doi: 10.1101/gr.209601.116. [Source of network/graph UMI deduplication methods and terminology (directional, adjacency) that HUMID's trie approach replaces.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of cD/cE SAM tag conventions for family size tracking and the duplex consensus calling pipeline that kam aims to replicate at the FASTQ level.]
