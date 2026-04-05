# Streaming and Memory-Bounded Molecule Assembly

## The Problem

At 10,000× coverage on a ctDNA panel, millions of reads need to be grouped into molecule families. HUMID's trie can stack overflow. fgbio loads everything into memory grouped by genomic coordinate (but requires alignment). kam-assemble needs a strategy that is both alignment-free and memory-bounded.

## Why This Is Hard Without Alignment

Alignment-based tools (fgbio) have a huge advantage: they group reads by genomic coordinate first, then by UMI within each coordinate. Since reads at different genomic positions can't be from the same molecule, this partitions the problem into small, independent groups.

Without alignment, you don't know which reads come from the same genomic region. Two reads with the same canonical UMI could be:
- From the same molecule (correct grouping)
- UMI collisions at different genomic positions (incorrect grouping)

The endpoint fingerprint helps distinguish these, but you still need to compare every read pair against every other read pair with the same UMI — which requires holding all reads with a given UMI in memory simultaneously.

## Memory Analysis

### Worst Case
With 1,024 possible 5bp UMIs and 10,000× coverage:
- ~10,000 read pairs per UMI bucket (assuming uniform distribution)
- Each read pair: ~300bp sequence + ~300bp quality + metadata ≈ 700 bytes
- Per UMI bucket: ~7MB
- All buckets simultaneously: ~7GB

This is manageable for a single sample on a modern machine but becomes problematic if:
- Multiple samples processed simultaneously
- Very deep sequencing (>50,000×)
- Whole-exome or whole-genome scale

### Best Case (With Endpoint Partitioning)
If we use the endpoint fingerprint to sub-partition UMI groups:
- Each (UMI, endpoint) group is a molecule family
- Typical family size: 5-20 read pairs
- Memory per family: <15KB
- Number of families in memory at once depends on processing order

## Strategy 1: Sort-Then-Group (Recommended for Panels)

### Approach
1. **Pass 1:** Stream all read pairs, extract (canonical_umi, endpoint_fingerprint) as a sort key, write to a temporary sorted file
2. **Pass 2:** Read the sorted file sequentially; all reads with the same sort key are adjacent and can be grouped with constant memory per group

### Implementation
```rust
// Phase 1: Extract keys and sort
struct SortableRead {
    canonical_umi: CanonicalUmiPair,    // primary sort key
    endpoint_fingerprint: u64,           // secondary sort key
    r1_seq: Vec<u8>,
    r2_seq: Vec<u8>,
    r1_qual: Vec<u8>,
    r2_qual: Vec<u8>,
    strand: Strand,
}

// External sort if data exceeds RAM
// Internal sort (Vec::sort) for panel-sized data
```

### Memory Bound
- Panel data (< ~50M read pairs): sort in memory, ~35GB RAM worst case but typically <8GB
- Larger data: external sort using temporary files, bounded to configurable RAM limit

### Trade-off
- Two passes over the data (or one pass + sort)
- Temporary disk usage proportional to input size
- Fully deterministic and reproducible

## Strategy 2: Hash-Partition (Best for Large Inputs)

### Approach
1. Hash the canonical UMI into N buckets (e.g., N=256)
2. Stream reads into bucket files on disk
3. Process each bucket independently (fits in memory)

### Memory Bound
Each bucket holds ~1/256th of the data. Process one bucket at a time.
- 10,000× panel: ~28MB per bucket
- 100,000× exome: ~2.8GB per bucket (still manageable)

### Implementation
```rust
fn partition_reads(input: &Path, n_buckets: usize) -> Vec<PathBuf> {
    let bucket_writers: Vec<BufWriter<File>> = (0..n_buckets)
        .map(|i| BufWriter::new(File::create(format!("bucket_{i}.tmp")).unwrap()))
        .collect();

    for read_pair in fastq_stream(input) {
        let parsed = parse_read_pair(&read_pair);
        let bucket = hash(&parsed.canonical_umi) % n_buckets;
        serialize_into(&mut bucket_writers[bucket], &parsed);
    }
    // Return bucket file paths
}
```

### Trade-off
- Single pass to partition, then per-bucket processing
- Disk I/O overhead (write all data twice)
- Embarrassingly parallelisable (process buckets with rayon)

## Strategy 3: Streaming with Bounded Buffer (For When Memory Is Truly Tight)

### Approach
Maintain a bounded LRU cache of active UMI groups. When a new read arrives:
1. Check if its canonical UMI is in the cache → add to existing group
2. If not in cache and cache is full → evict oldest group, emit its consensus, add new group
3. At end of stream → emit all remaining groups

### Risk
If reads from the same molecule are spread far apart in the input file, the group may be evicted and re-created, splitting the family. This produces incorrect results unless input is approximately sorted by UMI.

### When to Use
Only when input is already approximately sorted (e.g., by flowcell tile, which clusters reads from the same physical location). Not recommended as the default strategy.

## Recommendation for kam-assemble

### Default: Strategy 1 (Sort-Then-Group)
- For panel data: in-memory sort (fast, simple, ~8GB for typical panel at 10,000×)
- For larger data: external sort with configurable memory limit
- Deterministic, reproducible, Nextflow cache-friendly

### Configuration
```rust
pub enum AssemblyStrategy {
    /// Sort all reads in memory, then group. Fast for panels.
    InMemorySort,
    /// External sort with bounded memory. For large inputs.
    ExternalSort { max_memory_gb: usize },
    /// Hash-partition to disk, process buckets independently. For very large inputs.
    HashPartition { n_buckets: usize },
}
```

### Decision Logic
- Input < 10M read pairs → `InMemorySort`
- Input 10M–100M read pairs → `ExternalSort { max_memory_gb: 16 }`
- Input > 100M read pairs → `HashPartition { n_buckets: 256 }`

Auto-detect from input file size, or let user override.

## Interaction with Hamming Distance Clustering

Important: Hamming distance UMI clustering (grouping UMIs within edit distance 1) happens **after** the sort/partition step, within each UMI bucket. The sort key is the exact canonical UMI, but after loading a bucket, you check for nearby UMIs that should be merged.

For the hash-partition strategy, UMIs within Hamming distance 1 of each other may land in different buckets. Solution: use locality-sensitive hashing that maps similar UMIs to the same bucket, or process adjacent buckets together. For the sort strategy, nearby UMIs are adjacent in sorted order (since they differ by ≤1 character), so this is handled naturally.

## References

Langedijk, R. et al. (no formal publication date) HUMID: fast reference-free duplicate removal for FASTQ files with UMIs. GitHub repository. Available at: https://github.com/fengsong77/HUMID (Accessed: March 2026). [Source of the reference-free FASTQ deduplication approach and the stack overflow limitation in recursive trie traversal that motivates the memory-bounded strategies explored here.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of the genomic-coordinate-grouped in-memory approach critiqued as the alignment-dependent alternative; alignment requirement identified as the key difference from the streaming strategies here.]

Smith, T., Heger, A. and Sudbery, I. (2017) 'UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy', Genome Research, 27(3), pp. 491–499. doi: 10.1101/gr.209601.116. [Background on UMI grouping by genomic coordinate and the standard assumption that reads at different positions cannot share a molecule origin.]

Twist Biosciences (2024) Twist UMI Adapter System. Product documentation. Available at: https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system (Accessed: March 2026). [Source of the 5bp UMI with 1,024 possible values used in worst-case memory analysis (1,024 UMI buckets at 10,000× depth).]
