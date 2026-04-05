# K-mer Memory Strategies

## The Problem

Jellyfish uses ~30GB for WGS k-mer databases. kam's `MoleculeEvidence` struct is ~6× larger per entry than Jellyfish's simple integer count. Naive implementation would be untenable at scale.

## Memory Budget per K-mer

```
HashMap<u64, MoleculeEvidence>:
  8 bytes (key, packed k-mer) + 24 bytes (MoleculeEvidence) + ~16 bytes (hashmap overhead) = ~48 bytes/kmer

vs Jellyfish: ~8 bytes/kmer
```

### Scale implications:
- **Targeted panel (1Mb):** ~50M distinct k-mers × 48 bytes = ~2.4GB ✓
- **Exome (~50Mb):** Crossover point, needs different approach
- **WGS:** 3B k-mers × 48 bytes = 144GB ✗

## Solution 1: Allowlist Filtering (Primary for Targeted Mode)

Only count k-mers that appear in target regions. Pre-compute all k-mers from target sequences, use as allowlist during extraction.

```rust
let allowlist: HashSet<u64> = load_target_kmers(&targets, k);
for molecule in molecules {
    for kmer in molecule.kmers(k) {
        if allowlist.contains(&kmer.canonical()) {
            index.entry(kmer).or_default().add_evidence(&molecule);
        }
    }
}
```

Memory scales with panel size × k-mer diversity, NOT sequencing depth. 100,000× coverage on 2Mb panel uses <500MB.

**Critical:** Pad target capture zone by k-1 bases on each side to avoid boundary gaps in the graph.

## Solution 2: Bloom Filter Prefilter (For Larger Inputs)

Two-pass approach for exome or tumour-naive discovery:

- **Pass 1:** Stream all reads, insert every k-mer into a Bloom filter (~1.2GB for 1B k-mers at 1% FPR)
- **Pass 2:** Only insert into full HashMap if k-mer already in Bloom filter (filters singletons, which dominate raw reads)

Memory reduction: 10–50× because singleton error k-mers are eliminated.

## Solution 3: Compressed Representations

### Sorted array + binary search
After build phase (HashMap), serialize to `Vec<(u64, MoleculeEvidence)>` sorted. Binary search for queries. Removes ~16 bytes/entry hashmap overhead. Exactly 32 bytes × n_kmers.

### Bitpacked counters
Pack MoleculeEvidence tighter when max values are known:
```rust
struct PackedEvidence {
    counts: u64,    // n_molecules:20, n_duplex:20, n_simplex_fwd:12, n_simplex_rev:12
    error_probs: u32, // two f16s or quantised u16s
    _pad: u32,
}  // 16 bytes vs 24
```

### Count-min sketch (for de novo discovery)
Approximate counts with fixed bounded memory. Tradeoff: lose per-k-mer exact molecule provenance. Good for first-pass discovery followed by targeted validation.

## Solution 4: Memory-Mapped Disk-Backed Index

For WGS tumour-naive: build index on disk in sorted order, memory-map for querying. OS handles caching. Limited by disk, not RAM. Use `rocksdb` crate or custom sorted SSTable format.

## Recommended Architecture

Use a trait to abstract the storage backend:

```rust
pub trait KmerIndex {
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);
    fn contains(&self, kmer: u64) -> bool;
    fn molecule_count(&self, kmer: u64) -> u32;
}

// Concrete implementations:
// HashKmerIndex       — fast, memory-hungry, good for panels
// SortedKmerIndex     — slower insert, compact, good for queries after build
// PersistentKmerIndex — disk-backed, any scale, slower
// FilteredKmerIndex<T> — wraps any above, only stores target-region kmers
```

## Practical Default

For targeted ctDNA panels: `FilteredKmerIndex<HashKmerIndex>` — fast, memory-efficient, no configuration needed. Users doing larger work opt into other backends.

## References

Marçais, G. and Kingsford, C. (2011) 'A fast, lock-free approach for efficient parallel counting of occurrences of k-mers', Bioinformatics, 27(6), pp. 764–770. doi: 10.1093/bioinformatics/btr011. [Primary reference for Jellyfish. Source of the ~30GB WGS k-mer database figure and the ~8 bytes/k-mer memory footprint that is the baseline for comparison.]

Bloom, B.H. (1970) 'Space/time trade-offs in hash coding with allowable errors', Communications of the ACM, 13(7), pp. 422–426. doi: 10.1145/362686.362692. [Foundational reference for the Bloom filter prefilter strategy (Solution 2) used to eliminate singleton error k-mers in pass 1.]

Cormode, G. and Muthukrishnan, S. (2005) 'An improved data stream summary: the count-min sketch and its applications', Journal of Algorithms, 55(1), pp. 58–75. doi: 10.1016/j.jalgor.2003.12.001. [Foundational reference for the count-min sketch approximate counting approach described under Solution 3 for de novo discovery.]

RocksDB (2024) RocksDB: A persistent key-value store. Facebook/Meta Open Source. Available at: https://rocksdb.org (Accessed: March 2026). [Source tool cited for the memory-mapped disk-backed index (Solution 4) via the `rocksdb` Rust crate.]
