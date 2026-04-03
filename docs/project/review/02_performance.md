# Performance Review

## TL;DR

The code is written for correctness, not performance. This is the right priority for a first implementation, but the current algorithms have O(n^2) bottlenecks that will make this unusable on real sequencing data (millions of read pairs). The pipeline also has zero parallelism despite listing `rayon` and `dashmap` as key dependencies.

## Algorithmic Complexity Issues

### 1. O(n^2) UMI Grouping (Critical)

```rust
// kam-assemble/src/assembler.rs:158-173
for (idx, pair) in read_pairs.iter().enumerate() {
    let pos = unique_umis
        .iter()
        .position(|(u, _)| u == &pair.canonical_umi);  // O(n) scan
    match pos {
        Some(i) => { unique_umis[i].1 += 1; read_to_umi[idx] = i; }
        None => { unique_umis.push(...); read_to_umi[idx] = i; }
    }
}
```

For every read pair, this does a linear scan through all unique UMIs seen so far. With 10M read pairs and ~100K unique UMIs, that's 10M * 100K = 10^12 comparisons. This should be a `HashMap<CanonicalUmiPair, usize>` lookup, which would be O(1) amortized.

**Impact**: This alone makes the assembler unusable on real data. A 10M read pair sample would take hours instead of seconds.

### 2. O(n^2) Hamming Clustering

```rust
// kam-assemble/src/clustering.rs:96-102
for (cluster_idx, (seed_idx, _members)) in clusters.iter().enumerate() {
    let dist = umi_pair_hamming_distance(&pairs[*seed_idx].0, &pairs[idx].0);
    if dist <= max_distance {
        found = Some(cluster_idx);
        break;
    }
}
```

For each UMI pair, every existing cluster seed is checked. With 100K unique UMIs and O(K) clusters, this is O(N*K) per UMI. For 100K UMIs, that's potentially billions of Hamming distance calculations.

**Better approach**: Build a BK-tree or use locality-sensitive hashing (LSH) for Hamming distance queries. For a 10-base sequence with distance 1, you can enumerate all 30 possible 1-mutation neighbors and use a hash set lookup — O(1) per query.

### 3. O(n*m) Fingerprint Splitting

```rust
// kam-assemble/src/assembler.rs:234-236
let found = groups
    .iter_mut()
    .find(|(rep_fp, _)| fingerprints_compatible(*rep_fp, fp));
```

For each read in a cluster, every existing fingerprint group is scanned. In pathological cases (many UMI collisions at a hot locus), this could be slow. In practice, most clusters have 1-2 groups, so this is O(n) per cluster. But the greedy approach may also produce suboptimal groupings if fingerprint A is compatible with both B and C, but B and C are not compatible with each other.

### 4. O(n) `contains()` in Walk Cycle Detection

```rust
// kam-pathfind/src/walk.rs:128
if partial.contains(&next) {  // O(path_length)
```

Cycle detection in BFS path walking uses `Vec::contains()`, which is O(path_length) per check. For the default `max_path_length = 500`, this means 500 comparisons per successor per node. A `HashSet<u64>` alongside the path vec would make this O(1).

### 5. O(n) `observed_count()` in FilteredKmerIndex

```rust
// kam-index/src/allowlist.rs:128
pub fn observed_count(&self) -> usize {
    self.allowlist.iter().filter(|&&k| self.inner.contains(k)).count()
}
```

This scans the entire allowlist every time it's called. For a panel with 100K target k-mers, this is 100K hash lookups. If called per molecule (which the QC tracking might do), it becomes expensive. Cache the count or maintain it incrementally.

## Memory Issues

### 6. All Data in Memory

The entire pipeline loads all read pairs into a `Vec<ParsedReadPair>` before assembly:

```rust
let (read_pairs, parse_stats) = read_fastq_pairs(&args.r1, &args.r2, &parser_config)?;
```

For a 10M read pair sample with 150bp reads, this is roughly:
- Per read pair: 2 * 150 bytes (template) + 2 * 150 bytes (quality) + 2 * 5 bytes (UMI) + overhead ≈ 700 bytes
- Total: 10M * 700 = 7 GB

This should stream: parse, group by UMI on the fly (using a HashMap), then process families. The `docs/research/streaming_molecule_assembly.md` discusses this but it's not implemented.

### 7. read_header() Reads Entire File

```rust
// kam-core/src/serialize.rs:138-139
let mut buf = Vec::new();
reader.read_to_end(&mut buf)?;
```

`read_header()` reads the entire bincode file into memory just to deserialize the header. The header is maybe 20 bytes. For a 500MB molecule file, this allocates 500MB for no reason.

**Fix**: Read a fixed-size buffer (e.g., 1024 bytes) and deserialize from that.

### 8. Cloning in extract_and_index

```rust
// kam-index/src/extract.rs:86
let existing = index.get(canon).cloned();
let mut ev = existing.unwrap_or_default();
// ... modify ev ...
index.insert(canon, ev);
```

Every k-mer observation does a clone + insert. The `KmerIndex` trait doesn't have an `entry()` API, so this pattern is forced. For millions of k-mers across thousands of molecules, the repeated cloning and reinserting is wasteful. `HashKmerIndex` has `entry()` but the `extract` module uses the trait object.

## Missing Parallelism

### 9. Zero Use of Rayon

`rayon` is listed in CLAUDE.md as a key dependency, but it's not used anywhere in the code. The workspace `Cargo.toml` doesn't even include it. Zero `par_iter()`, zero `rayon::scope()`, zero thread pool usage.

Obvious parallelism opportunities:
- FASTQ parsing (read pairs are independent)
- Per-UMI-family consensus calling (families are independent)
- Per-target pathfinding (targets are independent)
- K-mer extraction from molecules (molecules are independent)

### 10. Zero Use of DashMap

`dashmap` is listed as a key dependency but is not used. The `HashKmerIndex` uses `std::collections::HashMap`. For concurrent k-mer extraction (once rayon is added), `DashMap` would be the right choice.

## Benchmarks That Don't Exist

The `docs/benchmarking/` directory is empty. There are no benchmarks, no performance tests, no profiling results. For a project that claims to replace Jellyfish (which processes millions of k-mers per second), this is a significant gap. Even a simple `cargo bench` with `criterion` would be valuable.

## What Would Be Needed for Production Scale

For a typical ctDNA panel sample (50M read pairs, 500-gene panel):

1. **Streaming FASTQ processing**: Sort-free hash-partition approach from the research docs
2. **HashMap for UMI grouping**: O(1) lookup instead of O(n) scan
3. **Rayon for parallelism**: At minimum, per-target pathfinding and per-family consensus
4. **Gzip support**: Every real FASTQ is compressed (use `flate2` or `niffler`)
5. **Memory-bounded k-mer index**: Use the `FilteredKmerIndex` from the allowlist module
6. **Streaming consensus**: Process families as they complete, don't buffer all molecules

## Positive Performance Aspects

- The `KmerIterator` sliding window is O(1) per k-mer after the first — this is correct
- The 2-bit encoding is space-efficient and allows direct comparison
- The de Bruijn graph construction uses prefix grouping for O(N*K) edge finding (where K is the average number of same-prefix k-mers), which is reasonable
- The `BufWriter` and `BufReader` wrappers in the serialization module are correct
- Endpoint fingerprinting uses `popcount` (hardware-accelerated on x86) for compatibility checking
