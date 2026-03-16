# Graph Building Challenges

## The Core Problem: Graph Connectivity

A de Bruijn graph has an edge between two k-mers if they overlap by k-1 bases. If you only load k-mers from target regions, you're building a **subgraph** of the full de Bruijn graph. The question is whether that subgraph is sufficient for path walking.

## Targeted Mode Problems

### Problem 1: Anchor K-mer Contamination

How targeted graph walking works:
1. Take reference sequence around variant of interest (~200bp)
2. Decompose into k-mers — anchors at start and end
3. Walk paths from start anchor to end anchor
4. Report which paths have evidence

**Issue:** K-mers at edges may be non-unique in the genome. Reads from other genomic locations bleed into the target graph through shared anchor k-mers, producing spurious variant paths.

**Fix:** Extend target capture zone so anchor k-mers are unique. For k=31, most genomic k-mers are unique. Target loading should flag non-unique anchors:

```rust
fn validate_anchor_uniqueness(target, index, k) -> AnchorValidation {
    // If anchor appears >100 times, it's in a repeat region
    // and path walking will produce garbage
}
```

### Problem 2: Missing Bridge Problem

When only loading target-region k-mers, the k-1 bases at window boundaries straddle inside/outside the target. Strict allowlists discard these boundary-straddling k-mers, creating gaps.

**Fix:** Pad target capture zone by k-1 bases on each side:

```rust
fn build_allowlist(target, k) -> HashSet<u64> {
    let padded = target.extend(k - 1);  // k-1 on each side
    padded.kmers(k).map(|km| km.canonical()).collect()
}
```

Without this, variants near target window edges will be missed.

### Problem 3: Complex Variants Spanning Target Windows

Large deletions, structural rearrangements, or tandem duplications can create junction k-mers drawing sequence from two different target entries. If those are in separate allowlists, the junction k-mer is never seen.

**Scope:** For SNVs and small indels (common ctDNA use case), this doesn't matter. For SVs, it matters a lot. The allowlist approach is well-suited to SNV/indel targeted detection and needs a different strategy for SV detection.

## De Novo Mode: The Deeper Problem

### The Connectivity Problem

A variant at 1% VAF in a 150bp read generates ~120 mutant k-mers forming a linear path. To find this path de novo, you need to traverse from wildtype to mutant k-mers at the divergence point. But if you've only loaded predefined region k-mers, you either already know where to look (targeted mode) or need all k-mers (memory problem).

### Approach 1: Chromosome-by-Chromosome Partitioning
Partition reads by likely chromosome using minimizers or lightweight sketch. Build graphs one chromosome at a time. ~23× memory reduction. Risk: sketch partitioning errors cause connectivity loss.

### Approach 2: Seed-and-Extend
1. First pass: find k-mers at unexpectedly low counts relative to neighbours (candidate variant k-mers)
2. For each seed, load local neighbourhood from reads
3. Build local subgraph and walk

Use count-min sketch for global first pass (fixed memory, approximate counts), find seeds, then do exact local graph building.

### Approach 3: Two-Level Indexing
Sparse global index tracking only k-mers above coverage threshold (≥5 molecule support). Gives graph skeleton for finding anchors. Where adjacent nodes have unexpectedly different coverage, do local re-query against raw reads for rare k-mers.

## K-mer Size Considerations

| k | Pros | Cons |
|---|------|------|
| 21 (small) | More connections, denser graph | Less unique, more collisions/spurious paths, more memory |
| 31 (default) | Good uniqueness for human genome, balanced | Standard choice for 150bp reads |
| 51 (large) | Very unique, fewer spurious paths | Needs longer reads (≥51bp), large indels invisible |

For 150bp Illumina reads with SNV/indel detection in ctDNA panels, **k=31 is the default**. Support k=21–41 range.

**Multi-k mode:** Run with multiple k values simultaneously, require consistency. A variant called at k=25 and k=31 is more reliable than one called at only one k. Costs proportionally more memory but fine for targeted panels with allowlist.

## Recommended Design

Separate graph implementations for targeted vs de novo:

```rust
pub trait VariantGraph {
    fn walk_paths(&self, start: u64, end: u64, max_length: usize) -> Vec<ScoredPath>;
}

struct TargetedGraph { ... }   // Pre-built from allowlist, fast, memory-bounded
struct DeNovoGraph { ... }     // Built lazily from raw reads, memory-intensive
```

Build targeted first. De novo uses same infrastructure but different construction strategy.

## References

Marçais, G. and Kingsford, C. (2011) 'A fast, lock-free approach for efficient parallel counting of occurrences of k-mers', Bioinformatics, 27(6), pp. 764–770. doi: 10.1093/bioinformatics/btr011. [Source for k-mer counting fundamentals and canonical k-mer representation; the `k=31` default for human-genome uniqueness is a community convention derived from experience with Jellyfish and similar tools.]

Bourgey, M. et al. (2019) 'km: find-mutation', GitHub repository and associated work at IRIC. Available at: https://github.com/iric-soft/km (Accessed: March 2026). [Source of the graph-walking approach over k-mer space to detect variant paths; the anchor k-mer and path-walking design pattern described here is adapted from km's approach.]

Iqbal, Z., Caccamo, M., Turner, I., Flicek, P. and McVean, G. (2012) 'De novo assembly and genotyping of variants using colored de Bruijn graphs', Nature Genetics, 44(2), pp. 226–232. doi: 10.1038/ng.1028. [Foundational reference for de Bruijn graph-based variant detection and multi-k consistency approaches.]

Cormen, T.H., Leiserson, C.E., Rivest, R.L. and Stein, C. (2009) Introduction to Algorithms. 3rd edn. Cambridge, MA: MIT Press. [General reference for count-min sketch and approximate data structure theory underlying the two-level indexing strategy.]
