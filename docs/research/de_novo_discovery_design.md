# De Novo Variant Discovery: Design Exploration

## Status

**Not being implemented in the first phase.** This document captures the design space for future implementation. The initial pipeline focuses on targeted detection (specific SNVs, indels, and regions).

## What De Novo Means in This Context

Targeted detection: "Is mutation X present at position Y?" — you provide a list of known variants or regions and kam checks for them.

De novo discovery: "What mutations exist in this sample?" — kam finds variants without prior knowledge of what to look for. This is critical for:
- Tumour-naive liquid biopsy (no prior tissue biopsy to define a variant set)
- Early cancer detection screening
- Discovering resistance mutations that emerged under treatment
- Clonal haematopoiesis of indeterminate potential (CHIP) detection

## Why De Novo Is Fundamentally Harder

### The Information Problem

In targeted mode, you know which k-mers to look for. The allowlist approach means you only index k-mers relevant to your targets — memory is bounded by panel size, not sequencing depth.

In de novo mode, you don't know which k-mers are interesting until you've seen them all. You need some form of global k-mer inventory before you can identify deviations from the expected (wildtype) signal. This brings back the memory problem that the allowlist solved.

### The Statistical Problem

In targeted mode, you're testing a specific hypothesis at a specific position. The multiple testing burden is small (number of targets).

In de novo mode, you're testing every position in the genome (or at least every position covered by your panel). At 0.1% VAF with a 2Mb panel, you're making ~2 million tests. The false positive rate must be extremely low per-test to keep the genome-wide rate acceptable. This requires either a very strong statistical model or a very conservative evidence threshold.

### The Graph Problem

In targeted mode, the de Bruijn graph is small (one subgraph per target, ~200bp each). Path walking is fast and bounded.

In de novo mode, you need to build a much larger graph and find all "interesting" branches — places where the graph forks from the expected path. This is closer to genome assembly than variant detection, and it inherits all the complexity of assembly: repeat regions, graph bubbles from heterozygous germline variants, systematic error hotspots creating false branches.

## Approach 1: Panel-Aware De Novo

If you have a panel of target regions (but not specific variants within those regions), you can do a focused de novo search:

1. Index all k-mers from reads that overlap your panel regions (similar to targeted, but you include ALL k-mers from those regions, not just k-mers from known variants)
2. Build a de Bruijn graph per panel region using the wildtype reference as the "expected" path
3. Find all branches in the graph that deviate from the expected path with sufficient molecule support
4. Score each branch using the MoleculeEvidence model

This is essentially targeted mode with "the target is the whole region, not a specific variant." The memory cost is higher (you're indexing all k-mers in the panel, not just those around known variants) but still bounded by panel size. For a 2Mb panel with k=31, that's ~2 million target k-mers plus variant k-mers observed — manageable.

```rust
enum TargetSpec {
    /// Specific variant: check if this exact variant is present
    SpecificVariant {
        target_id: String,
        ref_sequence: Vec<u8>,
        alt_sequence: Vec<u8>,
    },
    /// Region: find any variants within this region
    Region {
        target_id: String,
        ref_sequence: Vec<u8>,  // wildtype reference for the full region
        /// All k-mers from this region are indexed, not just variant-specific ones
    },
}
```

**This is the most practical near-term extension of targeted mode** and could be offered as `--discovery-mode panel` where the panel BED file defines regions to search within.

### Memory Estimate

- 2Mb panel, all k-mers: ~2M reference k-mers + observed variants
- At 48 bytes/k-mer: ~100MB for the reference k-mers, plus observed variant k-mers
- Total: likely <500MB — fine

### Limitations

- Still requires a panel definition (can't find variants outside defined regions)
- Requires a reference sequence for the panel regions to define "expected" paths
- This is NOT truly alignment-free — you're using a reference to define the expected graph. But you're not aligning reads to it; you're comparing k-mer paths.

## Approach 2: Reference-Guided Whole-Panel Graph

Build the complete de Bruijn graph for all reads in the sample, but prune it using a reference to identify the expected backbone:

1. **Phase 1 — Sketch:** Stream all reads through a count-min sketch (~1.2GB for 1B k-mers at 1% FPR). This gives approximate k-mer counts globally.
2. **Phase 2 — Identify candidates:** Find k-mers where the count-min sketch count is unexpectedly low relative to neighbouring k-mers (potential variant-containing regions) or where coverage suddenly changes.
3. **Phase 3 — Local assembly:** For each candidate region, do exact local graph building by re-streaming reads that contain the candidate k-mers. Build a local de Bruijn subgraph and walk it.
4. **Phase 4 — Score and call:** Use MoleculeEvidence from the local exact index to score variant paths.

This is a two-level approach: approximate global → exact local. Memory is bounded by the sketch size (fixed) plus the local subgraph size (small).

### Implementation Sketch

```rust
struct DeNovoDiscovery {
    /// Fixed-size approximate global k-mer counter
    global_sketch: CountMinSketch,

    /// Reference k-mer set for the panel (or genome) — defines "expected"
    reference_kmers: HashSet<u64>,

    /// Candidate regions identified by sketch analysis
    candidates: Vec<CandidateRegion>,
}

struct CandidateRegion {
    /// Anchor k-mers flanking the candidate variant
    start_anchor: u64,
    end_anchor: u64,

    /// Why this region was flagged
    reason: CandidateReason,

    /// Local exact k-mer index (built in phase 3)
    local_index: Option<HashKmerIndex>,
}

enum CandidateReason {
    /// k-mer count drop relative to neighbours (potential variant)
    CoverageDrop { expected: u32, observed: u32 },

    /// Novel k-mers present that aren't in the reference
    NovelKmers { count: u32 },

    /// Coverage asymmetry between strands
    StrandAsymmetry { fwd: u32, rev: u32 },
}
```

### Memory Estimate

- Count-min sketch: ~1.2GB (fixed, independent of input size)
- Reference k-mers for 2Mb panel: ~16MB
- Local subgraphs: ~10MB each, process one at a time
- Total peak: ~1.5GB — manageable

### Challenges

- **Sketch accuracy:** Count-min sketch can overestimate counts (never underestimate). This means some low-frequency variant k-mers may not be flagged as "unexpectedly low" because the sketch inflates their count. False negative rate depends on sketch size and hash functions.
- **Candidate explosion:** In noisy data, many regions will have apparent coverage drops due to random variation. Need a principled threshold for what constitutes "unexpectedly low."
- **Germline variants:** Common germline heterozygous variants will appear as graph branches. Need a mechanism to filter these (e.g., population frequency databases, matched normal sample).

## Approach 3: Streaming Graph Assembly

Build the graph incrementally as reads stream in, pruning low-evidence branches periodically:

1. Maintain a bounded-memory de Bruijn graph
2. As each consensus molecule is processed, add its k-mers to the graph
3. Periodically prune: remove branches with < N molecule support
4. At the end, all remaining branches are candidate variants

This is closest to how cortex_var and other assembly-based variant callers work.

### Challenges

- **Pruning threshold is critical:** too aggressive = miss real low-frequency variants. Too lenient = memory explodes.
- **Order-dependent:** different input orderings can produce different pruning results. Violates the determinism requirement for Nextflow caching.
- **Difficult to parallelize:** the graph is a shared mutable data structure.

**Not recommended for kam** due to determinism concerns and complexity.

## Approach 4: Two-Pass with Molecule-Aware Counting

Unique to kam because we have molecule-level identity:

1. **Pass 1:** Run the full molecule assembly (kam-assemble). Count distinct molecules per k-mer using a lightweight counter (just molecule count, not full MoleculeEvidence). This can use a simple HashMap<u64, u32> at ~12 bytes/k-mer.
2. **Pass 2:** Identify k-mers where molecule count is between a lower bound (e.g., 2 = not a singleton error) and an upper bound (e.g., < 50% of median coverage = not wildtype). These are candidate variant k-mers.
3. **Pass 3:** For candidate k-mers only, build the full MoleculeEvidence index. Construct local subgraphs and walk paths.

### Why Molecule Counting Helps

The key insight: most singleton k-mers (appearing in only 1 molecule) are sequencing errors. In raw read counts, error k-mers can appear many times (PCR duplicates of the same error), making them hard to distinguish from low-frequency variants. But in molecule counts, an error k-mer appears in exactly 1 molecule. A real variant at 0.1% VAF in 1,000 molecules appears in ~1 molecule — similar count, but the variant k-mer will also appear on both strands and in multiple overlapping k-mer positions, while the error won't.

### Memory Estimate for Pass 1

- All distinct k-mers in a 2Mb panel at 10,000×: ~50M k-mers
- At 12 bytes/k-mer (u64 key + u32 count): ~600MB
- After singleton removal: ~5M k-mers → 60MB
- This is the sweet spot: molecule-level deduplication before counting dramatically reduces the k-mer space

## Comparison of Approaches

| Approach | Memory | Passes | Complexity | Best For |
|----------|--------|--------|------------|----------|
| Panel-aware de novo | <500MB | 1 | Low | Known panel regions, unknown variants |
| Reference-guided sketch | ~1.5GB | 2-3 | Medium | Larger panels, unknown regions |
| Streaming assembly | Bounded but variable | 1 | High | Not recommended (determinism issues) |
| Two-pass molecule-aware | ~600MB peak | 2-3 | Medium | Best use of kam's duplex advantage |

## Recommended Implementation Path

### Near-term (extends targeted mode naturally)

**Panel-aware de novo (Approach 1):** Allow `TargetSpec::Region` alongside `TargetSpec::SpecificVariant`. Users provide a BED file of panel regions plus reference FASTA. kam indexes all k-mers in those regions and finds any deviations. This is a small extension of the existing targeted pipeline.

### Medium-term

**Two-pass molecule-aware (Approach 4):** Leverages kam's unique molecule-level counting. After kam-assemble produces consensus molecules, a lightweight first pass identifies candidate variant k-mers by molecule count, then a focused second pass builds full evidence for candidates only.

### Long-term

**Reference-guided sketch (Approach 2):** For whole-exome or larger panels where even Approach 4's first pass is too memory-intensive. The count-min sketch provides a fixed-memory global view.

## Interaction with Background Error Model

De novo discovery requires a background error model to distinguish true variants from systematic artefacts. In targeted mode, the background is estimated per-site from normal samples. In de novo mode, you need:

- A per-site background from a panel of normals (if available)
- Or a genome-wide background rate estimated from the data itself (e.g., the rate of non-reference k-mers across all positions)
- Known artefact signatures (OxoG, deamination) to filter

Without a good background model, de novo discovery at 0.1% VAF will have an unacceptable false positive rate. This is the primary reason to defer de novo to a later phase — the background model needs to be validated on targeted data first.

## References

Iqbal, Z., Caccamo, M., Turner, I., Flicek, P. and McVean, G. (2012) 'De novo assembly and genotyping of variants using colored de Bruijn graphs', Nature Genetics, 44(2), pp. 226–232. doi: 10.1038/ng.1028. [Foundational reference for de Bruijn graph-based variant detection without alignment. Source of the "coloured graph" concept where different samples contribute different colours, enabling variant discovery by graph topology.]

Narzisi, G., O'Rawe, J.A., Iossifov, I., Fang, H., Lee, Y.-h., Wang, Z., Wu, Y., Lyon, G.J., Wigler, M. and Schatz, M.C. (2014) 'Accurate de novo and transmitted indel detection in exome-capture data using microassembly', Nature Methods, 11(10), pp. 1033–1036. doi: 10.1038/nmeth.3069. [Source for the local assembly approach to variant detection: identify candidate regions, extract reads, build local graph, find variants. The "seed-and-extend" paradigm adapted here for k-mer space.]

Cormode, G. and Muthukrishnan, S. (2005) 'An improved data stream summary: the count-min sketch and its applications', Journal of Algorithms, 55(1), pp. 58–75. doi: 10.1016/j.jalgor.2003.12.001. [Foundational reference for the count-min sketch data structure used in Approach 2 for approximate global k-mer counting with bounded memory.]

Bloom, B.H. (1970) 'Space/time trade-offs in hash coding with allowable errors', Communications of the ACM, 13(7), pp. 422–426. doi: 10.1145/362686.362692. [Background reference for probabilistic data structures; the Bloom filter is an alternative to count-min sketch for the singleton elimination step in Approach 4.]

Marçais, G. and Kingsford, C. (2011) 'A fast, lock-free approach for efficient parallel counting of occurrences of k-mers', Bioinformatics, 27(6), pp. 764–770. doi: 10.1093/bioinformatics/btr011. [Source for the ~50M distinct k-mers estimate for a 2Mb panel at 10,000× coverage, derived from Jellyfish's characterisation of k-mer distributions in sequencing data.]
