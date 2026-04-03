# Architecture Review

## Workspace Layout

The crate layout is clean and follows the planned architecture exactly:

```
kam-core       → shared types (Molecule, ConsensusRead, MoleculeEvidence)
kam-assemble   → FASTQ parsing, UMI grouping, consensus calling
kam-index      → k-mer encoding, indexing, allowlist filtering
kam-pathfind   → de Bruijn graph, path walking, scoring
kam-call       → statistical calling, output formatting
kam            → CLI binary, pipeline orchestration
```

Each crate has a clear single responsibility. Dependencies flow one direction: `core ← assemble ← index ← pathfind ← call ← kam`. There are no circular dependencies, and `kam-core` is the only shared dependency. This is exactly right.

## What Works Well

### 1. Molecule as Atomic Unit

This is the core architectural insight and it's correctly implemented. `MoleculeEvidence` tracks `n_molecules`, `n_duplex`, `n_simplex_fwd`, `n_simplex_rev` — not read counts. This means a k-mer seen in 100 PCR duplicates from one molecule counts as 1, not 100. Every downstream stage (graph, scoring, calling) operates on these molecule counts. This is the fundamental advantage over Jellyfish.

### 2. Staged Pipeline with Zero-Copy Hot Path

The `kam run` command passes data between stages as in-memory Rust structs — no serialization overhead. The individual subcommands (`assemble`, `index`, `pathfind`, `call`) support serialized handoff via bincode for Nextflow compatibility. This dual-mode design is good: fast for single-machine runs, distributed for HPC.

### 3. QC at Every Stage

Each stage produces a typed QC struct (`AssemblyQc`, `IndexQc`, `PathfindQc`, `CallQc`) serialized to JSON. This supports the Nextflow validation-between-stages pattern. The QC structs have specific, useful fields (not just "passed" flags).

### 4. Chemistry Abstraction

`ReadStructure` is parameterized, not hardcoded. The Twist UMI duplex preset is a factory method, and other chemistries could be added as presets. This is the right abstraction level.

### 5. KmerIndex Trait

The `KmerIndex` trait abstracts over storage backends:

```rust
pub trait KmerIndex {
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);
    fn contains(&self, kmer: u64) -> bool;
    fn molecule_count(&self, kmer: u64) -> u32;
}
```

`HashKmerIndex` implements it, and `FilteredKmerIndex<T>` wraps any `KmerIndex` with allowlist filtering. This is composable and could support future backends (mmap'd, Bloom-filtered, etc.).

## Issues

### 1. Hardcoded Array Sizes Contradict ReadStructure Genericity

The `ReadStructure` says "umi_length" and "skip_length" are configurable. But every type that stores UMI or skip data uses fixed arrays:

```rust
pub struct CanonicalUmiPair {
    pub umi_a: [u8; 5],   // Hardcoded to 5
    pub umi_b: [u8; 5],   // Hardcoded to 5
}

pub struct ParsedReadPair {
    pub umi_r1: [u8; 5],   // Hardcoded
    pub umi_r2: [u8; 5],   // Hardcoded
    pub skip_r1: [u8; 2],  // Hardcoded
    pub skip_r2: [u8; 2],  // Hardcoded
    ...
}
```

If you tried to use `ReadStructure { umi_length: 8, skip_length: 3 }`, the parser would panic at the `try_into().expect()` calls. The abstraction promises genericity but the implementation doesn't deliver.

**Options**:
- Make it genuinely generic (use `Vec<u8>` or `SmallVec<[u8; 8]>`)
- Or make it honestly specific (remove `ReadStructure` configurability, document that only Twist 5M2S is supported)

The second option is better for now — don't promise what you can't deliver.

### 2. KmerIndex Trait is Too Narrow

The trait lacks several operations that the implementation needs:

- **No `entry()` API**: Forces the clone-get-modify-insert pattern in `extract.rs`
- **No iteration**: `HashKmerIndex` has `iter()` but the trait doesn't. This means code that needs to iterate must downcast
- **No `len()`**: Same issue — `HashKmerIndex` has it, the trait doesn't
- **Mutable-only `insert()`**: Takes `MoleculeEvidence` by value, forcing clones. An `update_or_insert()` method would be more efficient

The trait was designed for reads but the usage pattern requires writes. Adding `entry()` and `iter()` to the trait would clean up several call sites.

### 3. Molecule Type is Underused

The `Molecule` struct is defined in `kam-core` with an `evidence: Option<MoleculeEvidence>` field. But this field is always `None` after assembly — the evidence is accumulated in the k-mer index, not on the molecule. The molecule's evidence field is never populated anywhere in the pipeline.

This suggests an architectural mismatch: either the field should be removed, or the design should be updated so that per-molecule evidence *is* tracked (which would be useful for per-molecule QC).

### 4. The `run.rs` Integration is a Monolith

The `kam run` command at 437 lines is the largest single function in the codebase. It mixes:
- Configuration parsing
- Stage orchestration
- QC generation
- Output formatting
- Helper functions

This should be factored into a pipeline runner that takes stage functions as parameters, allowing individual stages to be swapped or skipped. The current structure makes it hard to run just stages 3-4 (pathfind + call) on pre-built data.

### 5. No Trait for Stage Execution

Each stage is a standalone function with different signatures:
- `read_fastq_pairs(path, path, config) -> (Vec, Stats)`
- `assemble_molecules(Vec, config) -> (Vec, Stats)`
- `extract_all(reads, k, &mut index)`
- `walk_paths(graph, start, end, config) -> Vec`
- `call_variant(id, ref, alt, ref_seq, alt_seq, config) -> VariantCall`

There's no common `Stage` trait or pipeline abstraction. This is fine for now but will become painful when adding:
- Resume/checkpoint support
- Parallel stage execution
- Stage-level profiling
- Plugin stages

### 6. Bincode Format is Fragile

The serialization format uses bincode with Serde derives. Bincode is not self-describing — if you add a field to `Molecule`, old bincode files become unreadable. The `FileHeader` has a `version: u32` field, but there's no migration logic. The docs acknowledge this gap.

**Recommendation**: Either use a self-describing format (MessagePack, CBOR) or implement explicit version-aware deserialization.

### 7. No Trait Object Safety for KmerIndex

`KmerIndex` is used as `&dyn KmerIndex` in several places (`extract_and_index`, `DeBruijnGraph::from_index`, `score_path`). This works because the trait methods don't use generics. But `FilteredKmerIndex<T>` requires a concrete `T`, so you can't wrap a `&dyn KmerIndex` in a `FilteredKmerIndex`. This limits composability.

## Data Flow

The data flow through the pipeline is clean:

```
FASTQ R1/R2
    → ParsedReadPair (parser)
    → Molecule (assembler, with consensus)
    → ConsensusReadInfo (adapter)
    → KmerIndex entries (extract)
    → DeBruijnGraph (graph builder)
    → GraphPath (walker)
    → ScoredPath (scorer)
    → VariantCall (caller)
    → TSV/VCF/JSON (output)
```

Each arrow is a clear transformation with no backtracking. The types at each stage are well-defined. The only awkward adapter is `molecules_to_consensus_reads()` in the index command, which extracts `ConsensusReadInfo` from `Molecule` — this should probably be a method on `Molecule`.

## Recommendations

1. **Honest about UMI sizes**: Remove `ReadStructure` configurability or make it real. Don't have both
2. **Enrich KmerIndex trait**: Add `entry()`, `iter()`, `len()` to the trait
3. **Factor run.rs**: Extract stage orchestration into a pipeline runner
4. **Decide on Molecule.evidence**: Either populate it or remove it
5. **Version bincode properly**: Add migration logic or switch to a self-describing format
