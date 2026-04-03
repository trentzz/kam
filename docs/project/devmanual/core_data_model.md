# Core Data Model

The atomic unit throughout the pipeline is the **Molecule**, not the Read. Everything flows from molecule identity.

## Primary Types

### Molecule

```rust
pub struct Molecule {
    pub id: u64,                           // hash of canonical UMI pair
    pub umi_fwd: [u8; 5],                  // forward strand UMI (Twist: 5bp)
    pub umi_rev: [u8; 5],                  // reverse strand UMI
    pub fwd_reads: Vec<ReadPair>,          // all reads from forward strand
    pub rev_reads: Vec<ReadPair>,          // all reads from reverse strand
    pub consensus_fwd: Option<ConsensusRead>,
    pub consensus_rev: Option<ConsensusRead>,
    pub duplex_consensus: Option<ConsensusRead>,  // both strands agree
}
```

### ConsensusRead

```rust
pub struct ConsensusRead {
    pub sequence: Vec<u8>,
    pub per_base_quality: Vec<f32>,           // probability of error (not just Phred)
    pub per_base_strand_support: Vec<(u8, u8)>, // (fwd_count, rev_count) per position
    pub family_size: (u8, u8),                // (fwd_reads, rev_reads)
}
```

### MoleculeEvidence (per k-mer)

```rust
pub struct MoleculeEvidence {
    pub n_molecules: u32,          // how many distinct molecules carry this kmer
    pub n_duplex: u32,             // of those, how many have duplex support
    pub n_simplex_fwd: u32,
    pub n_simplex_rev: u32,
    pub min_base_error_prob: f32,  // best quality observation
    pub mean_base_error_prob: f32,
}
```

**Critical distinction:** This counts **molecules, not reads**. A k-mer in 100 PCR duplicates of one molecule = n_molecules=1. A k-mer in 100 independent molecules = n_molecules=100.

### TwistReadPair (chemistry-specific input)

```rust
pub struct TwistReadPair {
    pub umi_r1: [u8; 5],
    pub umi_r2: [u8; 5],
    pub skip_r1: [u8; 2],       // should be consistent per adapter lot
    pub skip_r2: [u8; 2],
    pub template_r1: Vec<u8>,   // genomic sequence after skip
    pub template_r2: Vec<u8>,
    pub qual_r1: Vec<u8>,
    pub qual_r2: Vec<u8>,
    pub canonical_umi: CanonicalUmiPair,
    pub strand: Strand,
    pub endpoint_fingerprint: u64,  // hash of first+last 8bp of template
}
```

### CanonicalUmiPair

```rust
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CanonicalUmiPair {
    pub umi_a: [u8; 5],  // lexicographically smaller
    pub umi_b: [u8; 5],  // lexicographically larger
}
```

### TargetSpec

Two modes of targeted detection: specific known variants, or entire regions to scan.

```rust
pub enum TargetSpec {
    /// Check for a specific known variant (SNV or indel)
    SpecificVariant {
        id: String,                    // e.g., "TP53_R248W"
        ref_sequence: Vec<u8>,         // reference sequence around the variant
        alt_sequence: Vec<u8>,         // expected variant sequence
        chrom: Option<String>,         // genomic coordinate (optional, for VCF output)
        start_pos: Option<u64>,
        padding: usize,               // k-1 bases on each side (set by kam-index)
    },
    /// Scan a region for any variants (de novo within a defined region)
    Region {
        id: String,                    // e.g., "TP53_exon7"
        ref_sequence: Vec<u8>,         // wildtype reference for the full region
        chrom: Option<String>,
        start_pos: Option<u64>,
        padding: usize,
    },
}
```

`SpecificVariant` indexes only k-mers relevant to the known variant. `Region` indexes ALL k-mers in the region and finds any deviations from the reference path — more memory, but finds unexpected variants.

### ScoredPath

```rust
pub struct ScoredPath {
    pub sequence: Vec<u8>,             // the path sequence through the graph
    pub variant_type: VariantType,     // SNV, INS, DEL, MNV, COMPLEX
    pub molecule_evidence: MoleculeEvidence, // aggregate evidence across all k-mers in path
    pub weakest_kmer_evidence: MoleculeEvidence, // evidence at the least-supported k-mer
    pub path_length: usize,            // number of k-mers in the path
}

pub enum VariantType {
    Reference,
    Snv,
    Insertion,
    Deletion,
    Mnv,
    Complex,
}
```

### AnchorValidation

```rust
pub struct AnchorValidation {
    pub start_unique: bool,           // anchor k-mer at start appears < threshold times
    pub end_unique: bool,
    pub start_count: u32,             // actual count in index
    pub end_count: u32,
    pub warning: Option<String>,      // human-readable warning if non-unique
}
```

## Key Traits

### KmerIndex

```rust
pub trait KmerIndex {
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);
    fn contains(&self, kmer: u64) -> bool;
    fn molecule_count(&self, kmer: u64) -> u32;
}
```

### VariantGraph

```rust
pub trait VariantGraph {
    fn walk_paths(&self, start: u64, end: u64, max_length: usize) -> Vec<ScoredPath>;
    fn validate_anchors(&self, target: &TargetSequence) -> AnchorValidation;
}
```

### QueryMode

```rust
pub enum QueryMode {
    Targeted {
        targets: Vec<TargetSequence>,
    },
    Discovery {
        min_molecule_support: u32,
        min_duplex_fraction: f32,
        background_model: BackgroundModel,
    },
    Hybrid {
        targets: Vec<TargetSequence>,
        discovery_params: DiscoveryParams,
    },
}
```

### AssemblyStrategy (kam-assemble memory management)

```rust
pub enum AssemblyStrategy {
    /// Sort all reads in memory, then group. Fast for panels (<10M read pairs).
    InMemorySort,
    /// External sort with bounded memory. For medium inputs (10M-100M read pairs).
    ExternalSort { max_memory_gb: usize },
    /// Hash-partition to disk, process buckets independently. For large inputs (>100M read pairs).
    HashPartition { n_buckets: usize },
}
```

## Design Principles

1. **Internal representation is always maximally rich** — never simplified for interface reasons
2. **Integrated binary passes data as Rust structs** — no serialisation on the hot path
3. **Each crate testable independently** with synthetic data
4. **Statistical model centralised** in kam-call, tested independently
5. **Molecule is the atomic unit** — all counts, all grouping, all evidence is molecule-level

## Information Preservation Contract

At split points between crates, the key output formats (see `docs/research/output_format_specs.md` for full specs):

- **After molecule assembly:** Annotated FASTQ headers with SAM-style tags (MI, RX, FS, DS, SF, SR, FT, CQ, CP, SK). Universally FASTQ-compatible; extra tags ignored by tools that don't understand them.
- **After k-mer indexing:** Serialised KmerIndex (`bincode` binary) + optional Jellyfish-compatible `kmer\tcount` TSV for interop
- **After graph walking:** `Vec<ScoredPath>` with full provenance per path (molecule counts, duplex counts, weakest k-mer evidence)
- **After variant calling:** TSV (km-compatible), VCF (standard with custom fields), JSON (full evidence for reporting). Each stage also emits a QC JSON validated by the next stage.
