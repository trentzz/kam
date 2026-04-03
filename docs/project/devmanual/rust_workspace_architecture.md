# Rust Workspace Architecture

## Workspace Layout

```
kam/
├── Cargo.toml           (workspace)
├── kam-core/            (shared types: Molecule, KmerEvidence, etc.)
├── kam-assemble/        (molecule assembler — publishable standalone)
├── kam-index/           (k-mer indexer — publishable standalone)
├── kam-pathfind/        (de Bruijn graph + path walker)
├── kam-call/            (statistical caller)
└── kam/                 (the integrated binary that wires everything together)
```

## Why Workspace, Not Separate Tools

A workspace shares core data structures across crates while allowing independent publication. The internal representation is as rich as needed because workspace crates share internal types. When running the full pipeline via the `kam` binary, data passes between stages as Rust structs with zero serialisation cost — no information loss on the hot path.

## Crate Responsibilities

### kam-core
- Shared types: `Molecule`, `ConsensusRead`, `MoleculeEvidence`, `KmerIndex` trait
- Test data generator interface (`SyntheticDataset`)
- Error types (using `thiserror`)
- Chemistry presets (Twist, IDT, generic)
- **Most important design decision in the project.** Getting these types right means everything downstream builds on solid ground.

### kam-assemble (replaces HUMID)
- Molecule assembly from raw FASTQ
- Inline UMI extraction for Twist `5M2S+T` read structure
- Canonical UMI pairing and strand determination
- Hamming distance-based fuzzy UMI clustering (directional: high-count absorbs low-count)
- Endpoint fingerprinting for UMI collision detection (hash of first+last 8bp of template R1/R2; split families when fingerprints incompatible despite UMI match)
- Quality-weighted consensus calling per strand (weighted majority vote using Phred→error probability)
- Cross-strand duplex consensus (positions where fwd+rev agree get error probability ~10⁻⁷)
- UMI collision probability estimation per family: `P(collision) = 1 - (1 - 1/1024)^(n_molecules_at_position - 1)`
- Skip base validation and QC reporting
- **Memory strategy:** sort-then-group for panels (in-memory sort), external sort or hash-partition for larger inputs. See `docs/research/streaming_molecule_assembly.md`
- **Avoids HUMID's stack overflow:** iterative processing, no recursive trie traversal
- Outputs: annotated FASTQ with enriched headers (MI, RX, FS, DS, SF, SR, FT, CQ, CP, SK tags). See `docs/research/output_format_specs.md`
- **Highest reuse component** — anyone doing duplex UMI preprocessing benefits

### kam-index (replaces Jellyfish for this use case)
- Takes annotated FASTQ (from kam-assemble or standard FASTQ)
- Builds `KmerIndex` with `MoleculeEvidence` per k-mer (~48 bytes/entry vs Jellyfish's ~8 bytes)
- **Allowlist filtering from target sequences** — pre-compute target k-mers, only count those. Reduces memory from coverage-dependent to panel-size-dependent. 100,000× on 2Mb panel uses <500MB.
- **Target padding:** extend allowlist by k-1 bases on each side of target to avoid boundary gaps in downstream graph (see `docs/research/graph_building_challenges.md`, bridge problem)
- 2-bit k-mer encoding (A=00, C=01, G=10, T=11), packed into u64 for k≤31. For k>31 use u128.
- Canonical form: `min(kmer, reverse_complement(kmer))` but **store strand of observation separately** before canonicalising
- Can also emit Jellyfish-compatible `kmer\tcount` TSV for interop
- Multiple backend implementations via `KmerIndex` trait:
  - `HashKmerIndex` — fast, ~48 bytes/entry, good for panels
  - `SortedKmerIndex` — build in HashMap, compact to sorted Vec for queries (32 bytes/entry)
  - `FilteredKmerIndex<T>` — wraps any backend, only stores allowlisted k-mers
  - Future: `PersistentKmerIndex` (disk-backed via rocksdb), `BloomFilteredIndex` (two-pass singleton elimination for exome/WGS)
- See `docs/research/kmer_memory_strategies.md` for scaling analysis

### kam-pathfind (replaces km graph walking)
- De Bruijn graph construction from k-mer index (nodes = k-mers, edges = k-1 overlap)
- **Targeted mode:** anchored path walking between known flanking k-mers
  - Validate anchor k-mer uniqueness: if anchor appears >100× in index, flag as repeat region — path walking will produce spurious results. See `docs/research/graph_building_challenges.md`
  - Target capture zone must be padded by k-1 on each side (handled by kam-index allowlist)
  - Well-suited to SNV/indel detection; SVs spanning multiple targets need different strategy
- **De novo mode (Phase 2):** seed-and-extend with count-min sketch for global first pass, exact local graph building for candidate regions
- Path scoring using `MoleculeEvidence`: **how many independent duplex molecules support the variant path?** (not raw k-mer counts)
- `VariantGraph` trait abstracts targeted vs de novo implementations
- k configurable (default 31, support 21–41 for 150bp reads). Multi-k agreement scoring as future enhancement.

### kam-call (new — no existing equivalent)
- Statistical variant calling with confidence. See `docs/research/statistical_calling_models.md`
- **Phase 1:** Binomial model with molecule counts. VAF = k/M, 95% credible interval from Beta(k+1, M-k+1) posterior. Fast, simple, provides baseline.
- **Phase 2:** Beta-binomial model with per-site background from normal samples. Accounts for overdispersion from sequence context, damage artefacts (OxoG, deamination).
- **Duplex-aware scoring:** weight molecules by quality tier (duplex >> simplex ≥3 reads >> simplex 1-2 reads >> singleton)
- **Strand bias detection:** Fisher's exact test on 2×2 table (variant/wildtype × fwd/rev strand). Flag if p < 0.01.
- Per-site background error model: pre-computed from normal/germline samples as Beta(α, β) parameters per target position
- Targeted mode: test pre-specified positions against background
- De novo mode: find all supported variant paths above posterior threshold
- Output: TSV (primary, km-compatible), VCF (standard with custom INFO/FORMAT), JSON (full evidence for morning report). See `docs/research/output_format_specs.md`
- Filter labels: PASS, STRAND_BIAS, LOW_CONFIDENCE, LOW_DUPLEX, COLLISION_RISK

### kam (integrated binary)
- Wires all crates together with zero-copy data passing
- CLI via `clap` with subcommands
- Nextflow-friendly: structured QC JSON output per stage, deterministic for caching
- Chemistry presets (e.g., `--chemistry twist-umi-duplex`)

## Key Rust Crates (Dependencies)

- `needletail` — streaming FASTQ parsing (very fast)
- `dashmap` — concurrent HashMap for parallel k-mer indexing
- `rayon` — parallelism for molecule assembly and k-mer extraction
- `serde` + `bincode` — fast serialisation of k-mer index
- `clap` — CLI argument parsing
- `statrs` — beta distribution, binomial for statistical model
- `log` + `env_logger` — structured logging
- `thiserror` — error types
- `ahash` — fast hashing

## Build Order

1. **kam-core** — define data types (most important decision, ~1 day)
2. **kam-assemble** — molecule assembler (most self-contained, clearest correctness criteria)
3. **kam-index** — k-mer indexer (straightforward once consensus reads exist)
4. **kam-pathfind + kam-call** — together (tight interface between them). Start targeted only.
5. **kam** — integration binary (wire together, clean CLI)
