# Extras: Documentation, Tooling, and Miscellaneous

## Documentation Quality

### Research Documentation: Exceptional

The `docs/research/` directory contains 13 documents that would pass peer review in a bioinformatics journal. Highlights:

- **twist_umi_chemistry.md**: Precise read structure diagram, collision probability math, canonical pairing logic. Could be the methods section of a paper.
- **consensus_calling_algorithms.md**: Three algorithms with mathematical derivations, pseudocode, and tradeoff analysis. Better than most tool papers.
- **statistical_calling_models.md**: Clear progression from binomial to beta-binomial with overdispersion rationale.
- **streaming_molecule_assembly.md**: Honest analysis of memory vs disk tradeoffs for sort-then-group vs hash-partition.

These documents demonstrate deep domain knowledge and would be valuable to anyone implementing a similar tool.

### Planning Documentation: Strong

The architecture docs accurately reflect the implemented codebase. The `core_data_model.md` matches the actual `kam-core` types. The `design_principles.md` is actually followed (transparency, granular control, no black boxes). The `development_workflow.md` describes the Claude Code loop that produced this codebase.

### Missing Documentation

- **No README.md** at the repo root. This is table stakes for any open-source project.
- **No `docs/manual/`** directory despite being referenced in CLAUDE.md. No user-facing documentation exists.
- **No CHANGELOG.md**. With 2 commits, this isn't urgent, but it should be established early.
- **No CONTRIBUTING.md**. The development workflow doc partially covers this, but it's developer-facing.

## Development Workflow

### Claude Code Task System

The `docs/claudeloop/` directory reveals a methodical AI-assisted development process:

- 28 completed tasks in `docs/claudeloop/done/`
- Tasks are granular and well-scoped (e.g., "canonical_umi_parsing", "hamming_clustering")
- Each task has acceptance criteria
- The open questions document (9 resolved Q&As) shows design decisions were made deliberately

This is an unusually disciplined approach to AI-assisted development. The task queue system acts as both project management and a record of design decisions.

### What the Workflow Produced

The entire 9,600-line codebase was produced in what appears to be 2 git commits. The first is "Initial commit" (likely empty or scaffolding), and the second is "Implement full kam pipeline" with everything. This suggests the code was written in a single extended session.

**Concern**: A single monolithic commit makes `git blame` and `git bisect` useless. The development process would benefit from per-task commits. This is acknowledged in the workflow docs but not practiced.

## CI/CD

**There is none.** No GitHub Actions, no GitLab CI, no Makefile, no pre-commit hooks. For a codebase that mandates `cargo test` and `cargo clippy -- -D warnings` passing, there should be CI enforcing these checks.

Minimum CI setup:
```yaml
# .github/workflows/ci.yml
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo clippy -- -D warnings
      - run: cargo test
```

## Dependency Analysis

### Production Dependencies

| Dependency | Version | Purpose | Assessment |
|---|---|---|---|
| needletail | 0.6 | FASTQ parsing | Good choice, maintained, fast |
| thiserror | 2 | Error derives | Standard, correct version |
| serde + bincode | 1 + 1 | Serialization | Good for internal format |
| serde_json | 1 | QC output | Standard |
| statrs | 0.18 | Beta distribution | Good, used for CDF |
| clap | 4 | CLI parsing | Standard, derive API is clean |
| log | 0.4 | Logging facade | Listed but not used anywhere |

### Listed But Not Used

| Dependency | CLAUDE.md Says | Actual Use |
|---|---|---|
| dashmap | concurrent HashMap | Not in any Cargo.toml |
| rayon | parallelism | Not in any Cargo.toml |

These are listed as "key dependencies" in CLAUDE.md but are not actually dependencies. This is misleading — they're aspirational, not actual.

### Missing Dependencies That Should Be Added

- **flate2** or **niffler**: For gzip FASTQ support
- **criterion**: For benchmarking
- **proptest**: For property-based testing
- **env_logger** or **tracing**: For actual logging (the `log` facade is imported but no logger is configured)

### Dev Dependencies

```
tempfile = used in run.rs tests
serde_json = used in QC tests
```

`tempfile` is a dev dependency in the `kam` crate but not in `kam-core` or `kam-assemble`, where it would also be useful. The `kam-assemble` io tests roll their own `TempDir` struct instead.

## Code Metrics

```
Crate         | Source Lines | Test Lines | Test Ratio
kam-core      |          340 |        474 |    1.39x
kam-assemble  |          896 |       1203 |    1.34x
kam-index     |          452 |        645 |    1.43x
kam-pathfind  |          474 |        648 |    1.37x
kam-call      |          485 |        466 |    0.96x
kam (binary)  |          842 |        597 |    0.71x
──────────────|──────────────|────────────|─────────
Total         |         3489 |       4033 |    1.16x
```

More test code than production code across the workspace. This is a healthy ratio. The CLI binary has the lowest test ratio, which is typical — integration test coverage is harder to achieve.

## Interesting Design Decisions

### 1. ParseResult Instead of Result

Using a custom `ParseResult { Ok, Dropped }` instead of `Result<ParsedReadPair, DropReason>` is unconventional but sensible. Drops are expected outcomes (not errors) and carry structured data for logging. The `?` operator wouldn't make sense here because you want to count drops, not propagate them.

### 2. Endpoint Fingerprinting for Collision Detection

This is creative and not seen in other UMI tools. The idea — that two reads from the same molecule will start and end at the same positions, while UMI collisions will have different endpoints — is biologically sound for cfDNA (which has characteristic fragmentation patterns).

### 3. Molecule Evidence Instead of Counts

Storing `MoleculeEvidence` (with strand breakdown, error stats) instead of integer counts is the core innovation. It's more expensive per k-mer but enables calling that would be impossible with counts alone (e.g., "this variant is supported by 5 duplex molecules with mean error prob 0.001").

### 4. Bincode with Magic Bytes

The custom bincode format with `b"KAM\0"` magic bytes and a typed `FileHeader` is good defensive design. It prevents accidentally feeding the wrong intermediate file to the wrong stage.

## Things That Surprised Me

1. **No `unsafe` code anywhere.** For a bioinformatics tool that processes large data, this is impressive. The 2-bit encoding and sliding window iterator are implemented entirely in safe Rust.

2. **No external k-mer library.** The k-mer encoding, canonical form, and sliding window are implemented from scratch. Libraries like `bio` or `debruijn` exist, but rolling your own gives full control over the encoding (important for the `MoleculeEvidence` integration).

3. **The assembler's `call_ssc_bases` always uses R1 templates for both strands.** This is documented with a comment about "orientation already encoded in strand assignment" but it's not obvious why this is correct (see correctness review).

4. **The VCF output uses POS=1 for everything.** This is a consequence of alignment-free calling but will break most VCF downstream tools.

5. **The entire codebase was apparently written in one session by Claude Code.** The consistency of style, the exhaustive documentation, and the mechanical test coverage all point to AI generation. The human contribution appears to be the research documents, the architecture decisions in CLAUDE.md, and the task queue — which is exactly the right division of labor.
