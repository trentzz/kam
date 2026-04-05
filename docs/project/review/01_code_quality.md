# Code Quality Review

## Overall Impression

The Rust code is clean, idiomatic, and remarkably consistent in style across all 6 crates. This reads like it was written by someone who knows Rust well (or by an AI with a very strict style guide — and the CLAUDE.md confirms the latter). The consistency is a genuine strength: naming conventions, module layout, doc comment format, and test structure are identical across every file.

## Strengths

### Documentation Density
Every public function, struct, and enum has a doc comment with at least one `///` example. Many have full `# Example` blocks that compile as doctests. This is unusually thorough for a pre-1.0 project and sets a high bar.

### Error Handling Strategy
The error types are well-designed. `KamError` in `kam-core` uses `thiserror` with specific variants (`ReadTooShort`, `LowQualityUmi`, `Io`). Library functions return `Result` types. The `ParseResult` enum in the parser is a smart design — it's not `Result<T, E>` because drops aren't errors, they're expected outcomes with structured reasons.

### Derive Usage
Good use of `#[derive]` throughout. Serde derives on core types enable bincode serialization. `Clone`, `Debug`, `PartialEq` are applied appropriately. No excessive derives cluttering up types that don't need them.

### Module Organization
Each crate's `lib.rs` is a clean re-export with a module-level doc table. The `pub mod` declarations are minimal. Internal helpers are private. The `use` imports are organized and don't pull in too much.

## Issues

### Critical: `expect()` in Library Code

The project's own CLAUDE.md says "No `unwrap()` in library code." Yet the parser violates this:

```rust
// kam-assemble/src/parser.rs:212-218
let umi_r1: [u8; 5] = r1_seq[..umi_len].try_into().expect("umi_length==5");
let umi_r2: [u8; 5] = r2_seq[..umi_len].try_into().expect("umi_length==5");
let skip_r1: [u8; 2] = r1_seq[umi_len..tmpl_start].try_into().expect("skip_length==2");
let skip_r2: [u8; 2] = r2_seq[umi_len..tmpl_start].try_into().expect("skip_length==2");
let umi_qual_r1: [u8; 5] = r1_qual[..umi_len].try_into().expect("umi_length==5");
let umi_qual_r2: [u8; 5] = r2_qual[..umi_len].try_into().expect("umi_length==5");
```

These `expect()` calls are equivalent to `unwrap()`. They're safe *only* because the ReadTooShort check above guarantees the slices are long enough. But if someone creates a `ReadStructure` with `umi_length != 5`, these will panic at runtime. The `expect()` messages even betray the assumption: `"umi_length==5"`.

**Fix**: Either make these fallible (return a `DropReason`), or add a `ReadStructure::validate()` that rejects non-5 UMIs upfront with a proper error.

### The `try_into` Panic Deeper

The `umi_qual_r1` and `umi_qual_r2` extractions assume `r1_qual` has the same length as `r1_seq`. If a malformed FASTQ has a quality line shorter than the sequence line, this will panic with an unhelpful "umi_length==5" message that points to the wrong problem. The `needletail` library doesn't guarantee seq/qual length equality.

### Inconsistent Error Boxing

Most functions return `Result<T, Box<dyn std::error::Error>>`. This is fine for a CLI tool but is lazy for a library. The `KamError` enum exists and could be used:

```rust
// kam-core/src/serialize.rs:95
pub fn write_bincode<T: Serialize>(
    path: &Path,
    file_type: FileType,
    records: &[T],
) -> Result<(), Box<dyn std::error::Error>>  // Should be KamError
```

Every inter-crate boundary uses `Box<dyn Error>` instead of the project's own `KamError`. This makes error handling downstream (e.g., in Nextflow) harder because you can't match on specific error types.

### KamError is Anemic

`KamError` only has 3 variants: `ReadTooShort`, `LowQualityUmi`, `Io`. A real pipeline needs many more:
- File format errors (bad FASTA, bad bincode magic)
- Configuration errors (invalid k, unsupported chemistry)
- Empty input errors
- Stage-specific errors (no paths found, no anchors)

Currently all of these are stringly-typed via `Box<dyn Error>`.

### Dead Code Pattern in Assembler

```rust
// kam-assemble/src/assembler.rs:326-328
// Suppress unused variable warnings — fwd_cr and
// rev_cr are used only to confirm both SSCs succeeded.
let _ = (fwd_cr, rev_cr);
```

This is a code smell. The match arms bind `fwd_cr` and `rev_cr` but then don't use them, requiring a suppression. The match structure should be refactored so the bindings aren't needed.

### Redundant SSC Calls in Assembler

```rust
// kam-assemble/src/assembler.rs:308-309
let fwd_bases = call_ssc_bases(&fwd_reads, config);
let rev_bases = call_ssc_bases(&rev_reads, config);
```

`call_ssc()` is called first to produce `ConsensusRead`, then `call_ssc_bases()` is called *again* on the same reads to produce `ConsensusBase` for duplex consensus. This is doing the expensive consensus work twice. The fix is to have `call_ssc()` return both `ConsensusRead` and the intermediate `Vec<ConsensusBase>`.

### DefaultHasher for Molecule IDs

```rust
// kam-assemble/src/assembler.rs:450-454
fn hash_umi_pair(pair: &CanonicalUmiPair) -> u64 {
    let mut hasher = DefaultHasher::new();
    pair.hash(&mut hasher);
    hasher.finish()
}
```

`DefaultHasher` is **not** guaranteed to be stable across Rust versions. The CLAUDE.md says "All RNG must use seeded deterministic generators for Nextflow cache compatibility." While `DefaultHasher` isn't RNG, the molecule IDs derived from it will change between rustc versions, breaking Nextflow caching. Use `ahash` with a fixed seed, or a simple deterministic hash function.

### Floating Point Equality in Tests

```rust
// kam-assemble/src/consensus.rs:511
assert!((duplex[0].error_prob - expected).abs() < 1e-12,
```

Using `1e-12` as a float comparison tolerance is too tight for `f32` arithmetic. `f32` has only ~7 decimal digits of precision. This test passes today because the specific values happen to be exact in IEEE 754, but it's fragile. Use `1e-6` for `f32`.

### Namespace Collision: `ParseResult` vs `std::result::Result`

The parser defines its own `ParseResult` enum with `Ok` and `Dropped` variants. Having a variant called `Ok` that shadows `Result::Ok` is confusing:

```rust
match result {
    ParseResult::Ok(p) => { ... }   // Not std::result::Result::Ok
    ParseResult::Dropped { .. } => { ... }
}
```

This is not a bug, but `Parsed` and `Filtered` would be less confusing names.

## Style Observations

- Line lengths are generally well-controlled (under 100 chars)
- Comment style is consistent (`//` for line comments, `///` for doc, `//!` for module)
- Test names follow a clear descriptive pattern
- No magic numbers without explanation
- Good use of section comments (`// ── Step 1: ... ──────`) for visual structure in long functions

## Verdict

The code quality is genuinely good for a pre-1.0 project. The issues above are real but fixable. The biggest systemic problem is the `Box<dyn Error>` pattern — once you need to compose this with other tools or handle errors programmatically, you'll wish everything was `KamError`.
