# Task 003: Twist UMI FASTQ Parser

## Location
`kam-assemble/src/parser.rs`

## What to implement

A FASTQ read pair parser that extracts UMI, skip bases, and template from Twist UMI duplex data following the `5M2S+T` read structure. This is the first stage of molecule assembly.

**Key design decision:** The parser is k-agnostic. It does NOT enforce a minimum template length by default. Reads with very short templates are parsed successfully — downstream stages (k-mer extraction) decide what to do with them. An optional user-configurable minimum template length can cause reads to be rejected, but this is not the default.

## Interface

```rust
use kam_core::{CanonicalUmiPair, Strand, ReadStructure, KamError};

/// Configuration for the parser
#[derive(Debug, Clone)]
pub struct ParserConfig {
    pub read_structure: ReadStructure,
    pub min_template_length: Option<usize>,  // None = accept any length (k-agnostic)
    pub min_umi_quality: Option<u8>,          // None = accept any quality
}

/// Why a read pair was dropped by the parser
#[derive(Debug, Clone, Copy)]
pub enum DropReason {
    ReadTooShort,        // shorter than umi+skip (malformed, can't even extract UMI)
    TemplateTooShort,    // template below user-configured minimum
    LowUmiQuality,      // UMI base quality below threshold
}

#[derive(Debug, Clone)]
pub struct ParsedReadPair {
    pub umi_r1: [u8; 5],
    pub umi_r2: [u8; 5],
    pub skip_r1: [u8; 2],
    pub skip_r2: [u8; 2],
    pub template_r1: Vec<u8>,
    pub template_r2: Vec<u8>,
    pub qual_r1: Vec<u8>,      // quality for template bases only
    pub qual_r2: Vec<u8>,      // quality for template bases only
    pub umi_qual_r1: [u8; 5],  // quality for UMI bases
    pub umi_qual_r2: [u8; 5],
    pub canonical_umi: CanonicalUmiPair,
    pub strand: Strand,
}

/// Parse result: either a successfully parsed read pair, or a drop with reason.
/// Both variants carry enough info for logging.
pub enum ParseResult {
    Ok(ParsedReadPair),
    Dropped {
        reason: DropReason,
        detail: String,       // human-readable detail (e.g., "r1_template=12bp,min=20")
    },
}

pub fn parse_read_pair(
    r1_seq: &[u8],
    r1_qual: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    config: &ParserConfig,
) -> ParseResult;
```

## Parsing logic

1. If read shorter than `umi_length + skip_length` (can't even extract UMI): return `Dropped(ReadTooShort)`
2. Extract UMI from first 5bp of each read
3. Extract skip bases from positions 5-6
4. Template is everything from position 7 onward (may be empty or very short — that's fine)
5. Quality arrays split the same way
6. If `min_template_length` set and template shorter: return `Dropped(TemplateTooShort)` with detail
7. If `min_umi_quality` set and any UMI base below threshold: return `Dropped(LowUmiQuality)` with detail
8. Otherwise: determine canonical UMI pair and strand, return `Ok(ParsedReadPair)`

## Logging integration

The caller (not this function) is responsible for logging drops. `ParseResult::Dropped` carries the `reason` and `detail` strings so the caller can write them to the drop log if enabled. The parser itself does no IO.

A counter struct should be maintained by the caller:

```rust
pub struct ParseStats {
    pub n_processed: u64,
    pub n_passed: u64,
    pub n_read_too_short: u64,
    pub n_template_too_short: u64,
    pub n_low_umi_quality: u64,
}
```

Counters are always incremented (cheap). Detailed per-read logging is only written when the log sink is enabled. See `docs/planning/logging_architecture.md`.

## Tests required

1. Valid read pair with long template parses correctly — all fields at right positions
2. Read with 0bp template (exactly umi+skip length) parses OK when no min_template_length set
3. Read shorter than umi+skip (e.g., 4bp) returns `Dropped(ReadTooShort)`
4. Read with 10bp template returns `Dropped(TemplateTooShort)` when min_template_length=20
5. Read with 10bp template returns `Ok` when min_template_length is None
6. R1 and R2 of same length produce correct template lengths
7. UMI extraction from positions 0..5
8. Skip extraction from positions 5..7
9. Template extraction starts at position 7
10. Quality arrays match template length (not full read length)
11. Canonical UMI and strand determined correctly
12. Low UMI quality base returns `Dropped(LowUmiQuality)` when threshold set
13. Low UMI quality base returns `Ok` when no threshold set
14. Detail string in `Dropped` contains useful info (template length, which base failed, etc.)

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
- Doc comments with examples on `parse_read_pair`
- `DropReason` and `ParseStats` are well-documented
