# CHEM-CONFIG: Generalise Chemistry to Arbitrary UMI Lengths

**Status**: todo
**Priority**: critical
**Branch**: epic/CHEM-CONFIG

## Goal

Remove all hardcoded `[u8; 5]` UMI and `[u8; 2]` skip arrays from kam-core,
the parser, and downstream crates. Replace them with `Vec<u8>` driven by a
`ReadStructure` config. After this epic, kam works with any UMI length and any
skip length, enabling use on non-Twist chemistries (e.g. simplex 12 bp UMI).

## Motivation

kam-core currently hardcodes 5 bp UMI and 2 bp skip throughout. Every struct,
parser, and test assumes Twist UMI duplex chemistry. This prevents using kam on
other panels (e.g. Illumina/IDT UDI with 12 bp UMI, or simplex protocols with
no skip). Making chemistry configurable is a prerequisite for the public dataset
benchmarks (PUB-BENCH epic) and for any eventual production use beyond the
Twist panel.

## Design

### Core type changes (APPROVED)

`CanonicalUmiPair` and `Molecule` in `kam-core`:
- `umi_a: [u8; 5]` → `umi_a: Vec<u8>`
- `umi_b: [u8; 5]` → `umi_b: Vec<u8>`
- Skip bytes removed from core types; skip is a parser concern only.

These are APPROVED changes. Any commit touching `kam-core` types must contain
"APPROVED" in the commit message.

### ReadStructure config

```rust
pub struct ReadStructure {
    pub umi_length: usize,
    pub skip_length: usize,
}

impl ReadStructure {
    pub fn twist_umi_duplex() -> Self { ... }   // 5bp UMI, 2bp skip
    pub fn simplex_12bp() -> Self { ... }       // 12bp UMI, 0bp skip
    pub fn simplex_9bp() -> Self { ... }        // 9bp UMI, 0bp skip
}
```

The parser reads `umi_length` and `skip_length` from `ReadStructure` instead
of using compile-time constants.

### Hamming clustering

The Hamming distance function already operates on byte slices. No algorithmic
change needed. Remove any assertions that check `umi.len() == 5`.

### Consensus calling

Review for any hardcoded UMI length. Update if present.

## Child tasks

| ID | File | Status |
|----|------|--------|
| CHEM-001 | todo/chem_001_generalise_umi.md | todo |
| CHEM-002 | todo/chem_002_generalise_parser.md | todo |
| CHEM-003 | todo/chem_003_generalise_clustering.md | todo |
| CHEM-004 | todo/chem_004_read_structure_presets.md | todo |
| CHEM-005 | todo/chem_005_non_twist_integration.md | todo |

## Scope

- `kam-core/src/` — `CanonicalUmiPair`, `Molecule`: UMI fields to `Vec<u8>`
- `kam-assemble/src/parser.rs` — remove hardcoded lengths, use `ReadStructure`
- `kam-assemble/src/clustering.rs` — remove UMI-length assertions
- `kam-assemble/src/consensus.rs` — remove UMI-length assumptions if present
- `kam/src/config.rs` — add `ReadStructure` to pipeline config

## Out of scope

- Changing k-mer indexing (UMI is not part of k-mer encoding)
- Changing variant calling logic
- Paired-end read structure with non-symmetric UMIs
