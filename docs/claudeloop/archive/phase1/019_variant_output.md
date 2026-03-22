# Task 019: Variant Output Writers

## Location
`kam-call/src/output.rs`

## What to implement

Write variant calls to TSV, CSV, JSON, and VCF formats. See `docs/research/output_format_specs.md` for column definitions.

## Interface

```rust
use std::io::Write;
use crate::caller::{VariantCall, VariantType, VariantFilter};

/// Output format selection
#[derive(Debug, Clone, Copy)]
pub enum OutputFormat {
    Tsv,
    Csv,
    Json,
    Vcf,
}

/// Write variant calls to any supported format
pub fn write_variants(
    calls: &[VariantCall],
    format: OutputFormat,
    writer: &mut dyn Write,
) -> std::io::Result<()>;

/// Write TSV format (km-compatible)
pub fn write_tsv(calls: &[VariantCall], writer: &mut dyn Write) -> std::io::Result<()>;

/// Write CSV format
pub fn write_csv(calls: &[VariantCall], writer: &mut dyn Write) -> std::io::Result<()>;

/// Write JSON format (array of objects)
pub fn write_json(calls: &[VariantCall], writer: &mut dyn Write) -> std::io::Result<()>;

/// Write VCF format
pub fn write_vcf(calls: &[VariantCall], writer: &mut dyn Write) -> std::io::Result<()>;
```

## TSV columns (tab-separated, with header)

```
target_id	variant_type	ref_seq	alt_seq	vaf	vaf_ci_low	vaf_ci_high	n_molecules_ref	n_molecules_alt	n_duplex_alt	n_simplex_alt	strand_bias_p	confidence	filter
```

CSV is the same but comma-separated.

## JSON format

```json
[
  {
    "target_id": "TP53_exon7",
    "variant_type": "SNV",
    "ref_seq": "ACGT",
    "alt_seq": "ACTT",
    "vaf": 0.005,
    ...
  }
]
```

## VCF format

Minimal valid VCF 4.3 with custom INFO fields. Use target_id as CHROM, position 1 as POS (alignment-free, no real coordinates).

## Tests required

1. write_tsv produces correct header and one data line
2. write_csv uses commas instead of tabs
3. write_json produces valid JSON array
4. write_vcf produces valid VCF header and data line
5. Multiple variant calls produce multiple lines
6. Empty calls list → header only (TSV/CSV), empty array (JSON), header only (VCF)
7. Special characters in sequences don't break formatting
8. write_variants dispatches to correct format

## Definition of done

- `cargo test -p kam-call` passes
- `cargo clippy -p kam-call -- -D warnings` passes
- Doc comments
- Add `serde_json` to kam-call dependencies for JSON output
