# Feature: Interactive CLI Explorer (`kam explore`)

## Summary

A REPL-style interactive terminal interface for exploring kam's binary data files (bincode-serialized molecules, k-mer indices, variant calls). Reads directly from bincode without requiring full export. Inspired by the feel of gdb/claude code — type commands, get pretty-printed results, selectively export what you need.

## Motivation

Binary intermediate files are fast but opaque. Exporting the entire dataset to JSON/TSV just to find one molecule or check one k-mer is wasteful. The explorer lets you interactively query, search, filter, and selectively export from any kam data file.

## Invocation

```bash
kam explore molecules.bin          # explore molecule assembly output
kam explore kmer_index.bin         # explore k-mer index
kam explore variants.bin           # explore variant call results
```

Auto-detects file type from header magic bytes in the bincode file.

## Interface Design

### Look and Feel

Pretty-printed, colored terminal output. Prompt shows the file type and record count:

```
kam explore v0.1.0
Loaded: molecules.bin (612,847 molecules, 4,821,044 input reads)
Type 'help' for commands, 'quit' to exit.

molecules>
```

### Core Commands

#### Navigation and Inspection
```
molecules> summary
  Total molecules: 612,847
  Duplex families: 489,478 (79.8%)
  Simplex fwd:     73,927 (12.1%)
  Simplex rev:     49,442 (8.1%)
  Mean family size: 7.2
  Median family size: 6

molecules> head 5
  [table showing first 5 molecules with key fields]

molecules> show molecule:a3f8b2c1
  Molecule ID:    a3f8b2c1
  UMI pair:       ACGTA-TGCAT
  Family type:    DUPLEX
  Fwd reads:      5
  Rev reads:      3
  Family size:    8
  Consensus len:  143bp
  Mean quality:   0.0003
  Collision prob: 0.008
  Skip bases:     CT / CT
```

#### Searching and Filtering
```
molecules> search umi ACGTA-TGCAT
  Found 1 molecule matching UMI pair ACGTA-TGCAT
  [shows molecule detail]

molecules> filter family_type == DUPLEX AND family_size >= 10
  Matched 47,231 molecules (7.7%)
  [shows summary of matched set]

molecules> filter family_size < 3
  Matched 89,012 molecules (14.5%)
  [shows summary]

molecules> filter collision_prob > 0.05
  Matched 2,341 molecules (0.4%)
  [shows summary with warning about high collision risk]
```

#### K-mer Index Exploration
```
kmers> lookup ACGTACGTACGTACGTACGTACGTACGTACG
  K-mer:          ACGTACGTACGTACGTACGTACGTACGTACG
  Canonical:      ACGTACGTACGTACGTACGTACGTACGTACG
  Molecules:      891
  Duplex:         734
  Simplex fwd:    98
  Simplex rev:    59
  Min error prob: 0.000001
  Mean error prob: 0.0004

kmers> top 10 by n_molecules
  [table of 10 highest-coverage k-mers]

kmers> filter n_duplex == 0
  Matched 1,247 k-mers with no duplex support
```

#### Variant Exploration
```
variants> list
  [table of all called variants with VAF, confidence, filter status]

variants> show TP53_exon7:SNV
  Target:         TP53_exon7
  Type:           SNV
  Ref path:       ...ACGTACGT...
  Alt path:       ...ACGTTCGT...
  VAF:            0.0034 (95% CI: 0.0018-0.0062)
  Molecules ref:  891
  Molecules alt:  3
  Duplex alt:     2
  Strand bias p:  0.48
  Confidence:     0.997
  Filter:         PASS

variants> why-filtered KRAS_codon12:SNV
  Variant at KRAS_codon12 was filtered: STRAND_BIAS
  Reason: All 4 supporting molecules are forward-strand only (p=0.003)
  Strand bias threshold: p < 0.01
  To include this variant, rerun with --strand-bias-threshold 0.001
```

#### Selective Export
```
molecules> export filter family_type == DUPLEX --format tsv --output duplex_only.tsv
  Exported 489,478 molecules to duplex_only.tsv

molecules> export show molecule:a3f8b2c1 --format json
  {
    "id": "a3f8b2c1",
    "umi_pair": "ACGTA-TGCAT",
    ...
  }

molecules> export head 100 --format csv --output sample.csv
  Exported 100 molecules to sample.csv

kmers> export filter n_molecules > 500 --format tsv --output high_coverage_kmers.tsv
  Exported 12,847 k-mers to high_coverage_kmers.tsv
```

#### Help and Discovery
```
molecules> help
  Navigation:    summary, head N, tail N, show <id>
  Search:        search <field> <value>, filter <expression>
  Export:        export <command> --format <tsv|csv|json> [--output <file>]
  Info:          fields, stats <field>, histogram <field>
  Other:         help, quit, clear

molecules> fields
  id               u64      Molecule ID (hash of canonical UMI)
  umi_pair         String   Canonical UMI pair (dash-separated)
  family_type      Enum     DUPLEX, SIMPLEX_FWD, SIMPLEX_REV, SINGLETON
  family_size      u8       Total reads in family
  fwd_reads        u8       Forward strand read count
  rev_reads        u8       Reverse strand read count
  consensus_len    usize    Consensus sequence length (bp)
  mean_quality     f32      Mean per-base error probability
  collision_prob   f32      Estimated UMI collision probability
  skip_r1          String   Skip bases from R1
  skip_r2          String   Skip bases from R2

molecules> histogram family_size
  1-2   ████████░░░░░░░░░░░░  14.5%
  3-5   ████████████████░░░░  38.2%
  6-10  ████████████████████  41.1%
  11-20 ███░░░░░░░░░░░░░░░░░   5.8%
  21+   ░░░░░░░░░░░░░░░░░░░░   0.4%

molecules> stats family_size
  Count:   612,847
  Mean:    7.2
  Median:  6
  Std:     3.1
  Min:     1
  Max:     47
  P5:      2
  P95:     14
```

## Implementation Notes

### Rust Crates
- `ratatui` or `crossterm` — terminal UI / colored output
- `rustyline` — readline-style input with history and tab completion
- `comfy-table` — pretty table formatting
- `textplots` — terminal histograms/charts

### Architecture
- Deserialize bincode file into memory on load (fast — this is what bincode is good at)
- For very large files: memory-map and deserialize lazily (only materialize records when accessed)
- Filter expressions parsed into a simple AST, evaluated per-record
- Export streams results to file without materializing the full filtered set in memory

### Filter Expression Mini-Language
Simple, no need for a full DSL initially:
```
field op value [AND|OR field op value ...]
```
Operators: `==`, `!=`, `>`, `>=`, `<`, `<=`, `contains`
Values: integers, floats, strings (auto-detected)
Parentheses for grouping in later version.

## Priority

This is a quality-of-life feature, not on the critical path. Build after the core pipeline works. But design the bincode file format with exploration in mind:
- Include a header with file type, version, record count, field names
- Include summary statistics in the header so `summary` is instant
- Use a record-oriented layout so individual records can be deserialized without reading the whole file

## Relationship to Export Flags

The explorer subsumes the export functionality. But standalone export is still useful for scripting:

```bash
# Non-interactive export (for scripts/pipelines)
kam export molecules.bin --format tsv --output molecules.tsv
kam export molecules.bin --filter "family_type == DUPLEX" --format csv --output duplex.csv
kam export kmer_index.bin --format json --output kmers.json

# Interactive exploration
kam explore molecules.bin
```
