# Interactive Explorer

## Overview

`kam explore` opens an interactive terminal session for inspecting kam's binary data files.
Molecule assemblies, k-mer indices, and variant call results are stored as bincode. The explorer
reads these files directly and presents records in a human-readable format. You can browse,
filter, compute statistics, and export subsets without converting the entire file to TSV or JSON.

The interface is a REPL (read-eval-print loop). Type a command, press Enter, see the result.
Commands are short and composable. Tab completion and command history are supported.

---

## Quick start

```bash
kam explore results/molecules.bin
```

Output:

```
kam explore v0.3.0
Loaded: molecules.bin (612,847 molecules, 4,821,044 input reads)
Type 'help' for commands, 'quit' to exit.

molecules>
```

The file type is detected automatically from the bincode header. The prompt changes to reflect
the file type: `molecules>`, `kmers>`, or `variants>`.

---

## File types

| File | Produced by | What it contains |
|------|------------|-----------------|
| `molecules.bin` | `kam assemble` | Assembled molecules with UMI, family type, consensus length |
| `index.bin` | `kam index` | K-mer index with per-kmer molecule provenance |
| `variants.bin` | `kam call` | Variant calls with VAF, confidence, filter status |

Each of these is produced when you run `kam run` or the individual stage subcommands with
binary output enabled.

---

## Commands

### summary

Print high-level statistics about the loaded file. For molecules, this includes total count,
duplex fraction, and family size distribution. For variants, this includes call counts by
filter status and variant type.

```
molecules> summary
  Total molecules: 612,847
  Duplex families: 489,478 (79.8%)
  Simplex fwd:     73,927 (12.1%)
  Simplex rev:     49,442 (8.1%)
  Mean family size: 7.2
  Median family size: 6
```

### head

Show the first N records as a table. Defaults to 10.

```
molecules> head 5
  ID        UMI_PAIR       TYPE     FAM_SIZE  FWD  REV  CONS_LEN  MEAN_QUAL
  a3f8b2c1  ACGTA-TGCAT    DUPLEX   8         5    3    143       0.0003
  b1c2d3e4  GCTAT-ATAGC    DUPLEX   6         3    3    151       0.0002
  c4d5e6f7  TGCAT-ATGCA    SIM_FWD  3         3    0    148       0.0008
  d7e8f9a0  ATGCA-TGCAT    DUPLEX   12        7    5    139       0.0001
  e0f1a2b3  CGTAC-GTACG    SIM_REV  2         0    2    145       0.0012
```

### show

Display full detail for a single record. Use the record ID.

```
molecules> show a3f8b2c1
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

For variant records, `show` prints the full call detail including path sequences:

```
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
```

### filter

Select records matching an expression. The result summary shows how many records matched and
their aggregate statistics. Filtered records become the active set for subsequent `head`,
`stats`, and `export` commands.

```
molecules> filter family_type == DUPLEX AND family_size >= 10
  Matched 47,231 molecules (7.7%)

molecules> filter collision_prob > 0.05
  Matched 2,341 molecules (0.4%)

variants> filter filter == PASS AND vaf > 0.001
  Matched 12 variants
```

To clear the active filter and return to the full dataset:

```
molecules> filter clear
  Filter cleared. 612,847 records active.
```

### stats

Compute summary statistics for a numeric field.

```
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

When a filter is active, `stats` operates on the filtered set only:

```
molecules> filter family_type == DUPLEX
  Matched 489,478 molecules (79.8%)

molecules> stats family_size
  Count:   489,478
  Mean:    8.1
  Median:  7
  ...
```

### histogram

Draw a terminal histogram for a numeric field.

```
molecules> histogram family_size
  1-2   ████████░░░░░░░░░░░░  14.5%
  3-5   ████████████████░░░░  38.2%
  6-10  ████████████████████  41.1%
  11-20 ███░░░░░░░░░░░░░░░░░   5.8%
  21+   ░░░░░░░░░░░░░░░░░░░░   0.4%
```

### fields

List all available fields and their types for the current file type.

```
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
```

### export

Write the current view (all records, or the filtered set) to a file.

```
molecules> export --format tsv --output duplex_only.tsv
  Exported 489,478 molecules to duplex_only.tsv
```

Supported formats: `tsv`, `csv`, `json`. If no `--output` is given, the result prints to the
terminal (useful for small filtered sets).

```
molecules> filter family_size > 30
  Matched 214 molecules (0.03%)

molecules> export --format json
  [
    {"id": "f2a1b3c4", "umi_pair": "TGCAT-ACGTA", "family_type": "DUPLEX", ...},
    ...
  ]
```

### help

Print the command reference.

```
molecules> help
  Navigation:    summary, head N, show <id>
  Search:        filter <expression>, filter clear
  Analysis:      stats <field>, histogram <field>, fields
  Export:        export --format <tsv|csv|json> [--output <file>]
  Other:         help, quit, clear
```

### quit

Exit the explorer. Ctrl-D also works.

---

## Filter expression syntax

Expressions follow the pattern `field op value`, joined by `AND` or `OR`.

### Operators

| Operator | Meaning |
|----------|---------|
| `==` | Equals |
| `!=` | Not equals |
| `>` | Greater than |
| `>=` | Greater than or equal |
| `<` | Less than |
| `<=` | Less than or equal |
| `contains` | String contains substring |

### Values

- Integers: `42`, `1000`
- Floats: `0.05`, `1.5`
- Strings and enums: `DUPLEX`, `PASS`, `StrandBias`

Strings are unquoted. The parser auto-detects the type from the field definition.

### Compound expressions

```
filter family_type == DUPLEX AND family_size >= 5
filter family_type == SIMPLEX_FWD OR family_type == SIMPLEX_REV
filter vaf > 0.001 AND filter == PASS AND n_duplex_alt >= 1
```

Expressions are evaluated left to right. `AND` binds before `OR`. Parentheses are not yet
supported.

### Examples

Find high-collision-risk molecules:

```
molecules> filter collision_prob > 0.05
```

Find all PASS variants above 0.1% VAF:

```
variants> filter filter == PASS AND vaf > 0.001
```

Find simplex-only molecules with large families:

```
molecules> filter family_type != DUPLEX AND family_size >= 10
```

---

## Table output format

All tabular output (from `head`, `filter`, `export --format tsv`) uses fixed-width columns
with header labels. Column widths adapt to the data. Long sequences are truncated with `...`.

For molecule files, the default columns are: ID, UMI_PAIR, TYPE, FAM_SIZE, FWD, REV,
CONS_LEN, MEAN_QUAL.

For variant files, the default columns are: TARGET_ID, TYPE, VAF, N_REF, N_ALT, N_DUPLEX,
CONFIDENCE, FILTER.

---

## Tips

### Exploring molecules

Start with `summary` to understand the library. Check the duplex fraction: below 10% suggests
a library preparation issue. Use `histogram family_size` to see the family size distribution.

Filter for high-collision-risk molecules (`collision_prob > 0.05`) to find potential UMI
collisions. These molecules may contain mixed templates and could introduce artefacts.

Export singletons (`family_size == 1`) to quantify the proportion of the library that
contributes no consensus improvement.

### Exploring variant calls

Start with `summary` to see how many calls passed filters. Use
`filter filter == PASS AND n_duplex_alt == 0` to find PASS calls without duplex support.
These are your weakest true positives, worth manual review.

Use `stats vaf` on the PASS set to see the VAF distribution. A spike near 50% indicates
germline variants that should be filtered with `--max-vaf 0.35`.

### Finding outliers

Sort by a field using `filter` with a threshold. For example, to find the highest-VAF calls:

```
variants> filter vaf > 0.10
```

To find molecules with unusually poor consensus quality:

```
molecules> filter mean_quality > 0.01
```

### Combining with shell tools

Export a filtered set to TSV, then use standard tools for further analysis:

```
molecules> filter family_type == DUPLEX AND family_size >= 5
molecules> export --format tsv --output duplex_large.tsv
molecules> quit
```

```bash
# Count molecules per skip-base combination
cut -f10,11 duplex_large.tsv | sort | uniq -c | sort -rn | head
```
