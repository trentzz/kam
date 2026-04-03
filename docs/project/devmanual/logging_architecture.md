# Logging and Observability Architecture

## Design Principle

Every decision kam makes should be traceable. When a read is dropped, a family is split, a variant is filtered — the user should be able to find out exactly why. But this must be opt-in: zero-overhead when logging is off, full transparency when on.

## Log Output Structure

Each kam run produces a log directory (default: `logs/` or `--log-dir <path>`):

```
logs/
├── summary.json              # always written — high-level stats, pass/fail
├── reads_dropped.tsv         # reads rejected by parser, with reason per read
├── umi_clustering.tsv        # UMI merge decisions (which UMIs were grouped)
├── families.tsv              # per-family stats (size, strand counts, type, collision prob)
├── consensus_conflicts.tsv   # positions where consensus vote was ambiguous
├── kmers_skipped.tsv         # reads too short for chosen k, template issues
├── anchors_flagged.tsv       # non-unique anchors in graph walking
├── variants_filtered.tsv     # variants that didn't pass filters, with reason
└── full_trace.jsonl          # everything, structured, for programmatic analysis
```

## Log Levels / Granularity

Each log file is independently togglable. Configuration via CLI flags or a TOML config:

```toml
[logging]
log_dir = "logs/"

# Per-file toggles (all default to false except summary)
summary = true               # always on by default
reads_dropped = false         # every dropped read with reason
umi_clustering = false        # UMI merge decisions
families = false              # per-family statistics
consensus_conflicts = false   # ambiguous consensus positions
kmers_skipped = false         # reads too short for k
anchors_flagged = false       # non-unique anchor warnings
variants_filtered = false     # filtered variant details
full_trace = false            # everything (large, slow)
```

CLI equivalent:
```bash
kam assemble-molecules \
    --log-dir logs/ \
    --log reads_dropped \
    --log families \
    --log consensus_conflicts
```

`--log all` enables everything. `--log none` disables everything except summary.

## Read Drop Log Format

`reads_dropped.tsv` — one line per dropped read pair:

```
read_id    reason    detail    r1_len    r2_len    umi_r1    umi_r2    umi_qual_min
READ_001   template_too_short   r1_template=12bp,min=20   19   150   ACGTA   TGCAT   35
READ_002   low_umi_quality      r2_umi_pos3=Q8,threshold=Q20   150   150   ACGTA   NGCAT   8
READ_003   template_too_short   r2_template=5bp,min=20   150   12   TGCAT   ACGTA   32
READ_004   invalid_bases        r1_umi_has_N_at_pos2   150   150   ANGTA   TGCAT   2
```

Reasons enum:
- `template_too_short` — template length below user-configured minimum
- `low_umi_quality` — UMI base quality below threshold
- `invalid_bases` — N bases in UMI (configurable: drop or keep)
- `read_too_short` — entire read shorter than umi+skip (malformed)
- `unpaired` — R1 without matching R2 or vice versa

## Implementation: Zero-Cost When Off

Logging must not slow down the hot path when disabled. Pattern:

```rust
pub struct DropLog {
    writer: Option<BufWriter<File>>,
}

impl DropLog {
    pub fn log_drop(&mut self, read_id: &[u8], reason: DropReason, detail: &str) {
        if let Some(w) = &mut self.writer {
            // only pays IO cost when enabled
            writeln!(w, "{}\t{}\t{}",
                std::str::from_utf8(read_id).unwrap_or("?"),
                reason.as_str(),
                detail
            ).ok();
        }
    }
}
```

Each log sink is an `Option<BufWriter<File>>` — when `None`, the branch compiles to nothing meaningful. No allocation, no formatting, no syscalls.

## Counters Always On

Even when detailed logs are off, **counters** are always tracked (negligible cost):

```rust
pub struct ParseStats {
    pub n_read_pairs_processed: u64,
    pub n_template_too_short: u64,
    pub n_low_umi_quality: u64,
    pub n_invalid_bases: u64,
    pub n_read_too_short: u64,
    pub n_unpaired: u64,
    pub n_passed: u64,
}
```

These feed into `summary.json` (always written) and the QC JSON for Nextflow validation.

## Interaction with Nextflow

- `summary.json` is the QC output that the next Nextflow stage validates
- Detailed log files are published to `${params.outdir}/${sample_id}/logs/` via `publishDir`
- For production runs: only `summary` on (fast)
- For debugging/validation: enable specific logs as needed
- For development: `--log all` with `full_trace.jsonl`
