# Feature Completeness Review

## What's Implemented

The pipeline works end-to-end. You can run `kam run --r1 R1.fq --r2 R2.fq --targets targets.fa --output-dir out/` and get variant calls in TSV/CSV/JSON/VCF. Every stage produces QC JSON. The CLI is well-designed with sensible defaults and clear help text.

### Feature Inventory

| Feature | Status | Notes |
|---|---|---|
| FASTQ parsing (Twist 5M2S+T) | Done | Uncompressed only |
| UMI extraction | Done | Hardcoded to 5bp |
| Canonical UMI pairing | Done | Correct implementation |
| Hamming distance clustering | Done | O(n^2) but correct |
| Endpoint fingerprinting | Done | 64-bit fingerprint, 4-bit threshold |
| Quality-weighted SSC | Done | Correct algorithm |
| Duplex consensus | Done | PickBest and MaskAsN strategies |
| 2-bit k-mer encoding | Done | k=1..31, sliding window |
| Canonical k-mers | Done | min(fwd, revcomp) |
| Allowlist filtering | Done | Panel-bounded memory |
| Molecule-level indexing | Done | Key differentiator |
| De Bruijn graph construction | Done | Prefix-grouped edge finding |
| BFS path walking | Done | Cycle-safe, bounded |
| Path scoring | Done | min/mean molecules, duplex fraction |
| Anchor validation | Done | Configurable threshold |
| Beta posterior VAF | Done | With 95% credible intervals |
| Fisher strand bias | Done | Two-tailed exact test |
| Variant classification | Done | SNV, Insertion, Deletion, MNV, Complex |
| Multi-format output | Done | TSV, CSV, JSON, VCF |
| Per-stage QC JSON | Done | All 4 stages |
| Full pipeline (`kam run`) | Done | Zero-copy hot path |
| Individual subcommands | Done | assemble, index, pathfind, call |
| Bincode serialization | Done | With magic bytes and versioning |

## What's Missing

### Critical for Real Use

**1. Gzip/bgzf FASTQ Support**

Every real sequencing run produces compressed FASTQ. Without gzip support, users must decompress terabytes of data before running kam. `needletail` supports gzip natively — this is likely a one-line feature flag addition.

**2. Gzip FASTA Support**

Same issue for target sequence files, though these are small enough that manual decompression is tolerable.

**3. Progress Reporting**

The pipeline produces no progress output while running. For a job that might take hours on real data, users need feedback: read count, molecules assembled, targets processed, etc. A simple `eprintln!` every N records would suffice.

Currently, `run.rs` prints stage summaries *after* each stage completes. But during a long assembly stage, there's silence.

**4. Proper Exit Codes**

The `main.rs` uses `process::exit(1)` on error, but doesn't distinguish between different error types. A tool meant for pipeline integration should have distinct exit codes for: input error, runtime error, no variants found, etc.

**5. Thread Count Configuration**

The CLI accepts `--threads` but the value is never used anywhere. There's no thread pool, no rayon configuration, nothing. This is misleading.

### Important for Production

**6. Skip Base QC**

The research docs (`twist_skip_bases.md`) describe using skip bases as a QC signal — monotemplate bases should be consistent. The parser extracts skip bases but nothing checks them. A simple frequency analysis of skip bases would be a useful QC metric.

**7. Background Error Model**

The caller uses a fixed `background_error_rate = 1e-4` for all positions. The research docs describe per-site and per-trinucleotide background estimation. Without this, the caller will have inflated false positive rates at positions with higher systematic error.

**8. Multi-k Support**

The docs mention multi-k as a future feature. Currently, only one k value is used per run. Multi-k would improve sensitivity for variants near the ends of targets.

**9. Logging System**

The CLI has `--log-dir` and `--log` flags, and the docs describe an elaborate per-stage logging architecture with drop logs, UMI clustering logs, and family logs. None of this is implemented. The logging flags are accepted but ignored.

### Nice to Have

**10. `kam explore` Interactive Mode**

Described in `docs/features/interactive_explorer.md`. Not implemented.

**11. Collision Probability Reporting**

The assembler counts collisions but the collision probability formula from the output format spec is not computed or reported.

**12. Annotated FASTQ Output**

The output format spec describes custom FASTQ tags (MI, RX, FS, DS, etc.) for annotated FASTQ output from the assembly stage. Not implemented.

**13. Resume/Checkpoint**

No ability to resume a failed pipeline from a checkpoint. The bincode serialization exists for this purpose but there's no resume logic.

## CLI Design

The CLI is well-designed:

- Clear subcommands with descriptive help
- Long-form flags only (no cryptic `-x` short flags, except `-k` for k-mer size which is conventional)
- Sensible defaults (k=31, min_umi_quality=20, min_family_size=1)
- Repeatable flags (`--log umi --log family`)
- Output format as comma-separated string (`--output-format tsv,vcf`)

One issue: the `--chemistry` flag defaults to `"twist-umi-duplex"` but there's no validation that this is a known preset. Any string is accepted silently and the default Twist preset is always used.

## Output Formats

### TSV
Well-structured with all expected columns. Header is present. Matches the km-compatible format from the docs.

### CSV
Identical to TSV with commas instead of tabs. Works correctly.

### JSON
Pretty-printed array of objects. All fields present. Valid JSON confirmed by tests.

### VCF
VCF 4.3 with custom INFO fields (VAF, VAF_LO, VAF_HI, NREF, NALT, NDUPALT, NSIMALT, SBP, CONF). FILTER definitions for all filter types. Uses target_id as CHROM and position 1 (alignment-free, so no genomic coordinates).

**Issue**: VCF uses `POS=1` for all variants since there are no genomic coordinates. This is technically valid but some downstream VCF tools will choke on this. Consider using a non-standard format or at least documenting this clearly.

**Issue**: No sample column. Standard VCF has a FORMAT and sample column after INFO. The current output has only 8 columns (CHROM through INFO), which is valid for a sites-only VCF but may confuse tools expecting a sample column.

## Verdict

The feature set is complete enough to demonstrate the pipeline works. The missing features are well-documented in the research/planning docs and the deferred items are reasonable for a Phase 1 implementation. The critical gap is gzip support — without it, the tool literally cannot process real sequencing data.
