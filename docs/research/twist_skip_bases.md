# Twist UMI Adapter Skip/Spacer Bases

## What We Know

- Twist UMI duplex adapters use read structure `5M2S+T` — the `2S` is a 2bp skip/spacer
- The skip bases are a **fixed monotemplate** from the adapter ligation chemistry, not random and not genomic
- They are **consistent within an adapter lot**
- Cross-lot consistency is **not publicly confirmed**

## What's Not Publicly Available

- The actual dinucleotide identity is not published in Twist product pages, technical docs, white papers, or published papers using Twist adapters
- No bioinformatics tool (fgbio, fastp, UMI-tools, HUMID) hardcodes or validates Twist-specific skip bases
- fgbio's `S` (skip) operator simply discards the bases without inspecting their identity
- The adapter oligo specification sheet (Certificate of Analysis) shipped with the product may contain the full adapter sequence including spacer

## How to Determine Skip Bases Empirically

Parse positions 5-6 from a representative sample of reads in a run. Because the bases are fixed (not random), a single dominant dinucleotide will appear at overwhelming frequency (>95%). The remaining <5% are sequencing errors in the skip positions.

```rust
/// Count dinucleotide frequencies at skip positions across all reads
fn detect_skip_bases(reads: &[&[u8]], skip_start: usize) -> HashMap<[u8; 2], u64> {
    let mut counts = HashMap::new();
    for read in reads {
        if read.len() > skip_start + 1 {
            let skip: [u8; 2] = [read[skip_start], read[skip_start + 1]];
            *counts.entry(skip).or_default() += 1;
        }
    }
    counts
}
```

## kam's Approach

### Default: Auto-Detect and Report

1. During parsing, extract skip bases from every read
2. After processing (or during a quick first-pass sample), determine the dominant dinucleotide
3. Report skip base distribution in QC JSON:

```json
{
  "skip_base_distribution": {
    "CT": 4789234,
    "CC": 12847,
    "CA": 3201,
    "other": 8923
  },
  "dominant_skip": "CT",
  "dominant_skip_fraction": 0.995,
  "skip_base_consistent": true
}
```

### Optional: User-Supplied Validation

```bash
kam assemble-molecules --expected-skip-bases CT ...
```

When set:
- Reads where skip bases don't match expected are **flagged in the drop log** (if logging enabled)
- Reads are NOT hard-rejected by default (the skip mismatch is recorded but the read is still parsed)
- User can additionally set `--reject-skip-mismatch` to hard-reject mismatched reads
- Mismatch rate reported in QC

### QC Signal: Skip Consistency Within Families

After molecule assembly, check skip base consistency within each read family:
- All reads in a family should have the same skip bases (they're from the same adapter molecule)
- A family where some reads have different skip bases may indicate:
  - Adapter contamination
  - UMI collision (reads from different molecules grouped together)
  - Quality issues
- Report per-family skip consistency in the families log (when enabled)

## Configurable Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--expected-skip-bases` | None (auto-detect) | Expected 2bp skip sequence for validation |
| `--reject-skip-mismatch` | false | Hard-reject reads with unexpected skip bases |
| `--skip-quality-threshold` | None | Minimum Phred quality for skip bases (flag if below) |

## Recommendation for Users

1. Run once with `--log families` to see skip base distribution
2. If consistent (>99% one dinucleotide), supply it as `--expected-skip-bases` for future runs as an extra QC check
3. If using multiple adapter lots, check whether the skip bases are the same across lots
4. Contact Twist technical support for definitive lot-to-lot information

## References

Twist Biosciences (2024) Twist UMI Adapter System. Product documentation. Available at: https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system (Accessed: March 2026). [Primary source for the `5M2S+T` read structure, the fixed monotemplate skip/spacer description, and the statement that within-lot consistency is expected but cross-lot consistency is not publicly confirmed.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of the fgbio `S` (skip) operator behaviour — discards skip bases without inspecting their identity — cited as the standard tool comparison point.]
