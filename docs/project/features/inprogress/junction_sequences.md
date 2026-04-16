# Junction Sequences Input

## Status

In Progress.

## Summary

New FASTA input mode accepting raw junction sequences without coordinate
headers. Sequences observed in BAM/IGV can be provided directly to detect and
monitor fusions and junctions in ctDNA.

## Problem

The current fusion input required an exact coordinate-based header format
(`{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`). Users
who identify a fusion in IGV must know the exact genomic coordinates. Random
nucleotides at junctions and strand orientation concerns made it unclear how
to specify inputs.

These requirements created friction for the primary use case: a clinician or
researcher sees a chimeric read in a BAM file and wants to monitor that event
in serial cfDNA samples. They have the sequence but not necessarily the
coordinates.

## Solution

`--junction-sequences` accepts any FASTA file with any header format. Each
sequence is added to the k-mer allowlist and walked as a standalone target.
The total library depth is used as the VAF denominator.

The sequence itself encodes orientation and includes any inserted bases at the
junction. No special handling for strand or random nucleotides is needed.

Key implementation details:

1. Each FASTA entry is parsed with `needletail` and its k-mers are added to
   the allowlist alongside normal target and SV junction k-mers.
2. Each entry is treated as a standalone walk target. The reference path is
   the provided sequence; any molecule containing junction-spanning k-mers
   contributes to the alt count.
3. VAF denominator is the total library molecule count (not partner depth),
   because no partner locus coordinates are available.
4. Output uses the FASTA header as the `target_id`.

Config file equivalent: `[input] junction_sequences = "path/to/junctions.fa"`.

## Tests

- Unit test: allowlist extension includes junction k-mers
- Unit test: on-target molecule capture for junction-spanning reads
- Unit test: junction walk produces correct path with molecule counts
- Integration test: end-to-end run with `--junction-sequences` produces
  expected variant call

## Relevant Files

- `kam/src/cli.rs` — new `--junction-sequences` CLI flag
- `kam/src/config.rs` — `junction_sequences` field in config
- `kam/src/commands/run.rs` — loading and passing junction sequences
- `kam/src/commands/index.rs` — allowlist extension logic
