# BND Strand Orientation

## Status

Done.

## Summary

Support all 4 BND strand orientations (FF, FR, RF, RR) for coordinate-based
fusion input in `make_sv_targets.py` and VCF output.

## Problem

The implementation assumed forward-forward orientation when generating
junction sequences from coordinates. The VCF BND notation encodes four
possible strand orientations using bracket syntax:

| Orientation | VCF notation | Meaning |
|---|---|---|
| FF (forward-forward) | `N[chr:pos[` | Join forward A to forward B |
| FR (forward-reverse) | `N]chr:pos]` | Join forward A to reverse B |
| RF (reverse-forward) | `]chr:pos]N` | Join reverse A to forward B |
| RR (reverse-reverse) | `[chr:pos[N` | Join reverse A to reverse B |

## Solution

1. Added `--orientation` parameter to `make_sv_targets.py` accepting one of
   `FF`, `FR`, `RF`, `RR` per fusion entry.
2. Junction sequence generation reverse-complements partner B for FR,
   reverse-complements partner A for RF, and reverse-complements both for RR.
3. The VCF BND writer in `kam-call/src/output.rs` uses the correct bracket
   notation based on the detected orientation.

## Verification

- All four orientations produce correct junction sequences in unit tests.
- VCF BND output uses the correct bracket notation for each orientation.
- Existing fusion benchmarks pass with no regressions.

## Notes

When users provide junction sequences directly via `--junction-sequences`,
this issue does not arise. The observed sequence from the BAM already encodes
the correct orientation. The sequence is the forward-strand representation of
the chimeric read, with both partners in their observed orientations and any
inserted nucleotides included.

This is a key advantage of the `--junction-sequences` approach: it bypasses
the strand orientation problem entirely.
