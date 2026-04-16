# BND Strand Orientation

## Status

Todo.

## Summary

Support all 4 BND strand orientations (FF, FR, RF, RR) for coordinate-based
fusion input in `make_sv_targets.py` and VCF output.

## Problem

The current implementation assumes forward-forward orientation when generating
junction sequences from coordinates. The VCF BND notation encodes four
possible strand orientations using bracket syntax:

| Orientation | VCF notation | Meaning |
|---|---|---|
| FF (forward-forward) | `N[chr:pos[` | Join forward A to forward B |
| FR (forward-reverse) | `N]chr:pos]` | Join forward A to reverse B |
| RF (reverse-forward) | `]chr:pos]N` | Join reverse A to forward B |
| RR (reverse-reverse) | `[chr:pos[N` | Join reverse A to reverse B |

When auto-generating junction sequences from coordinates (via
`make_sv_targets.py`), the tool must reverse-complement the appropriate
partner segment based on the specified orientation. Currently, only FF is
handled.

## Solution (planned)

1. Add an `--orientation` parameter to `make_sv_targets.py` accepting one of
   `FF`, `FR`, `RF`, `RR` per fusion entry.
2. When generating the junction sequence, reverse-complement partner B for FR,
   reverse-complement partner A for RF, reverse-complement both for RR.
3. Update the VCF BND writer in `kam-call/src/output.rs` to use the correct
   bracket notation based on the detected orientation.

## Notes

When users provide junction sequences directly via `--junction-sequences`,
this issue does not arise. The observed sequence from the BAM already encodes
the correct orientation. The sequence IS the forward-strand representation of
the chimeric read, with both partners in their observed orientations and any
inserted nucleotides included.

This is a key advantage of the `--junction-sequences` approach: it bypasses
the strand orientation problem entirely.
