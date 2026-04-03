# NovelInsertion Classification Fix

**Status**: inprogress
**Epic**: SV-CLASSIFY

## Problem

The `is_tandem_duplication()` function in `kam-call/src/caller.rs` incorrectly
classifies genuine novel insertion sequences as tandem duplications. All
NovelInsertion benchmark variants are reported with SVTYPE=DUP instead of
SVTYPE=INS.

## Root Cause

Under investigation (SV-CLS-001).

## Fix

Pending investigation results.

## Verification

- All existing unit tests pass
- New regression tests cover benchmark-derived sequences
- Full NovelInsertion benchmark re-run shows correct classification
