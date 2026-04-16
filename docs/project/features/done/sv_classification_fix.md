# NovelInsertion Classification Fix

**Status**: done
**Epic**: SV-CLASSIFY

## Problem

The `is_tandem_duplication()` function in `kam-call/src/caller.rs` incorrectly
classified genuine novel insertion sequences as tandem duplications. All
NovelInsertion benchmark variants were reported with SVTYPE=DUP instead of
SVTYPE=INS.

## Root Cause

The tandem duplication check searched for the inserted sequence anywhere in the
reference window. Novel insertions with partial sequence similarity to the
reference triggered the check, even when the match was too short to indicate a
genuine tandem repeat.

## Fix

Tightened the `is_tandem_duplication()` function to require at least 80% of the
inserted sequence to match a contiguous reference segment. Partial matches below
this threshold are classified as NovelInsertion.

## Verification

- [x] All existing unit tests pass
- [x] New regression tests cover benchmark-derived sequences
- [x] Full NovelInsertion benchmark re-run shows correct classification
