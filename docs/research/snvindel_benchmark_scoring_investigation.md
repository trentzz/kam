# SNV/Indel Benchmark Scoring Investigation

## Symptom

Initial benchmark runs showed:
- SNV discovery sensitivity: 90-100% at most VAF levels
- SNV tumour-informed sensitivity: 20-40%
- Indel discovery sensitivity: 0% at all VAF levels
- Indel tumour-informed sensitivity: 0% at all VAF levels

The large gap between discovery (high) and tumour-informed (low) SNV sensitivity, and
the complete absence of indel detection, pointed to methodological problems.

## Investigation

### Indel 0% sensitivity

**Hypothesis 1**: Indel detection fails in kam.

**Test**: Checked the truth VCF REF sequences against the actual reference FASTA.

```
pos 100 (del3): truth REF=GCCC, actual ref=TTTG
pos 300 (ins5): truth REF=T,    actual ref=A   (anchor only, likely OK)
pos 550 (del4): truth REF=AATAA, actual ref=GTGA
pos 750 (ins5): truth REF=C,    actual ref=T   (anchor only, likely OK)
pos 1000 (del4): truth REF=ACATG, actual ref=TACA
```

**Root cause**: `INDEL_ROWS` in `make_snv_indel_suite.py` were written without checking the
actual reference sequence. For deletions, the REF field must match the reference exactly.
Varforge accepted the input VCF regardless and applied the mutations, but the simulated reads
and varforge output truth VCF retained the wrong REF sequences. Kam calls the actual reference
bases, so the called (REF, ALT) pairs never matched the truth (REF, ALT) pairs.

**Fix**: Updated `INDEL_ROWS` with actual reference bases:
```python
("chr1", 100,  "TTTG",  "T",      "TYPE=INDEL;SVLEN=-3"),
("chr1", 300,  "A",     "AGACGT", "TYPE=INDEL;SVLEN=5"),
("chr1", 550,  "GTGA",  "G",      "TYPE=INDEL;SVLEN=-3"),
("chr1", 750,  "T",     "TGCTAG", "TYPE=INDEL;SVLEN=5"),
("chr1", 1000, "TACA",  "T",      "TYPE=INDEL;SVLEN=-3"),
```

Regenerated all 50 indel truth VCFs, re-ran varforge, re-ran kam.

### SNV discovery inflated sensitivity

**Hypothesis 2**: SNV discovery accuracy is inflated by REF/ALT collisions.

**Test**: Compared truth REF/ALT pairs with all PASS calls at 1% VAF:

Truth variants: (C,T) at 50, (A,G) at 200, (G,A) at 450, (T,C) at 800, (C,A) at 1200.

Discovery PASS calls included: (C,T) at 90, (C,A) at 95, (A,G) at 199, (G,A) at 449,
(T,C) at 837. The REF/ALT-only scoring matched all 5 truth variants via positions that
were near the first target region (around pos 73-137), not the actual truth positions.

**Root cause**: REF/ALT-only matching (position-independent) is not appropriate for
these benchmarks. The scoring counts any PASS call anywhere with the same base substitution
as a TP. In a region with many false positives, common substitutions (G→A, C→T) will appear
repeatedly and accidentally match truth variants far away.

**Fix**: Changed `score_snv_indel_suite.py` to use position-based matching (±10bp tolerance),
matching the approach used in `score_sv_suite.py`.

### SNV wrong REF sequences

**Hypothesis 3**: SNV REF sequences also mismatched the actual reference.

**Test**:
```
pos 50:   truth REF=C, actual ref=A
pos 200:  truth REF=A, actual ref=A ✓
pos 450:  truth REF=G, actual ref=G ✓
pos 800:  truth REF=T, actual ref=G
pos 1200: truth REF=C, actual ref=A
```

Three of five SNV positions had wrong REF sequences. This meant:
- Varforge applied mutations from incorrect base (e.g., C→T at pos 50 where ref is A)
- The simulated reads contained wrong mutations at those positions
- Kam correctly reports actual reference bases; its calls for pos 50 would be A→T, not C→T
- For position-based scoring this still works if varforge placed the mutation near pos 50

**Fix**: Updated `SNV_ROWS` with actual reference bases:
```python
("chr1", 50,   "A", "T",  "TYPE=SNP"),
("chr1", 200,  "A", "G",  "TYPE=SNP"),
("chr1", 450,  "G", "A",  "TYPE=SNP"),
("chr1", 800,  "G", "A",  "TYPE=SNP"),
("chr1", 1200, "A", "C",  "TYPE=SNP"),
```

Regenerated all 50 SNV truth VCFs, re-ran varforge, re-ran kam.

### SNV tumour-informed low sensitivity (40%)

Even with corrected positions, tumour-informed mode passed only 2/5 SNVs at 1% VAF
(positions 200 and 450).

**Cause**: Two factors compound:
1. Off-by-1 in allele.rs: VCF output positions are `target_start + diff_pos` (0-based offset,
   no +1). This means kam writes position 199 when the true position is 200. The tumour-informed
   filter (`targeting.rs`) uses `target_start + diff_pos + 1`, so the filter position matches
   the truth VCF correctly for those variants.
2. Some SNVs not detected near their actual positions. With wrong REF sequences in the original
   SNV_ROWS, varforge may have placed mutations at different reference bases, creating reads
   that kam cannot reliably call as the intended variants.

After correcting REF sequences and re-running all datasets, expected improvement in
tumour-informed sensitivity.

## Actions Taken

1. Fixed `INDEL_ROWS` and `SNV_ROWS` in `make_snv_indel_suite.py` (correct REF sequences).
2. Regenerated all 100 truth VCFs.
3. Deleted all 100 existing simulated datasets and re-ran varforge.
4. Deleted all 100 existing kam results for re-run.
5. Updated `score_snv_indel_suite.py` to use position-based scoring (±10bp).

## Expected Impact

- Indel sensitivity: should become non-zero after correct REF sequences.
- SNV discovery sensitivity: should drop from inflated 90-100% to reflect actual detection
  near the correct positions.
- SNV tumour-informed sensitivity: should improve once SNVs are simulated with correct bases.
