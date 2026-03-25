# BENCH-SVN-002: NovelInsertion Varforge Configs and Truth VCFs

**Epic**: BENCH-SV-NEW (docs/claudetracking/overallplans/BENCH-SV-NEW.md)
**Priority**: high
**Depends on**: sv_exp_003_novel_insertion.md
**Status**: todo

## Goal

Generate 50 varforge configs and matching truth VCFs for NovelInsertion
benchmarking: 25 VAF levels × 2 replicates. Done looks like: all 50 configs
and truth VCFs exist in `docs/benchmarking/sv_new/configs/` and
`docs/benchmarking/sv_new/data/`, and varforge produces FASTQ for all 50
without error.

## Steps

1. Read the existing insertion suite generation script to understand the
   pattern (the existing SV suite may include an INS type).

2. Write `docs/benchmarking/sv_new/scripts/make_novins_suite.py`:
   - 25 VAF levels × 2 replicates.
   - Config naming: `novins_vaf{tag}_{rep}.yaml`.
   - Truth VCF naming: `truth_novins_vaf{tag}_{rep}.vcf`.
   - SV definition: 1 insertion of ≥60 bp of non-repetitive sequence on the
     2000 bp synthetic reference. The inserted sequence must not be a substring
     or near-copy of the reference window.

3. Generate the non-repetitive insertion sequence:
   - Use a fixed random byte string that has been verified (by inspection) to
     not repeat in the 2000 bp reference.
   - Hard-code this string in the generation script for reproducibility.

4. Run the script to generate all 50 configs and truth VCFs.

5. Run varforge for all 50 configs and verify all FASTQ pairs are produced.

## Notes

- The inserted sequence length (≥60 bp) must exceed the 50 bp NovelInsertion
  threshold in the classifier by a comfortable margin.
- If varforge uses a genome-wide reference for insertion, ensure the inserted
  sequence is also absent from the broader reference (use a synthetic sequence
  not found in any common genome).
- Store the canonical inserted sequence in
  `docs/benchmarking/sv_new/data/novins_insert_seq.txt` for reference.
