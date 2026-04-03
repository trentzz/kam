# BENCH-SV-001: New SV Types in Benchmark

**Epic**: BENCH-SV (`overallplans/BENCH-SV.md`)
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Add truth VCFs and varforge configs for four additional SV types: INS, large
DEL (>500 bp), complex inversion with flanking deletion (INVDEL), and
breakend/translocation (TRA).

Done when: each new type has a truth VCF in `docs/benchmarking/sv/data/` and
a varforge config in `docs/benchmarking/sv/configs/`.

## Steps

1. Add to `docs/benchmarking/sv/data/truth_svs_vaf010.vcf`:
   - INS: insertion of 100 bp de novo sequence at a new position
   - Large DEL: 600 bp deletion
   - INVDEL: inversion with 20 bp deletions at each flank
   - TRA: breakend record (inter-chromosomal requires chr2 in ref; use
     intra-chromosomal distant rearrangement instead if only chr1 is available)
2. Create separate per-type truth VCFs for isolated testing:
   - `truth_ins_vaf010.vcf`
   - `truth_largdel_vaf010.vcf`
   - `truth_invdel_vaf010.vcf`
3. Create varforge configs in `docs/benchmarking/sv/configs/`:
   - `sim_ins_vaf010.yaml`
   - `sim_largdel_vaf010.yaml`
   - `sim_invdel_vaf010.yaml`
4. Add result directories to `.gitignore`.

## Notes

- The synthetic reference (`docs/benchmarking/sv/data/ref.fa`) is 2000 bp.
  Large DEL at 600 bp is possible but leaves a short flanking region. Consider
  placing it near position 700 so flanks are ~700 bp and ~700 bp.
- INS uses a de novo sequence not present in the reference. Varforge must
  support a SEQ field in the VCF (or an INFO=SEQ tag) to embed the inserted
  sequence.
- For INVDEL, the VCF representation is non-standard. Check varforge docs for
  whether INVDEL is a supported SVTYPE or whether it must be modelled as two
  separate records.
