# Troubleshooting

Common problems and how to fix them.

---

## No variants found

**Symptom**: the output file is empty or contains only the header line.

**Checks**:

1. **Target FASTA is correct.** The target sequences must be the actual genomic sequences for
   the panel windows, not just coordinates. Check that the FASTA entries are non-empty and have
   the expected lengths (typically 100–200 bp per target).

2. **Target IDs match.** If you are using tumour-informed monitoring mode (`target_variants`),
   verify that the chromosome names in the VCF match the `target_id` format in the FASTA
   headers. For example, `chr17:7674220-7674320` in the FASTA must match `chr17` in the VCF
   CHROM column.

3. **K-mer size versus read length.** At k=31, reads shorter than 31 bp produce no k-mers. If
   your reads are trimmed or short-insert, either reduce `kmer_size` or remove the length filter.
   Check `index_qc.json`: if `n_kmers_observed` is 0 or very small, no reads are overlapping
   the targets.

4. **min_family_size too high.** If `min_family_size` is set to a high value (e.g. 5) but your
   sample has low read depth, most families will be discarded. Check `assemble_qc.json` field
   `n_families_below_min_size`. If this is large, reduce `min_family_size`.

5. **Anchors missing.** Check `pathfind_qc.json` field `no_anchor`. If this is high (e.g. more
   than 50% of targets), most targets have no molecule coverage at the anchor k-mers. This
   indicates either very low sequencing depth, high UMI quality filtering (check
   `n_low_umi_quality` in `assemble_qc.json`), or a mismatch between targets and the actual
   sequenced regions.

6. **All paths are reference only.** Check `pathfind_qc.json` field `ref_only`. If this is high,
   the walker is not finding alt paths. This is expected at very low VAF or if the variant
   detection window does not overlap the actual variant. Ensure the target FASTA windows are wide
   enough to capture the variant and flanking k-mers.

---

## All variants are filtered

**Symptom**: the output contains many rows but all have a non-PASS filter label.

**Checks**:

1. **LowConfidence**: too few supporting molecules relative to background error rate.
   - Increase depth or reduce `min_alt_molecules`.
   - If most calls have n_alt = 1, the sample may have very low VAF. Consider whether
     `min_alt_duplex_for_single` (duplex confirmation of single-molecule calls) is appropriate.

2. **StrandBias**: alt allele is predominantly on one strand.
   - Expected for artefacts from oxidative damage (8-oxoG), end-repair, or PCR.
   - For simplex protocols, inherit strand bias is unavoidable. Relax
     `strand_bias_threshold` (e.g. from 0.01 to 0.001) when running without duplex.
   - For duplex protocols, high strand bias after consensus calling usually indicates a real
     artefact. Reducing the threshold further is not recommended.

3. **LowDuplex**: if you have set `min_alt_duplex = 1`, ensure your sequencing depth is high
   enough to achieve meaningful duplex fractions. At 2M reads with Twist UMI duplex, only
   6–7% of variant-site molecules are duplex at 2% VAF. At least 20 alt molecules are needed
   for ~1–2 duplex alt molecules. Either increase depth or set `min_alt_duplex = 0`.

4. **HighVaf**: if most variants are being labelled HighVaf, your `max_vaf` setting is excluding
   them. Check whether the variants are truly somatic (expected VAF < 0.3) or germline (VAF ≈
   0.5). For somatic-only calling, `max_vaf = 0.35` is a reasonable cutoff.

5. **NotTargeted**: in monitoring mode, all PASS calls must match the `target_variants` VCF.
   Verify that the VCF uses the same coordinate system as the target FASTA. A 1-based vs 0-based
   mismatch will cause all calls to be labelled NotTargeted. Use `ti_position_tolerance = 1` to
   allow a 1 bp offset if needed.

---

## High false positive rate

**Symptom**: discovery mode produces many PASS calls that do not correspond to known variants.

**Options**:

1. **Enable tumour-informed monitoring mode.** If you have a matched tumour biopsy or known
   variant list, set `target_variants` in config. This suppresses all calls not in the truth
   set and produces near-zero false positives.

2. **Tighten the strand bias filter.** The default `strand_bias_threshold = 0.01` is permissive.
   Lowering it to 0.001 eliminates more artefacts at a small sensitivity cost.

3. **Require duplex confirmation.** Set `min_alt_duplex = 1` to require at least one duplex
   molecule at the variant site. This is most effective at sequencing depths ≥5M reads or with
   protocols that achieve ≥15% duplex fraction. At 2M reads, expect to lose 25–36% of true
   positives.

4. **Increase `min_alt_molecules`.** Setting `min_alt_molecules = 3` or higher reduces false
   positives from low-count noise at the cost of sensitivity for very-low-VAF variants.

5. **Exclude germline variants.** Set `max_vaf = 0.35` to suppress calls at germline VAF
   (≈ 0.5 for heterozygous, ≈ 1.0 for homozygous).

6. **Background biology.** The 35–72 false positives in discovery mode at 2M reads are real
   biological signal: germline variants and clonal haematopoiesis in the cfDNA. These cannot be
   eliminated by tightening statistical thresholds without also reducing sensitivity. Monitoring
   mode is the designed solution.

---

## Slow performance

**Symptom**: the pipeline takes longer than expected.

**Checks**:

1. **Single-threaded.** The current release is single-threaded. The `threads` config field is
   reserved for a future parallel implementation. Runtime at 2M reads is typically 20–35 seconds
   single-core.

2. **Large input files.** If the FASTQ files are very large (>5M read pairs), runtime scales
   roughly linearly. At 10M reads, expect 60–90 seconds.

3. **max_paths budget.** If you have complex targets (long regions, repetitive sequence), the DFS
   walker spends time exploring spurious paths up to the 100-path budget. This is the dominant
   cost for panels with many high-complexity targets. There is currently no configuration to
   reduce this; reducing `max_path_length` may help in some cases.

4. **Disk I/O.** The intermediate binary files (`molecules.bin`, `index.bin`, `paths.bin`) can be
   several hundred megabytes. Ensure the output directory is on a fast local disk (SSD preferred).

---

## UMI parsing errors

**Symptom**: high `n_read_too_short` or `n_low_umi_quality` in `assemble_qc.json`.

**Checks**:

1. **Wrong UMI length.** If `umi_length` is larger than the actual UMI in your reads, every read
   will fail the `ReadTooShort` filter. Verify the chemistry configuration by examining the first
   few reads in the FASTQ file and checking where the template starts.

2. **Wrong skip length.** A skip length longer than the actual spacer will include template bases
   in the skip region. This does not cause errors (the skip bases are extracted but not used
   downstream) but may slightly reduce template length.

3. **Low UMI quality.** If many reads fail `LowUmiQuality`, the UMI bases may have inherently
   lower quality than the template (common in some protocols). Try reducing `min_umi_quality`
   to 15 or even 10. UMI errors are corrected by the Hamming clustering step anyway; UMI quality
   filtering is a pre-filter for extreme cases.

4. **Simplex protocol with duplex config.** If your chemistry does not have a reverse-strand
   UMI (simplex), set `duplex = false` in the chemistry config. The canonical UMI pair
   computation still works for simplex reads but will produce no duplex molecules, which is
   correct behaviour.

---

## Unexpected variant types

**Symptom**: variant calls have unexpected types (e.g. LargeDeletion when expecting SNVs only).

**Checks**:

1. **SV junction k-mers in the allowlist.** If `sv_junctions` is set, k-mers from those
   sequences are in the allowlist. Some of those k-mers may match reads from unrelated variants
   and produce spurious SV calls. Verify the SV junction sequences are correct for your sample.

2. **SV thresholds.** The default `sv_min_alt_molecules = 1` is permissive. In discovery mode,
   raise it to 2 or more if you see unexpected SV calls.

3. **Fusion targets.** Fusion detection requires synthetic fusion sequences in `fusion_targets`.
   If this file contains sequences with broad k-mer matches (low-complexity sequences, repetitive
   regions), spurious fusion calls may appear. Verify the fusion target sequences are
   breakpoint-specific.
