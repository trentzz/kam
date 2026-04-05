# Structural Variant Detection in kam: Design

## Motivation

The current kam pipeline detects small variants (SNVs, indels up to ~30bp) within 100bp
target windows. Structural variants (SVs) are a distinct class: rearrangements of 50bp or
larger including deletions, insertions of foreign sequence, inversions, tandem
duplications, and translocations. Each type manifests differently in the de Bruijn graph
and requires targeted design decisions.

This document:
1. Defines which SV types to support in a first phase.
2. Analyses how each type manifests in the k-mer graph.
3. Identifies what changes are needed to the index, pathfind, and call stages.
4. Describes a validation strategy using varforge-generated synthetic data.

---

## SV type taxonomy

| Type | Description | Size range | Relevance to ctDNA |
|------|-------------|------------|-------------------|
| Deletion | Contiguous genomic sequence removed | 50bp – Mb | High: common in cancer |
| Insertion | Foreign sequence inserted at a site | Any | Medium: mobile elements, templated insertions |
| Inversion | Segment reversed in place | 50bp – Mb | Medium: diagnostic in specific cancers |
| Tandem duplication | Segment copied adjacently | 50bp – Mb | High: amplification events, gene fusions |
| Translocation | Breakpoint joining two non-adjacent loci | N/A | High: BCR-ABL, EML4-ALK, etc. |

**Phase 1 target**: deletions and tandem duplications. These are the most common
structural variants in solid tumours and have the most direct representation in the de
Bruijn graph. Inversions are next. Insertions and translocations are deferred to phase 2.

**Rationale for deferral**:
- Insertions of foreign sequence: the inserted k-mers are not in the target allowlist.
  Detection requires either a de novo allowlist extension or a breakpoint assembly
  approach. This is a significant architectural change.
- Translocations: chimeric reads span two different chromosomal targets. Detection requires
  either multi-target path walks or a separate breakpoint k-mer index. Deferred.

---

## How SVs manifest in the de Bruijn graph

### Deletions

A deletion of length ℓ removes a contiguous segment of the reference. For a target window
of length L with the deletion breakpoint at position d:

```
Reference path: [k-mer at 0]...[k-mer at d-1]--[k-mer at d]...[k-mer at d+ℓ-1]--[k-mer at d+ℓ]...[end]
Alt path:       [k-mer at 0]...[k-mer at d-1]--[k-mer at d+ℓ]...[end]
```

The alt path skips ℓ k-mers. It is shorter than the reference path by ℓ nodes. The end
anchor of the alt path is still the last k-mer of the reference window — as long as the
deletion is entirely within the window and the remaining sequence is ≥ k bases.

**Current handling**: the existing pathfind stage handles this via DFS. Small deletions
(ℓ ≤ 30bp) are detectable because the alt path finds a shortcut through the graph from a
node near d to a node near d+ℓ, and the alt path still terminates at the reference end
anchor.

**Problem for large deletions**: when ℓ is large (50bp+), the deleted segment may span
most of the 100bp target window. The remaining reference sequence on both sides of the
deletion may be less than k bases. The end anchor is the last k-mer of the reference
window, which begins at position L-k. If the deletion removes the sequence from d to d+ℓ
where d+ℓ > L-k, the alt path cannot reach the reference end anchor at all — the sequence
that forms the end anchor no longer exists in the alt path.

**Fix**: for large SVs, the target window must be designed to span the breakpoint. A 100bp
window centred on the breakpoint junction captures 50bp on each side. The alt path then
traverses from a k-mer on the left flank directly to a k-mer on the right flank, skipping
the deleted sequence entirely. The DFS finds this shortcut naturally. The end anchor is
within the right flank, which is still present.

Design implication: the target FASTA for SV detection must be constructed differently from
SNV targets. For a deletion, the target should be a 100bp window spanning the breakpoint,
not a window centred on the variant position.

### Tandem duplications

A tandem duplication copies a segment of length ℓ adjacent to the original. The alt path
traverses the duplicated region twice:

```
Reference: [flank-left]--[dup-region]--[flank-right]
Alt:       [flank-left]--[dup-region]--[dup-region]--[flank-right]
```

The junction k-mers (last ℓ bases of dup-region + first ℓ bases of dup-region, i.e., the
k-mers spanning the duplication junction) are unique to the alt path. These k-mers are not
present in the reference sequence and therefore are not in the current target allowlist.

**Problem**: the allowlist is built from the reference sequence. Junction k-mers from
the duplication are not in the reference, so they are filtered out and the alt path is
invisible to the DFS.

**Fix**: augment the allowlist with junction k-mers for each known SV. When a duplication
is known (from a prior tissue biopsy), its junction k-mers can be computed and added to the
allowlist. This is compatible with the tumour-informed monitoring approach already used for
SNVs.

For de novo discovery of duplications, the junction k-mers are unknown. Detection requires
either a de novo k-mer discovery pass (looking for k-mers with high molecule support that
are not in the reference) or a suffix-tree approach. This is a phase 2 problem.

### Inversions

A segment of length ℓ is reversed in place. The alt path traverses the segment in the
reverse-complement orientation:

```
Reference: [flank-left]--[segment: ATCG...]--[flank-right]
Alt:       [flank-left]--[segment: ...CGAT (rc)]--[flank-right]
```

The junction k-mers span the flank-to-inversion boundary. These k-mers contain
reverse-complement sequence that is not in the forward reference but IS the
reverse-complement of the reference. Because k-mer canonicalisation takes `min(kmer,
rc(kmer))`, the canonical form of a k-mer in the inverted segment is the same as the
canonical form of the corresponding k-mer on the opposite strand of the reference.

**Key insight**: if the k-mer extraction from consensus reads uses canonical form, and the
inversion creates the same canonical k-mers as the reference strand from another region,
then inversion k-mers will match allowlist entries from elsewhere in the target panel. The
DFS may traverse inversion paths using these canonically-equivalent k-mers.

**Problem**: the sequence reconstruction step (converting a k-mer path back to a DNA
sequence) uses raw k-mers, not canonical k-mers, for suffix/prefix overlap. Inversion
paths will contain raw k-mers that do not overlap correctly with flanking k-mers in the raw
graph. The DFS will not find a continuous path from start to end anchor through an inverted
segment unless the graph construction explicitly represents the inversion.

**Fix**: for known inversions, construct a target sequence that represents the alt
(inverted) allele directly, not the reference. The allowlist for the inversion target
includes the junction k-mers and the inverted segment k-mers. The DFS walks this alt
target and finds the reference path as the non-inverted path. This is symmetric: both
reference and alt targets are in the FASTA.

Alternatively, extend the graph construction to handle reverse-complement edges, allowing
the DFS to traverse segments in the reverse-complement direction. This is more general but
more complex.

---

## Required changes per stage

### Assemble stage

No changes required. The assemble stage operates on raw reads and produces molecule
consensus sequences regardless of variant type.

### Index stage

**Allowlist augmentation**: the current allowlist is built from the reference sequence of
each target. For SV detection:

- **Deletions**: the deletion junction k-mer (where the left flank meets the right flank
  post-deletion) must be added. For a known deletion, this is computable: concatenate the
  left ℓ bases before the deletion with the right ℓ bases after the deletion and extract
  k-mers from the junction.

- **Tandem duplications**: the duplication junction k-mer (where the duplicated segment
  meets itself) must be added. For a known duplication, concatenate the end of the
  duplicated segment with the start of the duplicated segment and extract k-mers.

- **Inversions**: the inversion junction k-mers must be added. Concatenate the left flank
  with the reverse-complement of the inverted segment to generate junction k-mers.

Implementation: accept an optional `--sv-junctions FASTA` file containing the synthetic
junction sequences for each known SV. The index stage adds these junction k-mers to the
allowlist in addition to the reference target k-mers. This keeps the allowlist change
contained and does not affect the hot path for standard SNV indexing.

### Pathfind stage

**Target FASTA for SVs**: the target FASTA must include a synthetic sequence representing
each SV junction. For a deletion:

```
>chr17:7674100-7674300_DEL_100bp
[50bp left flank][50bp right flank after deletion]
```

This 100bp target contains the breakpoint junction. The DFS walks from the left flank
anchor to the right flank anchor, finding both the reference path (which includes the
deleted sequence) and the deletion alt path (which skips it directly). The reference path
requires more k-mers than the target window can contain if the deletion is large; the DFS
may not find it. This is acceptable: for large deletions, we care about the alt path, not
the reference path.

**max_path_length**: for large deletions, the reference path may exceed `max_path_length`
(currently 150 k-mers). The reference path through a 200bp deletion window is 200-k+1 =
170 k-mers at k=31. This is within budget. Increase `max_path_length` to 200 for SV
targets, or make it configurable per target.

**WalkConfig extension**: add an optional per-target `WalkConfig` override (or a simple
`sv_mode: bool` flag that sets `max_path_length=200` and `max_paths=200`). The current
compiled-in defaults are sufficient for SNV targets; SV targets need slightly more headroom.

### Call stage

**Variant classification**: the current `VariantType` enum covers SNV, MNV, Insertion,
Deletion, Complex. Add:

- `LargeDeletion` — deletion ≥ 50bp
- `TandemDuplication` — alt path is longer by a segment identical to a flanking segment
- `Inversion` — alt path is reverse-complement of ref path segment

These are post-hoc labels: classification happens after calling, so no changes to the
calling logic itself. The confidence, strand bias, and VAF computations remain unchanged.

**VCF output**: the VCF writer must handle SV alleles. VCF 4.3 represents SVs as:
- `<DEL>` with `SVLEN=-ℓ` for deletions
- `<DUP>` with `SVLEN=ℓ` for duplications
- `<INV>` with `SVLEN=ℓ` for inversions

For exact sequence representation (small SVs, up to 100bp), the standard `REF`/`ALT`
allele columns are sufficient. For larger SVs where the allele is too long for the field,
use the symbolic allele notation.

---

## Target FASTA design for SVs

The target FASTA file for SV detection must be constructed separately from the SNV panel.
Two options:

**Option A: combined FASTA with SV targets appended**. Add SV junction targets to the
existing SNV panel FASTA. Each SV gets one target entry representing the breakpoint region.
The target ID encodes both the genomic coordinate and the SV type:

```
>chr17:7674100-7674300_DEL_100bp
```

The pipeline processes all targets uniformly; SV targets produce longer paths or show
deletion shortcuts.

**Option B: separate SV FASTA, separate pipeline run**. Run `kam run` once for SNV targets
and once for SV targets. Merge VCF outputs. Cleaner separation; allows different
`--kmer-size` or `--max-vaf` settings per run.

Recommendation: Option B for the first implementation. Combined targets can be a
future convenience feature.

---

## Validation with varforge

Varforge generates synthetic cfDNA data with ground-truth SVs. The validation workflow:

### 1. Generate synthetic SV data

For each SV type (deletion, duplication, inversion), generate a varforge dataset:

```yaml
# svs_deletion.yaml
sample:
  name: deletion_test
  coverage: 1000
  read_length: 150
  fragment_model: cfdna

umi:
  length: 5
  duplex: true
  chemistry: twist

mutations:
  random_count: 0
  svs:
    - type: deletion
      chrom: chr17
      start: 7674150
      end: 7674250
      vaf: 0.02
    - type: deletion
      chrom: chr3
      start: 178936091
      end: 178936200
      vaf: 0.05
```

Generate samples at VAF 0.5%, 1%, 2%, 5% for each SV type. Use the same target region
coordinates as the SNV benchmark panel where possible (to test that SV detection does not
disrupt existing SNV calling).

### 2. Construct SV junction targets

For each SV in the truth VCF, extract a 100bp window centred on the breakpoint junction.
Write to `sv_targets.fa`. Add junction k-mers to the index allowlist via `--sv-junctions`.

### 3. Run kam with SV targets

```bash
kam run \
  --r1 deletion_test_R1.fastq.gz \
  --r2 deletion_test_R2.fastq.gz \
  --targets sv_targets.fa \
  --output-dir results_sv/ \
  --target-variants truth_svs.vcf \
  --output-format vcf,tsv
```

### 4. Score against truth

The existing `run_titration_batch.py` and `score_variants.py` scripts work for SVs because
they match on `(CHROM, POS, REF, ALT)`. For SV targets with symbolic alleles, update the
scoring to also match on `SVTYPE` and `SVLEN` from the INFO field.

### 5. Sensitivity targets

Initial sensitivity targets for phase 1 SV detection (2M reads, k=31):

| SV type | VAF | Target sensitivity |
|---------|-----|-------------------|
| Deletion 50–200bp | 2% | ≥ 0.5 |
| Deletion 50–200bp | 1% | ≥ 0.3 |
| Tandem duplication 50–200bp | 2% | ≥ 0.4 |
| Inversion 50–200bp | 2% | ≥ 0.4 |

These are lower than SNV targets (SNV sensitivity at 2% is 0.75–0.80) because:
- End-anchor displacement is proportionally larger for SVs than for indels.
- Junction k-mers have lower multiplicity (only alt reads contribute, not ref reads).
- The allowlist augmentation may need tuning.

---

## Implementation order

1. **Target FASTA construction script** (`benchmarking/scripts/make_sv_targets.py`):
   reads a truth VCF with SVs, extracts 100bp junction windows from the reference genome,
   writes a FASTA file with one target per SV.

2. **Junction k-mer allowlist file** (`--sv-junctions`):
   extend `kam index` to accept an optional FASTA of junction sequences to add to the
   allowlist alongside the target reference k-mers. This is the minimum viable change to
   support deletion and duplication detection.

3. **WalkConfig: per-target `max_path_length` override**:
   allow the target FASTA to annotate a target with `#maxpath=200` in the header. Parse
   this in pathfind and override WalkConfig accordingly.

4. **VCF writer: symbolic SV alleles**:
   extend the VCF output to write `<DEL>`, `<DUP>`, `<INV>` with `SVLEN` when the variant
   length exceeds 50bp (or a configurable threshold).

5. **VariantType extension**: add `LargeDeletion`, `TandemDuplication`, `Inversion`.

6. **Benchmarking with varforge**: generate synthetic SV datasets, run the above pipeline,
   score, document.

---

## Known limitations at this design stage

**Soft end-anchor problem scales with SV size**. A 100bp deletion creates a 100-k-mer
displacement. The soft end-anchor fix (allowing termination within N k-mers of the
reference end) helps small indels but does not help large SVs where the alt path ends at a
completely different genomic coordinate. The junction-window design sidesteps this by
making the alt path fit within the window.

**De novo SV discovery is not addressed**. This design covers tumour-informed monitoring
of known SVs. De novo discovery of SVs (where the breakpoint coordinates are unknown)
requires a different approach: scanning for k-mers with high molecule support that are not
in the reference. This is phase 2.

**Translocations require multi-target path walking**. A translocation breakpoint connects
k-mers from two different chromosomal loci. The current DFS operates per target. Supporting
translocations requires either a chimeric target (synthetic sequence spanning the
junction) or a two-target walk. Phase 2.

**Very large SVs (>500bp) strain the 100bp window design**. For very large deletions, the
100bp junction window still works (only the flank k-mers matter). For very large
duplications (>500bp), the duplication junction sequence is far from both flanks and may
be hard to distinguish from reference at low VAF. Sensitivity will drop for large SVs.

---

## Open questions (DESIGN_QUESTION tags)

1. Should the SV junction allowlist be a separate file (`--sv-junctions`) or should the
   index stage automatically extract junction k-mers from known SV annotations embedded in
   the target FASTA header?

2. Should the pathfind stage receive per-target WalkConfig overrides via FASTA header
   annotation, or via a separate config file mapping target IDs to parameters?

3. For VCF output of large deletions, should the full deleted sequence be written in the
   REF field (as VCF spec allows), or should symbolic `<DEL>` notation be used? Symbolic
   is standard for SVs but makes scoring more complex.

4. How should the truth set for SV benchmarking be structured? The current
   `score_variants.py` matches on `(CHROM, POS, REF, ALT)`. SV VCFs use symbolic alleles
   and `SVTYPE`/`SVLEN` INFO fields. A new or extended scorer is needed.
