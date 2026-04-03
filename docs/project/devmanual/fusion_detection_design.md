# SV-EXP-004: Multi-Locus Fusion and Translocation Detection

## 1. Problem Statement

kam currently detects structural variants within a single target region. Each
target is a contiguous genomic sequence. The walker finds paths from the start
anchor to the end anchor of that one region. Junction k-mers span breakpoints
within the target window.

Fusions and translocations break this model. A fusion breakpoint joins two
distinct genomic loci, for example BCR (chr22) and ABL1 (chr9). A read
spanning the breakpoint contains k-mers from both loci. No single target
region contains the fusion junction sequence. The current pipeline has no way
to:

1. Define a junction that spans two targets.
2. Build a graph path that starts in one target and ends in another.
3. Attribute molecule evidence from two different loci to a single variant.
4. Emit VCF BND records with breakpoint coordinates on two chromosomes.

This document describes a design for fusion detection that works within kam's
existing alignment-free k-mer framework.


## 2. Proposed Approach

### Core idea: fusion targets

Introduce a new concept, the "fusion target". A fusion target is a synthetic
sequence formed by concatenating the breakpoint-adjacent segments of two
partner loci. It is not a real genomic region. It is the sequence that fusion-
bearing reads would produce around the breakpoint.

For a BCR-ABL1 fusion with a known breakpoint:

```
Partner A (BCR):   ...NNNNNNN[last 50bp]
Partner B (ABL1):  [first 50bp]NNNNNNN...
Fusion target:     [last 50bp of A][first 50bp of B]
```

The fusion target is 100 bp. It serves the same role as a normal target: the
walker finds paths from its start anchor to its end anchor. The reference path
through this synthetic target is the fusion sequence itself. An alt path (the
wildtype) would not exist in the graph because wild-type reads do not span
both loci.

This inverts the ref/alt logic compared to normal targets. For a normal
target, the reference path is the wild-type genome. For a fusion target, the
"reference" path is the fusion allele and wild-type molecules do not
contribute. The caller must handle this inversion.

### Why not cross-target walking?

An alternative is to modify the walker to start in one target's graph and end
in another target's graph. This requires:

- Merging two target graphs or building a global graph with cross-target edges.
- Defining anchor pairs across targets (which start anchor, which end anchor).
- Changing `walk_paths` to accept anchors from different targets.
- Handling the combinatorial explosion of possible start/end pairs.

This is invasive and complex. The synthetic fusion target approach reuses the
existing single-target pipeline unchanged. The walker, scorer, and caller all
operate on a single target as before. Only the input preparation and output
interpretation change.

### Why not junction-only detection?

Another alternative is to detect fusions purely from junction k-mer presence
in the index, without graph walking. Count how many molecules contain fusion
junction k-mers and call based on that count alone.

This is simpler but loses the path-level evidence that makes kam's calls
robust. Without a complete path, there is no sequence reconstruction, no
alignment to a reference, and no way to distinguish a real fusion from a
chimeric library artefact that happens to contain junction k-mers. The graph
walk provides the full alt sequence, enabling quality-weighted scoring and
breakpoint-precise calls.

Junction-only detection could serve as a fast pre-screen (see Section 7).


## 3. Data Flow

### 3.1 Input preparation

The user provides fusion definitions in a FASTA file via `--sv-junctions` (or
a new `--fusion-junctions` flag, see Section 5). Each fusion junction entry
has an ID encoding both partners:

```
>BCR_ABL1__chr22:23632500__chr9:130854000__fusion
ATGCGATCG...50bp_from_BCR...TGCAATGCC...50bp_from_ABL1...
```

The ID format encodes: `partnerA_partnerB__chromA:posA__chromB:posB__fusion`.
The `__fusion` suffix distinguishes fusion junctions from single-locus SV
junctions.

Alternatively, a separate fusion targets FASTA is provided via
`--fusion-targets`, containing the synthetic 100 bp fusion target sequences
described above. This is the preferred approach because it cleanly separates
fusion targets from single-locus SV junctions.

### 3.2 Allowlist construction

Fusion target k-mers are added to the allowlist alongside normal target
k-mers. This ensures that reads spanning the fusion breakpoint are identified
as on-target during the two-pass indexing step.

Critically, the fusion target shares k-mers with both partner loci. The first
50 bp overlap with partner A's normal target; the last 50 bp overlap with
partner B's normal target. Molecules from either locus that happen to contain
these shared k-mers are already indexed. Only the junction-spanning k-mers
(those that cross the concatenation point) are unique to the fusion target.
These are the k-mers that provide fusion-specific evidence.

### 3.3 Indexing

No changes to the indexing stage. The two-pass indexing already works:

1. Pass 1: any molecule with a k-mer in the allowlist is flagged as on-target.
   Molecules spanning the fusion breakpoint contain junction k-mers that are
   in the allowlist, so they are captured.
2. Pass 2: all k-mers from on-target molecules are indexed, including the
   fusion junction k-mers and any flanking k-mers.

### 3.4 Graph construction

The de Bruijn graph is built from all on-target raw k-mers. Fusion junction
k-mers are included because they were indexed in pass 2. The graph now
contains edges that connect partner A k-mers to partner B k-mers through the
junction.

### 3.5 Graph walking

The fusion target is walked like any normal target. The walker uses the first
k-mer of the fusion target as the start anchor and the last k-mer as the end
anchor.

- Start anchor: a k-mer from partner A (the 5' side of the fusion).
- End anchor: a k-mer from partner B (the 3' side of the fusion).

If fusion-bearing molecules are present, the graph contains a path from the
start anchor through the junction k-mers to the end anchor. This is the
fusion path.

If no fusion-bearing molecules are present, the junction k-mers are absent
from the graph. There is no path from start to end because partner A k-mers
are not connected to partner B k-mers. The walker returns an empty result,
which is the expected outcome for a negative sample.

### 3.6 Scoring

Path scoring works unchanged. Each k-mer in the fusion path is looked up in
the canonical index. Molecule evidence is aggregated. The weakest k-mer is
typically one of the junction-spanning k-mers (they are supported only by
fusion-bearing molecules, not by the bulk of wild-type molecules).

### 3.7 Calling

The caller needs a modified reference/alt interpretation for fusion targets.
For a normal target:

- Reference path = wild-type sequence.
- Alt path = variant sequence.

For a fusion target:

- There is no wild-type reference path (wild-type molecules do not span both
  loci).
- The "reference" for the fusion target is the fusion sequence itself.
- There is no "alt" in the traditional sense.

The caller must detect that this is a fusion target and handle it specially:

1. If the walker finds a path matching the fusion target sequence, a fusion is
   detected.
2. The VAF is estimated differently. The denominator is the total molecule
   count at the partner loci (from their normal targets), not the molecule
   count on the fusion target. This is because wild-type molecules do not
   produce k-mers on the fusion target, so using the fusion target's own
   depth as the denominator would overestimate VAF.
3. The numerator is the molecule count at the junction-spanning k-mers (the
   variant-specific k-mers of the fusion path).

### 3.8 VCF output

Fusions use BND notation (VCF 4.3 Section 5.4). Each fusion produces two
VCF records, one for each breakpoint partner:

```
chr22  23632500  bnd_BCR_ABL1_1  G  G]chr9:130854000]  .  PASS  SVTYPE=BND;MATEID=bnd_BCR_ABL1_2;...
chr9   130854000 bnd_BCR_ABL1_2  T  ]chr22:23632500]T  .  PASS  SVTYPE=BND;MATEID=bnd_BCR_ABL1_1;...
```

The existing `write_sv_vcf_record` function uses symbolic `<BND>` notation.
For proper BND output, the ALT field must contain the partner breakpoint
coordinate and the orientation brackets. This requires the writer to know both
partner coordinates, which it currently does not.


## 4. Changes Required

### 4.1 kam-index (allowlist.rs)

No code changes. Fusion target sequences are passed as additional targets to
`build_allowlist`. The existing function handles them correctly.

### 4.2 kam-pathfind (walk.rs, score.rs)

No code changes. The walker and scorer operate on a single target. Fusion
targets are walked and scored identically to normal targets.

### 4.3 kam-call

#### caller.rs: classify_variant

`classify_variant` currently infers the variant type from `ref_seq` and
`alt_seq` byte comparison. It cannot detect fusions because fusions are not
distinguishable from other variant types by sequence comparison alone.

Add a new parameter or a wrapper function:

```rust
pub fn classify_variant_with_context(
    ref_seq: &[u8],
    alt_seq: &[u8],
    target_id: &str,
) -> VariantType {
    if is_fusion_target(target_id) {
        return VariantType::Fusion;
    }
    classify_variant(ref_seq, alt_seq)
}
```

The `is_fusion_target` check parses the target ID for the `__fusion` suffix
(or checks a flag passed from the caller).

#### caller.rs: call_variant (fusion-specific VAF)

Add a `FusionContext` struct that carries the denominator molecule count from
the partner loci:

```rust
pub struct FusionContext {
    /// Total molecules at partner A's normal target.
    pub partner_a_depth: u32,
    /// Total molecules at partner B's normal target.
    pub partner_b_depth: u32,
}
```

When calling a fusion variant, use `(partner_a_depth + partner_b_depth) / 2`
as the denominator for VAF estimation instead of the reference path's molecule
count (which is zero or near-zero for a fusion target).

Add a new entry point:

```rust
pub fn call_fusion_variant(
    target_id: &str,
    fusion_evidence: &PathEvidence,
    fusion_context: &FusionContext,
    fusion_seq: &[u8],
    config: &CallerConfig,
) -> VariantCall { ... }
```

#### allele.rs: extract_minimal_allele

Fusions do not use `extract_minimal_allele`. BND records have a different
structure: the REF is a single anchor base, and the ALT encodes the partner
coordinate. Skip allele extraction for fusion types.

#### output.rs: write_sv_vcf_record

Replace the `<BND>` symbolic allele with proper BND notation. The function
needs the partner breakpoint coordinates, parsed from the fusion target ID.

Add a `FusionBreakpoint` struct:

```rust
pub struct FusionBreakpoint {
    pub chrom_a: String,
    pub pos_a: u64,
    pub chrom_b: String,
    pub pos_b: u64,
    pub id_prefix: String,
}
```

When `variant_type == Fusion`, emit two records:

```rust
fn write_fusion_vcf_records(
    call: &VariantCall,
    breakpoint: &FusionBreakpoint,
    filter_str: &str,
    base_info: &str,
    writer: &mut dyn Write,
) -> io::Result<()> {
    // Record 1: partner A
    writeln!(writer, "{}\t{}\t{}_1\t{}\t{}]{}:{}]\t.\t{}\t{};SVTYPE=BND;MATEID={}_2",
        breakpoint.chrom_a, breakpoint.pos_a,
        breakpoint.id_prefix,
        ref_base_a,
        ref_base_a, breakpoint.chrom_b, breakpoint.pos_b,
        filter_str, base_info, breakpoint.id_prefix)?;
    // Record 2: partner B
    writeln!(writer, "{}\t{}\t{}_2\t{}\t]{}:{}]{}\t.\t{}\t{};SVTYPE=BND;MATEID={}_1",
        breakpoint.chrom_b, breakpoint.pos_b,
        breakpoint.id_prefix,
        ref_base_b,
        breakpoint.chrom_a, breakpoint.pos_a, ref_base_b,
        filter_str, base_info, breakpoint.id_prefix)?;
    Ok(())
}
```

### 4.4 kam/src/commands/run.rs

This is where the largest changes occur.

#### Fusion target loading

Add logic to detect fusion targets in the targets FASTA (by ID suffix) or
load them from a separate `--fusion-targets` file. Track which target IDs are
fusion targets in a `HashSet<String>`.

#### Partner depth tracking

After pathfinding on normal targets, collect the reference path molecule depth
for each target. Store in a `HashMap<String, u32>`:

```rust
let mut target_depths: HashMap<String, u32> = HashMap::new();
for (target_id, scored_paths) in &all_scored {
    if let Some(ref_path) = scored_paths.iter().find(|p| p.is_reference) {
        target_depths.insert(
            target_id.clone(),
            ref_path.aggregate_evidence.mean_molecules.round() as u32,
        );
    }
}
```

#### Fusion calling

Process fusion targets after normal targets. For each fusion target:

1. Walk the fusion target graph (same as any target).
2. If a path is found, look up partner depths from `target_depths`.
3. Call using `call_fusion_variant` with the partner depth as the denominator.

#### Processing order

Normal targets must be processed first so that partner depths are available
when fusion targets are processed. Split the targets list into normal and
fusion, process normal first, then fusion.

### 4.5 Configuration (kam/src/config.rs, kam/src/cli.rs)

Add a new optional input field:

```rust
pub struct InputConfig {
    // ... existing fields ...
    pub fusion_targets: Option<PathBuf>,
}
```

And a CLI flag:

```
--fusion-targets <PATH>  FASTA of synthetic fusion target sequences
```


## 5. Configuration Format

### 5.1 Fusion target FASTA

Each entry is a synthetic sequence representing the fusion breakpoint region.
The FASTA header encodes both partners and their coordinates:

```
>BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion
ACGATCGATCG...50bp_BCR...TGCAATGCCAA...50bp_ABL1...
```

Header format:
```
>{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion
```

Fields:
- `name`: human-readable fusion name (e.g. `BCR_ABL1`).
- `chromA:startA-endA`: genomic coordinates of the partner A segment.
- `chromB:startB-endB`: genomic coordinates of the partner B segment.
- `__fusion`: suffix that marks this as a fusion target.

The sequence is the concatenation of partner A's breakpoint-adjacent segment
and partner B's breakpoint-adjacent segment. The breakpoint is at the exact
concatenation point.

### 5.2 config.toml

```toml
[input]
r1 = "sample_R1.fastq.gz"
r2 = "sample_R2.fastq.gz"
targets = "panel_targets.fa"
sv_junctions = "sv_junctions.fa"
fusion_targets = "fusion_targets.fa"    # new field
```

### 5.3 Tumour-informed mode

For tumour-informed monitoring, the `--target-variants` VCF can include BND
records. The targeting filter must recognise BND entries and match them against
fusion calls by partner coordinates rather than by single-locus position and
allele.

### 5.4 De novo fusion detection

De novo fusion detection (discovering fusions not specified in the panel) is
out of scope for the initial implementation. It requires enumerating all
possible target pairs and building synthetic fusion targets for each, which is
combinatorially expensive for large panels.

A future extension could use a two-step approach:

1. Pre-screen: for each molecule, check whether it contains k-mers from two
   different normal targets. Flag these as potential fusion-spanning molecules.
2. Targeted walk: for each flagged pair of targets, build a synthetic fusion
   target and run the standard pipeline.

This is deferred to a later phase.


## 6. Edge Cases and Limitations

### 6.1 Reciprocal fusions

A translocation produces two fusion junctions (e.g. BCR-ABL1 and ABL1-BCR).
Each requires its own fusion target entry. The user must provide both if both
are clinically relevant.

### 6.2 Multiple breakpoints per gene pair

The same gene pair can have different breakpoints (e.g. BCR exon 13 vs exon
14 fusions with ABL1). Each distinct breakpoint requires a separate fusion
target entry.

### 6.3 Fusion target length

The fusion target must be long enough for the walker to function. With k=31,
a 100 bp fusion target yields 70 k-mers, of which approximately 30 are
junction-spanning (unique to the fusion). Shorter targets reduce sensitivity.
Longer targets (150 bp, 75 bp per partner) provide more anchor context but
require reads to span further from the breakpoint.

Recommended length: 100 bp (50 bp per partner). This matches the typical
insert size of ctDNA fragments (150-170 bp), ensuring that most breakpoint-
spanning reads cover at least 50 bp on each side.

### 6.4 Shared k-mers between partners

If partners A and B share sequence near the breakpoint (e.g. microhomology at
the junction), some "junction" k-mers may already exist in wild-type
molecules. This inflates the apparent molecule count at those k-mers. The
variant-specific k-mer filtering (excluding anchor k-mers shared with the
reference) partially mitigates this, but only if the shared k-mers are
identified correctly.

Mitigation: when building the fusion target's reference k-mer set, include
k-mers from both partner A and partner B normal targets. Any k-mer present in
either partner's normal target is an anchor k-mer, not a fusion-specific
k-mer.

### 6.5 VAF denominator

The correct denominator for fusion VAF is not obvious. Options:

1. Mean depth across both partners: `(depth_A + depth_B) / 2`.
2. Minimum depth: `min(depth_A, depth_B)`.
3. Depth at the specific breakpoint region of each partner.

Option 1 is the simplest and most robust. Option 3 is more precise but
requires per-position depth tracking that kam does not currently provide.

Start with option 1. Revisit if validation shows systematic VAF bias.

### 6.6 No wild-type reference path

For fusion targets, the walker may find zero paths (no fusion) or one path
(fusion present). It will not find a wild-type reference path because
wild-type molecules do not span both loci. The `score_and_rank_paths` function
expects to identify a reference path. For fusion targets, all paths are "alt"
paths. The caller must handle this: skip the reference-path identification
step and treat the highest-evidence path as the fusion allele.

### 6.7 Chimeric artefacts

Library preparation chimeras can produce reads that span two loci without a
true genomic fusion. These appear as low-level fusion signal. The duplex
requirement provides protection: a chimeric artefact on one strand is unlikely
to be independently confirmed on the complementary strand. The existing duplex
filter applies without modification.

### 6.8 BND orientation

VCF BND notation encodes four orientations using bracket placement:

- `t]p]`: partner A forward, partner B forward.
- `t[p[`: partner A forward, partner B reverse.
- `]p]t`: partner A reverse, partner B forward.
- `[p[t`: partner A reverse, partner B reverse.

The fusion target ID must encode the orientation, or the writer must infer it
from the partner coordinate ranges (ascending = forward, descending =
reverse). For the initial implementation, assume forward-forward orientation
and extend later.


## 7. Implementation Order

### Phase 1: Fusion target input and allowlisting

1. Define the fusion target FASTA format and ID parsing function.
2. Add `--fusion-targets` CLI flag and config field.
3. Load fusion targets and add their k-mers to the allowlist.
4. Add `is_fusion_target` helper that checks the `__fusion` suffix.
5. Unit tests for ID parsing and allowlist construction.

Estimated scope: `kam/src/cli.rs`, `kam/src/config.rs`,
`kam/src/commands/run.rs`. Minimal risk, no existing behaviour changes.

### Phase 2: Fusion walking and scoring

1. Walk fusion targets through the existing graph (no walker changes needed).
2. Handle the "no reference path" case in `score_and_rank_paths`: when no path
   matches the fusion target sequence, treat the fusion target sequence as the
   reference for scoring purposes. The fusion path IS the reference.
3. Collect partner depths from normal target results.
4. Unit tests with synthetic fusion target sequences.

Estimated scope: `kam/src/commands/run.rs` (processing order split),
`kam-pathfind/src/score.rs` (minor: allow caller to specify reference
sequence explicitly).

### Phase 3: Fusion calling

1. Add `call_fusion_variant` to `kam-call/src/caller.rs`.
2. Add `classify_variant_with_context` or pass variant type override.
3. Implement fusion-specific VAF with partner depth denominator.
4. Unit tests for VAF estimation with fusion context.

Estimated scope: `kam-call/src/caller.rs`. Moderate complexity.

### Phase 4: VCF BND output

1. Add `FusionBreakpoint` struct and ID parser.
2. Implement `write_fusion_vcf_records` with proper BND notation.
3. Add MATEID header and paired record emission.
4. TSV/JSON output: include partner coordinates as extra fields.
5. Unit tests for VCF output format.

Estimated scope: `kam-call/src/output.rs`, `kam-call/src/allele.rs`. Moderate
complexity.

### Phase 5: Benchmarking with varforge

1. Add fusion targets to the varforge benchmark panel.
2. Generate synthetic fusion reads with known breakpoints.
3. Validate sensitivity and specificity across a VAF range (0.1% to 5%).
4. Compare junction k-mer depth vs path-based calling.

### Phase 6 (deferred): De novo fusion detection

1. Implement molecule-level multi-target overlap detection.
2. Build candidate fusion target pairs from overlap evidence.
3. Run the standard fusion pipeline on candidate pairs.

This phase is deferred until tumour-informed fusion detection is validated.
