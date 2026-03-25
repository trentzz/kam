# kam Pipeline: Stage 3 — Pathfind

## What this stage does

The pathfind stage builds a de Bruijn graph for each target region and enumerates all paths
through the graph from the start anchor to the end anchor. Each distinct path corresponds to a
candidate sequence: either the reference sequence or a variant.

The output is a list of scored paths, one per (target, sequence) pair. Each path carries
aggregate molecule evidence: how many molecules supported every k-mer in the path, how much
duplex support there was, and what the error rate was.

---

## De Bruijn graph

A de Bruijn graph for k-mers has one node per distinct k-mer and one directed edge per observed
(k-mer, next k-mer) overlap. Two consecutive k-mers `A` and `B` in a sequence have an edge
`A → B` if and only if the last `k-1` bases of `A` are the same as the first `k-1` bases of `B`.

For a target of length L with k-mer size k, the reference path through the graph has `L - k + 1`
nodes. At 100 bp target length and k=31, this is 70 nodes.

A variant (SNV, indel, SV) creates an alternative path through the graph:
- **SNV**: one node is replaced by a different node with the same anchor neighbours. The alt path
  has the same length as the ref path.
- **Deletion**: the alt path is shorter (missing k-mers at the deleted region).
- **Insertion**: the alt path is longer (extra k-mers at the inserted region).
- **SV**: the alt path may be dramatically shorter (large deletion), longer (tandem duplication),
  or structurally rearranged (inversion).

The graph is built from the k-mer index, not from all observed reads. Only k-mers in the index
(which are filtered to target overlaps only) appear as nodes.

### Graph construction

For each target, the graph is built by examining every k-mer in the index that overlaps the
target region. The graph adds a directed edge from k-mer A to k-mer B if the (k-1)-suffix of A
equals the (k-1)-prefix of B and both k-mers are present in the index.

The de Bruijn graph uses raw (non-canonical) k-mers internally. The suffix/prefix overlap
relationship only holds correctly for raw k-mers, not canonical ones. Canonicalisation is applied
when looking up evidence in the index.

---

## Anchor k-mers

The start anchor is the first k-mer of the target reference sequence. The end anchor is the last
k-mer. The DFS walk looks for paths from start anchor to end anchor.

### Anchor validation

Before walking, each anchor is validated against the index. An anchor is flagged as non-unique
when its molecule count in the index exceeds the uniqueness threshold
(`DEFAULT_ANCHOR_THRESHOLD = 100`).

A non-unique anchor indicates the k-mer appears in many molecules, typically because it falls in
a repetitive genomic region. Walking from a non-unique anchor can produce many spurious paths
because the graph has many outgoing edges from high-multiplicity nodes.

Anchor validation:
- Reports `start_unique` and `end_unique` flags.
- Emits a warning message if either anchor is non-unique.
- Does not stop the walk. It is informational only.

### Soft anchors (missing anchors)

If the start anchor or end anchor is not present in the index, the target has zero molecule
coverage at one end. The walk cannot proceed and the target is recorded as `no_paths`.

This is the most common cause of missed variants at very low read depth. At 200K reads, 15 ng,
2% VAF, approximately 30% of targets have a missing anchor.

---

## DFS path walk

`walk_paths` uses iterative depth-first search to enumerate all paths from start anchor to end
anchor.

### Algorithm

The DFS maintains:
- `stack`: pairs of (current node, next successor index to try). This replaces recursion.
- `current_path`: the single in-progress path, extended and retracted as DFS explores.
- `visited`: a `HashSet` of k-mers currently in `current_path`, for O(1) cycle detection.

At each DFS step:
1. Look at the top of the stack. Try the next unvisited successor of the current node.
2. If the successor is the end anchor: record the completed path, do not push to the DFS stack
   (avoid double-completion).
3. If the successor is not the end anchor and not yet visited and the path is within
   `max_path_length`: extend `current_path`, mark as visited, push a new frame.
4. If all successors have been tried: backtrack (pop from stack and `current_path`, remove from
   `visited`).

Cycle detection prevents infinite loops in graphs with repeat structure. A node already in the
current path is skipped as a successor.

### WalkConfig

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `max_path_length` | 150 | Maximum k-mers in a single path. Paths exceeding this are abandoned. |
| `max_paths` | 100 | Maximum paths to return. Walk stops once this many paths are found. |

`max_path_length` should be set to approximately `target_len / (k-1) + 50`. For a 100 bp target
with k=31, the reference path is 70 k-mers, so the default of 150 gives generous headroom for
large indels.

`max_paths` is the primary sensitivity bottleneck at complex targets. When the graph is dense
(many spurious k-mers from high depth, or repetitive sequence), the walk fills the budget of 100
paths with noise before reaching the true variant path.

### Path reconstruction

Once a path is found (a sequence of k-mers from start to end anchor), the corresponding DNA
sequence is reconstructed. The first k-mer contributes all k bases. Each subsequent k-mer
contributes only its last base (since consecutive k-mers overlap in k-1 bases). This is the
standard de Bruijn graph sequence reconstruction.

---

## SV fallback paths

For structural variants, the standard DFS walk may not find the alt path because:
- The deleted region means no k-mer chain spans the breakpoint in the standard reference-anchored
  walk.
- Inversion or duplication paths involve k-mers from distant genomic regions that are not
  present in the standard target allowlist.

When `sv_junctions` is provided, the pipeline constructs **synthetic fallback paths** directly
from the SV junction sequences, bypassing the DFS walk.

### DUP (tandem duplication) fallback

A synthetic TandemDuplication alt path is constructed from the junction sequence. The path
contains k-mers that span the duplication breakpoint, where the duplicated segment appears twice
in tandem. These k-mers are in the allowlist (added from `sv_junctions` at the index stage) and
are present in the index for samples that carry the duplication.

### InvDel fallback

An `InvDel` (inversion-with-deletion) alt path is constructed by locating the SV junction in
the index:
1. The inversion junction is searched in the index starting from the target's start anchor.
2. If found, a synthetic path is built that connects the start anchor to the inversion junction
   sequence and then to the target's end anchor.
3. The reconstructed sequence contains an inverted internal segment, with flanking deletions if
   the total length differs from the reference.

The junction can start before the nominal target window (a common structural scenario where the
inversion extends beyond the target region). The pipeline handles this by searching k-mers from
the junction sequence rather than requiring the junction to fall strictly within the target.

---

## Path scoring

After enumeration, each path is scored against the k-mer index to produce aggregate molecule
evidence.

For each k-mer in the path:
1. Compute the canonical form.
2. Look up the `MoleculeEvidence` from the index.
3. If the k-mer is absent from the index, treat it as zero evidence.

The aggregate statistics across all k-mers in the path:

| Statistic | Computation |
|-----------|-------------|
| `min_molecules` | Minimum `n_molecules` across all k-mers |
| `mean_molecules` | Mean `n_molecules` |
| `min_duplex` | Minimum `n_duplex` |
| `mean_duplex` | Mean `n_duplex` |
| `min_simplex_fwd` | Minimum `n_simplex_fwd` |
| `min_simplex_rev` | Minimum `n_simplex_rev` |
| `mean_error_prob` | Mean of `mean_base_error_prob` |

The minimum statistics are used for all downstream calling decisions, not the mean. The minimum
represents the weakest point in the path: the k-mer with the least support. A high minimum means
every k-mer in the path is well-covered; using the minimum is the conservative choice.

### Why minimum rather than sum

Using the sum of `n_simplex_fwd` across all ~70 k-mers would give a sum ~70× larger than the
actual molecule count. The Fisher's exact test for strand bias would then see 70× the apparent
read depth and flag genuine low-VAF variants as strand-biased. The minimum avoids this inflation.

### Reference path identification

After scoring, `score_and_rank_paths` identifies the reference path by comparing each path's
reconstructed sequence byte-for-byte against the target reference sequence. The first exact
match is marked `is_reference = true`.

### Variant-specific duplex

For each alt path, the variant-specific duplex count is computed separately from the overall
`min_duplex`:

1. Build a set of all canonical k-mers from the reference path.
2. For each k-mer in the alt path: if it is not in the reference k-mer set, it is
   variant-specific.
3. `min_variant_specific_duplex` = minimum `n_duplex` across all variant-specific k-mers.

The motivation: anchor k-mers at the ends of the path are shared between the reference and alt
paths. These anchor k-mers are covered by many reference molecules, so their `n_duplex` is high
(reflecting reference coverage, not variant coverage). Using `min_duplex` across all path k-mers
always gives a value at least as high as the anchor duplex count, making the LowDuplex filter
ineffective. The variant-specific count is calibrated to the actual duplex support at the variant
site.

### Path ranking

After scoring and identification, paths are sorted by `min_molecules` descending (strongest
evidence first), with ties broken by `mean_molecules`.

---

## Pathfind QC

The pathfind QC JSON (`pathfind_qc.json`) records per-target outcomes:

| Counter | Meaning |
|---------|---------|
| `n_targets` | Total targets in the FASTA |
| `no_anchor` | Targets where the start or end k-mer was absent from the index |
| `no_paths` | Targets where the walk found no paths at all |
| `ref_only` | Targets where only the reference path was found |
| `with_variants` | Targets where at least one non-reference path was found |

---

## Output file

Scored paths are serialised to a binary file in bincode format with file type tag
`FileType::ScoredPaths`. Each record is a `ScoredPathRecord`:

| Field | Type | Meaning |
|-------|------|---------|
| `target_id` | String | Target identifier |
| `sequence` | `Vec<u8>` | Reconstructed DNA sequence |
| `is_reference` | bool | True if this path matches the reference |
| `min_molecules` | u32 | Minimum molecule support |
| `mean_molecules` | f32 | Mean molecule support |
| `min_duplex` | u32 | Minimum duplex support |
| `mean_duplex` | f32 | Mean duplex support |
| `min_variant_specific_duplex` | u32 | Minimum duplex at variant-specific k-mers |
| `min_simplex_fwd` | u32 | Minimum forward simplex support |
| `min_simplex_rev` | u32 | Minimum reverse simplex support |
| `mean_error_prob` | f32 | Mean base error probability |

---

## Sensitivity limits

The two structural limits on sensitivity at this stage are:

**1. Indel end-anchor displacement.** For an indel of length ℓ, the variant path ends ℓ k-mers
before (deletion) or after (insertion) the reference path's end anchor. Reads supporting the
variant must span to this shifted position. At 2% VAF with 2M reads, 61–68% of indels fail
because variant-supporting reads do not cover the shifted anchor. This is documented in
`docs/research/depth_saturation_investigation.md`.

**2. max_paths budget.** At complex targets (high repetition, long targets, high depth), the DFS
exhausts the 100-path budget with spurious paths before finding the true variant. Increasing
`max_paths` would help but also increases runtime and produces more false positives.

Neither limit is fixed by the current implementation. The indel problem requires a relaxed
end-anchor approach; the budget problem requires reference-prioritised DFS or adaptive budgeting.
