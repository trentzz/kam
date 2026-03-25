# kam Pipeline: Stage 1 — Assemble

## What this stage does

The assemble stage converts raw paired-end FASTQ files into a set of molecules. Each molecule
represents one original DNA fragment as identified by its unique UMI pair. By the end of this
stage, all downstream analysis works in molecules, not reads.

The six steps run in sequence for every batch of read pairs:

1. Parse each R1/R2 pair: extract the UMI, skip, and template regions.
2. Group read pairs by canonical UMI.
3. Cluster UMIs within Hamming distance 1 to correct sequencing errors in the UMI itself.
4. Split clusters by endpoint fingerprint to detect UMI collisions.
5. Call single-strand consensus per strand family, then duplex consensus if both strands are
   present.
6. Emit one `Molecule` per fingerprint-compatible sub-group.

---

## Read structure and chemistry

### Twist UMI duplex (default)

The Twist UMI duplex panel uses the `5M2S+T` read structure on both R1 and R2:

```
R1: [NNNNN][NN][TTTTTT...T]
     0   4  5 6  7 onwards
     UMI  skip  template

R2: [NNNNN][NN][TTTTTT...T]
     0   4  5 6  7 onwards
     UMI  skip  template
```

- **UMI** (bases 0–4): 5 random bases. These identify the original DNA strand.
- **Skip** (bases 5–6): 2 monotemplate (fixed) bases. These are spacers introduced by the
  library prep. Their values can be used as a QC signal: deviations from the expected base at
  these positions indicate a chemistry or parsing problem.
- **Template** (bases 7+): the genomic sequence. This is the sequence that feeds into k-mer
  extraction and consensus calling.

### Configurable UMI and skip lengths

The UMI length and skip length are configurable via `config.toml` or CLI flags. This supports
other protocols:

```toml
[chemistry]
umi_length = 12   # e.g. a 12 bp UMI protocol
skip_length = 0   # no spacer region
duplex = false    # simplex (single-strand) chemistry
```

The `template_start()` offset is always `umi_length + skip_length`. For Twist: `5 + 2 = 7`.

The parser applies this split to every read pair before any filtering or grouping.

---

## Step 1: Parsing read pairs

`parse_read_pair` applies three filters in order. If any filter fires, the read pair is dropped
and the corresponding counter in `ParseStats` is incremented.

### Filter 1: ReadTooShort

Checks that both R1 and R2 are at least `umi_length + skip_length` bases long. If either read is
shorter, it is physically impossible to extract the full UMI and skip region. These reads are
almost always the result of a truncated FASTQ record.

Drop reason: `ReadTooShort`. Detail field: `"r1_len=N,r2_len=N,min=M"`.

### Filter 2: TemplateTooShort

Only applied when `min_template_length` is set. If the extracted template is shorter than the
configured minimum on either R1 or R2, the pair is dropped.

The default is no minimum. The parser is k-agnostic by design: a template of 0 bp passes this
filter when no minimum is set. Whether a k-mer can be extracted from the template is determined
later at the indexing stage.

Drop reason: `TemplateTooShort`. Detail field: `"r1_template=Nbp,r2_template=Nbp,min=Mbp"`.

### Filter 3: LowUmiQuality

Only applied when `min_umi_quality` is set. Checks every base quality in the UMI on both R1 and
R2. Quality bytes are Phred+33 encoded. If any UMI base has a Phred score below the threshold,
the pair is dropped.

The default is Phred 20. Set to 0 to disable this filter.

Drop reason: `LowUmiQuality`. Detail field: `"r1_umi_base=N,quality=Q,min=M"`.

### Canonical UMI and strand assignment

After all filters pass, the parser computes the canonical UMI pair and strand assignment.

The canonical pair is `min(umi_r1 + umi_r2, umi_r2 + umi_r1)` under lexicographic ordering:

- If `umi_r1 < umi_r2`: the canonical pair is `(umi_r1, umi_r2)` and R1 came from the forward
  strand.
- If `umi_r1 > umi_r2`: the canonical pair is `(umi_r2, umi_r1)` and R1 came from the reverse
  strand.
- If `umi_r1 == umi_r2`: the canonical pair is `(umi_r1, umi_r2)` and the strand is Forward by
  convention.

This canonical form is strand-agnostic: a forward strand read and the corresponding reverse
strand read from the same fragment both produce the same canonical UMI pair.

### Output: ParsedReadPair

A successfully parsed read pair produces a `ParsedReadPair` with:

| Field | Content |
|-------|---------|
| `umi_r1` | UMI bases from R1 (length = `umi_length`) |
| `umi_r2` | UMI bases from R2 (length = `umi_length`) |
| `skip_r1` | Skip bases from R1 (length = `skip_length`) |
| `skip_r2` | Skip bases from R2 (length = `skip_length`) |
| `template_r1` | Genomic sequence from R1 (bases `template_start` onwards) |
| `template_r2` | Genomic sequence from R2 (bases `template_start` onwards) |
| `qual_r1` | Base qualities matching `template_r1` |
| `qual_r2` | Base qualities matching `template_r2` |
| `umi_qual_r1` | Quality bytes for the R1 UMI |
| `umi_qual_r2` | Quality bytes for the R2 UMI |
| `canonical_umi` | Strand-agnostic UMI pair (`CanonicalUmiPair`) |
| `strand` | `Forward` or `Reverse` (which strand R1 came from) |

### ParseStats

| Counter | Meaning |
|---------|---------|
| `n_processed` | Total read pairs examined |
| `n_passed` | Read pairs that parsed successfully |
| `n_read_too_short` | Dropped for ReadTooShort |
| `n_template_too_short` | Dropped for TemplateTooShort |
| `n_low_umi_quality` | Dropped for LowUmiQuality |

---

## Step 2: UMI grouping (exact matching)

After parsing, all read pairs are grouped by canonical UMI using a `HashMap<CanonicalUmiPair, usize>`.
The hash lookup is O(1), so grouping N read pairs into U unique UMIs takes O(N) time.

For each unique canonical UMI, a read count is maintained. This count is used in Step 3 to direct
the clustering (higher-count UMIs act as seeds).

---

## Step 3: Hamming-distance UMI clustering (error correction)

With 5-base UMIs (4^5 = 1,024 possible values per strand), sequencing errors in the UMI are the
dominant source of UMI over-counting. A single-base sequencing error produces an apparent new UMI
that actually belongs to an existing molecule.

With longer UMIs (e.g. 12 bp), the collision space is larger but errors are still corrected by
this step.

### Algorithm: directional Hamming-1 clustering

The Hamming distance between two canonical UMI pairs is the total number of positions that differ
across both arms (maximum: 2 × `umi_length`).

Clustering uses a directional strategy:
1. Sort all unique UMIs by read count, descending.
2. For each UMI (in count order): if it has not yet been assigned to a cluster, start a new
   cluster with this UMI as seed. Scan all existing cluster seeds: if this UMI is within
   `max_hamming_distance` of a seed, absorb it into that cluster.
3. Once a UMI is absorbed, it cannot become a seed.

There is no transitive chaining. If A absorbs B and B is close to C, A only absorbs C if A and C
are directly within `max_hamming_distance`. This prevents long chains from collapsing into one
mega-cluster.

The default `max_hamming_distance` is 1.

### Hash-partition speedup

Naively, clustering U unique UMIs is O(U²). At 2M reads (350K unique molecules), this would be
too slow.

The optimisation partitions unique UMIs into 64 buckets based on the first 3 bases of `umi_a`.
There are 4^3 = 64 possible 3-base prefixes, so each bucket contains roughly 1/64 of all UMIs.
Clustering runs independently within each bucket, reducing the cost to approximately O(U²/64).

**Trade-off**: if the single-base mismatch falls in `umi_a[0..3]` (roughly 3/10 of all
Hamming-1 pairs), the two UMIs land in different buckets and are not merged. This ~30% miss rate
is acceptable: the unmerged UMIs remain as separate clusters, contributing independently to
variant evidence. The result is slightly more unique molecules (a conservative estimate), not a
loss of real variant signal.

At 2M reads, the speedup is ~64× with negligible impact on sensitivity.

---

## Step 4: UMI collision detection via endpoint fingerprinting

Two different DNA molecules can, by chance, share the same canonical UMI pair. With 5-base UMIs
(1,024 possible values per strand), the chance of a UMI collision in a cluster of 350K molecules
is non-negligible.

Endpoint fingerprinting detects collisions by checking whether reads in the same UMI cluster
actually came from the same genomic location.

### Fingerprint computation

A 64-bit fingerprint is computed from the template endpoints of a read pair. The fingerprint
packs 8 bases from each of 4 positions:

```
bits [63:48]  R1 first 8 bases
bits [47:32]  R1 last  8 bases
bits [31:16]  R2 first 8 bases
bits [15:0]   R2 last  8 bases
```

Each base is 2-bit encoded: A=00, C=01, G=10, T=11. N is treated as A.

Two fingerprints are compatible if the Hamming distance of their XOR is at most 4. This tolerates
approximately 2 base-level sequencing errors distributed across all 4 endpoints.

### Duplex-aware compatibility

For a genuine duplex pair, the forward and reverse reads produce complementary fingerprints
(upper and lower 32-bit halves are swapped). A direct XOR would exceed the threshold of 4 and
incorrectly split the pair. The compatibility check also accepts
`popcount(fp1 XOR rotate_left(fp2, 32)) <= 4`. This groups genuine forward/reverse pairs into
a single molecule.

### Splitting algorithm

Reads within a Hamming cluster are split into fingerprint-compatible sub-groups using a greedy
pass. Each read is added to the first existing group whose representative fingerprint is
compatible. If no compatible group exists, a new group starts.

Each extra sub-group (beyond the first) counts as one detected UMI collision in
`AssemblyStats::n_umi_collisions_detected`.

---

## Step 5: Consensus calling

For each fingerprint-compatible sub-group, reads are separated into forward-strand and
reverse-strand families. Consensus is called independently for each strand, then the two strands
are combined for duplex consensus.

### Single-strand consensus (SSC)

SSC collapses a family of reads from one strand into a per-base consensus using quality-weighted
majority vote.

At each position:
1. Convert each base's Phred quality score to error probability: `p = 10^(-Q/10)`.
2. Compute the vote weight for the called base: `w = 1 - p`.
3. Sum weights per nucleotide (A, C, G, T).
4. The consensus base is the nucleotide with the highest total weight (argmax).
5. The consensus error probability is `1 - winner_weight / total_weight`.
   - If all bases at a position are below `min_base_quality`, the position is masked:
     base = `N`, error probability = 1.0.
6. The output quality is capped at `max_consensus_quality` (default: Q60).

#### ConsensusConfig

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `min_base_quality` | 10 | Phred score below which a base is excluded from voting |
| `max_consensus_quality` | 60 | Maximum output Phred quality |

### Duplex consensus

Duplex consensus combines the forward-strand SSC and reverse-strand SSC into a final base call.

At each position:
- If both strands agree: `error_prob = fwd.error_prob * rev.error_prob`. The position is marked
  `is_duplex = true`.
- If both strands disagree:
  - `PickBest` (default): use the base from the strand with lower error probability.
    `is_duplex = false`.
  - `MaskAsN`: set base to `N`, error probability to 1.0. `is_duplex = false`.

### Duplex eligibility

A molecule qualifies for duplex consensus only when both the forward and reverse strand families
have at least `min_duplex_reads` reads. The default `min_duplex_reads = 1` means that a single
read on each strand is sufficient.

If only one strand is present, only SSC is produced. The molecule is classified as simplex.

### Family size filter

Before consensus is called, the total family size is checked against `min_family_size`. Families
below this threshold are dropped. The default `min_family_size = 1` accepts singletons.

Setting `min_family_size = 2` requires at least two reads per molecule. This eliminates true
singletons (one-read molecules), which have no error correction, at the cost of some sensitivity
at low read depth.

---

## Step 6: Molecule construction

Each surviving fingerprint sub-group becomes one `Molecule`:

| Field | Content |
|-------|---------|
| `id` | 64-bit stable hash of the canonical UMI pair |
| `umi_fwd` | The `umi_a` (lexicographically smaller) arm |
| `umi_rev` | The `umi_b` (lexicographically larger) arm |
| `consensus_fwd` | `ConsensusRead` from the forward strand, if forward reads were present |
| `consensus_rev` | `ConsensusRead` from the reverse strand, if reverse reads were present |
| `duplex_consensus` | `ConsensusRead` from duplex calling, if both strands qualified |
| `evidence` | `None` at this stage; populated by the index stage |

The `FamilyType` classification:

| Family type | Condition |
|-------------|-----------|
| `Duplex` | Both forward and reverse reads present |
| `SimplexFwd` | Forward reads only (≥2 reads) |
| `SimplexRev` | Reverse reads only (≥2 reads) |
| `Singleton` | Exactly one read total (either strand) |

Molecules are sorted by `id` (hash of the canonical UMI) for deterministic output order.

---

## AssemblyStats

All counters are accumulated during `assemble_molecules` and written to the assemble QC JSON.

| Counter | Meaning |
|---------|---------|
| `n_molecules` | Total molecules assembled after all filters |
| `n_duplex` | Molecules with both forward and reverse consensus |
| `n_simplex_fwd` | Molecules with forward reads only (≥2) |
| `n_simplex_rev` | Molecules with reverse reads only (≥2) |
| `n_singletons` | Molecules with exactly one read |
| `n_families_below_min_size` | Families dropped by `min_family_size` |
| `n_umi_collisions_detected` | Extra fingerprint splits (each = one detected collision) |

---

## AssemblerConfig summary

| Field | Default | Config key |
|-------|---------|-----------|
| `umi_length` | 5 | `[chemistry] umi_length` |
| `skip_length` | 2 | `[chemistry] skip_length` |
| `duplex` | true | `[chemistry] duplex` |
| `min_umi_quality` | 20 | `[chemistry] min_umi_quality` |
| `min_template_length` | (none) | `[chemistry] min_template_length` |
| `min_family_size` | 1 | `[assembly] min_family_size` |
| `max_hamming_distance` | 1 | Not currently exposed |
| `min_duplex_reads` | 1 | Not currently exposed |
| `consensus.min_base_quality` | 10 | Not currently exposed |
| `consensus.max_consensus_quality` | 60 | Not currently exposed |
| `disagreement_strategy` | PickBest | Not currently exposed |

---

## Output file

The assembled molecules are serialised to a binary file in bincode format with file type tag
`FileType::Molecules`. This file is the input to the index stage.

The assemble QC JSON is written to `assemble_qc.json` in the output directory.
