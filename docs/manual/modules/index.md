# kam Pipeline: Stage 2 — Index

## What this stage does

The index stage converts assembled molecules into a k-mer evidence index. For every consensus
read in every molecule, it extracts all k-mers of length `k`, canonicalises each k-mer, and
accumulates molecule-level evidence under the canonical key.

The output is a `HashKmerIndex` mapping each observed canonical k-mer to a `MoleculeEvidence`
record describing how many molecules carried it, in what strand configuration, and with what base
quality.

This index is the bridge between molecule assembly and graph-based variant finding. It preserves
the molecule identity that would be lost in a naive k-mer count.

---

## K-mer encoding

### 2-bit encoding

Each base is 2-bit encoded:

| Base | Bits |
|------|------|
| A | 00 |
| C | 01 |
| G | 10 |
| T | 11 |

A k-mer of length k is packed into a `u64` by shifting bits left as each base is added. For
k ≤ 31, the result fits in a 62-bit value (31 × 2 = 62 bits). k = 32 would require 64 bits
and would be ambiguous with the canonical form computation; the practical limit is k ≤ 31.

N bases and any unrecognised characters cause the k-mer to be skipped entirely. The
`KmerIterator` resets its window when it encounters an N, so a k-mer containing an N is never
emitted.

### Canonical form

The canonical form of a k-mer is `min(kmer, reverse_complement(kmer))` under the numerical
encoding. This makes the key strand-agnostic: the same genomic position yields the same canonical
k-mer whether it appears in a forward-strand read or a reverse-strand read.

The reverse complement of the 2-bit encoding is computed by:
1. Reversing the bit order of the k-mer within its 2k-bit field.
2. Flipping all bits (which maps A↔T and C↔G because A=00, T=11 are complements and C=01,
   G=10 are complements).

This works correctly for all k from 1 to 31.

---

## Target allowlist

Before extracting k-mers from consensus reads, the index stage builds an allowlist of canonical
k-mers from the target FASTA file. Only k-mers present in the allowlist are stored in the index.

This serves two purposes:

1. **Memory control**: at 2M reads with k=31, consensus reads produce millions of distinct
   k-mers. The vast majority are from off-target sequence (intronic DNA, contamination, library
   artefacts). Storing only k-mers that overlap the target panel regions keeps peak RSS under
   2 GB.

2. **Signal-to-noise**: k-mers outside the target windows cannot contribute to variant
   detection. Discarding them does not remove any useful information.

The allowlist is built by sliding a k-mer window across each target sequence in both the forward
and reverse-complement orientations, encoding each k-mer, computing the canonical form, and
inserting it into a `HashSet<u64>`.

### SV junction k-mers

When `sv_junctions` is specified in the config (or `--sv-junctions` at CLI), k-mers from those
sequences are added to the allowlist alongside normal target k-mers. This is required to detect
structural variants whose breakpoint k-mers are not present in the reference target sequences:

- **LargeDeletion**: the deletion junction sequence spans both sides of the deleted region.
- **TandemDuplication**: the duplication junction sequence contains k-mers from both the
  duplicated segment and its flanking reference.
- **Inversion**: inversion junction sequences contain the reverse-complement of the inverted
  segment.
- **InvDel**: inversion-with-deletion junction sequences combine both features.

Without adding SV junction k-mers to the allowlist, molecules supporting these structural
variants would not produce any k-mers in the index and the variants would be invisible.

### Fusion target k-mers

When `fusion_targets` is specified, k-mers from those synthetic sequences are also added to the
allowlist. Each fusion target encodes a junction between two distinct genomic loci. Their k-mers
must be in the allowlist so that fusion-spanning reads are captured during indexing. Fusion
detection is handled separately by the fusion calling step after the main pipeline.

---

## MoleculeEvidence

Each canonical k-mer in the index has one `MoleculeEvidence` record:

```rust
pub struct MoleculeEvidence {
    pub n_molecules: u32,
    pub n_duplex: u32,
    pub n_simplex_fwd: u32,
    pub n_simplex_rev: u32,
    pub min_base_error_prob: f32,
    pub mean_base_error_prob: f32,
}
```

| Field | Meaning |
|-------|---------|
| `n_molecules` | Total molecules (consensus reads) that contributed this k-mer |
| `n_duplex` | Of those, how many came from duplex consensus reads |
| `n_simplex_fwd` | How many from forward-strand simplex or singleton consensus |
| `n_simplex_rev` | How many from reverse-strand simplex consensus |
| `min_base_error_prob` | Minimum mean-across-k-bases error probability across all observations |
| `mean_base_error_prob` | Running mean of the per-k-mer mean error probability |

The counts are molecule counts, not read counts. A duplex molecule contributes one count to
`n_duplex`, not two (one per strand). This is the fundamental difference from a Jellyfish
k-mer count.

### FamilyType mapping

When indexing a consensus read, the family type determines which counter is incremented:

| FamilyType | Counter incremented |
|------------|---------------------|
| `Duplex` | `n_duplex` |
| `SimplexFwd` | `n_simplex_fwd` |
| `SimplexRev` | `n_simplex_rev` |
| `Singleton` | `n_simplex_fwd` (treated as simplex forward) |

---

## K-mer extraction from a consensus read

For each consensus read, the `KmerIterator` slides a window of length k across the sequence. At
each position `pos` from 0 to `len - k`:

1. Encode the k-mer starting at `pos` as a `u64`.
2. Check whether the k-mer is in the allowlist. If not, skip.
3. Compute the canonical form.
4. Look up or create the `MoleculeEvidence` entry for this canonical k-mer.
5. Increment `n_molecules` by 1.
6. Increment the appropriate strand counter.
7. Compute the mean error probability across the k bases of this window.
8. Update `min_base_error_prob` and `mean_base_error_prob` (running mean).

### Which consensus read is used per molecule

Each molecule has up to three consensus reads: `consensus_fwd`, `consensus_rev`, and
`duplex_consensus`. The indexing stage uses:

- `duplex_consensus` when present (highest-quality, error-corrected by both strands).
- `consensus_fwd` when present but no duplex.
- `consensus_rev` when present but no duplex.

In practice, a duplex molecule contributes only its duplex consensus read to the index. This
avoids double-counting: if both the duplex consensus and the individual strand consensuses were
indexed, the same molecule would contribute multiple times to `n_molecules`.

---

## HashKmerIndex

The index is backed by a `HashMap<u64, MoleculeEvidence>`. At approximately 48 bytes per entry
(the `u64` key plus the evidence struct), 882K k-mers observed at 2M reads occupies roughly
42 MB — well within memory budget.

Key operations:

| Method | Complexity | Use |
|--------|------------|-----|
| `get(kmer)` | O(1) average | Evidence lookup |
| `insert(kmer, ev)` | O(1) average | Evidence update |
| `contains(kmer)` | O(1) average | Presence check |
| `molecule_count(kmer)` | O(1) average | Quick count query |
| `entry(kmer)` | O(1) average | Mutable update without double lookup |
| `iter()` | O(n) | Serialisation, export |
| `len()` | O(1) | Size query |

The `with_capacity(n)` constructor pre-allocates for `n` entries to avoid rehashing during bulk
insertion. This is used when the allowlist size is known in advance.

---

## Index QC

After extraction, the index QC JSON (`index_qc.json`) records:

| Field | Meaning |
|-------|---------|
| `n_kmers_observed` | Total distinct canonical k-mers in the index |
| `n_molecules_indexed` | Total consensus reads processed |
| `kmer_size` | The k value used |
| `n_targets` | Number of targets in the FASTA |

---

## Configuration

| Parameter | Default | Config key |
|-----------|---------|-----------|
| `kmer_size` (k) | 31 | `[indexing] kmer_size` |
| `sv_junctions` | (none) | `[input] sv_junctions` |
| `fusion_targets` | (none) | `[input] fusion_targets` |

The allowlist is always built from the target FASTA. There is no option to disable it.

---

## Output file

The index is serialised to a binary file in bincode format with file type tag `FileType::KmerIndex`.
This file is the input to the pathfind stage.
