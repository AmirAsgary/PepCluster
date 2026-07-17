# PepCluster — documentation

BLOSUM62-aware **anchor-residue clustering** for immunopeptides, with a fast
Rust backend. Groups peptides by the similarity of their MHC-I anchor residues
(by default the first 3 + last 3 amino acids, weighting P2 and PΩ 2×).

- Repo: <https://github.com/AmirAsgary/PepCluster>
- PyPI: `pip install pepcluster` (prebuilt wheels — no Rust toolchain needed)

---

## 1. Install

```bash
pip install pepcluster            # or: pip install -U pepcluster
pepcluster --version
```

`pepcluster.HAS_RUST` is `True` when the compiled backend is active (the
default from a wheel). A pure-Python fallback runs otherwise, with **identical**
results.

---

## 2. Quick start

```bash
pepcluster -i peptides.fasta -o out -t 0.6
```

Writes cluster assignments and per-cluster FASTA files under `out/`.

Python:

```python
import pepcluster
stats = pepcluster.cluster_fasta("peptides.fasta", "out", threshold=0.6)
print(stats["n_clusters"], "clusters")
```

---

## 3. Command-line options

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | *required* | Input FASTA |
| `-o, --outdir` | `anchor_clusters` | Output directory |
| `-t, --threshold` | `0.6` | Similarity threshold in `[0, 1]` (higher = stricter) |
| `--min-cluster-size` | `2` | Min members for a per-cluster FASTA |
| `--n-front` | `3` | N-terminal anchor length |
| `--n-back` | `3` | C-terminal anchor length |
| `--anchors` | `"2;3"` | Which positions are binding anchors (see §5) |
| `--anchor-weight` | `2.0` | Weight of anchor positions (all others `1.0`) |
| `--refinement` | off | Lloyd-style refinement after greedy clustering (see §6) |
| `--iterations` | `3` | Max refinement passes |
| `--refine-cap` | `32` | Cap on reassignment comparisons per anchor (`<=0` = no cap) |
| `--no-merge` | off | Skip the refinement merge step |
| `--fast-medoid` | off | O(N) medoid instead of exact O(k²) |
| `--merge-cap` | `0` | Cap on merge comparisons per cluster (`0` = no cap) |
| `--backend` | `auto` | `auto` \| `rust` \| `python` |
| `-q, --quiet` | — | Suppress progress output |

**Threshold guide:** `0.8` strict · `0.6` moderate (recommended) · `0.4` relaxed.

---

## 4. Output files

```
out/
├── clusters.tsv          # cluster_id, representative_anchor, representative_peptide, header, sequence, anchor  (one row per peptide)
├── cluster_summary.tsv   # cluster_id, representative_anchor, representative_peptide, size  (sorted by size)
├── summary.txt           # run statistics
└── fasta/
    ├── cluster_0.fasta   # per-cluster FASTA, ready for MSA (>= --min-cluster-size)
    └── SHORT_peptides.fasta  # peptides too short to form an anchor (if any)
```

- **`representative_peptide`** — the cluster's central member (least average
  distance to the rest); a good label/seed. Always a real member sequence.
- **`representative_anchor`** — the anchor that originally seeded the cluster.

---

## 5. Choosing the anchors (`--anchors`)

The anchor is the first `--n-front` + last `--n-back` residues. `--anchors`
selects which positions *inside the anchor* are binding anchors: they get
`--anchor-weight` and define the coarse-alphabet blocking.

Format **`"FRONT;BACK"`** — each side a comma-separated list of **1-based**
indices into that side's residues; either side may be empty.

| `--anchors` | Positions (in the 6-mer) | Meaning |
|-------------|--------------------------|---------|
| `"2;3"` *(default)* | 2, 6 | P2 + PΩ (MHC-I) |
| `"2;2,3"` | 2, 5, 6 | P2 + last two C-terminal residues |
| `"1,2;3"` | 1, 2, 6 | first two N-terminal + PΩ |
| `";3"` | 6 | C-terminal anchor only |
| `"2;"` | 2 | N-terminal anchor only |

```bash
# P2 + last two residues, anchors weighted 3x
pepcluster -i peptides.fasta -o out --anchors "2;2,3" --anchor-weight 3

# a 2+2 anchor, anchor at the 2nd of each side
pepcluster -i peptides.fasta -o out --n-front 2 --n-back 2 --anchors "2;2"
```

Up to 8 anchor positions are supported.

---

## 6. Refinement and making it fast

`--refinement` runs a Lloyd-style pass (medoid update → reassignment → merge).
Two of its sub-steps are quadratic by default, which is fine at moderate
thresholds but blows up at the extremes on large data:

| Sub-step | Default cost | Blows up when… | Fast flag |
|----------|--------------|----------------|-----------|
| medoid update | O(k²) / cluster | a few clusters are **huge** (low threshold) | `--fast-medoid` → O(N) |
| reassignment | O(cap) / anchor | *(already capped)* | `--refine-cap` |
| merge | O(clusters × neighbours) | **many** clusters (high threshold) | `--merge-cap N` (or `--no-merge`) |

The fast paths are **opt-in**; defaults reproduce earlier versions exactly. They
are near-lossless (they change results only at cluster boundaries).

**Large-dataset recipe** — bound every sub-step:

```bash
pepcluster -i peptides.fasta -o out -t 0.6 --refinement \
    --fast-medoid --refine-cap 32 --merge-cap 32
```

Rule of thumb: add `--fast-medoid` at **low** thresholds (giant clusters) and
`--merge-cap 32` at **high** thresholds (many clusters).

---

## 7. Python API

```python
import pepcluster

# End-to-end (FASTA -> TSV + per-cluster FASTA); returns a stats dict
pepcluster.cluster_fasta("peptides.fasta", "out", threshold=0.6,
                         anchors="2;3", refinement=True,
                         fast_medoid=True, merge_cap=32)

# Low-level: cluster a dict of unique anchors -> frequency
mapping, n_cmp, n_early = pepcluster.cluster_anchors(
    {"YLLAGV": 3, "YMLAGV": 1}, 0.6, anchor_positions=[1, 5], anchor_weight=2.0)

# Refine an existing clustering
refined, stats = pepcluster.refine_clusters(
    counts, mapping, 0.6, iterations=3, cap=32, merge=True,
    anchor_positions=[1, 5], anchor_weight=2.0, fast_medoid=False, merge_cap=0)

# Central (medoid) anchor per cluster
reps = pepcluster.cluster_representatives(counts, mapping)

# Parse an --anchors spec to 0-based positions
pepcluster.parse_anchors("2;3", 3, 3)   # -> [1, 5]
```

---

## 8. How it works (brief)

1. Reduce each peptide to its anchor (first `n_front` + last `n_back` residues).
2. Deduplicate by exact anchor; cluster *unique* anchors weighted by frequency.
3. Similarity = BLOSUM62 normalized to `B(a,b)/√(B(a,a)·B(b,b))`, summed over
   positions; anchor positions carry `--anchor-weight`.
4. Block anchors by the coarse 10-letter alphabet at the anchor positions, so
   only plausibly-similar anchors are compared (heavy positions checked first,
   with early termination).
5. Greedy centroid clustering within each block; optional refinement.

The Rust and pure-Python backends compute in f64 with identical deterministic
orderings, so their outputs are **bit-identical**.

---

## 9. Notes

- **Single-threaded** — extra cores don't speed it up (allocate `-c 2` on SLURM).
- **Memory** scales linearly with the number of unique anchors.
- **License:** MIT.
