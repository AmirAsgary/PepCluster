# PepCluster

**BLOSUM62-aware anchor-residue clustering for immunopeptides — with a fast Rust backend.**

PepCluster groups peptides by the similarity of their MHC-I **anchor residues**
(the first 3 + last 3 amino acids) using a BLOSUM62-normalized similarity metric
with double weight on the primary anchor positions P2 and PΩ. It is built for
large immunopeptidomics datasets: a Rust extension does the heavy lifting
(10–100× faster than pure Python), and a pure-Python fallback keeps it working
everywhere.

Unlike general-purpose sequence tools (e.g. MMseqs2), PepCluster distinguishes
anchor from non-anchor positions, which is what actually drives MHC-I binding
specificity — producing biologically interpretable clusters on short (8–14 aa)
peptides where full-length estimators break down.

---

## Install

```bash
pip install pepcluster
```

Prebuilt wheels are published for Linux, macOS, and Windows, so **no Rust
toolchain is required** for end users. If you install on a platform without a
wheel, pip builds from source (needs a Rust compiler — see
[Building from source](#building-from-source)).

---

## Quick start

### Command line

```bash
pepcluster -i examples/peptides.fasta -o out -t 0.6
```

This writes cluster assignments and per-cluster FASTA files under `out/`
(see [Output files](#output-files)).

### Python

```python
import pepcluster

# End-to-end: FASTA in → TSV + per-cluster FASTA out
stats = pepcluster.cluster_fasta("peptides.fasta", "out", threshold=0.6)
print(stats["n_clusters"], "clusters")

# Low-level: cluster a dict of unique 6-mer anchors → frequency
mapping, n_cmp, n_early = pepcluster.cluster_anchors(
    {"YLLAGV": 3, "YMLAGV": 1, "GYAWTK": 2}, 0.6)
# mapping: {anchor -> representative anchor}

# Optional Lloyd-style refinement on top of the greedy result
refined, refine_stats = pepcluster.refine_clusters(
    {"YLLAGV": 3, "YMLAGV": 1}, mapping, 0.6, iterations=3)
```

`pepcluster.HAS_RUST` tells you whether the compiled backend is active
(`cluster_anchors` / `refine_clusters` automatically use Rust when available and
fall back to identical pure-Python implementations otherwise).

---

## CLI options

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | *required* | Input FASTA file |
| `-o, --outdir` | `anchor_clusters` | Output directory |
| `-t, --threshold` | `0.6` | BLOSUM similarity threshold (0.0–1.0) |
| `--min-cluster-size` | `2` | Min members for a per-cluster FASTA |
| `--n-front` | `3` | N-terminal anchor length |
| `--n-back` | `3` | C-terminal anchor length |
| `--anchors` | `"2;3"` | Which positions are binding anchors, as `FRONT;BACK` (1-based per side). See [Choosing the anchors](#choosing-the-anchors) |
| `--anchor-weight` | `2.0` | Weight given to anchor positions (all others are `1.0`) |
| `--refinement` | off | Apply Lloyd-style refinement after greedy clustering |
| `--iterations` | `3` | Max refinement passes (with `--refinement`) |
| `--refine-cap` | `32` | Max centroid comparisons per anchor in refinement reassignment (`<=0` = no cap). Lower = faster |
| `--no-merge` | off | Skip the refinement centroid-merge step (much faster on many-cluster data) |
| `--fast-medoid` | off | O(N) medoid instead of exact O(k²) (much faster when a few clusters are huge) |
| `--merge-cap` | `0` | Cap candidate centroids per cluster in the merge step (`0` = no cap; e.g. `32` on many-cluster data) |
| `--backend` | `auto` | `auto` \| `rust` \| `python` |
| `-q, --quiet` | — | Suppress progress output |

### Making refinement fast

Refinement (`--refinement`) has three sub-steps, and **two of them are quadratic
by default** — which is fine at moderate thresholds but blows up at the extremes:

| Sub-step | Cost (default) | Blows up when… | Fast flag |
|----------|----------------|----------------|-----------|
| medoid update | O(k²) per cluster | a few clusters are **huge** (low threshold) | `--fast-medoid` → O(N) |
| reassignment | O(cap) per anchor | *(already capped)* | `--refine-cap` |
| merge | O(clusters × neighbours) | there are **many** clusters (high threshold) | `--merge-cap N` (or `--no-merge`) |

- **`--refine-cap N`** bounds candidate centroids per anchor in reassignment
  (own-block-first, largest-cluster-first). Default `32` is near-lossless.
- **`--fast-medoid`** replaces the exact all-pairs medoid with an O(N)
  per-position decomposition. On a single cluster of 30k anchors it is ~470×
  faster (and the gap grows quadratically with cluster size).
- **`--merge-cap N`** bounds candidate centroids per cluster in the merge step
  (~130× faster with ~1M clusters). **`--no-merge`** skips merging entirely.

For a huge dataset, turn all three on so every sub-step is bounded:

```bash
pepcluster -i peptides.fasta -o out -t 0.6 --refinement \
    --fast-medoid --refine-cap 32 --merge-cap 32
```

Defaults are unchanged (exact medoid, uncapped merge) for backward
compatibility — the fast paths are opt-in.

**Threshold guide:**

| Value | Effect |
|-------|--------|
| 0.8 | Strict — mostly exact matches + very conservative substitutions |
| 0.6 | Moderate — allows 1–2 conservative substitutions (recommended) |
| 0.4 | Relaxed — broader groups for exploratory analysis |

---

## Choosing the anchors

The anchor is the peptide's first `--n-front` + last `--n-back` residues.
`--anchors` picks which positions *inside that anchor* are the **binding
anchors**: they get `--anchor-weight` (default 2×) in the similarity score and
they define the coarse-alphabet blocking.

The format is **`"FRONT;BACK"`**, where each side is a comma-separated list of
**1-based** indices into that side's residues. Either side may be empty.

The default **`"2;3"`** means *the 2nd of the first 3 residues* (**P2**) and
*the 3rd of the last 3* (**PΩ**) — the classic MHC-I anchors.

| `--anchors` | Anchor positions (1-based, in the 6-mer) | Meaning |
|-------------|------------------------------------------|---------|
| `"2;3"` *(default)* | 2, 6 | P2 + PΩ — MHC-I |
| `"2;2,3"` | 2, 5, 6 | P2 plus the last two C-terminal residues |
| `"1,2;3"` | 1, 2, 6 | first two N-terminal residues plus PΩ |
| `";3"` | 6 | C-terminal anchor only |
| `"2;"` | 2 | N-terminal anchor only |

```bash
# emphasise P2 and the last two residues, and weight anchors 3x
pepcluster -i peptides.fasta -o out --anchors "2;2,3" --anchor-weight 3
```

Anchor positions are indexed **relative to `--n-front` / `--n-back`**, so the
two work together — e.g. a 2+2 anchor with anchors at the 2nd of each side:

```bash
pepcluster -i peptides.fasta -o out --n-front 2 --n-back 2 --anchors "2;2"
```

Up to 8 anchor positions are supported (they form the blocking key).

---

## Output files

```
out/
├── clusters.tsv            # cluster_id, representative_anchor, representative_peptide, header, sequence, anchor (every peptide)
├── cluster_summary.tsv     # cluster_id, representative_anchor, representative_peptide, size (sorted by size)
├── summary.txt             # run statistics
└── fasta/
    ├── cluster_0.fasta     # per-cluster FASTA, ready for MSA (>= --min-cluster-size members)
    ├── cluster_1.fasta
    └── SHORT_peptides.fasta # peptides too short to form an anchor (if any)
```

**`representative_peptide`** is the cluster's *central* member — the peptide
with the least average distance (highest weighted average similarity) to every
other peptide in the cluster. It is computed in linear time and is always a
real member sequence, so it's a good label or seed for each cluster.
`representative_anchor` is the anchor that originally seeded the cluster.

---

## How it works

1. **Anchor extraction.** Each peptide is reduced to its anchor: the first
   `--n-front` (3) and last `--n-back` (3) amino acids. Peptides shorter than
   that are set aside in `SHORT_peptides.fasta`.
2. **Deduplicate.** Peptides are grouped by their exact anchor, so clustering
   operates on *unique* anchors weighted by frequency.
3. **Similarity metric.** Two anchors are compared position-by-position with a
   BLOSUM62 score normalized to `sim(a,b) = B(a,b) / sqrt(B(a,a)·B(b,b))`.
   The **anchor positions** (`--anchors`, default P2 and PΩ) carry
   `--anchor-weight` (default 2×); every other position carries 1×. The score
   is a weighted mean in `[−…, 1]`.
4. **Blocking.** Unique anchors are bucketed by a reduced 10-letter alphabet at
   the anchor positions (10ᵏ bins for k anchors — 100 by default), so only
   plausibly-similar anchors are ever compared. High-weight positions are
   checked first with early termination.
5. **Greedy clustering.** Within each block, anchors are processed
   most-frequent-first; each joins the first centroid above `threshold` or
   becomes a new centroid.
6. **Optional refinement** (`--refinement`). A Lloyd-style pass iterates:
   medoid update → cross-block reassignment → centroid merging, until stable.

The Rust backend (`pepcluster._core`) and the pure-Python reference
(`pepcluster.clustering`) implement identical logic — both compute in f64 with
identical, deterministic orderings — and produce **bit-identical** cluster
assignments; the test suite asserts this parity across anchor lengths,
positions and weights.

---

## Performance

| Dataset | Python | Rust |
|---------|--------|------|
| 7K peptides | <1 s | <1 s |
| 2.5M peptides | ~3 min | ~15 s |

Speed comes from anchor deduplication, coarse-alphabet blocking, and weighted
early-termination in the similarity check.

---

## Building from source

Requires a [Rust toolchain](https://rustup.rs) and
[maturin](https://www.maturin.rs).

```bash
# one-time
pip install maturin

# build + install into the current environment (editable-ish)
maturin develop --release

# or build a wheel
maturin build --release      # wheel lands in target/wheels/
```

The project uses maturin's mixed layout: Rust lives in `src/lib.rs`
(compiled to `pepcluster._core`), Python in `python/pepcluster/`.

Run the tests with:

```bash
pip install pytest
pytest
```

---

## Releasing (maintainers)

Wheels are built for Linux / macOS / Windows by
[`.github/workflows/CI.yml`](.github/workflows/CI.yml) and published to PyPI on
version tags via [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/)
(OpenID Connect — no API token or stored secret).

1. One-time: on <https://pypi.org/manage/account/publishing/> add a pending
   publisher — Owner `AmirAsgary`, Repository `PepCluster`, Workflow `CI.yml`
   (leave the environment blank).
2. Bump the version in `pyproject.toml` **and** `Cargo.toml`.
3. Tag and push:

   ```bash
   git tag v0.1.0
   git push origin v0.1.0
   ```

The `release` job then builds all wheels + an sdist and uploads them to PyPI.

---

## License

[MIT](LICENSE) © 2026 Amir Asgary

---

## Citation

If you use PepCluster in your research, please cite this repository:

> Asgary, A. *PepCluster: BLOSUM62-aware anchor-residue clustering for
> immunopeptides.* https://github.com/AmirAsgary/PepCluster
