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
| `--obg-block-search` | off | Search neighbouring blocks (fewer, tighter clusters; see §8.5) |
| `--obg-max-probes` | `0` | With OBG, max blocks searched per anchor (`0` = all) |
| `--obg-min-block-upper-bound` | `0.0` | With OBG, raise the eligibility cut to `max(threshold, X)` |
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

## 8. Methodology

This section describes the algorithm in full: the similarity model, the blocking
scheme, the greedy clustering, the representative selection, the optional
refinement, and the engineering that keeps the two backends bit-identical.

### 8.1 Motivation

For MHC-I peptides, binding specificity is dominated by a small number of
**anchor residues** — canonically position 2 (P2) and the C-terminus (PΩ) — while
the middle of the peptide is far more variable. General sequence-clustering
tools treat every position equally and, on short peptides (8–14 aa), their
internal identity estimators degrade. PepCluster instead reduces each peptide to
just its anchor region and scores similarity with **position-specific weights**,
so clusters are organised by the residues that actually drive binding.

The method has three design goals: (i) be *biologically meaningful* (anchor-aware,
substitution-aware via BLOSUM62); (ii) *scale* to tens of millions of peptides;
and (iii) be *deterministic and reproducible* across machines and backends.

### 8.2 Anchor extraction and deduplication

Each peptide of length `≥ n_front + n_back` is reduced to its **anchor**: the
first `n_front` and last `n_back` residues concatenated (default 3 + 3 = a
6-mer). Peptides shorter than the anchor are set aside (`SHORT_peptides.fasta`)
and not clustered.

Many distinct peptides share the same anchor. PepCluster **deduplicates**: it
clusters the set of *unique anchors*, each weighted by its **frequency** `f` =
the number of peptides carrying it. All downstream steps operate on unique
anchors with frequencies, which is what makes 11M peptides tractable (typically
a few million unique anchors, and the heavy steps depend on that smaller count).

### 8.3 The similarity metric

**Per-residue similarity.** Two amino acids `x, y` are compared with a
BLOSUM62 score normalised to a cosine-like `[-…, 1]` scale:

```
normsim(x, y) = B(x, y) / sqrt( B(x, x) · B(y, y) )
```

where `B` is the standard BLOSUM62 matrix. By construction `normsim(x, x) = 1`
(the maximum), and conservative substitutions score high while dissimilar pairs
score low or negative. This value is precomputed for all 20×20 residue pairs.

**Position weighting.** Each of the `L = n_front + n_back` anchor positions has a
raw weight: `anchor_weight` (default 2.0) if it is a chosen anchor position
(`--anchors`), else 1.0. Raw weights are **normalised to sum to 1**:

```
w_p = raw_p / Σ_q raw_q ,      Σ_p w_p = 1
```

**Anchor similarity.** The similarity of two anchors `a, b` (both length `L`) is
the weighted sum over positions:

```
S(a, b) = Σ_{p=0}^{L-1}  w_p · normsim(a[p], b[p])
```

Key properties used throughout:

- **Additivity over positions** — `S` is a plain sum of per-position terms. This
  is the single fact that enables the O(N) medoid (§8.6) and the early-exit prune
  (§8.5).
- **Range and maximum** — because `normsim ≤ 1` with equality only at identity,
  and `Σ w_p = 1`, the maximum possible score is `S(a, a) = 1`.
- **Symmetry** — `S(a, b) = S(b, a)`.

The `threshold` (`-t`) is compared directly against `S`: two anchors join the
same cluster when `S ≥ threshold`. Higher threshold ⇒ stricter ⇒ more, smaller
clusters.

### 8.4 Reduced-alphabet blocking

Comparing every anchor against every other is O(U²) and infeasible. PepCluster
partitions anchors into **blocks** so only plausibly-similar anchors are ever
compared.

**Coarse alphabet.** The 20 amino acids are mapped to 10 physicochemical groups:

```
AST · VILM · FYW · DE · KR · NQ · G · H · C · P
```

**Block key.** An anchor's block key is the tuple of coarse groups **at its
anchor positions only** (the heavily-weighted, binding-relevant positions). With
the default two anchor positions this gives up to `10 × 10 = 100` blocks; with
`k` anchor positions, up to `10^k`. Internally the key is packed one byte per
anchor position into a 64-bit integer (so ≤ 8 anchor positions are supported).

**Rationale and trade-off.** Two anchors can only be similar if their *anchor*
residues are at least in the same physicochemical class, so anchors in the same
block are the only credible cluster-mates for the greedy pass. This is a
heuristic: it can *over-split* — two genuinely similar anchors whose anchor
residues fall in different coarse groups land in different blocks and are not
compared during the greedy pass. Refinement (§8.7) mitigates this by also
looking at **neighbouring blocks**.

**Neighbouring blocks.** Two blocks are *neighbours* if they share the same
coarse group in **at least one** anchor position. Refinement's reassignment and
merge steps search a centroid's own block plus its neighbours, allowing anchors
to cross the block boundaries the greedy pass could not.

### 8.5 Greedy centroid clustering

The core (`--refinement` off) is a single-pass, order-dependent greedy pass:

1. Sort unique anchors by **frequency descending** (stable; ties keep input
   order). Frequent anchors become centroids first, which tends to place the
   most representative sequences at cluster centres.
2. Partition anchors into blocks by their block key.
3. Within each block, process anchors in the sorted order. For each anchor:
   - compare it to the block's existing centroids;
   - assign it to the **first** centroid with `S ≥ threshold`;
   - if none qualifies, it becomes a **new centroid**.

The result is a mapping `anchor → centroid anchor`. Only within-block
comparisons happen, which is what makes this pass fast.

**Early termination (admissible prune).** Positions are examined
**heaviest-weight first** (`check_pos`). Let `remaining_after[k]` be the sum of
weights of the positions not yet examined after step `k`. Because each remaining
position contributes at most `w_p · 1`, the quantity

```
upper_bound = partial_score + remaining_after[k+1]
```

is a true upper bound on the final `S`. If `upper_bound < threshold`, the pair
can never reach the threshold and the comparison is abandoned early (returns
−1). This never discards a real match (it is *admissible*), and on strict
thresholds it prunes the vast majority of comparisons after only one or two
positions — e.g. at `t = 0.8` on 11M peptides, ~99.8% of comparisons early-exit.

**Upper-bound-guided multi-probe search (`--obg-block-search`).** Blocking is a
heuristic and can over-split (§8.4). With OBG on, an anchor also searches
*neighbouring* blocks, restricted by an admissible **block-level upper bound**.
For two blocks `K1, K2`, the best achievable similarity between any anchor in
one and any anchor in the other is

```
UB(K1, K2) = Σ_{non-anchor p} w_p · 1  +  Σ_{anchor p} w_p · maxgroupsim(g1_p, g2_p)
```

where `maxgroupsim(g1, g2)` is the maximum normalized BLOSUM62 similarity between
any residue of coarse group `g1` and any of `g2` (1.0 on the diagonal), and
non-anchor positions are unconstrained by the block key so they contribute their
full weight. Because `UB` is an upper bound, any block with `UB < threshold`
provably contains no match and is skipped. Eligible blocks are ranked by `UB`
(own block first, highest bound next) and searched greedily; `--obg-max-probes`
caps how many are searched and `--obg-min-block-upper-bound` raises the
eligibility cut to `max(threshold, X)`. The result is fewer, tighter clusters
(the greedy recovers cross-block matches it would otherwise miss) at bounded
extra cost. Off by default.

### 8.6 Cluster representatives (medoid)

Every cluster reports a **representative peptide** — the central member with the
least average distance to the rest. Because distance is anchor-based, all
peptides sharing an anchor are equidistant, so the representative is the peptide
carrying the cluster's **medoid anchor**: the member `a` maximising its total
frequency-weighted similarity to the whole cluster,

```
Sigma(a) = Σ_j f_j · S(a, a_j)          (sum over all members j)
representative = argmax_a Sigma(a)
```

**O(N) computation.** Evaluating `Sigma` for every member pairwise is O(k²).
Exploiting additivity (§8.3), regroup the double sum by position. For each
position `p` build a per-residue aggregate over the whole cluster:

```
S_p[x] = Σ_j f_j · w_p · normsim(x, a_j[p])
```

then

```
Sigma(a) = Σ_p S_p[a[p]]
```

Building the `S_p` tables touches each member **once per position**, and scoring
each member is one lookup per position — no pairwise loop. This is **O(k · L)**
per cluster instead of O(k²), so representatives are computed in linear time even
for clusters of millions of anchors. Ties are broken deterministically (higher
`Sigma`, then higher frequency, then lexicographically smaller anchor).

### 8.7 Optional refinement (Lloyd-style)

`--refinement` runs up to `--iterations` passes; each pass has three sub-steps
and the pass loop stops as soon as a pass changes nothing.

**(1) Medoid update.** Each cluster's centroid is replaced by its member with
the highest frequency-weighted average similarity to the *other* members:

```
avg(i) = [ Σ_{j≠i} f_j · S(i, j) ] / [ Σ_{j≠i} f_j ]        centroid ← argmax_i avg(i)
```

- *Exact* (default): the pairwise form above, **O(k²)** per cluster.
- *Fast* (`--fast-medoid`): the same per-position decomposition as §8.6. With
  `Sigma(i) = Σ_p S_p[a_i[p]]`, `S(i,i) = 1`, and `F = Σ_j f_j`,

  ```
  avg(i) = ( Sigma(i) − f_i ) / ( F − f_i )
  ```

  giving **O(k · L)**. This is what rescues low thresholds, where a handful of
  clusters can each hold millions of anchors (k² is otherwise ~10¹²⁺ per pass).

**(2) Cross-block reassignment.** Each anchor may move to a better centroid. Its
candidate centroids are those in its own block plus neighbouring blocks (§8.4),
examined **own-block-first, largest-cluster-first**, and the anchor moves to the
highest-scoring candidate `≥ threshold` that beats its current centroid. The scan
is bounded by **`--refine-cap`** comparisons per anchor (default 32): since the
best match is almost always in the same block and near the front of the ordered
list, the cap is near-lossless while making this step **O(U · cap)**.

**(3) Centroid merge.** Whole clusters that drifted apart in the greedy pass are
recombined: centroids are processed largest-first, and a *smaller* neighbouring
centroid is absorbed into a larger one when `S ≥ threshold`. By default this
scans all neighbouring centroids (**O(C × neighbours)**); **`--merge-cap N`**
bounds it to the `N` most-promising candidates per centroid (own-block-first,
largest-first), making it **O(C · N)**. `--no-merge` skips the step entirely.

**Convergence.** Each sub-step only accepts a change when it strictly improves
the local objective (higher medoid score, strictly-better reassignment target,
merge only of a smaller into a strictly-not-smaller cluster), and the pass loop
terminates when a full pass produces no change or after `--iterations` passes.

**Which flags for which regime.** Refinement's two quadratic sub-steps blow up at
opposite ends of the threshold range:

| Regime | Cluster shape | Bottleneck sub-step | Fix |
|--------|---------------|---------------------|-----|
| Low `-t` | a few **giant** clusters | medoid O(k²) | `--fast-medoid` |
| Moderate `-t` | balanced | — | (defaults are fine) |
| High `-t` | **many** small clusters | merge O(C×neighbours) | `--merge-cap` / `--no-merge` |

### 8.8 Determinism and backend parity

The Rust extension and the pure-Python fallback are engineered to produce
**bit-identical** results, which matters for reproducible science:

- **f64 everywhere.** Both backends accumulate similarities in IEEE double
  precision (the Rust tables are `f64`, not `f32`), so the arithmetic matches.
- **Fixed summation order.** Positions are summed in a fixed `check_pos` order;
  per-cluster aggregations iterate members in a canonical **anchor-sorted** order
  and residues in ascending byte order. Since floating-point addition is not
  associative, pinning the order is required for exact agreement.
- **Deterministic tie-breaking.** Frequency sorts break ties by anchor string;
  argmaxes keep the first candidate under a strict `>` over the sorted order.
  Results therefore do not depend on hash-map iteration order and are stable
  run-to-run.

This parity is asserted by the test suite across thresholds, anchor lengths
(4/6/8), anchor positions, weights, and every refinement flag combination.

### 8.9 Complexity summary

Let `U` = unique anchors, `L` = anchor length, `C` = number of clusters,
`k` = anchors in a cluster.

| Stage | Cost |
|-------|------|
| Anchor extraction + dedup | O(total peptides · L) |
| Greedy clustering | O(within-block comparisons · L), heavily pruned by early exit |
| Representatives | O(U · L) |
| Refinement — medoid (exact / fast) | O(Σ k²) · L  /  **O(U · L)** |
| Refinement — reassignment | O(passes · U · cap · L) |
| Refinement — merge (uncapped / capped) | O(passes · C × neighbours)  /  **O(passes · C · merge_cap)** |

Memory is O(U · L) for the anchors plus O(U) bookkeeping.

### 8.10 Parameter guidance

- **`-t / --threshold`** is the main knob: `0.8` strict, `0.6` moderate
  (recommended), `0.4` relaxed. It trades cluster count against cluster size.
- **`--anchors` / `--anchor-weight`** encode the biology: which positions drive
  the grouping, and how strongly. Default `2;3` @ 2× = P2 + PΩ for MHC-I.
- **`--refinement`** improves cluster assignments at extra cost; on large data
  pair it with the fast flags from §8.7 so every sub-step is bounded.

---

## 9. Notes

- **Single-threaded** — extra cores don't speed it up (allocate `-c 2` on SLURM).
- **Memory** scales linearly with the number of unique anchors.
- **License:** MIT.
