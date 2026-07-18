# Improvement 1: Upper-Bound-Guided Multi-Probe Block Search

## Where it goes

Replace:

```text
Create blocks
→ Compare only anchors in the same block
→ Greedy clustering
```

with:

```text
Create blocks
→ Calculate upper bounds between block pairs
→ Rank eligible neighboring blocks
→ Search selected blocks
→ Greedy clustering
```

This improvement combines two related ideas:

1. **Upper-bound-guided search** determines which block pairs could theoretically contain anchors above the similarity threshold.
2. **Multi-probe search** optionally limits the search to the most promising eligible blocks.

For every pair of blocks, calculate the maximum similarity that any pair of anchors from those blocks could achieve. Blocks whose upper bound is below the clustering threshold are safely ignored.

Eligible blocks are ranked by upper bound and searched from highest to lowest.

Suggested flags:

```text
--obg-block-search
--obg-max-probes <integer>
--obg-min-block-upper-bound <float>
```

Behavior:

```text
--obg-block-search disabled
    Use the current same-block-only search.

--obg-block-search enabled with no additional limits
    Search every block whose upper bound is at least the clustering threshold.

--obg-max-probes K
    Search at most the top K eligible blocks.

--obg-min-block-upper-bound X
    Search only eligible blocks with upper bound ≥ X.
```

Both limits may be active together. In that case, searching stops when either limit is reached.

---

# Improvement 2: Central Region k-mer Refinement

After the initial anchor-based clustering, PepCluster can optionally use central-residue similarity during the refinement and peptide-reassignment stage. It should not affect the initial blocking or greedy anchor clustering, preserving most of the speed advantage.

Because peptides with identical anchors may have different middle regions, the initial clustering may still operate on unique anchors, but central-region refinement must operate at the peptide level.

When both peptides contain central-region k-mers, their final similarity is:

$$
S_{\mathrm{final}}
==================

(1-w_{\mathrm{central}})S_{\mathrm{anchor}}
+
w_{\mathrm{central}}S_{\mathrm{central}}.
$$

The default central-region weight is:

$$
w_{\mathrm{central}} = 0.2,
$$

giving anchor similarity 80% of the final score. This weight is user-adjustable.

If either peptide is too short to produce central-region k-mers, the central term is ignored:

$$
S_{\mathrm{final}} = S_{\mathrm{anchor}}.
$$

By default, middle regions are represented using position-binned multiset k-mers. This representation preserves both repeated k-mers and their approximate relative positions without requiring sequence alignment.

The default k-mer length is two residues, but it is adjustable.

## Central similarity calculations

### Multiset k-mers only

Use weighted Jaccard over k-mer counts:

$$
J
=

\frac{
\sum_k \min\left(c_A(k),c_B(k)\right)
}{
\sum_k \max\left(c_A(k),c_B(k)\right)
}.
$$

### Position-binned k-mers only

Treat every $(\mathrm{bin},\mathrm{k\text{-}mer})$ pair as a binary feature and use ordinary set Jaccard.

### Position-binned multiset k-mers

Use weighted Jaccard over the counts of each $(\mathrm{bin},\mathrm{k\text{-}mer})$ feature.

## Relative position binning

Bins should represent relative positions within the middle region.

Example with three bins and `k=2`:

```text
9-mer:
middle length = 3
two k-mers
first k-mer  → bin 0
second k-mer → bin 2

15-mer:
longer middle region
k-mers distributed across bins 0, 1 and 2
```

Thus, early middle-region k-mers from peptides of different lengths can still match each other, as can late middle-region k-mers.

## Adjacent-bin matching

To avoid overly strict hard-bin matching, adjacent bins may match with reduced weight:

```text
same bin       weight 1.0
adjacent bin   weight 0.5
distant bin    weight 0.0
```

Each k-mer occurrence should be matched at most once. Same-bin matches should be assigned first, followed by adjacent-bin matches.

Suggested flags:

```text
--central-region-refinement
--crr-weight 0.2
--crr-kmer-size 2
--crr-bins 3
--crr-multiset
--crr-position-binning
--crr-adjacent-bin-weight 0.5
```

Recommended defaults:

```text
--central-region-refinement     disabled
--crr-weight                    0.2
--crr-kmer-size                 2
--crr-bins                      3
--crr-multiset                  enabled
--crr-position-binning          enabled
--crr-adjacent-bin-weight       0.5
```

---

# Improvement 3: Cluster-Profile-Aware Merging

During the refinement stage, cluster merging can consider both centroid-anchor similarity and central-region similarity across all peptides in both clusters.

Anchor similarity is calculated between the two cluster centroids. Central similarity is calculated using aggregated position-binned k-mer profiles built from all peptides in each cluster. K-mer counts are normalized so that larger clusters do not dominate the comparison.

For cluster $C$, the normalized frequency of feature $k$ is:

$$
f_C(k)
======

\frac{
\operatorname{count}_C(k)
}{
\sum_j \operatorname{count}_C(j)
}.
$$

The two cluster profiles are compared using weighted Jaccard:

$$
S_{\mathrm{cluster\text{-}kmer}}
================================

\frac{
\sum_k
\min\left(f_{C_1}(k),f_{C_2}(k)\right)
}{
\sum_k
\max\left(f_{C_1}(k),f_{C_2}(k)\right)
}.
$$

The final merge score is:

$$
S_{\mathrm{merge}}
==================

(1-w_{\mathrm{merge}})
S_{\mathrm{centroid\text{-}anchor}}
+
w_{\mathrm{merge}}
S_{\mathrm{cluster\text{-}kmer}}.
$$

The default cluster-profile weight is:

$$
w_{\mathrm{merge}} = 0.2.
$$

Clusters are merged only when:

$$
S_{\mathrm{merge}}
\geq
t_{\mathrm{merge}}.
$$

Because one aggregated profile is compared per cluster, this remains much faster than comparing every peptide in one cluster against every peptide in the other.

Cluster-profile merging should reuse the central-region k-mer settings:

```text
--crr-kmer-size
--crr-bins
--crr-multiset
--crr-position-binning
--crr-adjacent-bin-weight
```

Suggested flags:

```text
--cluster-profile-merge-weight 0.2
--cluster-profile-merge-threshold 0.6
--no-cluster-profile-merge
```

Recommended behavior:

```text
Refinement disabled
    No cluster merging is performed.

Refinement enabled, but central-region refinement disabled
    Merge clusters using centroid-anchor similarity only.

Refinement and central-region refinement enabled
    Cluster-profile-aware merging is enabled by default.
    The merge score combines centroid-anchor similarity with
    aggregated cluster k-mer-profile similarity.

--no-cluster-profile-merge
    Disable cluster-profile-aware merging and use
    centroid-anchor similarity only.
```

A separate activation flag is unnecessary because cluster-profile-aware merging is automatically used whenever both refinement and central-region refinement are active.

Cluster profiles should be updated incrementally after each merge rather than recomputed from all member peptides.

---

# Improvement 4: SIMD, Bit Packing and Multithreading

PepCluster can be accelerated without changing its clustering logic or similarity definitions.

## Bit packing

Each amino acid can be encoded using five bits. A short anchor can therefore be stored inside one machine word instead of as a string.

This reduces:

* memory usage;
* hashing overhead;
* cache misses;
* residue-access overhead.

Bit packing should be enabled automatically when the anchor length fits within the selected packed representation.

## SIMD

SIMD instructions can evaluate several residue or anchor comparisons simultaneously.

Instead of comparing one anchor pair at a time, the implementation may compare multiple candidate anchors in parallel using CPU vector registers.

SIMD should be enabled automatically when supported by the CPU, with a scalar fallback for unsupported systems.

## Multithreading

Independent work can be distributed across CPU cores, including:

* processing different blocks;
* candidate comparisons;
* cluster-profile construction;
* some refinement operations.

Deterministic output must be preserved through:

* stable centroid ordering;
* deterministic tie-breaking;
* stable merge ordering;
* avoiding unordered parallel updates to shared cluster state.

Suggested flags:

```text
--bit-packing
--no-bit-packing

--simd
--no-simd

--threads <integer>

--deterministic
--no-deterministic
```

Recommended defaults:

```text
bit packing      enabled automatically when supported
SIMD             enabled automatically when supported
threads          1
deterministic    enabled
```

Thread behavior:

```text
--threads 1
    Single-threaded execution.

--threads 0
    Use all available CPU cores.

--threads N
    Use exactly N worker threads.
```

Bit packing and SIMD should not change clustering results. Multithreading should also produce identical results when deterministic mode is enabled.
