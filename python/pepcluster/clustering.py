"""
clustering.py — BLOSUM-aware anchor clustering for immunopeptides.

Clusters peptides by the similarity of their anchor residues (first 3 + last 3
aa) using a BLOSUM62-normalized similarity metric with position-specific
weighting (P2 and PΩ weighted 2× for MHC-I relevance).

This module contains:
  * the pure-Python reference implementation of the clustering
    (`cluster_anchors_py`) and refinement (`refine_clusters_py`);
  * FASTA parsing / anchor extraction helpers;
  * `cluster_fasta()` — the full end-to-end pipeline (read FASTA → cluster →
    optionally refine → write outputs), callable as a library or from the CLI.

The Rust backend (`pepcluster._core`) exposes drop-in equivalents of the two
core functions and is preferred automatically when available.
"""

import math
import time
from collections import defaultdict
from pathlib import Path

# ============================================================================
# BLOSUM62 — standard 20×20
# ============================================================================
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"

BLOSUM62_FLAT = [
#    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,
]

N_AA = len(AA_ORDER)
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

# Precompute normalized similarity: sim(a,b) = B(a,b) / sqrt(B(a,a) * B(b,b))
_self = {aa: BLOSUM62_FLAT[i * N_AA + i] for i, aa in enumerate(AA_ORDER)}
SIM = {}
for i, a in enumerate(AA_ORDER):
    for j, b in enumerate(AA_ORDER):
        denom = math.sqrt(_self[a] * _self[b])
        SIM[(a, b)] = BLOSUM62_FLAT[i * N_AA + j] / denom if denom > 0 else 0.0

# ============================================================================
# Reduced alphabet for blocking (10 groups)
# ============================================================================
COARSE = {}
for gid, aas in enumerate(["AST", "VILM", "FYW", "DE", "KR", "NQ",
                            "G", "H", "C", "P"]):
    for aa in aas:
        COARSE[aa] = gid

# ============================================================================
# Position weights: [P1, P2, P3, PΩ-2, PΩ-1, PΩ]
#   P2 and PΩ get 2× weight (primary MHC-I anchors)
# ============================================================================
_RAW_W = [1.0, 2.0, 1.0, 1.0, 1.0, 2.0]
_WSUM = sum(_RAW_W)
WEIGHTS = [w / _WSUM for w in _RAW_W]


# ============================================================================
# Precompute fast lookup arrays
# ============================================================================
# 6 flat arrays of size 128×128, one per position, with weight baked in.
# Indexed by ord(aa1)*128 + ord(aa2) → weighted normalized similarity.
WSIM = []
for _k in range(6):
    _arr = [0.0] * (128 * 128)
    for (_a, _b), _v in SIM.items():
        _arr[ord(_a) * 128 + ord(_b)] = WEIGHTS[_k] * _v
    WSIM.append(_arr)

# Precompute max possible score from remaining positions (for early exit).
_MAX_REMAINING = [0.0] * 7  # index 6 = 0.0 sentinel
for _k in range(5, -1, -1):
    _MAX_REMAINING[_k] = _MAX_REMAINING[_k + 1] + WEIGHTS[_k]

# Check order for positions: check high-weight positions (1, 5) first
# so early termination kicks in sooner.
_CHECK_ORDER = [1, 5, 0, 2, 3, 4]
_REMAINING_AFTER = [0.0] * 7
for _i in range(5, -1, -1):
    _REMAINING_AFTER[_i] = _REMAINING_AFTER[_i + 1] + WEIGHTS[_CHECK_ORDER[_i]]

# Reorder WSIM to match _CHECK_ORDER
_WSIM_ORD = [WSIM[_CHECK_ORDER[i]] for i in range(6)]

# Coarse-group lookup array (indexed by ord)
_COARSE_ORD = [-1] * 128
for _aa, _gid in COARSE.items():
    _COARSE_ORD[ord(_aa)] = _gid


# ============================================================================
# Core logic
# ============================================================================

def extract_anchor(seq, nf=3, nb=3):
    """Return first ``nf`` + last ``nb`` residues, or None if too short."""
    if len(seq) < nf + nb:
        return None
    return seq[:nf] + seq[-nb:]


def _to_ords(anchor):
    """Convert 6-char anchor to tuple of ord values."""
    return (ord(anchor[0]), ord(anchor[1]), ord(anchor[2]),
            ord(anchor[3]), ord(anchor[4]), ord(anchor[5]))


def anchor_sim_fast(a, b, threshold):
    """
    Weighted BLOSUM62-normalized similarity with early termination.
    Checks high-weight positions first; bails if remaining positions
    can't push score above threshold. Returns score or -1.0 on early exit.

    ``a`` and ``b`` are ordinal tuples (see :func:`_to_ords`).
    """
    s = 0.0
    # Unrolled loop over _CHECK_ORDER = [1, 5, 0, 2, 3, 4]
    # Position index 1 (P2, weight 2×)
    s += _WSIM_ORD[0][a[1] * 128 + b[1]]
    if s + _REMAINING_AFTER[1] < threshold:
        return -1.0
    # Position index 5 (PΩ, weight 2×)
    s += _WSIM_ORD[1][a[5] * 128 + b[5]]
    if s + _REMAINING_AFTER[2] < threshold:
        return -1.0
    # Position index 0 (P1)
    s += _WSIM_ORD[2][a[0] * 128 + b[0]]
    if s + _REMAINING_AFTER[3] < threshold:
        return -1.0
    # Position index 2 (P3)
    s += _WSIM_ORD[3][a[2] * 128 + b[2]]
    # Position index 3 (PΩ-2)
    s += _WSIM_ORD[4][a[3] * 128 + b[3]]
    # Position index 4 (PΩ-1)
    s += _WSIM_ORD[5][a[4] * 128 + b[4]]
    return s


def block_key_fast(ords):
    """Reduced-alphabet key at P2 (idx 1) and PΩ (idx 5)."""
    return (_COARSE_ORD[ords[1]], _COARSE_ORD[ords[5]])


def parse_fasta(path):
    """Yield ``(header_id, full_sequence)`` from a FASTA file."""
    hdr = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(parts)
                hdr = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(parts)


def cluster_anchors_py(anchor_counts, threshold):
    """
    Greedy centroid clustering on unique anchors (pure-Python reference
    implementation; the Rust backend ``pepcluster._core.cluster_anchors`` is a
    drop-in equivalent).

    1. Convert anchors to ordinal tuples
    2. Partition into blocks by coarse(P2)+coarse(PΩ)
    3. Within each block, process anchors most-frequent-first
    4. Assign to first centroid above threshold, or become new centroid

    Args:
        anchor_counts: dict[str, int] — unique 6-mer anchor → frequency
        threshold:     float          — minimum similarity to join a cluster

    Returns:
        (mapping, n_comparisons, n_early_exits)
        mapping: dict[str, str] — anchor → centroid anchor
    """
    # Convert to ordinals; sort by frequency descending
    items = sorted(anchor_counts.items(), key=lambda x: -x[1])
    str_to_ords = {}
    blocks = defaultdict(list)
    for anchor_str, cnt in items:
        ords = _to_ords(anchor_str)
        str_to_ords[anchor_str] = ords
        blocks[block_key_fast(ords)].append(anchor_str)

    mapping = {}
    n_cmp = 0
    n_early = 0

    for bk, anchors in blocks.items():
        centroids_ords = []    # list of ord-tuples for centroids
        centroids_str = []     # matching string anchors
        for anchor_str in anchors:
            a_ords = str_to_ords[anchor_str]
            matched = False
            for ci in range(len(centroids_ords)):
                n_cmp += 1
                score = anchor_sim_fast(a_ords, centroids_ords[ci], threshold)
                if score < 0:
                    n_early += 1
                    continue
                if score >= threshold:
                    mapping[anchor_str] = centroids_str[ci]
                    matched = True
                    break
            if not matched:
                centroids_ords.append(a_ords)
                centroids_str.append(anchor_str)
                mapping[anchor_str] = anchor_str

    return mapping, n_cmp, n_early


# ============================================================================
# Optional refinement (Lloyd-style, opt-in via --refinement)
# ============================================================================

def refine_clusters_py(anchor_counts, mapping, threshold,
                       iterations=3, verbose=False):
    """
    Lloyd-style refinement on top of greedy clustering output (pure-Python
    reference; the Rust backend ``pepcluster._core.refine_clusters`` is a
    drop-in equivalent).

    Each pass performs three sub-steps:
      1. Medoid update     — replace each centroid with the member that
                              maximises frequency-weighted mean similarity
                              to the cluster's other members.
      2. Cross-block reassign — for each anchor, find the best centroid
                              above ``threshold`` across its own block plus
                              neighbouring blocks (blocks sharing at least
                              one coarse-alphabet key at P2 or PΩ). Move
                              if strictly more similar than current.
      3. Centroid merge    — if two centroids satisfy sim >= threshold,
                              absorb the smaller cluster into the larger.

    Stops early when no change occurs in a pass.

    Args:
        anchor_counts: dict[str, int]  unique anchor -> peptide frequency
        mapping:       dict[str, str]  unique anchor -> centroid
        threshold:     float           same threshold used during clustering
        iterations:    int             max refinement passes
        verbose:       bool            print per-pass stats

    Returns:
        (refined_mapping, stats_dict)
    """
    str_to_ords = {a: _to_ords(a) for a in anchor_counts}
    cur_mapping = dict(mapping)  # do not mutate caller's dict

    def cluster_members(m):
        out = defaultdict(list)
        for a, c in m.items():
            out[c].append(a)
        return out

    initial_clusters = len(set(cur_mapping.values()))
    total_medoid_changes = 0
    total_reassignments = 0
    total_merges = 0
    passes_run = 0

    for pass_idx in range(iterations):
        passes_run = pass_idx + 1
        changed = False

        # ── 1. Medoid update ───────────────────────────────────────────
        members = cluster_members(cur_mapping)
        centroid_remap = {}  # old centroid -> new centroid
        pass_medoid_changes = 0
        for centroid, mems in members.items():
            if len(mems) == 1:
                centroid_remap[centroid] = centroid
                continue
            best_member = centroid
            best_score = -1e18
            for i, mi in enumerate(mems):
                ords_i = str_to_ords[mi]
                s = 0.0
                w = 0.0
                for j, mj in enumerate(mems):
                    if i == j:
                        continue
                    sim = anchor_sim_fast(ords_i, str_to_ords[mj], -1.0)
                    fj = anchor_counts[mj]
                    s += fj * sim
                    w += fj
                avg = s / w if w > 0 else 0.0
                if avg > best_score:
                    best_score = avg
                    best_member = mi
            centroid_remap[centroid] = best_member
            if best_member != centroid:
                pass_medoid_changes += 1
                changed = True
        total_medoid_changes += pass_medoid_changes
        cur_mapping = {a: centroid_remap[c] for a, c in cur_mapping.items()}

        # ── 2. Cross-block reassignment ────────────────────────────────
        centroids_in_block = defaultdict(list)
        for c in set(cur_mapping.values()):
            centroids_in_block[block_key_fast(str_to_ords[c])].append(c)

        all_blocks = list(centroids_in_block.keys())
        neighbours_of = {}
        for bk in all_blocks:
            p, q = bk
            neighbours_of[bk] = [b for b in all_blocks
                                 if b[0] == p or b[1] == q]

        pass_reassigns = 0
        for a in anchor_counts:
            a_ords = str_to_ords[a]
            bk = block_key_fast(a_ords)
            cur_c = cur_mapping[a]
            cur_score = anchor_sim_fast(a_ords, str_to_ords[cur_c], -1.0)
            best_c, best_score = cur_c, cur_score
            for nb in neighbours_of.get(bk, [bk]):
                for c in centroids_in_block[nb]:
                    if c == cur_c:
                        continue
                    sc = anchor_sim_fast(a_ords, str_to_ords[c], threshold)
                    if sc < 0:
                        continue
                    if sc > best_score:
                        best_score = sc
                        best_c = c
            if best_c != cur_c and best_score >= threshold:
                cur_mapping[a] = best_c
                pass_reassigns += 1
                changed = True
        total_reassignments += pass_reassigns

        # ── 3. Centroid merge ──────────────────────────────────────────
        members = cluster_members(cur_mapping)
        centroid_freq = {c: sum(anchor_counts[m] for m in mems)
                         for c, mems in members.items()}
        sorted_cents = sorted(centroid_freq.keys(),
                              key=lambda c: -centroid_freq[c])
        centroids_in_block = defaultdict(list)
        block_of = {}
        for c in sorted_cents:
            bk = block_key_fast(str_to_ords[c])
            centroids_in_block[bk].append(c)
            block_of[c] = bk

        merge_map = {}    # absorbed centroid -> absorbing centroid
        absorbed = set()
        for c1 in sorted_cents:
            if c1 in absorbed:
                continue
            p1, q1 = block_of[c1]
            for bk2 in (b for b in centroids_in_block
                        if b[0] == p1 or b[1] == q1):
                for c2 in centroids_in_block[bk2]:
                    if c2 == c1 or c2 in absorbed:
                        continue
                    if centroid_freq[c2] >= centroid_freq[c1]:
                        continue  # only smaller absorbed into larger
                    sc = anchor_sim_fast(str_to_ords[c1], str_to_ords[c2],
                                         threshold)
                    if sc >= threshold:
                        merge_map[c2] = c1
                        absorbed.add(c2)
                        total_merges += 1
                        changed = True

        if merge_map:
            def resolve(c):
                while c in merge_map:
                    c = merge_map[c]
                return c
            cur_mapping = {a: resolve(c) for a, c in cur_mapping.items()}

        if verbose:
            n_now = len(set(cur_mapping.values()))
            print(f"      pass {pass_idx + 1}: "
                  f"medoid_changes={pass_medoid_changes}, "
                  f"reassignments={pass_reassigns}, "
                  f"merges={len(merge_map)}, "
                  f"clusters_now={n_now:,}", flush=True)

        if not changed:
            break

    stats = {
        "passes":           passes_run,
        "medoid_changes":   total_medoid_changes,
        "reassignments":    total_reassignments,
        "merges":           total_merges,
        "initial_clusters": initial_clusters,
        "final_clusters":   len(set(cur_mapping.values())),
    }
    return cur_mapping, stats


# ============================================================================
# Cluster representatives (medoid anchors)
# ============================================================================

def cluster_representatives(anchor_counts, mapping):
    """
    Find the representative (medoid) anchor of every cluster.

    The representative is the member anchor with the **least average distance**
    to all peptides in its cluster — equivalently, the highest frequency-
    weighted average similarity. Any peptide carrying this anchor is a valid
    "central" representative of the cluster.

    Fast by construction: the similarity is an additive sum over the 6 anchor
    positions, so instead of the naive O(k^2) all-pairs medoid we aggregate,
    per position, a peptide-weighted amino-acid frequency over the cluster and
    score each member in O(1) per position. Total cost is O(total_anchors)
    plus a tiny per-cluster constant — no pairwise loop.

    For a peptide with anchor ``a``, its total similarity to the whole cluster
    is ``Sigma(a) = sum_p sum_j count_j * wsim_p(a_p, a_j_p)``. Every peptide
    sharing anchor ``a`` has the same average distance, so the medoid is simply
    ``argmax_a Sigma(a)`` over the cluster's member anchors.

    Args:
        anchor_counts: dict[str, int] — unique anchor → peptide frequency
        mapping:       dict[str, str] — anchor → centroid (clustering output)

    Returns:
        dict[str, str] — centroid anchor → medoid (representative) anchor
    """
    groups = defaultdict(list)
    for a, c in mapping.items():
        groups[c].append(a)

    reps = {}
    for centroid, anchors in groups.items():
        if len(anchors) == 1:
            reps[centroid] = anchors[0]
            continue

        # Peptide-weighted amino-acid frequency at each of the 6 positions,
        # stored as (ord(aa), weight) lists for fast inner loops.
        freq = [defaultdict(float) for _ in range(6)]
        for a in anchors:
            cnt = anchor_counts[a]
            f0, f1, f2, f3, f4, f5 = freq
            f0[a[0]] += cnt
            f1[a[1]] += cnt
            f2[a[2]] += cnt
            f3[a[3]] += cnt
            f4[a[4]] += cnt
            f5[a[5]] += cnt

        # Per-position, weighted similarity score for each residue present:
        #   sp[p][x] = sum_aa freq_p[aa] * WSIM[p][ord(x)*128 + ord(aa)]
        sp = []
        for p in range(6):
            wp = WSIM[p]
            present = [(ord(aa), w) for aa, w in freq[p].items()]
            tbl = {}
            for res in freq[p]:
                base = ord(res) * 128
                s = 0.0
                for oaa, w in present:
                    s += w * wp[base + oaa]
                tbl[res] = s
            sp.append(tbl)

        sp0, sp1, sp2, sp3, sp4, sp5 = sp
        # Deterministic argmax: highest Sigma, then highest count, then the
        # lexicographically smallest anchor (via sorted iteration + strict >).
        best_anchor = None
        best_key = None
        for a in sorted(anchors):
            sigma = (sp0[a[0]] + sp1[a[1]] + sp2[a[2]]
                     + sp3[a[3]] + sp4[a[4]] + sp5[a[5]])
            key = (sigma, anchor_counts[a])
            if best_key is None or key > best_key:
                best_key = key
                best_anchor = a
        reps[centroid] = best_anchor

    return reps


# ============================================================================
# Backend resolution
# ============================================================================

def _resolve_backend(backend="auto"):
    """
    Return ``(cluster_fn, refine_fn, name)`` for the requested backend.

    ``backend`` may be "auto" (Rust if available, else Python), "rust", or
    "python". Raises ImportError if "rust" is requested but unavailable.
    """
    if backend in ("auto", "rust"):
        try:
            from ._core import cluster_anchors, refine_clusters
            return cluster_anchors, refine_clusters, "rust"
        except ImportError:
            if backend == "rust":
                raise
    return cluster_anchors_py, refine_clusters_py, "python"


# ============================================================================
# End-to-end pipeline
# ============================================================================

def cluster_fasta(input, outdir, threshold=0.6, min_cluster_size=2,
                  n_front=3, n_back=3, refinement=False, iterations=3,
                  backend="auto", verbose=True):
    """
    Run the full anchor-clustering pipeline on a FASTA file.

    Reads peptides, clusters them by anchor similarity, optionally refines,
    assigns peptides to clusters, and writes:
        <outdir>/clusters.tsv          cluster assignment for every peptide
        <outdir>/cluster_summary.tsv   cluster stats sorted by size
        <outdir>/fasta/cluster_*.fasta per-cluster FASTA (>= min_cluster_size)
        <outdir>/summary.txt           run statistics

    Args:
        input:            path to input FASTA
        outdir:           output directory (created if missing)
        threshold:        BLOSUM similarity threshold in [0, 1] (default 0.6)
        min_cluster_size: min members for a per-cluster FASTA (default 2)
        n_front:          N-terminal anchor length (default 3)
        n_back:           C-terminal anchor length (default 3)
        refinement:       apply Lloyd-style refinement after greedy clustering
        iterations:       max refinement passes (only if refinement=True)
        backend:          "auto" | "rust" | "python"
        verbose:          print progress to stdout

    Returns:
        dict of run statistics (also written to summary.txt).
    """
    cluster_fn, refine_fn, backend_name = _resolve_backend(backend)

    outdir = Path(outdir)
    fastadir = outdir / "fasta"
    outdir.mkdir(parents=True, exist_ok=True)
    fastadir.mkdir(exist_ok=True)

    def log(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    t0 = time.time()

    # ── Step 1: read & extract ────────────────────────────────────────
    log("[1/4] Reading FASTA and extracting anchors …", flush=True)
    peptides = []
    short = []
    acounts = defaultdict(int)
    anchor_to_pep = {}  # anchor -> first (header, sequence) seen, for representatives

    for hdr, seq in parse_fasta(input):
        seq = seq.upper()
        anc = extract_anchor(seq, n_front, n_back)
        if anc is None:
            short.append((hdr, seq))
        else:
            peptides.append((hdr, seq, anc))
            acounts[anc] += 1
            if anc not in anchor_to_pep:
                anchor_to_pep[anc] = (hdr, seq)

    n_total = len(peptides) + len(short)
    n_unique = len(acounts)
    log(f"      {n_total:>12,}  total peptides")
    log(f"      {len(peptides):>12,}  valid (>={n_front + n_back} aa)")
    log(f"      {len(short):>12,}  too short")
    log(f"      {n_unique:>12,}  unique anchors")

    # ── Step 2: cluster unique anchors ────────────────────────────────
    log(f"\n[2/4] Clustering unique anchors (threshold {threshold}, "
        f"backend: {backend_name}) …", flush=True)
    mapping, n_cmp, n_early = cluster_fn(dict(acounts), threshold)
    n_clusters = len(set(mapping.values()))
    log(f"      {n_clusters:>12,}  clusters")
    log(f"      {n_cmp:>12,}  pairwise comparisons (within blocks)")
    early_pct = 100.0 * n_early / n_cmp if n_cmp else 0
    log(f"      {n_early:>12,}  early-terminated ({early_pct:.1f}%)")

    # ── Step 2.5: optional refinement ─────────────────────────────────
    refine_stats = None
    if refinement:
        log(f"\n[2.5/4] Refining clusters (max {iterations} passes, "
            f"backend: {backend_name}) …", flush=True)
        mapping, refine_stats = refine_fn(
            dict(acounts), mapping, threshold, iterations)
        log(f"      passes run:      {refine_stats['passes']}")
        log(f"      medoid changes:  {refine_stats['medoid_changes']:,}")
        log(f"      reassignments:   {refine_stats['reassignments']:,}")
        log(f"      merges:          {refine_stats['merges']:,}")
        log(f"      clusters: {refine_stats['initial_clusters']:,} "
            f"→ {refine_stats['final_clusters']:,}")
        n_clusters = refine_stats["final_clusters"]

    # ── Step 3: assign peptides → clusters + pick representatives ──────
    log("\n[3/4] Assigning peptides and finding representatives …", flush=True)
    clusters = defaultdict(list)
    for hdr, seq, anc in peptides:
        clusters[mapping[anc]].append((hdr, seq, anc))

    ranked = sorted(clusters.items(), key=lambda kv: (-len(kv[1]), kv[0]))

    # Central (medoid) peptide per cluster: the member with the least average
    # distance to the whole cluster. O(N), see cluster_representatives().
    medoid_anchor = cluster_representatives(dict(acounts), mapping)
    rep_seq_of = {ctr: anchor_to_pep[medoid_anchor[ctr]][1] for ctr, _ in ranked}

    # ── Step 4: write outputs ─────────────────────────────────────────
    log("[4/4] Writing outputs …", flush=True)

    # 4a. cluster_summary.tsv
    with open(outdir / "cluster_summary.tsv", "w") as f:
        f.write("cluster_id\trepresentative_anchor\trepresentative_peptide"
                "\tsize\n")
        for idx, (ctr, members) in enumerate(ranked):
            f.write(f"cluster_{idx}\t{ctr}\t{rep_seq_of[ctr]}\t{len(members)}\n")

    # 4b. clusters.tsv (full mapping)
    with open(outdir / "clusters.tsv", "w") as f:
        f.write("cluster_id\trepresentative_anchor\trepresentative_peptide"
                "\tpeptide_header\tsequence\tanchor\n")
        for idx, (ctr, members) in enumerate(ranked):
            cname = f"cluster_{idx}"
            rep_seq = rep_seq_of[ctr]
            for hdr, seq, anc in members:
                f.write(f"{cname}\t{ctr}\t{rep_seq}\t{hdr}\t{seq}\t{anc}\n")

    # 4c. per-cluster FASTA
    n_fasta = 0
    for idx, (ctr, members) in enumerate(ranked):
        if len(members) < min_cluster_size:
            continue
        with open(fastadir / f"cluster_{idx}.fasta", "w") as f:
            for hdr, seq, _ in members:
                f.write(f">{hdr}\n{seq}\n")
        n_fasta += 1

    # 4d. short peptides
    if short:
        with open(fastadir / "SHORT_peptides.fasta", "w") as f:
            for hdr, seq in short:
                f.write(f">{hdr}\n{seq}\n")

    # ── Summary ───────────────────────────────────────────────────────
    elapsed = time.time() - t0
    sizes = [len(m) for _, m in ranked]

    report = [
        "=" * 62,
        "  ANCHOR CLUSTERING SUMMARY",
        "=" * 62,
        f"  Input:               {input}",
        f"  Total peptides:      {n_total:,}",
        f"  Valid peptides:      {len(peptides):,}",
        f"  Too short:           {len(short):,}",
        f"  Unique anchors:      {n_unique:,}",
        f"  Clusters:            {n_clusters:,}",
        f"  FASTA files:         {n_fasta:,}  (>={min_cluster_size} members)",
        "",
        f"  Threshold:           {threshold}",
        f"  Backend:             {backend_name}",
        f"  Matrix:              BLOSUM62 (normalized)",
        f"  Weights:             P1={_RAW_W[0]:.0f}  P2={_RAW_W[1]:.0f}  "
        f"P3={_RAW_W[2]:.0f}  PO-2={_RAW_W[3]:.0f}  PO-1={_RAW_W[4]:.0f}  "
        f"PO={_RAW_W[5]:.0f}",
        f"  Blocking:            coarse alphabet at P2 + PO  (10x10 = 100 bins)",
        f"  Comparisons:         {n_cmp:,}",
        f"  Early-terminated:    {n_early:,}  ({early_pct:.1f}%)",
        "",
    ]

    if refine_stats is not None:
        report += [
            "  Refinement:          ENABLED",
            f"    Backend:           {backend_name}",
            f"    Max iterations:    {iterations}",
            f"    Passes run:        {refine_stats['passes']}",
            f"    Medoid changes:    {refine_stats['medoid_changes']:,}",
            f"    Reassignments:     {refine_stats['reassignments']:,}",
            f"    Merges:            {refine_stats['merges']:,}",
            f"    Clusters before:   {refine_stats['initial_clusters']:,}",
            f"    Clusters after:    {refine_stats['final_clusters']:,}",
            "",
        ]
    else:
        report += [
            "  Refinement:          disabled (pass refinement=True to enable)",
            "",
        ]

    report += [
        "  Cluster size distribution:",
        f"    singletons:        {sum(1 for s in sizes if s == 1):,}",
        f"    2-10:              {sum(1 for s in sizes if 2 <= s <= 10):,}",
        f"    11-100:            {sum(1 for s in sizes if 11 <= s <= 100):,}",
        f"    101-1000:          {sum(1 for s in sizes if 101 <= s <= 1000):,}",
        f"    >1000:             {sum(1 for s in sizes if s > 1000):,}",
        f"    largest:           {max(sizes):,}" if sizes else "",
        "",
        f"  Elapsed:             {elapsed:.1f} s",
        "=" * 62,
    ]

    with open(outdir / "summary.txt", "w") as f:
        for line in report:
            f.write(line + "\n")

    if verbose:
        print()
        for line in report:
            print(line)

    stats = {
        "input":            str(input),
        "outdir":           str(outdir),
        "total_peptides":   n_total,
        "valid_peptides":   len(peptides),
        "too_short":        len(short),
        "unique_anchors":   n_unique,
        "n_clusters":       n_clusters,
        "n_fasta":          n_fasta,
        "threshold":        threshold,
        "backend":          backend_name,
        "comparisons":      n_cmp,
        "early_terminated": n_early,
        "refinement":       refine_stats,
        "elapsed_sec":      elapsed,
    }
    return stats
