"""
clustering.py — BLOSUM-aware anchor clustering for immunopeptides.

Clusters peptides by the similarity of their anchor residues (the first
``n_front`` + last ``n_back`` amino acids) using a BLOSUM62-normalized
similarity metric with position-specific weighting.

Which positions inside the anchor count as *binding anchors* is configurable
(``--anchors`` / ``anchor_positions``). Anchor positions carry ``anchor_weight``
(default 2x) in the similarity and define the coarse-alphabet blocking. The
default — the 2nd of the first 3 residues and the 3rd of the last 3, i.e. P2
and PΩ — is the classic MHC-I choice.

This module contains:
  * the pure-Python reference implementation of the clustering
    (`cluster_anchors_py`) and refinement (`refine_clusters_py`);
  * FASTA parsing / anchor extraction helpers;
  * `cluster_representatives()` — the central (medoid) anchor of each cluster;
  * `cluster_fasta()` — the full end-to-end pipeline, callable as a library or
    from the CLI.

The Rust backend (`pepcluster._core`) exposes drop-in equivalents of the two
core functions and is preferred automatically when available. Both backends
compute in f64 and use identical, deterministic orderings, so their results are
bit-identical.
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

# Coarse-group lookup array indexed by ord(); 255 = unknown residue.
# (255, not -1, so the packed block keys match the Rust backend byte-for-byte.)
_COARSE_ORD = [255] * 128
for _aa, _gid in COARSE.items():
    _COARSE_ORD[ord(_aa)] = _gid

# Defaults: 6-residue anchor (3 front + 3 back); anchors at P2 and PΩ.
DEFAULT_N_FRONT = 3
DEFAULT_N_BACK = 3
DEFAULT_ANCHOR_POSITIONS = (1, 5)
DEFAULT_ANCHOR_WEIGHT = 2.0
DEFAULT_ANCHOR_SPEC = "2;3"

# At most 8 anchor positions, so a block key packs into one byte each.
MAX_ANCHOR_POSITIONS = 8


# ============================================================================
# Per-configuration similarity / blocking tables
# ============================================================================

class Tables:
    """Precomputed weighted-similarity and blocking tables for one
    (anchor length, anchor positions, anchor weight) configuration."""

    __slots__ = ("alen", "anchor_positions", "anchor_weight", "weights",
                 "wsim_by_pos", "wsim_ord", "check_pos", "remaining_after")

    def __init__(self, alen, anchor_positions, anchor_weight):
        self.alen = alen
        self.anchor_positions = tuple(anchor_positions)
        self.anchor_weight = anchor_weight

        # Raw weights: anchor_weight at anchor positions, 1.0 elsewhere.
        raw = [1.0] * alen
        for p in anchor_positions:
            raw[p] = anchor_weight
        wsum = sum(raw)
        self.weights = [w / wsum for w in raw]

        # Weighted 128x128 similarity table per position (weight baked in).
        self.wsim_by_pos = []
        for p in range(alen):
            arr = [0.0] * (128 * 128)
            wp = self.weights[p]
            for (a, b), v in SIM.items():
                arr[ord(a) * 128 + ord(b)] = wp * v
            self.wsim_by_pos.append(arr)

        # Check order: heaviest positions first (ties by position ascending),
        # so early termination kicks in sooner. Must match the Rust backend.
        self.check_pos = sorted(range(alen), key=lambda p: (-raw[p], p))
        self.wsim_ord = [self.wsim_by_pos[p] for p in self.check_pos]

        # remaining_after[k] = max score still obtainable from check_pos[k:]
        self.remaining_after = [0.0] * (alen + 1)
        for k in range(alen - 1, -1, -1):
            self.remaining_after[k] = (self.remaining_after[k + 1]
                                       + self.weights[self.check_pos[k]])


def check_positions(anchor_positions, alen):
    """Validate and canonicalise (sort, dedup) anchor positions."""
    v = sorted(set(int(p) for p in anchor_positions))
    if not v:
        raise ValueError("anchor_positions must contain at least one position")
    if len(v) > MAX_ANCHOR_POSITIONS:
        raise ValueError(
            f"at most {MAX_ANCHOR_POSITIONS} anchor positions are supported, "
            f"got {len(v)}")
    if v[0] < 0 or v[-1] >= alen:
        raise ValueError(
            f"anchor position {v[-1] if v[-1] >= alen else v[0]} is out of "
            f"range for anchor length {alen}")
    return v


def build_tables(alen=DEFAULT_N_FRONT + DEFAULT_N_BACK,
                 anchor_positions=DEFAULT_ANCHOR_POSITIONS,
                 anchor_weight=DEFAULT_ANCHOR_WEIGHT):
    """Build :class:`Tables` for a given anchor configuration."""
    return Tables(alen, check_positions(anchor_positions, alen), anchor_weight)


DEFAULT_TABLES = build_tables()


def parse_anchors(spec=DEFAULT_ANCHOR_SPEC, n_front=DEFAULT_N_FRONT,
                  n_back=DEFAULT_N_BACK):
    """
    Parse an ``--anchors`` spec into 0-based positions within the anchor.

    Format: ``"<front>;<back>"``, where each side is a comma-separated list of
    **1-based** indices into that side's residues. Either side may be empty.

    The default ``"2;3"`` means: the 2nd of the first ``n_front`` residues (P2)
    and the 3rd of the last ``n_back`` residues (PΩ) — i.e. positions ``[1, 5]``
    of a 6-residue anchor.

    Examples (with n_front=3, n_back=3):
        "2;3"    -> [1, 5]      P2 and PΩ  (default, MHC-I)
        "2;2,3"  -> [1, 4, 5]   P2 plus the last two C-terminal residues
        "1,2;3"  -> [0, 1, 5]   first two N-terminal residues plus PΩ
        ";3"     -> [5]         C-terminal anchor only
        "2;"     -> [1]         N-terminal anchor only
    """
    if spec is None:
        spec = DEFAULT_ANCHOR_SPEC
    text = str(spec).strip()
    parts = text.split(";")
    if len(parts) != 2:
        raise ValueError(
            f"--anchors must look like 'FRONT;BACK' (e.g. '2;3'), got {spec!r}")

    positions = []
    for part, size, offset, name in ((parts[0], n_front, 0, "front"),
                                     (parts[1], n_back, n_front, "back")):
        part = part.strip()
        if not part:
            continue
        for tok in part.split(","):
            tok = tok.strip()
            if not tok:
                continue
            try:
                idx = int(tok)
            except ValueError:
                raise ValueError(
                    f"--anchors: {tok!r} is not an integer") from None
            if not 1 <= idx <= size:
                raise ValueError(
                    f"--anchors: {name} index {idx} is out of range 1..{size} "
                    f"(n-{name} = {size})")
            positions.append(offset + idx - 1)

    return check_positions(positions, n_front + n_back)


# ============================================================================
# Core primitives
# ============================================================================

def extract_anchor(seq, nf=DEFAULT_N_FRONT, nb=DEFAULT_N_BACK):
    """Return first ``nf`` + last ``nb`` residues, or None if too short."""
    if len(seq) < nf + nb:
        return None
    return seq[:nf] + seq[-nb:]


def _to_ords(anchor):
    """Convert an anchor string to a tuple of ord values."""
    return tuple(ord(c) for c in anchor)


def anchor_sim_fast(a, b, threshold, tables=DEFAULT_TABLES):
    """
    Weighted BLOSUM62-normalized similarity with early termination.

    Positions are checked heaviest-first; the scan bails out as soon as the
    remaining positions cannot lift the score to ``threshold``, returning -1.0.
    Pass ``threshold=-1.0`` to always get the full score.

    ``a`` and ``b`` are ordinal tuples (see :func:`_to_ords`).
    """
    s = 0.0
    n = tables.alen
    check_pos = tables.check_pos
    wsim_ord = tables.wsim_ord
    remaining_after = tables.remaining_after
    for k in range(n):
        p = check_pos[k]
        s += wsim_ord[k][a[p] * 128 + b[p]]
        if k + 1 < n and s + remaining_after[k + 1] < threshold:
            return -1.0
    return s


def block_key_fast(ords, tables=DEFAULT_TABLES):
    """
    Coarse-alphabet block key: one byte per anchor position, first anchor
    position in the most-significant byte (integer order == tuple order).
    """
    key = 0
    for p in tables.anchor_positions:
        key = (key << 8) | _COARSE_ORD[ords[p]]
    return key


def _is_neighbour(k1, k2, n_ap):
    """True if two block keys agree on at least one anchor position."""
    for i in range(n_ap):
        sh = 8 * (n_ap - 1 - i)
        if ((k1 >> sh) & 0xFF) == ((k2 >> sh) & 0xFF):
            return True
    return False


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


def _tables_for(anchor_counts, anchor_positions, anchor_weight):
    """Infer anchor length from the data and build the matching tables."""
    alen = None
    for a in anchor_counts:
        if alen is None:
            alen = len(a)
        elif len(a) != alen:
            raise ValueError(
                f"all anchors must have the same length; found {alen} and "
                f"{len(a)}")
    if not alen:
        raise ValueError("anchors must be non-empty")
    return build_tables(alen, anchor_positions, anchor_weight)


# ============================================================================
# Greedy clustering
# ============================================================================

def cluster_anchors_py(anchor_counts, threshold,
                       anchor_positions=DEFAULT_ANCHOR_POSITIONS,
                       anchor_weight=DEFAULT_ANCHOR_WEIGHT):
    """
    Greedy centroid clustering on unique anchors (pure-Python reference; the
    Rust backend ``pepcluster._core.cluster_anchors`` is a drop-in equivalent).

    1. Block anchors by the coarse alphabet at the anchor positions
    2. Within each block, process anchors most-frequent-first
    3. Assign to the first centroid above threshold, or become a new centroid

    Args:
        anchor_counts:    dict[str, int] — unique anchor → peptide frequency
        threshold:        float          — minimum similarity to join a cluster
        anchor_positions: 0-based anchor positions (weighted up + blocked on)
        anchor_weight:    weight of anchor positions (others are 1.0)

    Returns:
        (mapping, n_comparisons, n_early_exits)
        mapping: dict[str, str] — anchor → centroid anchor
    """
    if not anchor_counts:
        return {}, 0, 0
    tables = _tables_for(anchor_counts, anchor_positions, anchor_weight)

    # Sort by frequency descending (stable: ties keep insertion order).
    items = sorted(anchor_counts.items(), key=lambda x: -x[1])
    str_to_ords = {}
    blocks = defaultdict(list)
    for anchor_str, _cnt in items:
        ords = _to_ords(anchor_str)
        str_to_ords[anchor_str] = ords
        blocks[block_key_fast(ords, tables)].append(anchor_str)

    mapping = {}
    n_cmp = 0
    n_early = 0

    for _bk, anchors in blocks.items():
        centroids_ords = []
        centroids_str = []
        for anchor_str in anchors:
            a_ords = str_to_ords[anchor_str]
            matched = False
            for ci in range(len(centroids_ords)):
                n_cmp += 1
                score = anchor_sim_fast(a_ords, centroids_ords[ci], threshold,
                                        tables)
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
# Optional refinement (Lloyd-style)
# ============================================================================

def refine_clusters_py(anchor_counts, mapping, threshold, iterations=3,
                       cap=32, merge=True,
                       anchor_positions=DEFAULT_ANCHOR_POSITIONS,
                       anchor_weight=DEFAULT_ANCHOR_WEIGHT, verbose=False):
    """
    Lloyd-style refinement on top of greedy clustering output (pure-Python
    reference; the Rust backend ``pepcluster._core.refine_clusters`` is a
    drop-in equivalent).

    Each pass performs three sub-steps:
      1. Medoid update     — replace each centroid with the member that
                              maximises frequency-weighted mean similarity
                              to the cluster's other members.
      2. Cross-block reassign — for each anchor, find the best centroid above
                              ``threshold`` across its own block plus
                              neighbouring blocks (blocks agreeing on at least
                              one anchor position), bounded by ``cap``
                              comparisons.
      3. Centroid merge    — if two centroids satisfy sim >= threshold,
                              absorb the smaller cluster into the larger.
                              Skipped when ``merge`` is False.

    Stops early when no change occurs in a pass.

    Args:
        anchor_counts:    dict[str, int]  unique anchor -> peptide frequency
        mapping:          dict[str, str]  unique anchor -> centroid
        threshold:        float           same threshold used during clustering
        iterations:       int             max refinement passes
        cap:              int             max centroid comparisons per anchor in
                                          the reassignment step (candidate cap;
                                          own-block-first, largest-cluster-first).
                                          <= 0 means no cap.
        merge:            bool            run the centroid-merge sub-step
        anchor_positions: 0-based anchor positions
        anchor_weight:    weight of anchor positions
        verbose:          bool            print per-pass stats

    Returns:
        (refined_mapping, stats_dict)
    """
    if not anchor_counts:
        return {}, {"passes": 0, "medoid_changes": 0, "reassignments": 0,
                    "merges": 0, "initial_clusters": 0, "final_clusters": 0}
    tables = _tables_for(anchor_counts, anchor_positions, anchor_weight)
    n_ap = len(tables.anchor_positions)

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
            # Canonical member order (anchor ascending) so the medoid argmax
            # and its ties resolve identically to the Rust backend.
            mems = sorted(mems)
            best_member = centroid
            best_score = -1e18
            for i, mi in enumerate(mems):
                ords_i = str_to_ords[mi]
                s = 0.0
                w = 0.0
                for j, mj in enumerate(mems):
                    if i == j:
                        continue
                    sim = anchor_sim_fast(ords_i, str_to_ords[mj], -1.0, tables)
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

        # ── 2. Cross-block reassignment (candidate-capped) ─────────────
        members = cluster_members(cur_mapping)
        centroid_freq = {c: sum(anchor_counts[m] for m in mems)
                         for c, mems in members.items()}

        centroids_in_block = defaultdict(list)
        for c in centroid_freq:
            centroids_in_block[block_key_fast(str_to_ords[c], tables)].append(c)
        # Larger clusters first, then anchor asc — deterministic order that
        # matches the Rust backend exactly (required once we cap the scan).
        for v in centroids_in_block.values():
            v.sort(key=lambda c: (-centroid_freq[c], c))

        sorted_blocks = sorted(centroids_in_block.keys())
        neighbours_of = {}
        for bk in centroids_in_block:
            neighbours_of[bk] = [bk] + [b for b in sorted_blocks
                                        if b != bk and _is_neighbour(bk, b, n_ap)]

        pass_reassigns = 0
        for a in anchor_counts:
            a_ords = str_to_ords[a]
            bk = block_key_fast(a_ords, tables)
            cur_c = cur_mapping[a]
            cur_score = anchor_sim_fast(a_ords, str_to_ords[cur_c], -1.0, tables)
            best_c, best_score = cur_c, cur_score
            examined = 0
            stop = False
            for nb in neighbours_of.get(bk, ()):
                if stop:
                    break
                for c in centroids_in_block.get(nb, ()):
                    if c == cur_c:
                        continue
                    if cap > 0 and examined >= cap:
                        stop = True
                        break
                    sc = anchor_sim_fast(a_ords, str_to_ords[c], threshold,
                                         tables)
                    examined += 1
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

        # ── 3. Centroid merge (optional; skipped when merge=False) ─────
        merge_map = {}    # absorbed centroid -> absorbing centroid
        if merge:
            members = cluster_members(cur_mapping)
            centroid_freq = {c: sum(anchor_counts[m] for m in mems)
                             for c, mems in members.items()}
            sorted_cents = sorted(centroid_freq.keys(),
                                  key=lambda c: (-centroid_freq[c], c))
            centroids_in_block = defaultdict(list)
            block_of = {}
            for c in sorted_cents:
                bk = block_key_fast(str_to_ords[c], tables)
                centroids_in_block[bk].append(c)
                block_of[c] = bk

            absorbed = set()
            for c1 in sorted_cents:
                if c1 in absorbed:
                    continue
                bk1 = block_of[c1]
                for bk2 in (b for b in centroids_in_block
                            if _is_neighbour(bk1, b, n_ap)):
                    for c2 in centroids_in_block[bk2]:
                        if c2 == c1 or c2 in absorbed:
                            continue
                        if centroid_freq[c2] >= centroid_freq[c1]:
                            continue  # only smaller absorbed into larger
                        sc = anchor_sim_fast(str_to_ords[c1], str_to_ords[c2],
                                             threshold, tables)
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

def cluster_representatives(anchor_counts, mapping,
                            anchor_positions=DEFAULT_ANCHOR_POSITIONS,
                            anchor_weight=DEFAULT_ANCHOR_WEIGHT):
    """
    Find the representative (medoid) anchor of every cluster.

    The representative is the member anchor with the **least average distance**
    to all peptides in its cluster — equivalently, the highest frequency-
    weighted average similarity. Any peptide carrying this anchor is a valid
    "central" representative of the cluster.

    Fast by construction: the similarity is an additive sum over the anchor
    positions, so instead of the naive O(k^2) all-pairs medoid we aggregate,
    per position, a peptide-weighted amino-acid frequency over the cluster and
    score each member in O(1) per position. Total cost is O(total_anchors)
    plus a tiny per-cluster constant — no pairwise loop.

    For a peptide with anchor ``a``, its total similarity to the whole cluster
    is ``Sigma(a) = sum_p sum_j count_j * wsim_p(a_p, a_j_p)``. Every peptide
    sharing anchor ``a`` has the same average distance, so the medoid is simply
    ``argmax_a Sigma(a)`` over the cluster's member anchors.

    Returns:
        dict[str, str] — centroid anchor → medoid (representative) anchor
    """
    if not anchor_counts:
        return {}
    tables = _tables_for(anchor_counts, anchor_positions, anchor_weight)
    alen = tables.alen
    wsim_by_pos = tables.wsim_by_pos

    groups = defaultdict(list)
    for a, c in mapping.items():
        groups[c].append(a)

    reps = {}
    for centroid, anchors in groups.items():
        if len(anchors) == 1:
            reps[centroid] = anchors[0]
            continue

        # Peptide-weighted amino-acid frequency at each anchor-string position.
        freq = [defaultdict(float) for _ in range(alen)]
        for a in anchors:
            cnt = anchor_counts[a]
            for p in range(alen):
                freq[p][a[p]] += cnt

        # Per-position, weighted score for each residue that occurs there:
        #   sp[p][x] = sum_aa freq_p[aa] * wsim_by_pos[p][ord(x)*128 + ord(aa)]
        sp = []
        for p in range(alen):
            wp = wsim_by_pos[p]
            present = [(ord(aa), w) for aa, w in freq[p].items()]
            tbl = {}
            for res in freq[p]:
                base = ord(res) * 128
                s = 0.0
                for oaa, w in present:
                    s += w * wp[base + oaa]
                tbl[res] = s
            sp.append(tbl)

        # Deterministic argmax: highest Sigma, then highest count, then the
        # lexicographically smallest anchor (via sorted iteration + strict >).
        best_anchor = None
        best_key = None
        for a in sorted(anchors):
            sigma = 0.0
            for p in range(alen):
                sigma += sp[p][a[p]]
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
                  n_front=DEFAULT_N_FRONT, n_back=DEFAULT_N_BACK,
                  anchors=DEFAULT_ANCHOR_SPEC,
                  anchor_weight=DEFAULT_ANCHOR_WEIGHT,
                  refinement=False, iterations=3,
                  refine_cap=32, merge=True, backend="auto", verbose=True):
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
        anchors:          which positions are binding anchors, as "FRONT;BACK"
                          with 1-based indices per side (default "2;3" = P2 and
                          PΩ). See :func:`parse_anchors`.
        anchor_weight:    weight given to anchor positions (default 2.0; all
                          other positions are 1.0)
        refinement:       apply Lloyd-style refinement after greedy clustering
        iterations:       max refinement passes (only if refinement=True)
        refine_cap:       max centroid comparisons per anchor during refinement
                          reassignment (default 32; <= 0 = no cap)
        merge:            run the refinement centroid-merge step (default True)
        backend:          "auto" | "rust" | "python"
        verbose:          print progress to stdout

    Returns:
        dict of run statistics (also written to summary.txt).
    """
    cluster_fn, refine_fn, backend_name = _resolve_backend(backend)

    if n_front < 0 or n_back < 0 or n_front + n_back < 1:
        raise ValueError("n_front + n_back must be at least 1")
    anchor_positions = parse_anchors(anchors, n_front, n_back)
    alen = n_front + n_back

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
    anchor_label = "+".join(str(p + 1) for p in anchor_positions)
    log(f"      {n_total:>12,}  total peptides")
    log(f"      {len(peptides):>12,}  valid (>={alen} aa)")
    log(f"      {len(short):>12,}  too short")
    log(f"      {n_unique:>12,}  unique anchors ({alen}-mer: "
        f"{n_front} front + {n_back} back)")
    log(f"      {'':>12}  anchor positions {anchor_label} "
        f"(weight {anchor_weight}x)")

    if not acounts:
        raise ValueError(
            f"no peptide is long enough for a {alen}-residue anchor "
            f"(n_front={n_front}, n_back={n_back})")

    # ── Step 2: cluster unique anchors ────────────────────────────────
    log(f"\n[2/4] Clustering unique anchors (threshold {threshold}, "
        f"backend: {backend_name}) …", flush=True)
    mapping, n_cmp, n_early = cluster_fn(dict(acounts), threshold,
                                         anchor_positions, anchor_weight)
    n_clusters = len(set(mapping.values()))
    log(f"      {n_clusters:>12,}  clusters")
    log(f"      {n_cmp:>12,}  pairwise comparisons (within blocks)")
    early_pct = 100.0 * n_early / n_cmp if n_cmp else 0
    log(f"      {n_early:>12,}  early-terminated ({early_pct:.1f}%)")

    # ── Step 2.5: optional refinement ─────────────────────────────────
    refine_stats = None
    if refinement:
        log(f"\n[2.5/4] Refining clusters (max {iterations} passes, "
            f"cap {refine_cap}, merge {'on' if merge else 'off'}, "
            f"backend: {backend_name}) …", flush=True)
        mapping, refine_stats = refine_fn(
            dict(acounts), mapping, threshold, iterations, refine_cap, merge,
            anchor_positions, anchor_weight)
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
    medoid_anchor = cluster_representatives(dict(acounts), mapping,
                                            anchor_positions, anchor_weight)
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
    n_blocks = 10 ** len(anchor_positions)

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
        f"  Anchor:              {alen}-mer  ({n_front} front + {n_back} back)",
        f"  Anchor spec:         {anchors}  -> positions {anchor_label}",
        f"  Anchor weight:       {anchor_weight}x  (other positions 1x)",
        f"  Blocking:            coarse alphabet at positions {anchor_label} "
        f"({n_blocks:,} bins)",
        f"  Comparisons:         {n_cmp:,}",
        f"  Early-terminated:    {n_early:,}  ({early_pct:.1f}%)",
        "",
    ]

    if refine_stats is not None:
        report += [
            "  Refinement:          ENABLED",
            f"    Backend:           {backend_name}",
            f"    Max iterations:    {iterations}",
            f"    Candidate cap:     {refine_cap}",
            f"    Merge step:        {'on' if merge else 'off'}",
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
        "n_front":          n_front,
        "n_back":           n_back,
        "anchors":          anchors,
        "anchor_positions": list(anchor_positions),
        "anchor_weight":    anchor_weight,
        "comparisons":      n_cmp,
        "early_terminated": n_early,
        "refinement":       refine_stats,
        "elapsed_sec":      elapsed,
    }
    return stats
