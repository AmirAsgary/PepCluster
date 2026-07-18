"""Type stubs for the compiled Rust extension ``pepcluster._core``."""

from typing import Dict, Sequence, Tuple

def cluster_anchors(
    anchor_counts: Dict[str, int],
    threshold: float,
    anchor_positions: Sequence[int] = ...,
    anchor_weight: float = 2.0,
    obg_block_search: bool = False,
    obg_max_probes: int = 0,
    obg_min_block_upper_bound: float = 0.0,
) -> Tuple[Dict[str, str], int, int]:
    """Greedy centroid clustering of unique anchors.

    ``anchor_positions`` are 0-based positions within the anchor that count as
    binding anchors: they carry ``anchor_weight`` in the similarity and define
    the coarse-alphabet blocking (default ``[1, 5]`` = P2 and PΩ of a 6-mer).

    ``obg_block_search`` enables Upper-Bound-Guided multi-probe search: each
    anchor also considers centroids in neighbouring blocks whose upper bound
    reaches the threshold (ranked by bound, own block first). ``obg_max_probes``
    caps blocks searched per anchor (incl. own; ``<=0`` = all).
    ``obg_min_block_upper_bound`` raises the eligibility cut to
    ``max(threshold, this)``.

    Returns ``(mapping, n_comparisons, n_early_exits)`` where ``mapping`` maps
    each anchor to its centroid anchor.
    """
    ...

def refine_clusters(
    anchor_counts: Dict[str, int],
    mapping: Dict[str, str],
    threshold: float,
    iterations: int,
    cap: int = 32,
    merge: bool = True,
    anchor_positions: Sequence[int] = ...,
    anchor_weight: float = 2.0,
    fast_medoid: bool = False,
    merge_cap: int = 0,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Lloyd-style refinement of an existing clustering.

    ``cap`` bounds the centroid comparisons per anchor in the reassignment step
    (candidates examined own-block-first, largest-cluster-first; ``<= 0`` means
    no cap). ``merge`` toggles the centroid-merge sub-step. ``fast_medoid`` uses
    the O(N) medoid decomposition instead of the exact O(k^2) medoid.
    ``merge_cap`` bounds candidate centroids examined per centroid in the merge
    step (``<= 0`` = no cap). ``anchor_positions`` / ``anchor_weight`` must match
    those used for clustering.

    Returns ``(refined_mapping, stats)``. ``stats`` has keys: ``passes``,
    ``medoid_changes``, ``reassignments``, ``merges``, ``initial_clusters``,
    ``final_clusters``.
    """
    ...
