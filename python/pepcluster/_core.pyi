"""Type stubs for the compiled Rust extension ``pepcluster._core``."""

from typing import Dict, Tuple

def cluster_anchors(
    anchor_counts: Dict[str, int],
    threshold: float,
) -> Tuple[Dict[str, str], int, int]:
    """Greedy centroid clustering of unique 6-mer anchors.

    Returns ``(mapping, n_comparisons, n_early_exits)`` where ``mapping`` maps
    each anchor to its centroid anchor.
    """
    ...

def refine_clusters(
    anchor_counts: Dict[str, int],
    mapping: Dict[str, str],
    threshold: float,
    iterations: int,
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Lloyd-style refinement of an existing clustering.

    Returns ``(refined_mapping, stats)``. ``stats`` has keys: ``passes``,
    ``medoid_changes``, ``reassignments``, ``merges``, ``initial_clusters``,
    ``final_clusters``.
    """
    ...
