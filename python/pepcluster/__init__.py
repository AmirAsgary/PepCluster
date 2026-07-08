"""
pepcluster — BLOSUM62-aware anchor-residue clustering for immunopeptides.

Groups peptides by the similarity of their MHC-I anchor residues (first 3 +
last 3 amino acids) using a BLOSUM62-normalized metric with double weight on
positions 2 and Ω. A fast Rust backend is used automatically when available,
with a pure-Python fallback.

Typical usage::

    import pepcluster

    # end-to-end: FASTA in, TSV/FASTA cluster files out
    pepcluster.cluster_fasta("peptides.fasta", "out", threshold=0.6)

    # low-level: cluster a dict of unique 6-mer anchors -> frequency
    mapping, n_cmp, n_early = pepcluster.cluster_anchors(
        {"YLAFLV": 3, "YMAFLV": 1}, 0.6)
"""

from .clustering import (
    cluster_fasta,
    cluster_representatives,
    cluster_anchors_py,
    refine_clusters_py,
    extract_anchor,
    parse_fasta,
    anchor_sim_fast,
)

__version__ = "0.1.1"

# Prefer the compiled Rust backend; fall back to the pure-Python reference.
try:
    from ._core import cluster_anchors, refine_clusters  # type: ignore
    HAS_RUST = True
except ImportError:  # pragma: no cover - exercised only when unbuilt
    from .clustering import (
        cluster_anchors_py as cluster_anchors,
        refine_clusters_py as refine_clusters,
    )
    HAS_RUST = False

__all__ = [
    "__version__",
    "HAS_RUST",
    "cluster_fasta",
    "cluster_representatives",
    "cluster_anchors",
    "refine_clusters",
    "cluster_anchors_py",
    "refine_clusters_py",
    "extract_anchor",
    "parse_fasta",
    "anchor_sim_fast",
]
