"""
Command-line interface for pepcluster.

Installed as the ``pepcluster`` console script and runnable via
``python -m pepcluster``.
"""

import argparse

from . import __version__
from .clustering import cluster_fasta


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="pepcluster",
        description="BLOSUM62-aware anchor clustering for immunopeptides")
    ap.add_argument("-i", "--input", required=True, help="Input FASTA")
    ap.add_argument("-o", "--outdir", default="anchor_clusters",
                    help="Output dir (default: anchor_clusters)")
    ap.add_argument("-t", "--threshold", type=float, default=0.6,
                    help="Similarity threshold in [0, 1] (default: 0.6)")
    ap.add_argument("--min-cluster-size", type=int, default=2,
                    help="Min members for per-cluster FASTA (default: 2)")
    ap.add_argument("--n-front", type=int, default=3,
                    help="N-terminal anchor length (default: 3)")
    ap.add_argument("--n-back", type=int, default=3,
                    help="C-terminal anchor length (default: 3)")
    ap.add_argument("--refinement", action="store_true",
                    help="Apply Lloyd-style refinement after greedy "
                         "clustering (off by default)")
    ap.add_argument("--iterations", type=int, default=3,
                    help="Max refinement passes (default: 3; only used "
                         "with --refinement)")
    ap.add_argument("--refine-cap", type=int, default=32,
                    help="Max centroid comparisons per anchor during "
                         "refinement reassignment (default: 32; <=0 = no cap). "
                         "Lower is faster; only used with --refinement")
    ap.add_argument("--no-merge", action="store_true",
                    help="Skip the centroid-merge step during refinement "
                         "(faster; only used with --refinement)")
    ap.add_argument("--backend", choices=["auto", "rust", "python"],
                    default="auto",
                    help="Clustering backend (default: auto — Rust if built, "
                         "else Python)")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="Suppress progress output")
    ap.add_argument("--version", action="version",
                    version=f"pepcluster {__version__}")
    args = ap.parse_args(argv)

    cluster_fasta(
        input=args.input,
        outdir=args.outdir,
        threshold=args.threshold,
        min_cluster_size=args.min_cluster_size,
        n_front=args.n_front,
        n_back=args.n_back,
        refinement=args.refinement,
        iterations=args.iterations,
        refine_cap=args.refine_cap,
        merge=not args.no_merge,
        backend=args.backend,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
