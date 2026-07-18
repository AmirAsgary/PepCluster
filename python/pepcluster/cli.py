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
    ap.add_argument("--anchors", default="2;3",
                    help="Which positions are binding anchors, as 'FRONT;BACK' "
                         "with 1-based indices into each side (default: '2;3' "
                         "= 2nd of the first 3 residues (P2) and 3rd of the "
                         "last 3 (POmega)). Each side takes a comma-separated "
                         "list and may be empty, e.g. '2;2,3', '1,2;3', ';3'. "
                         "Anchor positions get --anchor-weight and define the "
                         "blocking")
    ap.add_argument("--anchor-weight", type=float, default=2.0,
                    help="Weight given to anchor positions; all other "
                         "positions have weight 1.0 (default: 2.0)")
    ap.add_argument("--obg-block-search", action="store_true",
                    help="Upper-Bound-Guided multi-probe block search: also "
                         "compare each anchor against centroids in neighbouring "
                         "blocks whose upper bound reaches the threshold, not "
                         "just its own block. Gives fewer, tighter clusters at "
                         "some extra cost (default: off = same-block only)")
    ap.add_argument("--obg-max-probes", type=int, default=0,
                    help="With --obg-block-search, search at most this many "
                         "blocks per anchor (including its own block); "
                         "0 = unlimited (default: 0)")
    ap.add_argument("--obg-min-block-upper-bound", type=float, default=0.0,
                    help="With --obg-block-search, only search blocks whose "
                         "upper bound is at least this; the effective cut is "
                         "max(threshold, this) (default: 0.0)")
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
    ap.add_argument("--fast-medoid", action="store_true",
                    help="Use the O(N) medoid decomposition instead of the "
                         "exact O(k^2) medoid during refinement. Much faster "
                         "when a few clusters are very large (e.g. low "
                         "thresholds). Only used with --refinement")
    ap.add_argument("--merge-cap", type=int, default=0,
                    help="Max candidate centroids examined per centroid in the "
                         "refinement merge step (default: 0 = no cap). A value "
                         "like 32 makes merge fast when there are many clusters "
                         "(e.g. high thresholds). Only used with --refinement")
    ap.add_argument("--central-region-profiling", action="store_true",
                    help="Build a central-region k-mer profile per cluster and "
                         "(unless --no-cluster-profile-merge) merge clusters by "
                         "a combined anchor + central-profile score. Needs "
                         "--refinement to take effect")
    ap.add_argument("--crr-kmer-size", type=int, default=2,
                    help="Central-region k-mer length (default: 2)")
    ap.add_argument("--crr-bins", type=int, default=3,
                    help="Number of relative-position bins (default: 3)")
    ap.add_argument("--crr-adjacent-bin-smoothing", type=float, default=0.5,
                    help="Adjacent-bin smoothing weight (default: 0.5)")
    ap.add_argument("--no-cluster-profile-merge", action="store_true",
                    help="With --central-region-profiling, keep the anchor-only "
                         "merge instead of the combined-score merge")
    ap.add_argument("--cluster-profile-merge-weight", type=float, default=0.2,
                    help="Weight of the central-profile term in the merge score "
                         "(default: 0.2)")
    ap.add_argument("--cluster-profile-merge-threshold", type=float, default=0.6,
                    help="Combined-score threshold to merge two clusters "
                         "(default: 0.6)")
    ap.add_argument("--threads", type=int, default=1,
                    help="Worker threads for the Rust backend's greedy "
                         "clustering and refinement (1 = serial, default; "
                         "0 = all cores; N = exactly N). Results are identical "
                         "regardless of thread count. The Python backend is "
                         "always serial")
    ap.add_argument("--backend", choices=["auto", "rust", "python"],
                    default="auto",
                    help="Clustering backend (default: auto — Rust if built, "
                         "else Python)")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="Suppress progress output")
    ap.add_argument("--version", action="version",
                    version=f"pepcluster {__version__}")
    args = ap.parse_args(argv)

    try:
        cluster_fasta(
            input=args.input,
            outdir=args.outdir,
            threshold=args.threshold,
            min_cluster_size=args.min_cluster_size,
            n_front=args.n_front,
            n_back=args.n_back,
            anchors=args.anchors,
            anchor_weight=args.anchor_weight,
            obg_block_search=args.obg_block_search,
            obg_max_probes=args.obg_max_probes,
            obg_min_block_upper_bound=args.obg_min_block_upper_bound,
            refinement=args.refinement,
            iterations=args.iterations,
            refine_cap=args.refine_cap,
            merge=not args.no_merge,
            fast_medoid=args.fast_medoid,
            merge_cap=args.merge_cap,
            central_region_profiling=args.central_region_profiling,
            crr_kmer_size=args.crr_kmer_size,
            crr_bins=args.crr_bins,
            crr_smoothing=args.crr_adjacent_bin_smoothing,
            cluster_profile_merge=not args.no_cluster_profile_merge,
            merge_weight=args.cluster_profile_merge_weight,
            merge_threshold=args.cluster_profile_merge_threshold,
            threads=args.threads,
            backend=args.backend,
            verbose=not args.quiet,
        )
    except ValueError as exc:
        ap.error(str(exc))


if __name__ == "__main__":
    main()
