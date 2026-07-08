"""Tests for pepcluster: end-to-end pipeline + Rust/Python backend parity."""

from collections import defaultdict
from pathlib import Path

import pytest

import pepcluster
from pepcluster.clustering import (
    parse_fasta,
    extract_anchor,
    cluster_anchors_py,
    refine_clusters_py,
)

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "peptides.fasta"
THRESHOLD = 0.6


def _anchor_counts():
    counts = defaultdict(int)
    for _hdr, seq in parse_fasta(EXAMPLE):
        anc = extract_anchor(seq.upper())
        if anc is not None:
            counts[anc] += 1
    return dict(counts)


def test_example_fasta_present():
    assert EXAMPLE.exists(), f"missing example FASTA at {EXAMPLE}"


def test_cluster_anchors_py_basic():
    counts = _anchor_counts()
    assert counts, "no valid anchors extracted"
    mapping, n_cmp, n_early = cluster_anchors_py(counts, THRESHOLD)
    # Every anchor is mapped, and centroids are themselves.
    assert set(mapping) == set(counts)
    for centroid in set(mapping.values()):
        assert mapping[centroid] == centroid
    # The three SLL...AGV peptides share one anchor → collapse to one centroid.
    assert mapping["SLLAGV"] == "SLLAGV"


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
def test_rust_python_cluster_parity():
    counts = _anchor_counts()
    py_map, _, _ = cluster_anchors_py(counts, THRESHOLD)
    rs_map, _, _ = pepcluster.cluster_anchors(counts, THRESHOLD)
    assert rs_map == py_map


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
def test_rust_python_refine_parity():
    counts = _anchor_counts()
    base_map, _, _ = cluster_anchors_py(counts, THRESHOLD)
    py_map, py_stats = refine_clusters_py(counts, base_map, THRESHOLD, 3)
    rs_map, rs_stats = pepcluster.refine_clusters(counts, base_map, THRESHOLD, 3)
    assert rs_map == py_map
    assert rs_stats["final_clusters"] == py_stats["final_clusters"]


def test_cluster_fasta_writes_outputs(tmp_path):
    stats = pepcluster.cluster_fasta(
        str(EXAMPLE), str(tmp_path), threshold=THRESHOLD, verbose=False)
    assert (tmp_path / "clusters.tsv").exists()
    assert (tmp_path / "cluster_summary.tsv").exists()
    assert (tmp_path / "summary.txt").exists()
    # One peptide ("ACDE") is too short for a 6-mer anchor.
    assert stats["too_short"] == 1
    assert stats["valid_peptides"] == 18
    assert stats["n_clusters"] >= 1
