"""Tests for pepcluster: end-to-end pipeline + Rust/Python backend parity."""

from collections import defaultdict
from pathlib import Path

import pytest

import pepcluster
from pepcluster.clustering import (
    parse_fasta,
    extract_anchor,
    anchor_sim_fast,
    cluster_anchors_py,
    refine_clusters_py,
    cluster_representatives,
    _to_ords,
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


def _random_counts(n_unique, seed):
    import random
    rng = random.Random(seed)
    aa = "ARNDCQEGHILKMFPSTWYV"
    s = set()
    while len(s) < n_unique:
        s.add("".join(rng.choice(aa) for _ in range(6)))
    return {a: rng.randint(1, 5) for a in s}


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
@pytest.mark.parametrize("thr", [0.3, 0.6])
@pytest.mark.parametrize("cap,merge", [(32, True), (8, True), (0, True),
                                       (32, False), (0, False)])
def test_rust_python_refine_parity_random(thr, cap, merge):
    """The capped/merge-optional refinement must be bit-identical across the
    Rust and Python backends on random data (needs matching candidate order)."""
    counts = _random_counts(3000, seed=1)
    base, _, _ = cluster_anchors_py(counts, thr)
    py_map, py_stats = refine_clusters_py(counts, base, thr, 3, cap, merge)
    rs_map, rs_stats = pepcluster.refine_clusters(counts, base, thr, 3, cap, merge)
    assert rs_map == py_map
    assert dict(rs_stats) == dict(py_stats)


def test_no_merge_and_cap_semantics():
    """--no-merge yields zero merges (and >= clusters); a tiny cap still runs."""
    counts = _random_counts(2000, seed=7)
    base, _, _ = cluster_anchors_py(counts, 0.5)

    _m, s_merge = refine_clusters_py(counts, base, 0.5, 3, 32, True)
    _m, s_nomerge = refine_clusters_py(counts, base, 0.5, 3, 32, False)
    assert s_nomerge["merges"] == 0
    assert s_nomerge["final_clusters"] >= s_merge["final_clusters"]

    # An unbounded cap examines at least as much as a tiny cap (never fewer
    # reassignments than the capped run misses out on is not guaranteed, but
    # both must complete and stay valid mappings).
    m_cap, _ = refine_clusters_py(counts, base, 0.5, 1, 4, True)
    assert set(m_cap) == set(counts)


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


# ── Cluster representatives ────────────────────────────────────────────────

def _brute_force_medoid(counts, mapping):
    """Reference O(k^2) peptide-level medoid per cluster: for each member the
    total frequency-weighted similarity to the whole cluster, argmax = medoid.
    Returns {centroid: max_sigma} so we can compare optimal scores directly."""
    groups = defaultdict(list)
    for a, c in mapping.items():
        groups[c].append(a)
    best_sigma = {}
    for centroid, anchors in groups.items():
        ords = {a: _to_ords(a) for a in anchors}
        best = None
        for a in anchors:
            sigma = sum(counts[b] * anchor_sim_fast(ords[a], ords[b], -1.0)
                        for b in anchors)
            best = sigma if best is None else max(best, sigma)
        best_sigma[centroid] = best
    return best_sigma


def test_cluster_representatives_are_optimal():
    counts = _anchor_counts()
    mapping, _, _ = cluster_anchors_py(counts, THRESHOLD)
    reps = cluster_representatives(counts, mapping)
    best_sigma = _brute_force_medoid(counts, mapping)

    # Every cluster has a representative that is one of its own members.
    members = defaultdict(set)
    for a, c in mapping.items():
        members[c].add(a)
    for centroid, rep in reps.items():
        assert rep in members[centroid]

    # The fast decomposition must pick a genuine medoid: the representative's
    # total similarity equals the brute-force maximum.
    ords = {a: _to_ords(a) for a in counts}
    for centroid, rep in reps.items():
        sigma_rep = sum(counts[b] * anchor_sim_fast(ords[rep], ords[b], -1.0)
                        for b in members[centroid])
        assert sigma_rep == pytest.approx(best_sigma[centroid], rel=1e-5)


def test_output_has_representative_peptide_column(tmp_path):
    pepcluster.cluster_fasta(str(EXAMPLE), str(tmp_path),
                             threshold=THRESHOLD, verbose=False)
    header = (tmp_path / "clusters.tsv").read_text().splitlines()[0].split("\t")
    assert "representative_peptide" in header
    summ = (tmp_path / "cluster_summary.tsv").read_text().splitlines()
    assert "representative_peptide" in summ[0].split("\t")

    # Every representative_peptide must be a real member sequence of its cluster.
    rep_col = header.index("representative_peptide")
    seq_col = header.index("sequence")
    cid_col = header.index("cluster_id")
    seqs_by_cluster = defaultdict(set)
    reps_by_cluster = {}
    for line in (tmp_path / "clusters.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        seqs_by_cluster[f[cid_col]].add(f[seq_col])
        reps_by_cluster[f[cid_col]] = f[rep_col]
    for cid, rep in reps_by_cluster.items():
        assert rep in seqs_by_cluster[cid]
