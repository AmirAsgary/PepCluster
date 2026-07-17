"""Tests for pepcluster: pipeline, representatives, configurable anchors, and
exact Rust/Python backend parity."""

import random
from collections import defaultdict
from pathlib import Path

import pytest

import pepcluster
from pepcluster.clustering import (
    parse_fasta,
    parse_anchors,
    extract_anchor,
    anchor_sim_fast,
    build_tables,
    cluster_anchors_py,
    refine_clusters_py,
    cluster_representatives,
    _to_ords,
)

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "peptides.fasta"
THRESHOLD = 0.6
DEFAULT_AP = [1, 5]


def _anchor_counts(n_front=3, n_back=3):
    counts = defaultdict(int)
    for _hdr, seq in parse_fasta(EXAMPLE):
        anc = extract_anchor(seq.upper(), n_front, n_back)
        if anc is not None:
            counts[anc] += 1
    return dict(counts)


def _random_counts(n_unique, seed, alen=6):
    # sorted(), not set-iteration: Python's string hash is randomized per
    # process, so iterating a set would make the data non-reproducible.
    rng = random.Random(seed)
    aa = "ARNDCQEGHILKMFPSTWYV"
    s = set()
    while len(s) < n_unique:
        s.add("".join(rng.choice(aa) for _ in range(alen)))
    return {a: rng.randint(1, 5) for a in sorted(s)}


def test_example_fasta_present():
    assert EXAMPLE.exists(), f"missing example FASTA at {EXAMPLE}"


def test_cluster_anchors_py_basic():
    counts = _anchor_counts()
    assert counts, "no valid anchors extracted"
    mapping, _n_cmp, _n_early = cluster_anchors_py(counts, THRESHOLD)
    assert set(mapping) == set(counts)
    for centroid in set(mapping.values()):
        assert mapping[centroid] == centroid
    assert mapping["SLLAGV"] == "SLLAGV"


# ── --anchors parsing ──────────────────────────────────────────────────────

@pytest.mark.parametrize("spec,expected", [
    ("2;3",   [1, 5]),      # default: P2 and POmega
    ("2;2,3", [1, 4, 5]),
    ("1,2;3", [0, 1, 5]),
    (";3",    [5]),
    ("2;",    [1]),
    (" 2 ; 3 ", [1, 5]),
    ("3,3;1", [2, 3]),      # dedup + sort
])
def test_parse_anchors(spec, expected):
    assert parse_anchors(spec, 3, 3) == expected


def test_parse_anchors_default_matches_p2_pomega():
    assert parse_anchors() == list(pepcluster.DEFAULT_ANCHOR_POSITIONS)


@pytest.mark.parametrize("spec", ["4;3", "2;4", "0;3", "2", "2;3;4", "x;3", ";"])
def test_parse_anchors_rejects_bad_specs(spec):
    with pytest.raises(ValueError):
        parse_anchors(spec, 3, 3)


def test_parse_anchors_respects_n_front_n_back():
    # With a 2+2 anchor, "2;2" = 2nd of front, 2nd of back -> positions 1 and 3.
    assert parse_anchors("2;2", 2, 2) == [1, 3]
    with pytest.raises(ValueError):
        parse_anchors("3;1", 2, 2)  # front index 3 out of range for n_front=2


# ── Anchor weighting actually changes the metric ───────────────────────────

def test_anchor_positions_change_similarity_weighting():
    a, b = _to_ords("AAAAAA"), _to_ords("AAAAAV")  # differ only at position 5
    t_last = build_tables(6, [5], 2.0)   # position 5 is a weighted anchor
    t_first = build_tables(6, [0], 2.0)  # position 5 is an ordinary position
    # Mismatching a heavily-weighted position must cost more.
    assert anchor_sim_fast(a, b, -1.0, t_last) < anchor_sim_fast(a, b, -1.0, t_first)


# ── Rust/Python parity ─────────────────────────────────────────────────────

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


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
@pytest.mark.parametrize("thr", [0.3, 0.6])
@pytest.mark.parametrize("cap,merge", [(32, True), (8, True), (0, True),
                                       (32, False), (0, False)])
def test_rust_python_refine_parity_random(thr, cap, merge):
    """Capped / merge-optional refinement must be bit-identical across backends."""
    counts = _random_counts(3000, seed=1)
    base, _, _ = cluster_anchors_py(counts, thr)
    py_map, py_stats = refine_clusters_py(counts, base, thr, 3, cap, merge)
    rs_map, rs_stats = pepcluster.refine_clusters(counts, base, thr, 3, cap, merge)
    assert rs_map == py_map
    assert dict(rs_stats) == dict(py_stats)


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
@pytest.mark.parametrize("alen,ap", [
    (6, [1, 5]),      # default
    (6, [1, 4, 5]),   # "2;2,3"
    (6, [5]),         # ";3"  — C-terminal anchor only
    (6, [0, 1, 5]),   # "1,2;3"
    (4, [1, 3]),      # n_front=2, n_back=2 -> "2;2"
    (8, [1, 7]),      # n_front=4, n_back=4 -> "2;4"
])
def test_rust_python_parity_custom_anchors(alen, ap):
    """Greedy + refinement stay bit-identical for any anchor length/positions."""
    counts = _random_counts(2000, seed=5, alen=alen)
    py_map, py_c, py_e = cluster_anchors_py(counts, 0.5, ap, 2.0)
    rs_map, rs_c, rs_e = pepcluster.cluster_anchors(counts, 0.5, ap, 2.0)
    assert rs_map == py_map and (rs_c, rs_e) == (py_c, py_e)

    py_r, py_s = refine_clusters_py(counts, py_map, 0.5, 2, 32, True, ap, 2.0)
    rs_r, rs_s = pepcluster.refine_clusters(counts, rs_map, 0.5, 2, 32, True, ap, 2.0)
    assert rs_r == py_r
    assert dict(rs_s) == dict(py_s)


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
def test_anchor_weight_is_honoured_by_both_backends():
    counts = _random_counts(1500, seed=9)
    for w in (1.0, 3.0):
        py_map, _, _ = cluster_anchors_py(counts, 0.5, DEFAULT_AP, w)
        rs_map, _, _ = pepcluster.cluster_anchors(counts, 0.5, DEFAULT_AP, w)
        assert rs_map == py_map
    # A different weight must actually change the clustering.
    m1, _, _ = pepcluster.cluster_anchors(counts, 0.5, DEFAULT_AP, 1.0)
    m2, _, _ = pepcluster.cluster_anchors(counts, 0.5, DEFAULT_AP, 4.0)
    assert m1 != m2


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
def test_rust_rejects_bad_anchor_positions():
    counts = _random_counts(50, seed=3)
    with pytest.raises(ValueError):
        pepcluster.cluster_anchors(counts, 0.6, [6], 2.0)   # out of range for 6-mer
    with pytest.raises(ValueError):
        pepcluster.cluster_anchors(counts, 0.6, [], 2.0)    # empty


def test_no_merge_and_cap_semantics():
    counts = _random_counts(2000, seed=7)
    base, _, _ = cluster_anchors_py(counts, 0.5)

    _m, s_merge = refine_clusters_py(counts, base, 0.5, 3, 32, True)
    _m, s_nomerge = refine_clusters_py(counts, base, 0.5, 3, 32, False)
    assert s_nomerge["merges"] == 0
    assert s_nomerge["final_clusters"] >= s_merge["final_clusters"]

    m_cap, _ = refine_clusters_py(counts, base, 0.5, 1, 4, True)
    assert set(m_cap) == set(counts)


# ── fast_medoid (O(N)) + merge_cap ─────────────────────────────────────────

@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
@pytest.mark.parametrize("thr", [0.3, 0.6])
@pytest.mark.parametrize("fast_medoid", [False, True])
@pytest.mark.parametrize("merge_cap", [0, 8, 32])
def test_rust_python_parity_fast_medoid_and_merge_cap(thr, fast_medoid, merge_cap):
    """New flags must stay bit-identical across the Rust and Python backends."""
    counts = _random_counts(3000, seed=2)
    base, _, _ = cluster_anchors_py(counts, thr)
    py_map, py_s = refine_clusters_py(counts, base, thr, 3, 32, True,
                                      DEFAULT_AP, 2.0, fast_medoid, merge_cap)
    rs_map, rs_s = pepcluster.refine_clusters(counts, base, thr, 3, 32, True,
                                              DEFAULT_AP, 2.0, fast_medoid, merge_cap)
    assert rs_map == py_map
    assert dict(rs_s) == dict(py_s)


@pytest.mark.skipif(not pepcluster.HAS_RUST, reason="Rust backend not built")
@pytest.mark.parametrize("alen,ap", [(6, [1, 5]), (4, [1, 3]), (8, [1, 7])])
def test_fast_medoid_merge_cap_parity_custom_anchors(alen, ap):
    counts = _random_counts(1500, seed=6, alen=alen)
    base, _, _ = cluster_anchors_py(counts, 0.4, ap, 2.0)
    py_map, py_s = refine_clusters_py(counts, base, 0.4, 2, 32, True, ap, 2.0,
                                      True, 16)
    rs_map, rs_s = pepcluster.refine_clusters(counts, base, 0.4, 2, 32, True, ap,
                                              2.0, True, 16)
    assert rs_map == py_map and dict(rs_s) == dict(py_s)


def test_merge_cap_reduces_or_equals_merges():
    """A merge cap can only find fewer-or-equal merges than the uncapped scan."""
    counts = _random_counts(2500, seed=4)
    base, _, _ = cluster_anchors_py(counts, 0.35)  # low thr -> mergeable clusters
    _m, s_uncapped = refine_clusters_py(counts, base, 0.35, 1, 32, True,
                                        DEFAULT_AP, 2.0, False, 0)
    _m, s_capped = refine_clusters_py(counts, base, 0.35, 1, 32, True,
                                      DEFAULT_AP, 2.0, False, 4)
    assert s_capped["merges"] <= s_uncapped["merges"]
    assert s_capped["final_clusters"] >= s_uncapped["final_clusters"]


def test_fast_medoid_matches_exact_when_unambiguous():
    """On a hand-built cluster with a clear centre, fast and exact medoids agree."""
    # B is central: closer to both A and C than they are to each other.
    counts = {"YLLAGV": 10, "YMLAGV": 6, "YVLAGV": 6}
    base = {a: "YLLAGV" for a in counts}
    m_exact, _ = refine_clusters_py(counts, base, 0.6, 1, 0, False,
                                    DEFAULT_AP, 2.0, False, 0)
    m_fast, _ = refine_clusters_py(counts, base, 0.6, 1, 0, False,
                                   DEFAULT_AP, 2.0, True, 0)
    assert set(m_exact.values()) == set(m_fast.values())


def test_cluster_fasta_fast_medoid_and_merge_cap(tmp_path):
    stats = pepcluster.cluster_fasta(
        str(EXAMPLE), str(tmp_path), threshold=0.6, refinement=True,
        fast_medoid=True, merge_cap=16, verbose=False)
    assert stats["n_clusters"] >= 1
    assert (tmp_path / "clusters.tsv").exists()


# ── Cluster representatives ────────────────────────────────────────────────

def _brute_force_best_sigma(counts, mapping, tables):
    """Reference O(k^2) peptide-level medoid score per cluster."""
    groups = defaultdict(list)
    for a, c in mapping.items():
        groups[c].append(a)
    best = {}
    for centroid, anchors in groups.items():
        ords = {a: _to_ords(a) for a in anchors}
        best[centroid] = max(
            sum(counts[b] * anchor_sim_fast(ords[a], ords[b], -1.0, tables)
                for b in anchors)
            for a in anchors)
    return best


def test_cluster_representatives_are_optimal():
    counts = _anchor_counts()
    mapping, _, _ = cluster_anchors_py(counts, THRESHOLD)
    reps = cluster_representatives(counts, mapping)
    tables = build_tables(6, DEFAULT_AP, 2.0)
    best_sigma = _brute_force_best_sigma(counts, mapping, tables)

    members = defaultdict(set)
    for a, c in mapping.items():
        members[c].add(a)
    for centroid, rep in reps.items():
        assert rep in members[centroid]

    ords = {a: _to_ords(a) for a in counts}
    for centroid, rep in reps.items():
        sigma_rep = sum(counts[b] * anchor_sim_fast(ords[rep], ords[b], -1.0,
                                                    tables)
                        for b in members[centroid])
        assert sigma_rep == pytest.approx(best_sigma[centroid], rel=1e-5)


def test_cluster_representatives_optimal_with_custom_anchors():
    ap, w = [0, 1, 5], 3.0
    counts = _random_counts(400, seed=11)
    mapping, _, _ = cluster_anchors_py(counts, 0.5, ap, w)
    reps = cluster_representatives(counts, mapping, ap, w)
    tables = build_tables(6, ap, w)
    best_sigma = _brute_force_best_sigma(counts, mapping, tables)

    members = defaultdict(set)
    for a, c in mapping.items():
        members[c].add(a)
    ords = {a: _to_ords(a) for a in counts}
    for centroid, rep in reps.items():
        assert rep in members[centroid]
        sigma_rep = sum(counts[b] * anchor_sim_fast(ords[rep], ords[b], -1.0,
                                                    tables)
                        for b in members[centroid])
        assert sigma_rep == pytest.approx(best_sigma[centroid], rel=1e-5)


# ── End-to-end pipeline ────────────────────────────────────────────────────

def test_cluster_fasta_writes_outputs(tmp_path):
    stats = pepcluster.cluster_fasta(
        str(EXAMPLE), str(tmp_path), threshold=THRESHOLD, verbose=False)
    assert (tmp_path / "clusters.tsv").exists()
    assert (tmp_path / "cluster_summary.tsv").exists()
    assert (tmp_path / "summary.txt").exists()
    assert stats["too_short"] == 1
    assert stats["valid_peptides"] == 18
    assert stats["n_clusters"] >= 1
    assert stats["anchor_positions"] == [1, 5]


def test_output_has_representative_peptide_column(tmp_path):
    pepcluster.cluster_fasta(str(EXAMPLE), str(tmp_path),
                             threshold=THRESHOLD, verbose=False)
    header = (tmp_path / "clusters.tsv").read_text().splitlines()[0].split("\t")
    assert "representative_peptide" in header
    summ = (tmp_path / "cluster_summary.tsv").read_text().splitlines()
    assert "representative_peptide" in summ[0].split("\t")

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


def test_cluster_fasta_custom_anchors(tmp_path):
    """--anchors selects different positions and changes the clustering."""
    a = pepcluster.cluster_fasta(str(EXAMPLE), str(tmp_path / "a"),
                                 threshold=0.7, anchors="2;3", verbose=False)
    b = pepcluster.cluster_fasta(str(EXAMPLE), str(tmp_path / "b"),
                                 threshold=0.7, anchors="1,2;2,3", verbose=False)
    assert a["anchor_positions"] == [1, 5]
    assert b["anchor_positions"] == [0, 1, 4, 5]
    assert a["n_clusters"] != b["n_clusters"]


def test_cluster_fasta_custom_n_front_n_back(tmp_path):
    """Regression: non-3/3 anchor lengths used to crash with a KeyError."""
    stats = pepcluster.cluster_fasta(
        str(EXAMPLE), str(tmp_path), threshold=0.6,
        n_front=2, n_back=2, anchors="2;2", verbose=False)
    assert stats["anchor_positions"] == [1, 3]
    assert stats["n_clusters"] >= 1
    anchors = [line.split("\t")[-1]
               for line in (tmp_path / "clusters.tsv").read_text()
               .splitlines()[1:]]
    assert anchors and all(len(a) == 4 for a in anchors)


def test_cluster_fasta_rejects_bad_anchors(tmp_path):
    with pytest.raises(ValueError):
        pepcluster.cluster_fasta(str(EXAMPLE), str(tmp_path),
                                 anchors="4;3", verbose=False)
