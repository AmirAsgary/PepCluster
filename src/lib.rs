//! pepcluster._core — Rust backend for BLOSUM62-aware anchor clustering.
//!
//! Exposes two functions to Python:
//!   * `cluster_anchors(anchor_counts, threshold)` — greedy centroid clustering
//!   * `refine_clusters(anchor_counts, mapping, threshold, iterations)` — Lloyd-style refinement
//!
//! Both operate on the 6-residue anchor (peptide first 3 + last 3 amino acids)
//! using a BLOSUM62-normalized similarity with P2 and PΩ weighted 2×.

use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::{HashMap, HashSet};

// ============================================================================
// BLOSUM62 matrix (20×20, standard order ARNDCQEGHILKMFPSTWYV)
// ============================================================================
const AA_ORDER: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

#[rustfmt::skip]
const BLOSUM62: [i8; 400] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
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
];

// ============================================================================
// Reduced alphabet (10 groups) for blocking
// ============================================================================
const COARSE_GROUPS: &[&[u8]] = &[
    b"AST", b"VILM", b"FYW", b"DE", b"KR", b"NQ",
    b"G", b"H", b"C", b"P",
];

// Position weights: P1=1, P2=2, P3=1, PΩ-2=1, PΩ-1=1, PΩ=2  (sum=8)
// Check order: [1, 5, 0, 2, 3, 4] — high-weight positions first for early exit
const CHECK_ORDER: [usize; 6] = [1, 5, 0, 2, 3, 4];
const RAW_WEIGHTS: [f32; 6] = [1.0, 2.0, 1.0, 1.0, 1.0, 2.0];
const WEIGHT_SUM: f32 = 8.0;

/// Precomputed similarity tables and blocking tables
struct SimTables {
    /// wsim[pos][a * 128 + b] = normalized BLOSUM62 * weight for position `pos`
    /// Positions stored in CHECK_ORDER: [1, 5, 0, 2, 3, 4]
    wsim: [[f32; 128 * 128]; 6],
    /// remaining_after[i] = max possible score from positions CHECK_ORDER[i+1..]
    remaining_after: [f32; 7],
    /// coarse group id for each ASCII byte (0..9), 255 = unknown
    coarse: [u8; 128],
}

impl SimTables {
    fn new() -> Self {
        // Build AA → index mapping
        let mut aa_idx = [255u8; 128];
        for (i, &aa) in AA_ORDER.iter().enumerate() {
            aa_idx[aa as usize] = i as u8;
        }

        // Self-scores for normalization
        let mut self_score = [0.0f32; 20];
        for i in 0..20 {
            self_score[i] = BLOSUM62[i * 20 + i] as f32;
        }

        // Normalized similarity: sim(a,b) = B(a,b) / sqrt(B(a,a) * B(b,b))
        let mut norm_sim = [[0.0f32; 20]; 20];
        for i in 0..20 {
            for j in 0..20 {
                let denom = (self_score[i] * self_score[j]).sqrt();
                if denom > 0.0 {
                    norm_sim[i][j] = BLOSUM62[i * 20 + j] as f32 / denom;
                }
            }
        }

        // Build weighted similarity tables in CHECK_ORDER
        let mut wsim = [[0.0f32; 128 * 128]; 6];
        for (check_idx, &pos) in CHECK_ORDER.iter().enumerate() {
            let w = RAW_WEIGHTS[pos] / WEIGHT_SUM;
            for &a in AA_ORDER.iter() {
                let ai = aa_idx[a as usize] as usize;
                for &b in AA_ORDER.iter() {
                    let bi = aa_idx[b as usize] as usize;
                    wsim[check_idx][(a as usize) * 128 + (b as usize)] =
                        w * norm_sim[ai][bi];
                }
            }
        }

        // remaining_after[i] = sum of weights for CHECK_ORDER positions i+1..6
        let mut remaining_after = [0.0f32; 7];
        for i in (0..6).rev() {
            let pos = CHECK_ORDER[i];
            remaining_after[i] = remaining_after[i + 1] + RAW_WEIGHTS[pos] / WEIGHT_SUM;
        }

        // Coarse group lookup
        let mut coarse = [255u8; 128];
        for (gid, &group) in COARSE_GROUPS.iter().enumerate() {
            for &aa in group {
                coarse[aa as usize] = gid as u8;
            }
        }

        SimTables { wsim, remaining_after, coarse }
    }

    /// Compute weighted BLOSUM62-normalized similarity with early termination.
    /// Returns the score if >= threshold, otherwise returns -1.0.
    #[inline]
    fn anchor_sim(&self, a: &[u8; 6], b: &[u8; 6], threshold: f32) -> f32 {
        let mut s: f32 = 0.0;

        // Check position 1 (P2, weight 2×)
        s += self.wsim[0][(a[1] as usize) * 128 + (b[1] as usize)];
        if s + self.remaining_after[1] < threshold { return -1.0; }

        // Check position 5 (PΩ, weight 2×)
        s += self.wsim[1][(a[5] as usize) * 128 + (b[5] as usize)];
        if s + self.remaining_after[2] < threshold { return -1.0; }

        // Check position 0 (P1)
        s += self.wsim[2][(a[0] as usize) * 128 + (b[0] as usize)];
        if s + self.remaining_after[3] < threshold { return -1.0; }

        // Remaining positions — no early exit (cost not worth it)
        s += self.wsim[3][(a[2] as usize) * 128 + (b[2] as usize)];
        s += self.wsim[4][(a[3] as usize) * 128 + (b[3] as usize)];
        s += self.wsim[5][(a[4] as usize) * 128 + (b[4] as usize)];

        s
    }

    /// Full similarity computation without early termination.
    /// Used by the refinement pass where we need the actual score (e.g.
    /// for medoid scoring and for comparing against the current centroid).
    #[inline]
    fn anchor_sim_full(&self, a: &[u8; 6], b: &[u8; 6]) -> f32 {
        let mut s: f32 = 0.0;
        s += self.wsim[0][(a[1] as usize) * 128 + (b[1] as usize)];
        s += self.wsim[1][(a[5] as usize) * 128 + (b[5] as usize)];
        s += self.wsim[2][(a[0] as usize) * 128 + (b[0] as usize)];
        s += self.wsim[3][(a[2] as usize) * 128 + (b[2] as usize)];
        s += self.wsim[4][(a[3] as usize) * 128 + (b[3] as usize)];
        s += self.wsim[5][(a[4] as usize) * 128 + (b[4] as usize)];
        s
    }

    #[inline]
    fn block_key(&self, a: &[u8; 6]) -> (u8, u8) {
        (self.coarse[a[1] as usize], self.coarse[a[5] as usize])
    }
}

/// Convert a 6-char string slice to a fixed-size byte array
#[inline]
fn to_anchor(s: &str) -> Option<[u8; 6]> {
    let bytes = s.as_bytes();
    if bytes.len() != 6 { return None; }
    Some([bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5]])
}

/// Core clustering function exposed to Python.
///
/// Args:
///     anchor_counts: dict[str, int] — unique 6-mer anchors with their frequency
///     threshold: float — minimum similarity to join a cluster
///
/// Returns:
///     (mapping, n_comparisons, n_early_exits)
///     mapping: dict[str, str] — anchor → centroid anchor
#[pyfunction]
fn cluster_anchors(
    py: Python<'_>,
    anchor_counts: &Bound<'_, PyDict>,
    threshold: f64,
) -> PyResult<(PyObject, u64, u64)> {
    let threshold = threshold as f32;
    let tables = SimTables::new();

    // Collect and parse anchors, sort by frequency descending
    let mut items: Vec<(String, [u8; 6], u64)> = Vec::new();
    for (key, val) in anchor_counts.iter() {
        let anchor_str: String = key.extract()?;
        let count: u64 = val.extract()?;
        if let Some(arr) = to_anchor(&anchor_str) {
            items.push((anchor_str, arr, count));
        }
    }
    items.sort_by(|a, b| b.2.cmp(&a.2));

    // Partition into blocks by coarse(P2) + coarse(PΩ)
    let mut blocks: HashMap<(u8, u8), Vec<usize>> = HashMap::new();
    for (idx, (_, arr, _)) in items.iter().enumerate() {
        let bk = tables.block_key(arr);
        blocks.entry(bk).or_default().push(idx);
    }

    // Greedy centroid clustering within each block
    let mut centroid_of: Vec<usize> = vec![0; items.len()]; // index → centroid index
    let mut n_cmp: u64 = 0;
    let mut n_early: u64 = 0;

    for (_, members) in blocks.iter() {
        let mut centroids: Vec<usize> = Vec::new();

        for &idx in members.iter() {
            let anchor = &items[idx].1;
            let mut matched = false;

            for &ci in centroids.iter() {
                n_cmp += 1;
                let score = tables.anchor_sim(anchor, &items[ci].1, threshold);
                if score < 0.0 {
                    n_early += 1;
                    continue;
                }
                if score >= threshold {
                    centroid_of[idx] = ci;
                    matched = true;
                    break;
                }
            }

            if !matched {
                centroids.push(idx);
                centroid_of[idx] = idx;
            }
        }
    }

    // Build result dict: anchor_str → centroid_str
    let result = PyDict::new(py);
    for (idx, (anchor_str, _, _)) in items.iter().enumerate() {
        let ci = centroid_of[idx];
        result.set_item(anchor_str, &items[ci].0)?;
    }

    Ok((result.into(), n_cmp, n_early))
}

/// Python module definition
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cluster_anchors, m)?)?;
    m.add_function(wrap_pyfunction!(refine_clusters, m)?)?;
    Ok(())
}

// ============================================================================
// Optional refinement (Lloyd-style; mirrors the Python refine_clusters)
// ============================================================================

/// Refinement pass on top of greedy clustering output.
///
/// Each pass performs three sub-steps:
///   1. Medoid update     — replace each centroid with the member that
///                          maximises frequency-weighted mean similarity
///                          to the cluster's other members.
///   2. Cross-block reassign — for each anchor, find the best centroid
///                          above `threshold` across its own block plus
///                          neighbouring blocks (blocks sharing at least
///                          one coarse-alphabet key at P2 or PΩ). Move
///                          if strictly more similar than current.
///   3. Centroid merge    — if two centroids satisfy sim >= threshold,
///                          absorb the smaller cluster into the larger.
///
/// Stops early when no change occurs in a pass.
///
/// Args:
///     anchor_counts: dict[str, int]  — unique anchor → frequency
///     mapping:       dict[str, str]  — unique anchor → centroid
///     threshold:     float           — same as initial clustering
///     iterations:    int             — max passes
///
/// Returns:
///     (refined_mapping, stats_dict)
///     stats_dict keys: passes, medoid_changes, reassignments, merges,
///                       initial_clusters, final_clusters
#[pyfunction]
fn refine_clusters(
    py: Python<'_>,
    anchor_counts: &Bound<'_, PyDict>,
    mapping: &Bound<'_, PyDict>,
    threshold: f64,
    iterations: usize,
) -> PyResult<(PyObject, PyObject)> {
    let threshold = threshold as f32;
    let tables = SimTables::new();

    // ── Parse anchor_counts into parallel vectors ─────────────────────
    let mut names: Vec<String> = Vec::new();
    let mut anchors: Vec<[u8; 6]> = Vec::new();
    let mut freqs: Vec<u64> = Vec::new();
    let mut str_to_idx: HashMap<String, usize> = HashMap::new();

    for (key, val) in anchor_counts.iter() {
        let s: String = key.extract()?;
        let cnt: u64 = val.extract()?;
        if let Some(arr) = to_anchor(&s) {
            str_to_idx.insert(s.clone(), names.len());
            names.push(s);
            anchors.push(arr);
            freqs.push(cnt);
        }
    }
    let n = names.len();

    // ── Parse mapping into cur_centroid[i] ────────────────────────────
    // Default: every anchor is its own centroid (only matters if mapping is
    // missing an entry, which shouldn't happen in normal use).
    let mut cur_centroid: Vec<usize> = (0..n).collect();
    for (key, val) in mapping.iter() {
        let a: String = key.extract()?;
        let c: String = val.extract()?;
        if let (Some(&ai), Some(&ci)) = (str_to_idx.get(&a), str_to_idx.get(&c)) {
            cur_centroid[ai] = ci;
        }
    }

    let initial_clusters: u64 = {
        let s: HashSet<usize> = cur_centroid.iter().copied().collect();
        s.len() as u64
    };

    let mut total_medoid_changes: u64 = 0;
    let mut total_reassignments:  u64 = 0;
    let mut total_merges:         u64 = 0;
    let mut passes_run:           u64 = 0;

    for pass_idx in 0..iterations {
        passes_run = (pass_idx + 1) as u64;
        let mut changed = false;

        // ── 1. Medoid update ──────────────────────────────────────────
        let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..n {
            members.entry(cur_centroid[i]).or_default().push(i);
        }

        let mut centroid_remap: HashMap<usize, usize> = HashMap::with_capacity(members.len());
        for (&centroid, mems) in members.iter() {
            if mems.len() == 1 {
                centroid_remap.insert(centroid, centroid);
                continue;
            }
            let mut best_member = centroid;
            let mut best_score: f32 = f32::NEG_INFINITY;
            for &i in mems.iter() {
                let mut s: f32 = 0.0;
                let mut w: f32 = 0.0;
                for &j in mems.iter() {
                    if i == j { continue; }
                    let sim = tables.anchor_sim_full(&anchors[i], &anchors[j]);
                    let fj  = freqs[j] as f32;
                    s += fj * sim;
                    w += fj;
                }
                let avg = if w > 0.0 { s / w } else { 0.0 };
                if avg > best_score {
                    best_score  = avg;
                    best_member = i;
                }
            }
            centroid_remap.insert(centroid, best_member);
            if best_member != centroid {
                total_medoid_changes += 1;
                changed = true;
            }
        }
        for i in 0..n {
            cur_centroid[i] = centroid_remap[&cur_centroid[i]];
        }

        // ── 2. Cross-block reassignment ───────────────────────────────
        let current_centroids: HashSet<usize> =
            cur_centroid.iter().copied().collect();
        let mut centroids_in_block: HashMap<(u8, u8), Vec<usize>> = HashMap::new();
        for &c in current_centroids.iter() {
            centroids_in_block
                .entry(tables.block_key(&anchors[c]))
                .or_default()
                .push(c);
        }

        // Precompute neighbour-block lists
        let all_blocks: Vec<(u8, u8)> = centroids_in_block.keys().copied().collect();
        let mut neighbours_of: HashMap<(u8, u8), Vec<(u8, u8)>> = HashMap::new();
        for &bk in all_blocks.iter() {
            let (p, q) = bk;
            let nb: Vec<(u8, u8)> = all_blocks.iter()
                .filter(|b| b.0 == p || b.1 == q)
                .copied()
                .collect();
            neighbours_of.insert(bk, nb);
        }

        let fallback_nb = vec![]; // empty fallback if anchor's block has no centroids
        let mut pass_reassigns: u64 = 0;
        for i in 0..n {
            let bk    = tables.block_key(&anchors[i]);
            let cur_c = cur_centroid[i];
            let cur_score =
                tables.anchor_sim_full(&anchors[i], &anchors[cur_c]);
            let mut best_c     = cur_c;
            let mut best_score = cur_score;

            let nbs = neighbours_of.get(&bk).unwrap_or(&fallback_nb);
            for &nb in nbs.iter() {
                if let Some(centroids) = centroids_in_block.get(&nb) {
                    for &c in centroids.iter() {
                        if c == cur_c { continue; }
                        let sc = tables.anchor_sim(
                            &anchors[i], &anchors[c], threshold);
                        if sc < 0.0 { continue; }
                        if sc > best_score {
                            best_score = sc;
                            best_c     = c;
                        }
                    }
                }
            }
            if best_c != cur_c && best_score >= threshold {
                cur_centroid[i] = best_c;
                pass_reassigns += 1;
                changed = true;
            }
        }
        total_reassignments += pass_reassigns;

        // ── 3. Centroid merge ─────────────────────────────────────────
        let mut cluster_freq: HashMap<usize, u64> = HashMap::new();
        for i in 0..n {
            *cluster_freq.entry(cur_centroid[i]).or_insert(0) += freqs[i];
        }
        let mut sorted_cents: Vec<usize> =
            cluster_freq.keys().copied().collect();
        sorted_cents.sort_by(|&a, &b| {
            cluster_freq[&b].cmp(&cluster_freq[&a])
                .then(a.cmp(&b))  // deterministic tiebreak
        });

        let mut centroids_in_block: HashMap<(u8, u8), Vec<usize>> = HashMap::new();
        let mut block_of: HashMap<usize, (u8, u8)> = HashMap::new();
        for &c in sorted_cents.iter() {
            let bk = tables.block_key(&anchors[c]);
            centroids_in_block.entry(bk).or_default().push(c);
            block_of.insert(c, bk);
        }

        let mut merge_map: HashMap<usize, usize> = HashMap::new();
        let mut absorbed:  HashSet<usize>       = HashSet::new();

        for &c1 in sorted_cents.iter() {
            if absorbed.contains(&c1) { continue; }
            let (p1, q1) = block_of[&c1];
            let nb_keys: Vec<(u8, u8)> = centroids_in_block.keys()
                .filter(|b| b.0 == p1 || b.1 == q1)
                .copied()
                .collect();
            for bk2 in &nb_keys {
                if let Some(centroids2) = centroids_in_block.get(bk2) {
                    for &c2 in centroids2 {
                        if c2 == c1 || absorbed.contains(&c2) { continue; }
                        if cluster_freq[&c2] >= cluster_freq[&c1] {
                            continue;  // only smaller absorbed into larger
                        }
                        let sc = tables.anchor_sim(
                            &anchors[c1], &anchors[c2], threshold);
                        if sc >= threshold {
                            merge_map.insert(c2, c1);
                            absorbed.insert(c2);
                            total_merges += 1;
                            changed = true;
                        }
                    }
                }
            }
        }

        if !merge_map.is_empty() {
            // Transitive resolution: follow merge chains to the final target
            for i in 0..n {
                let mut c = cur_centroid[i];
                while let Some(&next) = merge_map.get(&c) {
                    c = next;
                }
                cur_centroid[i] = c;
            }
        }

        if !changed { break; }
    }

    // ── Build result dict (anchor_str → centroid_str) ─────────────────
    let result = PyDict::new(py);
    for i in 0..n {
        result.set_item(&names[i], &names[cur_centroid[i]])?;
    }

    let final_clusters: u64 = {
        let s: HashSet<usize> = cur_centroid.iter().copied().collect();
        s.len() as u64
    };

    // ── Build stats dict ──────────────────────────────────────────────
    let stats = PyDict::new(py);
    stats.set_item("passes",           passes_run)?;
    stats.set_item("medoid_changes",   total_medoid_changes)?;
    stats.set_item("reassignments",    total_reassignments)?;
    stats.set_item("merges",           total_merges)?;
    stats.set_item("initial_clusters", initial_clusters)?;
    stats.set_item("final_clusters",   final_clusters)?;

    Ok((result.into(), stats.into()))
}
