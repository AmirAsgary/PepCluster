//! pepcluster._core — Rust backend for BLOSUM62-aware anchor clustering.
//!
//! Exposes two functions to Python:
//!   * `cluster_anchors(anchor_counts, threshold, anchor_positions, anchor_weight)`
//!   * `refine_clusters(anchor_counts, mapping, threshold, iterations, cap, merge,
//!                      anchor_positions, anchor_weight)`
//!
//! Anchors are the peptide's first `n_front` + last `n_back` residues (any
//! length). `anchor_positions` selects which positions *within the anchor* are
//! the binding anchors: they carry `anchor_weight` (default 2x) in the weighted
//! BLOSUM62 similarity and define the coarse-alphabet blocking. The default,
//! `[1, 5]` on a 6-residue anchor, is P2 and PΩ for MHC-I.

use pyo3::exceptions::PyValueError;
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

/// At most 8 anchor positions, so the block key packs into a u64 (one byte each).
const MAX_ANCHOR_POSITIONS: usize = 8;

/// Precomputed similarity + blocking tables for one (anchor length, anchor
/// positions, anchor weight) configuration.
struct SimTables {
    /// Anchor length (n_front + n_back).
    alen: usize,
    /// Positions in check order: highest weight first, so early termination
    /// kicks in sooner. Ties broken by position ascending (deterministic).
    check_pos: Vec<usize>,
    /// wsim_by_pos[p][a * 128 + b] = normalized BLOSUM62 * weight for position
    /// `p`. f64 (not f32) so the arithmetic is bit-identical to Python. Indexed
    /// by anchor position so the O(N) medoid decomposition can use it directly.
    wsim_by_pos: Vec<Vec<f64>>,
    /// remaining_after[k] = max score still obtainable from check_pos[k..].
    remaining_after: Vec<f64>,
    /// coarse group id for each ASCII byte (0..9), 255 = unknown
    coarse: [u8; 128],
    /// Anchor positions (sorted, deduped) — weighted up and used for blocking.
    anchor_positions: Vec<usize>,
    /// Normalized per-position weight (sums to 1 over all positions).
    weights: Vec<f64>,
    /// maxgroupsim[g1][g2] = max normalized BLOSUM62 similarity between any
    /// residue of coarse group g1 and any residue of coarse group g2. Used to
    /// upper-bound the similarity achievable between two blocks (OBG search).
    maxgroupsim: [[f64; 10]; 10],
    /// Σ of normalized weights over the anchor positions (the blocked ones).
    anchor_wsum: f64,
}

impl SimTables {
    fn new(alen: usize, anchor_positions: &[usize], anchor_weight: f64) -> Self {
        // Build AA → index mapping
        let mut aa_idx = [255u8; 128];
        for (i, &aa) in AA_ORDER.iter().enumerate() {
            aa_idx[aa as usize] = i as u8;
        }

        // Self-scores for normalization
        let mut self_score = [0.0f64; 20];
        for i in 0..20 {
            self_score[i] = BLOSUM62[i * 20 + i] as f64;
        }

        // Normalized similarity: sim(a,b) = B(a,b) / sqrt(B(a,a) * B(b,b))
        let mut norm_sim = [[0.0f64; 20]; 20];
        for i in 0..20 {
            for j in 0..20 {
                let denom = (self_score[i] * self_score[j]).sqrt();
                if denom > 0.0 {
                    norm_sim[i][j] = BLOSUM62[i * 20 + j] as f64 / denom;
                }
            }
        }

        // Raw weights: anchor_weight at anchor positions, 1.0 everywhere else.
        let mut raw = vec![1.0f64; alen];
        for &p in anchor_positions {
            raw[p] = anchor_weight;
        }
        let wsum: f64 = raw.iter().sum();

        // Check order: heaviest positions first, ties by position ascending.
        let mut check_pos: Vec<usize> = (0..alen).collect();
        check_pos.sort_by(|&a, &b| {
            raw[b].partial_cmp(&raw[a]).unwrap().then(a.cmp(&b))
        });

        // Weighted similarity tables, one per position (weight baked in).
        let mut wsim_by_pos: Vec<Vec<f64>> = Vec::with_capacity(alen);
        for pos in 0..alen {
            let w = raw[pos] / wsum;
            let mut t = vec![0.0f64; 128 * 128];
            for &a in AA_ORDER.iter() {
                let ai = aa_idx[a as usize] as usize;
                for &b in AA_ORDER.iter() {
                    let bi = aa_idx[b as usize] as usize;
                    t[(a as usize) * 128 + (b as usize)] = w * norm_sim[ai][bi];
                }
            }
            wsim_by_pos.push(t);
        }

        // remaining_after[k] = sum of normalized weights of check_pos[k..]
        let mut remaining_after = vec![0.0f64; alen + 1];
        for k in (0..alen).rev() {
            remaining_after[k] = remaining_after[k + 1] + raw[check_pos[k]] / wsum;
        }

        // Coarse group lookup
        let mut coarse = [255u8; 128];
        for (gid, &group) in COARSE_GROUPS.iter().enumerate() {
            for &aa in group {
                coarse[aa as usize] = gid as u8;
            }
        }

        // Max normalized similarity between any residues of two coarse groups.
        // The diagonal is 1.0 (a group can pair a residue with itself).
        let mut maxgroupsim = [[0.0f64; 10]; 10];
        for (g1, &grp1) in COARSE_GROUPS.iter().enumerate() {
            for (g2, &grp2) in COARSE_GROUPS.iter().enumerate() {
                // Start at -inf, not 0: some group pairs have a negative maximum
                // similarity, and clamping to 0 would inflate the bound (and
                // diverge from the Python backend's plain max()).
                let mut m = f64::NEG_INFINITY;
                for &x in grp1 {
                    let xi = aa_idx[x as usize] as usize;
                    for &y in grp2 {
                        let yi = aa_idx[y as usize] as usize;
                        if norm_sim[xi][yi] > m { m = norm_sim[xi][yi]; }
                    }
                }
                maxgroupsim[g1][g2] = m;
            }
        }

        let weights: Vec<f64> = (0..alen).map(|p| raw[p] / wsum).collect();
        let anchor_wsum: f64 = anchor_positions.iter().map(|&p| weights[p]).sum();

        SimTables {
            alen,
            check_pos,
            wsim_by_pos,
            remaining_after,
            coarse,
            anchor_positions: anchor_positions.to_vec(),
            weights,
            maxgroupsim,
            anchor_wsum,
        }
    }

    /// Max similarity between coarse groups g1 and g2. Same group (or any
    /// unknown residue, coded 255) can reach the identity maximum of 1.0.
    #[inline]
    fn maxgroupsim_of(&self, g1: u8, g2: u8) -> f64 {
        if g1 == g2 || g1 >= 10 || g2 >= 10 {
            1.0
        } else {
            self.maxgroupsim[g1 as usize][g2 as usize]
        }
    }

    /// Upper bound on the anchor similarity achievable between any anchor in
    /// block `k1` and any anchor in block `k2`. Non-anchor positions are
    /// unconstrained by the block key, so they contribute their full weight
    /// (max similarity 1). Admissible: never below the true achievable maximum.
    #[inline]
    fn block_upper_bound(&self, k1: u64, k2: u64) -> f64 {
        let m = self.anchor_positions.len();
        let mut ub = 1.0 - self.anchor_wsum; // all non-anchor positions at max
        for (i, &p) in self.anchor_positions.iter().enumerate() {
            let sh = 8 * (m - 1 - i);
            let g1 = ((k1 >> sh) & 0xFF) as u8;
            let g2 = ((k2 >> sh) & 0xFF) as u8;
            ub += self.weights[p] * self.maxgroupsim_of(g1, g2);
        }
        ub
    }

    /// Weighted BLOSUM62-normalized similarity with early termination.
    /// Returns the score, or -1.0 if it cannot possibly reach `threshold`.
    /// Positions are summed in check order (heaviest first) — same as Python.
    #[inline]
    fn anchor_sim(&self, a: &[u8], b: &[u8], threshold: f64) -> f64 {
        let mut s: f64 = 0.0;
        let n = self.alen;
        for k in 0..n {
            let p = self.check_pos[k];
            s += self.wsim_by_pos[p][(a[p] as usize) * 128 + (b[p] as usize)];
            if k + 1 < n && s + self.remaining_after[k + 1] < threshold {
                return -1.0;
            }
        }
        s
    }

    /// Full similarity, no early termination (medoid scoring / current centroid).
    #[inline]
    fn anchor_sim_full(&self, a: &[u8], b: &[u8]) -> f64 {
        let mut s: f64 = 0.0;
        for k in 0..self.alen {
            let p = self.check_pos[k];
            s += self.wsim_by_pos[p][(a[p] as usize) * 128 + (b[p] as usize)];
        }
        s
    }

    /// Coarse-alphabet block key: one byte per anchor position, first anchor
    /// position in the most-significant byte (so integer order == tuple order).
    #[inline]
    fn block_key(&self, a: &[u8]) -> u64 {
        let mut key: u64 = 0;
        for &p in self.anchor_positions.iter() {
            key = (key << 8) | (self.coarse[a[p] as usize] as u64);
        }
        key
    }

    /// Two blocks are neighbours if they agree on at least one anchor position.
    #[inline]
    fn is_neighbour(&self, k1: u64, k2: u64) -> bool {
        let n = self.anchor_positions.len();
        for i in 0..n {
            let sh = 8 * (n - 1 - i);
            if ((k1 >> sh) & 0xFF) == ((k2 >> sh) & 0xFF) {
                return true;
            }
        }
        false
    }
}

/// Row `i` of the flat anchor buffer.
#[inline]
fn arow(flat: &[u8], i: usize, alen: usize) -> &[u8] {
    &flat[i * alen..(i + 1) * alen]
}

/// Exact medoid of a cluster: the member with the highest frequency-weighted
/// average similarity to the others. O(k^2). `mems` must be anchor-sorted so
/// the strict-`>` argmax breaks ties toward the smallest anchor (matches Python).
fn medoid_exact(tables: &SimTables, flat: &[u8], freqs: &[u64],
                mems: &[usize], centroid: usize) -> usize {
    let alen = tables.alen;
    let mut best_member = centroid;
    let mut best_score: f64 = f64::NEG_INFINITY;
    for &i in mems.iter() {
        let mut s: f64 = 0.0;
        let mut w: f64 = 0.0;
        for &j in mems.iter() {
            if i == j { continue; }
            let sim = tables.anchor_sim_full(arow(flat, i, alen), arow(flat, j, alen));
            let fj = freqs[j] as f64;
            s += fj * sim;
            w += fj;
        }
        let avg = if w > 0.0 { s / w } else { 0.0 };
        if avg > best_score { best_score = avg; best_member = i; }
    }
    best_member
}

/// O(k) medoid via the additive per-position decomposition. Since
/// `sim(i,j) = Σ_p wsim_p(a_i[p], a_j[p])`, the total similarity of member `i`
/// to the whole cluster is `Sigma(i) = Σ_p S_p[a_i[p]]`, where
/// `S_p[x] = Σ_j f_j·wsim_p(x, a_j[p])`. With `sim(i,i)=1` and `F` the cluster
/// frequency, `avg(i) = (Sigma(i) − f_i) / (F − f_i)`. No pairwise loop.
///
/// `mems` must be anchor-sorted, and present residues are summed in ascending
/// byte order, so the result is bit-identical to the Python backend.
fn medoid_fast(tables: &SimTables, flat: &[u8], freqs: &[u64],
               mems: &[usize], centroid: usize) -> usize {
    let alen = tables.alen;

    // Total cluster frequency (exact integer sum, order-independent).
    let mut total_f: u64 = 0;
    for &i in mems.iter() { total_f += freqs[i]; }
    let total_ff = total_f as f64;

    // Per position: S_p[x] for each residue x that occurs there.
    let mut s_tables: Vec<[f64; 128]> = vec![[0.0f64; 128]; alen];
    for p in 0..alen {
        let mut freq_res = [0u64; 128];
        for &i in mems.iter() {
            freq_res[arow(flat, i, alen)[p] as usize] += freqs[i];
        }
        let present: Vec<usize> = (0..128usize).filter(|&r| freq_res[r] > 0).collect();
        let wp = &tables.wsim_by_pos[p];
        for &x in present.iter() {
            let base = x * 128;
            let mut sx = 0.0f64;
            for &aa in present.iter() {
                sx += (freq_res[aa] as f64) * wp[base + aa];
            }
            s_tables[p][x] = sx;
        }
    }

    let mut best_member = centroid;
    let mut best_score: f64 = f64::NEG_INFINITY;
    for &i in mems.iter() {
        let anc = arow(flat, i, alen);
        let mut sigma = 0.0f64;
        for p in 0..alen {
            sigma += s_tables[p][anc[p] as usize];
        }
        let fi = freqs[i] as f64;
        let denom = total_ff - fi;
        let avg = if denom > 0.0 { (sigma - fi) / denom } else { 0.0 };
        if avg > best_score { best_score = avg; best_member = i; }
    }
    best_member
}

/// Validate + canonicalise (sort, dedup) the anchor positions.
fn check_positions(ap: &[usize], alen: usize) -> PyResult<Vec<usize>> {
    let mut v = ap.to_vec();
    v.sort_unstable();
    v.dedup();
    if v.is_empty() {
        return Err(PyValueError::new_err(
            "anchor_positions must contain at least one position",
        ));
    }
    if v.len() > MAX_ANCHOR_POSITIONS {
        return Err(PyValueError::new_err(format!(
            "at most {} anchor positions are supported, got {}",
            MAX_ANCHOR_POSITIONS,
            v.len()
        )));
    }
    if let Some(&mx) = v.last() {
        if mx >= alen {
            return Err(PyValueError::new_err(format!(
                "anchor position {} is out of range for anchor length {}",
                mx, alen
            )));
        }
    }
    Ok(v)
}

/// Parse the anchor_counts dict into (names, flat anchor bytes, counts, alen),
/// validating that every anchor has the same length.
fn parse_anchor_counts(
    anchor_counts: &Bound<'_, PyDict>,
) -> PyResult<(Vec<String>, Vec<u8>, Vec<u64>, usize)> {
    let mut names: Vec<String> = Vec::with_capacity(anchor_counts.len());
    let mut counts: Vec<u64> = Vec::with_capacity(anchor_counts.len());
    for (key, val) in anchor_counts.iter() {
        names.push(key.extract()?);
        counts.push(val.extract()?);
    }
    if names.is_empty() {
        return Ok((names, Vec::new(), counts, 0));
    }
    let alen = names[0].as_bytes().len();
    if alen == 0 {
        return Err(PyValueError::new_err("anchors must be non-empty"));
    }
    let mut flat: Vec<u8> = Vec::with_capacity(names.len() * alen);
    for nm in names.iter() {
        let b = nm.as_bytes();
        if b.len() != alen {
            return Err(PyValueError::new_err(format!(
                "all anchors must have the same length; found {} and {}",
                alen,
                b.len()
            )));
        }
        flat.extend_from_slice(b);
    }
    Ok((names, flat, counts, alen))
}

/// Greedy centroid clustering of unique anchors.
///
/// Args:
///     anchor_counts:    dict[str, int] — unique anchors → peptide frequency
///     threshold:        float — minimum similarity to join a cluster
///     anchor_positions: list[int] — 0-based positions within the anchor that
///                       are binding anchors (weighted up + used for blocking)
///     anchor_weight:    float — weight of anchor positions (others are 1.0)
///     obg_block_search: bool — search neighbouring blocks whose upper bound is
///                       at least the threshold, not just the anchor's own block
///     obg_max_probes:   int — max blocks searched per anchor (incl. own block);
///                       <= 0 means unlimited
///     obg_min_block_upper_bound: float — only search blocks whose upper bound
///                       is at least this (the effective cut is
///                       max(threshold, this))
///
/// Returns:
///     (mapping, n_comparisons, n_early_exits)
#[pyfunction]
#[pyo3(signature = (anchor_counts, threshold, anchor_positions=vec![1, 5], anchor_weight=2.0,
                    obg_block_search=false, obg_max_probes=0, obg_min_block_upper_bound=0.0))]
fn cluster_anchors(
    py: Python<'_>,
    anchor_counts: &Bound<'_, PyDict>,
    threshold: f64,
    anchor_positions: Vec<usize>,
    anchor_weight: f64,
    obg_block_search: bool,
    obg_max_probes: i64,
    obg_min_block_upper_bound: f64,
) -> PyResult<(PyObject, u64, u64)> {
    let (names_in, flat_in, counts_in, alen) = parse_anchor_counts(anchor_counts)?;
    if names_in.is_empty() {
        return Ok((PyDict::new(py).into(), 0, 0));
    }
    let ap = check_positions(&anchor_positions, alen)?;
    let tables = SimTables::new(alen, &ap, anchor_weight);
    let cutoff = threshold.max(obg_min_block_upper_bound);

    // Sort by frequency descending (stable: ties keep dict insertion order,
    // matching the Python backend).
    let mut order: Vec<usize> = (0..names_in.len()).collect();
    order.sort_by(|&a, &b| counts_in[b].cmp(&counts_in[a]));

    let n = order.len();
    let mut names: Vec<String> = Vec::with_capacity(n);
    let mut flat: Vec<u8> = Vec::with_capacity(n * alen);
    for &oi in order.iter() {
        names.push(names_in[oi].clone());
        flat.extend_from_slice(arow(&flat_in, oi, alen));
    }

    // Greedy centroid clustering. Anchors are processed most-frequent-first;
    // each joins the first centroid at or above `threshold`, else becomes a new
    // centroid. Centroids are indexed by block key. When OBG search is on, an
    // anchor also considers centroids in neighbouring blocks whose upper bound
    // reaches the threshold, ranked by upper bound (own block always first).
    let mut centroids_by_block: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut centroid_of: Vec<usize> = vec![0; n];
    let mut n_cmp: u64 = 0;
    let mut n_early: u64 = 0;

    for i in 0..n {
        let ki = tables.block_key(arow(&flat, i, alen));

        // Neighbour blocks to probe, beyond the anchor's own block.
        let mut probe_blocks: Vec<u64> = Vec::new();
        if obg_block_search {
            let mut cand: Vec<(f64, u64)> = Vec::new();
            for &kb in centroids_by_block.keys() {
                if kb == ki { continue; }
                let ub = tables.block_upper_bound(ki, kb);
                if ub >= cutoff {
                    cand.push((ub, kb));
                }
            }
            // Highest upper bound first; ties by block key for determinism.
            cand.sort_by(|a, b| {
                b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1))
            });
            if obg_max_probes > 0 {
                // Own block counts as one probe.
                let keep = (obg_max_probes - 1).max(0) as usize;
                cand.truncate(keep);
            }
            probe_blocks.extend(cand.into_iter().map(|(_, kb)| kb));
        }

        // Search own block first, then the ranked neighbour blocks.
        let mut matched: Option<usize> = None;
        if let Some(cs) = centroids_by_block.get(&ki) {
            for &ci in cs.iter() {
                n_cmp += 1;
                let s = tables.anchor_sim(
                    arow(&flat, i, alen), arow(&flat, ci, alen), threshold);
                if s < 0.0 { n_early += 1; continue; }
                if s >= threshold { matched = Some(ci); break; }
            }
        }
        if matched.is_none() {
            'nb: for &kb in probe_blocks.iter() {
                if let Some(cs) = centroids_by_block.get(&kb) {
                    for &ci in cs.iter() {
                        n_cmp += 1;
                        let s = tables.anchor_sim(
                            arow(&flat, i, alen), arow(&flat, ci, alen), threshold);
                        if s < 0.0 { n_early += 1; continue; }
                        if s >= threshold { matched = Some(ci); break 'nb; }
                    }
                }
            }
        }

        match matched {
            Some(ci) => centroid_of[i] = ci,
            None => {
                centroids_by_block.entry(ki).or_default().push(i);
                centroid_of[i] = i;
            }
        }
    }

    // Build result dict: anchor_str → centroid_str
    let result = PyDict::new(py);
    for i in 0..n {
        result.set_item(&names[i], &names[centroid_of[i]])?;
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
// Optional refinement (Lloyd-style; mirrors the Python refine_clusters_py)
// ============================================================================

/// Refinement pass on top of greedy clustering output.
///
/// Each pass performs three sub-steps:
///   1. Medoid update     — replace each centroid with the member that
///                          maximises frequency-weighted mean similarity
///                          to the cluster's other members.
///   2. Cross-block reassign — for each anchor, find the best centroid above
///                          `threshold` across its own block plus neighbouring
///                          blocks (blocks agreeing on at least one anchor
///                          position). Bounded by `cap` comparisons per anchor.
///   3. Centroid merge    — if two centroids satisfy sim >= threshold, absorb
///                          the smaller cluster into the larger. Skipped when
///                          `merge` is false.
///
/// Stops early when no change occurs in a pass.
///
/// Args:
///     anchor_counts:    dict[str, int]  — unique anchor → frequency
///     mapping:          dict[str, str]  — unique anchor → centroid
///     threshold:        float           — same as initial clustering
///     iterations:       int             — max passes
///     cap:              int             — max centroid comparisons per anchor
///                                         in reassignment. <= 0 means no cap.
///     merge:            bool            — run the centroid-merge sub-step
///     anchor_positions: list[int]       — 0-based anchor positions
///     anchor_weight:    float           — weight of anchor positions
///     fast_medoid:      bool            — use the O(N) per-position medoid
///                                         decomposition instead of the exact
///                                         O(k^2) all-pairs medoid
///     merge_cap:        int             — max candidate centroids examined per
///                                         centroid in the merge step (own-block
///                                         first, largest first). <= 0 = no cap.
///
/// Returns:
///     (refined_mapping, stats_dict)
#[pyfunction]
#[pyo3(signature = (anchor_counts, mapping, threshold, iterations, cap=32, merge=true,
                    anchor_positions=vec![1, 5], anchor_weight=2.0,
                    fast_medoid=false, merge_cap=0))]
fn refine_clusters(
    py: Python<'_>,
    anchor_counts: &Bound<'_, PyDict>,
    mapping: &Bound<'_, PyDict>,
    threshold: f64,
    iterations: usize,
    cap: i64,
    merge: bool,
    anchor_positions: Vec<usize>,
    anchor_weight: f64,
    fast_medoid: bool,
    merge_cap: i64,
) -> PyResult<(PyObject, PyObject)> {
    let (names, flat, freqs, alen) = parse_anchor_counts(anchor_counts)?;
    let n = names.len();
    if n == 0 {
        let stats = PyDict::new(py);
        for k in ["passes", "medoid_changes", "reassignments", "merges",
                  "initial_clusters", "final_clusters"] {
            stats.set_item(k, 0u64)?;
        }
        return Ok((PyDict::new(py).into(), stats.into()));
    }
    let ap = check_positions(&anchor_positions, alen)?;
    let tables = SimTables::new(alen, &ap, anchor_weight);

    let mut str_to_idx: HashMap<&str, usize> = HashMap::with_capacity(n);
    for (i, nm) in names.iter().enumerate() {
        str_to_idx.insert(nm.as_str(), i);
    }

    // ── Parse mapping into cur_centroid[i] ────────────────────────────
    // Default: every anchor is its own centroid (only matters if mapping is
    // missing an entry, which shouldn't happen in normal use).
    let mut cur_centroid: Vec<usize> = (0..n).collect();
    for (key, val) in mapping.iter() {
        let a: String = key.extract()?;
        let c: String = val.extract()?;
        if let (Some(&ai), Some(&ci)) =
            (str_to_idx.get(a.as_str()), str_to_idx.get(c.as_str()))
        {
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
        // Canonical member order (anchor ascending) so the medoid argmax and
        // its ties resolve deterministically — identical across runs and the
        // Python backend, independent of hash iteration order.
        for v in members.values_mut() {
            v.sort_by(|&a, &b| names[a].cmp(&names[b]));
        }

        let mut centroid_remap: HashMap<usize, usize> = HashMap::with_capacity(members.len());
        for (&centroid, mems) in members.iter() {
            if mems.len() == 1 {
                centroid_remap.insert(centroid, centroid);
                continue;
            }
            let best_member = if fast_medoid {
                medoid_fast(&tables, &flat, &freqs, mems, centroid)
            } else {
                medoid_exact(&tables, &flat, &freqs, mems, centroid)
            };
            centroid_remap.insert(centroid, best_member);
            if best_member != centroid {
                total_medoid_changes += 1;
                changed = true;
            }
        }
        for i in 0..n {
            cur_centroid[i] = centroid_remap[&cur_centroid[i]];
        }

        // ── 2. Cross-block reassignment (candidate-capped) ────────────
        let current_centroids: HashSet<usize> =
            cur_centroid.iter().copied().collect();

        // Cluster frequency drives the "most promising first" candidate order.
        let mut cluster_freq2: HashMap<usize, u64> = HashMap::new();
        for i in 0..n {
            *cluster_freq2.entry(cur_centroid[i]).or_insert(0) += freqs[i];
        }

        let mut centroids_in_block: HashMap<u64, Vec<usize>> = HashMap::new();
        for &c in current_centroids.iter() {
            centroids_in_block
                .entry(tables.block_key(arow(&flat, c, alen)))
                .or_default()
                .push(c);
        }
        // Within each block, examine larger clusters first (then anchor asc for
        // a deterministic order that matches the Python backend exactly).
        for v in centroids_in_block.values_mut() {
            v.sort_by(|&a, &b| {
                cluster_freq2[&b].cmp(&cluster_freq2[&a])
                    .then_with(|| names[a].cmp(&names[b]))
            });
        }

        // Neighbour blocks, own block first, then the rest in sorted order.
        let mut sorted_blocks: Vec<u64> =
            centroids_in_block.keys().copied().collect();
        sorted_blocks.sort_unstable();
        let mut neighbours_of: HashMap<u64, Vec<u64>> = HashMap::new();
        for &bk in centroids_in_block.keys() {
            let mut nb: Vec<u64> = vec![bk];
            for &b in sorted_blocks.iter() {
                if b != bk && tables.is_neighbour(bk, b) {
                    nb.push(b);
                }
            }
            neighbours_of.insert(bk, nb);
        }

        let fallback_nb = vec![]; // empty fallback if anchor's block has no centroids
        let mut pass_reassigns: u64 = 0;
        for i in 0..n {
            let bk    = tables.block_key(arow(&flat, i, alen));
            let cur_c = cur_centroid[i];
            let cur_score = tables.anchor_sim_full(
                arow(&flat, i, alen), arow(&flat, cur_c, alen));
            let mut best_c     = cur_c;
            let mut best_score = cur_score;

            let nbs = neighbours_of.get(&bk).unwrap_or(&fallback_nb);
            let mut examined: i64 = 0;
            'scan: for &nb in nbs.iter() {
                if let Some(centroids) = centroids_in_block.get(&nb) {
                    for &c in centroids.iter() {
                        if c == cur_c { continue; }
                        if cap > 0 && examined >= cap { break 'scan; }
                        let sc = tables.anchor_sim(
                            arow(&flat, i, alen), arow(&flat, c, alen), threshold);
                        examined += 1;
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

        // ── 3. Centroid merge (optional; skipped when merge=false) ────
        if merge {
            let mut cluster_freq: HashMap<usize, u64> = HashMap::new();
            for i in 0..n {
                *cluster_freq.entry(cur_centroid[i]).or_insert(0) += freqs[i];
            }
            let mut sorted_cents: Vec<usize> =
                cluster_freq.keys().copied().collect();
            sorted_cents.sort_by(|&a, &b| {
                cluster_freq[&b].cmp(&cluster_freq[&a])
                    .then_with(|| names[a].cmp(&names[b]))  // anchor asc (matches Python)
            });

            // Blocks hold centroids in sorted_cents order (freq desc, anchor
            // asc) already, so each block's candidate list is largest-first.
            let mut centroids_in_block: HashMap<u64, Vec<usize>> = HashMap::new();
            let mut block_of: HashMap<usize, u64> = HashMap::new();
            for &c in sorted_cents.iter() {
                let bk = tables.block_key(arow(&flat, c, alen));
                centroids_in_block.entry(bk).or_default().push(c);
                block_of.insert(c, bk);
            }

            // Own block first, then neighbour blocks in sorted order — the same
            // deterministic traversal the reassignment uses, so the merge cap
            // is bit-identical across backends (and unchanged when uncapped).
            let mut sorted_blocks: Vec<u64> =
                centroids_in_block.keys().copied().collect();
            sorted_blocks.sort_unstable();
            let mut neighbours_of: HashMap<u64, Vec<u64>> = HashMap::new();
            for &bk in centroids_in_block.keys() {
                let mut nb: Vec<u64> = vec![bk];
                for &b in sorted_blocks.iter() {
                    if b != bk && tables.is_neighbour(bk, b) {
                        nb.push(b);
                    }
                }
                neighbours_of.insert(bk, nb);
            }

            let mut merge_map: HashMap<usize, usize> = HashMap::new();
            let mut absorbed:  HashSet<usize>       = HashSet::new();

            for &c1 in sorted_cents.iter() {
                if absorbed.contains(&c1) { continue; }
                let bk1 = block_of[&c1];
                let mut examined: i64 = 0;
                'mscan: for &bk2 in neighbours_of[&bk1].iter() {
                    if let Some(centroids2) = centroids_in_block.get(&bk2) {
                        for &c2 in centroids2 {
                            if c2 == c1 { continue; }
                            if merge_cap > 0 && examined >= merge_cap { break 'mscan; }
                            examined += 1;
                            if absorbed.contains(&c2) { continue; }
                            if cluster_freq[&c2] >= cluster_freq[&c1] {
                                continue;  // only smaller absorbed into larger
                            }
                            let sc = tables.anchor_sim(
                                arow(&flat, c1, alen), arow(&flat, c2, alen), threshold);
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
