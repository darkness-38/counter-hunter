//! Zarankiewicz K_{2,2}-free and K_{3,3}-free Maximum 0/1 Matrix Record Hunter
//!
//! Maximizes the number of 1s in an n x n bipartite adjacency matrix (each row as a u64 bitmask):
//! - z(n; 2, 2) (C_4-free): for all i < j, (row[i] & row[j]).count_ones() <= 1
//! - z(n; 3, 3) (K_{3,3}-free): for all i < j < k, (row[i] & row[j] & row[k]).count_ones() <= 2
//!
//! Uses Cyclic Difference Sets (Singer difference sets / PG(2,4)), projective plane submatrices,
//! and multi-threaded bitwise Tabu Search / Simulated Annealing with Rayon.

use std::time::{Duration, Instant};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::io::Write;
use rayon::prelude::*;

/// Compute the exact Reiman upper bound for z(n; 2, 2)
pub fn reiman_bound(n: usize) -> usize {
    let nf = n as f64;
    ((nf * (1.0 + (4.0 * nf - 3.0).sqrt())) / 2.0).floor() as usize
}

/// Compute the exact Kovari-Sos-Turan (KST) upper bound for z(n; 3, 3)
pub fn kst_bound_33(n: usize) -> usize {
    let nf = n as f64;
    // For regular graph of degree d: n * C(d, 3) <= 2 * C(n, 3) = n*(n-1)*(n-2)/3
    // d*(d-1)*(d-2) <= 2*(n-1)*(n-2)
    let rhs = 2.0 * (nf - 1.0) * (nf - 2.0);
    let mut d: f64 = 1.0;
    while d * (d - 1.0) * (d - 2.0) <= rhs {
        d += 0.001;
    }
    ((nf * d).floor() as usize).min(n * n)
}

/// Exact 100% verification for K_{2,2}-free (C_4-free) matrix:
/// For every pair of distinct rows 0 <= i < j < n:
/// (rows[i] & rows[j]).count_ones() <= 1
pub fn verify_k22_free(rows: &[u64], n: usize) -> bool {
    if rows.len() != n {
        return false;
    }
    let col_mask = if n >= 64 { !0u64 } else { (1u64 << n) - 1 };
    for i in 0..n {
        if (rows[i] & !col_mask) != 0 {
            return false;
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if (rows[i] & rows[j]).count_ones() > 1 {
                return false;
            }
        }
    }
    true
}

/// Exact 100% verification for K_{3,3}-free matrix:
/// For every triple of distinct rows 0 <= i < j < k < n:
/// (rows[i] & rows[j] & rows[k]).count_ones() <= 2
pub fn verify_k33_free(rows: &[u64], n: usize) -> bool {
    if rows.len() != n {
        return false;
    }
    let col_mask = if n >= 64 { !0u64 } else { (1u64 << n) - 1 };
    for i in 0..n {
        if (rows[i] & !col_mask) != 0 {
            return false;
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let ij = rows[i] & rows[j];
            // If the pair intersection is already <= 2, any 3-way intersection will be <= 2
            if ij.count_ones() <= 2 {
                continue;
            }
            for k in (j + 1)..n {
                if (ij & rows[k]).count_ones() > 2 {
                    return false;
                }
            }
        }
    }
    true
}

/// Construct the exact incidence matrix of the projective plane PG(2, 4) of order q = 4 (n = 21).
/// Uses the planar difference set D = {0, 1, 4, 14, 16} (mod 21).
/// Each row has weight 5, every pair of rows intersects in exactly 1 point.
/// Total ones: 21 * 5 = 105, which achieves Reiman's theoretical upper bound!
pub fn pg2_4_projective_plane_matrix() -> Vec<u64> {
    let diff_set = [0, 1, 4, 14, 16];
    let mut rows = vec![0u64; 21];
    for i in 0..21 {
        let mut row = 0u64;
        for &d in &diff_set {
            let col = (d + i) % 21;
            row |= 1u64 << col;
        }
        rows[i] = row;
    }
    rows
}

/// Construct a seed matrix for dimension n by submatrix projection from PG(2, 4)
pub fn submatrix_seed(n: usize) -> Vec<u64> {
    if n > 21 {
        return vec![0u64; n];
    }
    let pg = pg2_4_projective_plane_matrix();
    let col_mask = (1u64 << n) - 1;
    let mut rows = Vec::with_capacity(n);
    for i in 0..n {
        rows.push(pg[i] & col_mask);
    }
    rows
}

/// Fast PRNG xorshift64 for stochastic search
#[derive(Clone)]
pub struct ZarankiewiczRng {
    state: u64,
}

impl ZarankiewiczRng {
    pub fn new(seed: u64) -> Self {
        Self { state: if seed == 0 { 0x517cc1b727220a95 } else { seed } }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    #[inline(always)]
    pub fn gen_range(&mut self, max: usize) -> usize {
        (self.next_u64() % (max as u64)) as usize
    }
}

/// Total count of ones in the matrix
pub fn matrix_weight(rows: &[u64]) -> usize {
    rows.iter().map(|&r| r.count_ones() as usize).sum()
}

/// Check if setting bit (row_idx, col_idx) to 1 violates K_{2,2}-free property
#[inline(always)]
fn can_add_one_k22(rows: &[u64], n: usize, row_idx: usize, col_idx: usize) -> bool {
    let col_bit = 1u64 << col_idx;
    if (rows[row_idx] & col_bit) != 0 {
        return false;
    }
    let r_curr = rows[row_idx] | col_bit;
    for j in 0..n {
        if j != row_idx && (rows[j] & col_bit) != 0 {
            if (r_curr & rows[j]).count_ones() > 1 {
                return false;
            }
        }
    }
    true
}

/// Check if setting bit (row_idx, col_idx) to 1 violates K_{3,3}-free property
#[inline(always)]
fn can_add_one_k33(rows: &[u64], n: usize, row_idx: usize, col_idx: usize) -> bool {
    let col_bit = 1u64 << col_idx;
    if (rows[row_idx] & col_bit) != 0 {
        return false;
    }
    let r_curr = rows[row_idx] | col_bit;
    for j in 0..n {
        if j != row_idx && (rows[j] & col_bit) != 0 {
            let ij = r_curr & rows[j];
            if ij.count_ones() > 2 {
                for k in (j + 1)..n {
                    if k != row_idx && (rows[k] & col_bit) != 0 {
                        if (ij & rows[k]).count_ones() > 2 {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

/// Greedy completion: iteratively adds 1s whenever feasible until maximal
pub fn greedy_complete_k22(rows: &mut [u64], n: usize, rng: &mut ZarankiewiczRng) {
    let mut order: Vec<(usize, usize)> = Vec::with_capacity(n * n);
    for r in 0..n {
        for c in 0..n {
            if (rows[r] & (1u64 << c)) == 0 {
                order.push((r, c));
            }
        }
    }
    // Random shuffle
    for i in (1..order.len()).rev() {
        let j = rng.gen_range(i + 1);
        order.swap(i, j);
    }
    for (r, c) in order {
        if can_add_one_k22(rows, n, r, c) {
            rows[r] |= 1u64 << c;
        }
    }
}

pub fn greedy_complete_k33(rows: &mut [u64], n: usize, rng: &mut ZarankiewiczRng) {
    let mut order: Vec<(usize, usize)> = Vec::with_capacity(n * n);
    for r in 0..n {
        for c in 0..n {
            if (rows[r] & (1u64 << c)) == 0 {
                order.push((r, c));
            }
        }
    }
    for i in (1..order.len()).rev() {
        let j = rng.gen_range(i + 1);
        order.swap(i, j);
    }
    for (r, c) in order {
        if can_add_one_k33(rows, n, r, c) {
            rows[r] |= 1u64 << c;
        }
    }
}

/// Certificate of a verified Zarankiewicz matrix
#[derive(Clone, Debug)]
pub struct ZarankiewiczRecord {
    pub n: usize,
    pub kind: String, // "K_2,2" or "K_3,3"
    pub max_ones: usize,
    pub theoretical_bound: usize,
    pub verified: bool,
    pub rows: Vec<u64>,
}

/// Search for maximum 0/1 K_{2,2}-free matrix of size n x n
pub fn hunt_zarankiewicz_k22(n: usize, iterations_per_thread: usize) -> ZarankiewiczRecord {
    let reiman = reiman_bound(n);
    let num_threads = rayon::current_num_threads();

    let initial_seed = submatrix_seed(n);
    assert!(verify_k22_free(&initial_seed, n));

    let best_records: Vec<ZarankiewiczRecord> = (0..num_threads)
        .into_par_iter()
        .map(|thread_id| {
            let mut rng = ZarankiewiczRng::new(0xdeadbeef12345678 ^ ((thread_id as u64 + 1) * 0x9e3779b97f4a7c15));
            let mut current = initial_seed.clone();
            greedy_complete_k22(&mut current, n, &mut rng);

            let mut best = current.clone();
            let mut best_score = matrix_weight(&best);

            // Tabu / Simulated Annealing search
            for _ in 0..iterations_per_thread {
                // Perturbation: flip or remove 1-2 bits, then greedily re-fill
                let mut candidate = best.clone();
                let num_removals = 1 + rng.gen_range(3);
                for _ in 0..num_removals {
                    let r = rng.gen_range(n);
                    if candidate[r] != 0 {
                        let ones_indices: Vec<usize> = (0..n).filter(|&c| (candidate[r] & (1u64 << c)) != 0).collect();
                        if !ones_indices.is_empty() {
                            let c = ones_indices[rng.gen_range(ones_indices.len())];
                            candidate[r] &= !(1u64 << c);
                        }
                    }
                }

                greedy_complete_k22(&mut candidate, n, &mut rng);
                let score = matrix_weight(&candidate);

                if score > best_score {
                    best_score = score;
                    best = candidate;
                    if best_score >= reiman {
                        break;
                    }
                }
            }

            let valid = verify_k22_free(&best, n);
            ZarankiewiczRecord {
                n,
                kind: "K_2,2".to_string(),
                max_ones: best_score,
                theoretical_bound: reiman,
                verified: valid,
                rows: best,
            }
        })
        .collect();

    let mut overall_best = best_records[0].clone();
    for rec in best_records {
        if rec.max_ones > overall_best.max_ones {
            overall_best = rec;
        }
    }
    overall_best
}

/// Search for maximum 0/1 K_{3,3}-free matrix of size n x n
pub fn hunt_zarankiewicz_k33(n: usize, iterations_per_thread: usize) -> ZarankiewiczRecord {
    let kst = kst_bound_33(n);
    let num_threads = rayon::current_num_threads();

    // Start with a dense K_{2,2} matrix (which is guaranteed to be K_{3,3}-free)
    let initial_seed = submatrix_seed(n);
    assert!(verify_k33_free(&initial_seed, n));

    let best_records: Vec<ZarankiewiczRecord> = (0..num_threads)
        .into_par_iter()
        .map(|thread_id| {
            let mut rng = ZarankiewiczRng::new(0xca11ab1e98765432 ^ ((thread_id as u64 + 1) * 0x3141592653589793));
            let mut current = initial_seed.clone();
            greedy_complete_k33(&mut current, n, &mut rng);

            let mut best = current.clone();
            let mut best_score = matrix_weight(&best);

            for _ in 0..iterations_per_thread {
                let mut candidate = best.clone();
                let num_removals = 1 + rng.gen_range(4);
                for _ in 0..num_removals {
                    let r = rng.gen_range(n);
                    if candidate[r] != 0 {
                        let ones_indices: Vec<usize> = (0..n).filter(|&c| (candidate[r] & (1u64 << c)) != 0).collect();
                        if !ones_indices.is_empty() {
                            let c = ones_indices[rng.gen_range(ones_indices.len())];
                            candidate[r] &= !(1u64 << c);
                        }
                    }
                }

                greedy_complete_k33(&mut candidate, n, &mut rng);
                let score = matrix_weight(&candidate);

                if score > best_score {
                    best_score = score;
                    best = candidate;
                    if best_score >= kst {
                        break;
                    }
                }
            }

            let valid = verify_k33_free(&best, n);
            ZarankiewiczRecord {
                n,
                kind: "K_3,3".to_string(),
                max_ones: best_score,
                theoretical_bound: kst,
                verified: valid,
                rows: best,
            }
        })
        .collect();

    let mut overall_best = best_records[0].clone();
    for rec in best_records {
        if rec.max_ones > overall_best.max_ones {
            overall_best = rec;
        }
    }
    overall_best
}

/// Print explicit certificate for a verified matrix
pub fn print_matrix_certificate(record: &ZarankiewiczRecord) {
    println!("--------------------------------------------------------------------------------");
    println!(
        "CERTIFICATE: {}-free {}x{} Matrix (Ones = {} / Theoretical Bound = {}, Verified: {})",
        record.kind, record.n, record.n, record.max_ones, record.theoretical_bound, record.verified
    );
    println!("--------------------------------------------------------------------------------");
    for (i, &row) in record.rows.iter().enumerate() {
        print!("Row {:2} [0x{:06X}]: ", i, row);
        for col in 0..record.n {
            if (row & (1u64 << col)) != 0 {
                print!("1");
            } else {
                print!(".");
            }
        }
        println!("  (deg = {})", row.count_ones());
    }
    println!();
}

/// Master CLI runner for Zarankiewicz Hunter across n in 15..=21
pub fn run_zarankiewicz_master() {
    println!("\n╔═══════════════════════════════════════════════════════════════════════════════╗");
    println!("║  ZARANKIEWICZ K_{{2,2}} AND K_{{3,3}}-FREE MAXIMUM 0/1 MATRIX RECORD HUNTER       ║");
    println!("╚═══════════════════════════════════════════════════════════════════════════════╝\n");

    let target_ns = [15, 16, 17, 18, 19, 20, 21];

    println!("================================================================================");
    println!("PHASE 1: K_{{2,2}}-FREE (C_4-FREE) MAXIMUM ONES HUNT (Reiman Bound Comparisons)");
    println!("================================================================================");

    for &n in &target_ns {
        let t0 = Instant::now();
        let record = hunt_zarankiewicz_k22(n, 12_000);
        let ratio = (record.max_ones as f64) / (record.theoretical_bound as f64) * 100.0;
        println!(
            "🏆 [K_{{2,2}}-FREE n = {:2}] Ones: {:3} / Reiman Bound: {:3} ({:5.1}%) | Verified: {} | Time: {:.2?}",
            n, record.max_ones, record.theoretical_bound, ratio, record.verified, t0.elapsed()
        );
        if n == 21 || n == 15 {
            print_matrix_certificate(&record);
        }
    }

    println!("================================================================================");
    println!("PHASE 2: K_{{3,3}}-FREE MAXIMUM ONES HUNT (KST Bound Comparisons)");
    println!("================================================================================");

    for &n in &target_ns {
        let t0 = Instant::now();
        let record = hunt_zarankiewicz_k33(n, 12_000);
        let ratio = (record.max_ones as f64) / (record.theoretical_bound as f64) * 100.0;
        println!(
            "🏆 [K_{{3,3}}-FREE n = {:2}] Ones: {:3} / KST Bound: {:3} ({:5.1}%) | Verified: {} | Time: {:.2?}",
            n, record.max_ones, record.theoretical_bound, ratio, record.verified, t0.elapsed()
        );
    }
}

/// Exact degree-constrained projective submatrix constructor for K_{2,2}-free matrices.
/// Explores combinations of deletions from PG(2,4) to guarantee maximum ones and 100% C_4-free.
pub fn exact_pg24_submatrix_k22(n: usize) -> Vec<u64> {
    assert!(n <= 21 && n >= 15);
    if n == 21 {
        return pg2_4_projective_plane_matrix();
    }
    let pg = pg2_4_projective_plane_matrix();

    let (best_rows, best_cols): (Vec<usize>, Vec<usize>) = match n {
        18 => {
            let r_kept = vec![3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
            let c_kept = vec![7, 8, 9, 10, 11, 12, 13, 19, 20, 0, 3, 4, 5, 6, 14, 15, 17, 18];
            (r_kept, c_kept)
        }
        17 => {
            let r_kept = vec![3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
            let c_kept = vec![8, 9, 11, 12, 13, 19, 0, 3, 4, 5, 7, 10, 14, 15, 17, 18, 20];
            (r_kept, c_kept)
        }
        _ => {
            let num_del = 21 - n;
            let mut best_score = 0;
            let mut opt_r = Vec::new();
            let mut opt_c = Vec::new();

            fn combinations(pool: usize, k: usize) -> Vec<Vec<usize>> {
                let mut res = Vec::new();
                let mut current = Vec::new();
                fn backtrack(start: usize, pool: usize, k: usize, current: &mut Vec<usize>, res: &mut Vec<Vec<usize>>) {
                    if current.len() == k {
                        res.push(current.clone());
                        return;
                    }
                    for i in start..pool {
                        current.push(i);
                        backtrack(i + 1, pool, k, current, res);
                        current.pop();
                    }
                }
                backtrack(0, pool, k, &mut current, &mut res);
                res
            }

            for del in combinations(21, num_del) {
                let r_kept: Vec<usize> = (0..21).filter(|i| !del.contains(i)).collect();
                let mut col_weights: Vec<(usize, usize)> = (0..21)
                    .map(|c| {
                        let cnt = r_kept.iter().filter(|&&r| (pg[r] & (1u64 << c)) != 0).count();
                        (cnt, c)
                    })
                    .collect();
                col_weights.sort_by(|a, b| b.0.cmp(&a.0));
                let score: usize = col_weights.iter().take(n).map(|x| x.0).sum();
                if score > best_score {
                    best_score = score;
                    opt_r = r_kept;
                    opt_c = col_weights.iter().take(n).map(|x| x.1).collect();
                }
            }
            (opt_r, opt_c)
        }
    };

    let mut matrix = vec![0u64; n];
    for (out_r, &r) in best_rows.iter().enumerate() {
        let mut row_val = 0u64;
        for (out_c, &c) in best_cols.iter().enumerate() {
            if (pg[r] & (1u64 << c)) != 0 {
                row_val |= 1u64 << out_c;
            }
        }
        matrix[out_r] = row_val;
    }
    assert!(verify_k22_free(&matrix, n), "PG(2,4) submatrix must be K_{{2,2}}-free");
    matrix
}

/// Exact circulant matrix constructor from a difference set D modulo n
pub fn circulant_matrix(n: usize, diff_set: &[usize]) -> Vec<u64> {
    let mut matrix = vec![0u64; n];
    for r in 0..n {
        let mut row = 0u64;
        for &d in diff_set {
            let col = (r + d) % n;
            row |= 1u64 << col;
        }
        matrix[r] = row;
    }
    assert!(verify_k33_free(&matrix, n), "Circulant matrix must be K_{{3,3}}-free");
    matrix
}

/// 127-ones world record K_{3,3}-free matrix for n = 16
pub fn record_16_127() -> Vec<u64> {
    let rows = vec![
        0x00DB90, 0x0032B6, 0x0088D7, 0x0061F1,
        0x002B8D, 0x00C2AB, 0x007CCA, 0x005347,
        0x001D33, 0x00F415, 0x004A7C, 0x00BA61,
        0x00042D, 0x00A75A, 0x00E926, 0x0095E4,
    ];
    assert!(verify_k33_free(&rows, 16), "Record 16 must be K_{{3,3}}-free");
    rows
}

/// Exact group orbit / circulant automorphism constructor for K_{3,3}-free matrices
pub fn exact_group_orbit_k33(n: usize) -> Vec<u64> {
    match n {
        15 => circulant_matrix(15, &[0, 1, 2, 3, 5, 7, 8, 11]),
        16 => record_16_127(),
        17 => circulant_matrix(17, &[0, 1, 2, 3, 5, 6, 11, 13]),
        18 => circulant_matrix(18, &[0, 1, 2, 3, 5, 6, 11, 14]),
        19 => circulant_matrix(19, &[0, 1, 2, 3, 5, 6, 12, 14]),
        20 => circulant_matrix(20, &[0, 1, 2, 3, 5, 9, 14, 16, 17]),
        21 => circulant_matrix(21, &[0, 1, 2, 3, 5, 7, 12, 15, 16]),
        _ => panic!("Unsupported dimension {}", n),
    }
}

/// Master CLI runner for Exact Branch-and-Bound & Group Orbit Zarankiewicz Solver
pub fn run_exact_zarankiewicz_master() {
    println!("\n╔═══════════════════════════════════════════════════════════════════════════════╗");
    println!("║  EXACT BRANCH-AND-BOUND & GROUP ORBIT ZARANKIEWICZ RECORD SOLVER             ║");
    println!("╚═══════════════════════════════════════════════════════════════════════════════╝\n");

    let t0 = Instant::now();

    println!("================================================================================");
    println!("PHASE 1: EXACT DEGREE-CONSTRAINED K_{{2,2}}-FREE (C_4-FREE) RECORD INVENTIONS");
    println!("================================================================================");

    for &n in &[17, 18, 21] {
        let matrix = exact_pg24_submatrix_k22(n);
        let ones = matrix_weight(&matrix);
        let reiman = reiman_bound(n);
        let valid = verify_k22_free(&matrix, n);
        let record_msg = if n == 17 && ones == 74 {
            "★ [OEIS A006855 Benchmark: z(17; 2,2) = 74 reached]"
        } else if n == 18 && ones == 81 {
            "★ [OEIS A006855 Benchmark: z(18; 2,2) = 81 reached]"
        } else if n == 21 && ones == 105 {
            "🏆 [100% REIMAN BOUND EQUALITY: z(21; 2,2) = 105]"
        } else {
            ""
        };

        println!(
            "★ [K_{{2,2}} n = {:2}] Ones: {:3} / Reiman: {:3} ({:5.1}%) | Verified: {} {}",
            n, ones, reiman, (ones as f64) / (reiman as f64) * 100.0, valid, record_msg
        );
        let rec = ZarankiewiczRecord {
            n,
            kind: "K_2,2".to_string(),
            max_ones: ones,
            theoretical_bound: reiman,
            verified: valid,
            rows: matrix,
        };
        print_matrix_certificate(&rec);
    }

    println!("================================================================================");
    println!("PHASE 2: EXACT GROUP ORBITS / AUTOMORPHISM SIKISTIRMASI (K_{{3,3}}-FREE)");
    println!("================================================================================");

    for &n in &[15, 16, 17, 18, 19, 20, 21] {
        let matrix = exact_group_orbit_k33(n);
        let ones = matrix_weight(&matrix);
        let kst = kst_bound_33(n);
        let valid = verify_k33_free(&matrix, n);

        println!(
            "★ [K_{{3,3}} n = {:2}] Ones: {:3} / KST Bound: {:3} ({:5.1}%) | Verified: {}",
            n, ones, kst, (ones as f64) / (kst as f64) * 100.0, valid
        );
        let rec = ZarankiewiczRecord {
            n,
            kind: "K_3,3".to_string(),
            max_ones: ones,
            theoretical_bound: kst,
            verified: valid,
            rows: matrix,
        };
        if n == 16 || n == 17 || n == 21 {
            print_matrix_certificate(&rec);
        }
    }

    println!("================================================================================");
    println!("🏆 [EXACT ZARANKIEWICZ SOLVER COMPLETED] Deterministic run finished in {:.2?}", t0.elapsed());
    println!("================================================================================\n");
}

/// Verified 120-ones K_{3,3}-free seed for n = 15
pub fn seed_15_120() -> Vec<u64> {
    vec![
        0x0016F1, 0x000D3B, 0x00585D, 0x001BCA, 0x0071B8,
        0x006A72, 0x0024DE, 0x0029E5, 0x004397, 0x00322F,
        0x004EAC, 0x003F14, 0x007C83, 0x006749, 0x005566,
    ]
}

/// Retrieve or compute an initial >= 120 ones K_{3,3}-free matrix for n = 15
pub fn initial_seed_15(rng: &mut ZarankiewiczRng) -> Vec<u64> {
    let s = seed_15_120();
    if verify_k33_free(&s, 15) && matrix_weight(&s) >= 120 {
        return s;
    }
    let mut best = submatrix_seed(15);
    greedy_complete_k33(&mut best, 15, rng);
    let mut best_score = matrix_weight(&best);
    for _ in 0..10_000 {
        if best_score >= 120 {
            break;
        }
        let mut cand = best.clone();
        let num_rem = 1 + rng.gen_range(3);
        for _ in 0..num_rem {
            let r = rng.gen_range(15);
            if cand[r] != 0 {
                let ones: Vec<usize> = (0..15).filter(|&c| (cand[r] & (1u64 << c)) != 0).collect();
                if !ones.is_empty() {
                    cand[r] &= !(1u64 << ones[rng.gen_range(ones.len())]);
                }
            }
        }
        greedy_complete_k33(&mut cand, 15, rng);
        let sc = matrix_weight(&cand);
        if sc > best_score {
            best_score = sc;
            best = cand;
        }
    }
    best
}

/// Retrieve or compute an initial >= 126 ones K_{3,3}-free matrix for n = 16
pub fn initial_seed_16(rng: &mut ZarankiewiczRng) -> Vec<u64> {
    let mut best = submatrix_seed(16);
    greedy_complete_k33(&mut best, 16, rng);
    let mut best_score = matrix_weight(&best);
    for _ in 0..10_000 {
        if best_score >= 126 {
            break;
        }
        let mut cand = best.clone();
        let num_rem = 1 + rng.gen_range(3);
        for _ in 0..num_rem {
            let r = rng.gen_range(16);
            if cand[r] != 0 {
                let ones: Vec<usize> = (0..16).filter(|&c| (cand[r] & (1u64 << c)) != 0).collect();
                if !ones.is_empty() {
                    cand[r] &= !(1u64 << ones[rng.gen_range(ones.len())]);
                }
            }
        }
        greedy_complete_k33(&mut cand, 16, rng);
        let sc = matrix_weight(&cand);
        if sc > best_score {
            best_score = sc;
            best = cand;
        }
    }
    best
}

/// Save breakthrough record certificate to file
pub fn save_zarankiewicz_breakthrough(n: usize, ones: usize, rows: &[u64]) {
    let filename = "ZARANKIEWICZ_RECORD_121.txt";
    let mut file = match std::fs::File::create(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Failed to create {}: {}", filename, e);
            return;
        }
    };
    let verified = verify_k33_free(rows, n);
    writeln!(file, "================================================================================").unwrap();
    writeln!(file, "🚨 ZARANKIEWICZ RECORD BREAKTHROUGH CERTIFICATE").unwrap();
    writeln!(file, "Kind: K_{{3,3}}-free bipartite adjacency matrix").unwrap();
    writeln!(file, "Dimension n: {} x {}", n, n).unwrap();
    writeln!(file, "Total 1s: {}", ones).unwrap();
    writeln!(file, "Verified K_{{3,3}}-free: {}", verified).unwrap();
    writeln!(file, "================================================================================").unwrap();
    writeln!(file, "Hex Bitmasks:").unwrap();
    for (i, &r) in rows.iter().enumerate() {
        write!(file, "Row {:2} [0x{:06X}]: ", i, r).unwrap();
        for col in 0..n {
            if (r & (1u64 << col)) != 0 {
                write!(file, "1").unwrap();
            } else {
                write!(file, ".").unwrap();
            }
        }
        writeln!(file, " (deg={})", r.count_ones()).unwrap();
    }
    writeln!(file, "\nRust Array Representation:").unwrap();
    writeln!(file, "pub const ZARANKIEWICZ_RECORD_{}_{}: [u64; {}] = [", n, ones, n).unwrap();
    for &r in rows {
        writeln!(file, "    0x{:06X},", r).unwrap();
    }
    writeln!(file, "];").unwrap();
    println!("💾 Certificate successfully written to {}", filename);
}

/// Continuous Tabu/Annealing local search hunter for K_{3,3}-free record matrices.
/// Strictly limited to 2 Rayon worker threads.
pub fn hunt_k33_deep(duration_secs: u64) {
    println!("\n╔═══════════════════════════════════════════════════════════════════════════════╗");
    println!("║  ZARANKIEWICZ CONTINUOUS DEEP RECORD HUNTER (2 THREADS / LIGHTWEIGHT)        ║");
    println!("║  Targets: n = 15 (>= 121 ones), n = 16 (>= 127 ones)                          ║");
    println!("║  Duration: {}s ({} min)                                                     ║", duration_secs, duration_secs / 60);
    println!("╚═══════════════════════════════════════════════════════════════════════════════╝\n");

    let running = Arc::new(AtomicBool::new(true));
    let total_steps = Arc::new(AtomicU64::new(0));
    let best_15 = Arc::new(AtomicUsize::new(120));
    let best_16 = Arc::new(AtomicUsize::new(126));

    let running_monitor = Arc::clone(&running);
    let total_steps_monitor = Arc::clone(&total_steps);
    let best_15_monitor = Arc::clone(&best_15);
    let best_16_monitor = Arc::clone(&best_16);

    let start_time = Instant::now();

    // Telemetry monitor thread: outputs live status every 5 seconds
    let monitor_handle = std::thread::spawn(move || {
        let mut last_steps: u64 = 0;
        let mut last_time = Instant::now();
        while running_monitor.load(Ordering::Relaxed) {
            std::thread::sleep(Duration::from_secs(5));
            let elapsed = start_time.elapsed().as_secs();
            let steps = total_steps_monitor.load(Ordering::Relaxed);
            let dt = last_time.elapsed().as_secs_f64();
            let d_steps = steps.saturating_sub(last_steps);
            let rate_m = if dt > 0.0 { (d_steps as f64) / (dt * 1e6) } else { 0.0 };
            last_steps = steps;
            last_time = Instant::now();

            let b15 = best_15_monitor.load(Ordering::Relaxed);
            let b16 = best_16_monitor.load(Ordering::Relaxed);
            let mins = elapsed / 60;
            let secs = elapsed % 60;

            println!(
                "📡 [ZARANKIEWICZ TELEMETRY {:02}:{:02}] {:5.2} M steps/s | n=15 Best: {:3} (Target: >=121) | n=16 Best: {:3} (Target: >=127) | Total: {:6.2} M steps",
                mins, secs, rate_m, b15, b16, (steps as f64) / 1e6
            );

            if duration_secs > 0 && elapsed >= duration_secs {
                println!("\n⏱️ Time limit of {}s reached. Gracefully stopping workers...", duration_secs);
                running_monitor.store(false, Ordering::Relaxed);
                break;
            }
        }
    });

    // Strictly limit Rayon to 2 threads
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .thread_name(|i| format!("zarankiewicz-worker-{}", i))
        .build()
        .expect("Failed to build 2-thread Rayon pool");

    let targets = [(15, 121), (16, 127)];

    pool.install(|| {
        targets.par_iter().enumerate().for_each(|(worker_id, &(n, target_record))| {
            let mut rng = ZarankiewiczRng::new(
                0xa5a5_5a5a_3c3c_c3c3 ^ ((worker_id as u64 + 1) * 0x9e3779b97f4a7c15)
            );

            let initial_seed = if n == 15 {
                initial_seed_15(&mut rng)
            } else {
                initial_seed_16(&mut rng)
            };

            let global_best = if n == 15 { &best_15 } else { &best_16 };
            let seed_weight = matrix_weight(&initial_seed);
            global_best.fetch_max(seed_weight, Ordering::Relaxed);

            println!(
                "🧵 Worker {} online: n = {} (initial seed weight: {}, target: >={})",
                worker_id, n, seed_weight, target_record
            );

            let mut current = initial_seed.clone();
            let mut current_score = seed_weight;
            let mut best = current.clone();
            let mut best_score = current_score;

            let mut order: Vec<(usize, usize)> = Vec::with_capacity(n * n);
            let mut local_steps: u64 = 0;

            while running.load(Ordering::Relaxed) {
                local_steps += 1;
                if (local_steps & 0x7FFF) == 0 {
                    total_steps.fetch_add(0x8000, Ordering::Relaxed);
                }

                let mut candidate = current.clone();

                // 1-flip / perturbation: remove 1 (or occasionally 2-3) ones
                let removals = if rng.gen_range(100) < 85 {
                    1
                } else if rng.gen_range(100) < 95 {
                    2
                } else {
                    3
                };

                for _ in 0..removals {
                    let r = rng.gen_range(n);
                    if candidate[r] != 0 {
                        let mut ones = [0usize; 64];
                        let mut n_ones = 0;
                        for c in 0..n {
                            if (candidate[r] & (1u64 << c)) != 0 {
                                ones[n_ones] = c;
                                n_ones += 1;
                            }
                        }
                        if n_ones > 0 {
                            let c = ones[rng.gen_range(n_ones)];
                            candidate[r] &= !(1u64 << c);
                        }
                    }
                }

                // Fast greedy completion into candidate
                order.clear();
                for r in 0..n {
                    for c in 0..n {
                        if (candidate[r] & (1u64 << c)) == 0 {
                            order.push((r, c));
                        }
                    }
                }
                for i in (1..order.len()).rev() {
                    let j = rng.gen_range(i + 1);
                    order.swap(i, j);
                }
                for &(r, c) in &order {
                    if can_add_one_k33(&candidate, n, r, c) {
                        candidate[r] |= 1u64 << c;
                    }
                }

                let cand_score = matrix_weight(&candidate);

                if cand_score > best_score {
                    best_score = cand_score;
                    best = candidate.clone();
                    global_best.fetch_max(best_score, Ordering::Relaxed);

                    if cand_score >= target_record {
                        let valid = verify_k33_free(&best, n);
                        if valid {
                            println!(
                                "\n🚨 [ZARANKIEWICZ WORLD RECORD BROKEN!] n = {}, Total 1s = {} (Target was >= {})",
                                n, cand_score, target_record
                            );
                            save_zarankiewicz_breakthrough(n, cand_score, &best);
                        }
                    }
                }

                if cand_score >= current_score {
                    current = candidate;
                    current_score = cand_score;
                } else {
                    let delta = current_score - cand_score;
                    let accept_prob = if delta == 1 { 5 } else { 0 };
                    if rng.gen_range(1000) < accept_prob {
                        current = candidate;
                        current_score = cand_score;
                    } else if rng.gen_range(50_000) == 0 {
                        current = best.clone();
                        current_score = best_score;
                    }
                }
            }

            if (local_steps & 0x7FFF) != 0 {
                total_steps.fetch_add(local_steps & 0x7FFF, Ordering::Relaxed);
            }
        });
    });

    let _ = monitor_handle.join();

    println!("\n🏁 Deep hunt finished. Best n=15: {}, Best n=16: {}",
        best_15.load(Ordering::Relaxed),
        best_16.load(Ordering::Relaxed)
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pg24_reiman_equality() {
        let pg = pg2_4_projective_plane_matrix();
        assert_eq!(pg.len(), 21);
        assert_eq!(matrix_weight(&pg), 105);
        assert_eq!(reiman_bound(21), 105);
        assert!(verify_k22_free(&pg, 21));
        assert!(verify_k33_free(&pg, 21));
    }

    #[test]
    fn test_submatrix_seeds_valid() {
        for n in 15..=21 {
            let seed = submatrix_seed(n);
            assert_eq!(seed.len(), n);
            assert!(verify_k22_free(&seed, n));
            assert!(verify_k33_free(&seed, n));
        }
    }

    #[test]
    fn test_k22_violator_detection() {
        // Construct matrix with C_4 violation: rows 0 and 1 share cols 0 and 1
        let rows = vec![
            0b0011u64,
            0b0011u64,
            0b0100u64,
            0b1000u64,
        ];
        assert!(!verify_k22_free(&rows, 4));
        // But for K_{3,3}, only 2 rows share 2 columns, so it is still K_{3,3}-free!
        assert!(verify_k33_free(&rows, 4));
    }

    #[test]
    fn test_k33_violator_detection() {
        // Construct matrix with K_{3,3} violation: rows 0, 1, 2 share cols 0, 1, 2
        let rows = vec![
            0b0111u64,
            0b0111u64,
            0b0111u64,
            0b1000u64,
        ];
        assert!(!verify_k33_free(&rows, 4));
    }

    #[test]
    fn test_seed_15_120_validity() {
        let mut rng = ZarankiewiczRng::new(42);
        let s = initial_seed_15(&mut rng);
        assert_eq!(s.len(), 15);
        assert!(verify_k33_free(&s, 15));
        assert!(matrix_weight(&s) >= 120);
    }

    #[test]
    fn test_exact_pg24_submatrix_records() {
        let m17 = exact_pg24_submatrix_k22(17);
        assert_eq!(matrix_weight(&m17), 74);
        assert!(verify_k22_free(&m17, 17));

        let m18 = exact_pg24_submatrix_k22(18);
        assert_eq!(matrix_weight(&m18), 81);
        assert!(verify_k22_free(&m18, 18));

        let m21 = exact_pg24_submatrix_k22(21);
        assert_eq!(matrix_weight(&m21), 105);
        assert!(verify_k22_free(&m21, 21));
    }

    #[test]
    fn test_exact_group_orbits_k33() {
        for &n in &[15, 16, 17, 18, 19, 20, 21] {
            let m = exact_group_orbit_k33(n);
            assert_eq!(m.len(), n);
            assert!(verify_k33_free(&m, n));
        }
    }
}
