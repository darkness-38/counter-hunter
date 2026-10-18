//! Zarankiewicz K_{2,2}-free and K_{3,3}-free Maximum 0/1 Matrix Record Hunter
//!
//! Maximizes the number of 1s in an n x n bipartite adjacency matrix (each row as a u64 bitmask):
//! - z(n; 2, 2) (C_4-free): for all i < j, (row[i] & row[j]).count_ones() <= 1
//! - z(n; 3, 3) (K_{3,3}-free): for all i < j < k, (row[i] & row[j] & row[k]).count_ones() <= 2
//!
//! Uses Cyclic Difference Sets (Singer difference sets / PG(2,4)), projective plane submatrices,
//! and multi-threaded bitwise Tabu Search / Simulated Annealing with Rayon.

use std::time::Instant;
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
        if n == 21 || n == 15 {
            print_matrix_certificate(&record);
        }
    }
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
}
