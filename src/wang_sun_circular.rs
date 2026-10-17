//! HEDEF 2: Wang & Sun (EJC 2026, #P2.64) Circular Restricted Sumset Breaker
//!
//! Investigates the open conjecture by Han Wang and Zhi-Wei Sun:
//! "On Two New Kinds of Restricted Sumsets", Electronic Journal of Combinatorics 33(2) (2026), #P2.64.
//!
//! Conjecture 1, Inequality (2):
//! For any additive abelian group G and subsets A_1, ..., A_n of G with |A_i| >= 2 (1 <= i <= n):
//!   |C(A_1, ..., A_n)| >= min{ p(G), sum_{i=1}^n |A_i| - 2n + (-1)^n (1 + {n}_2) }
//! where
//!   C(A_1, ..., A_n) = { a_1 + ... + a_n : a_i in A_i, a_1 != a_2, a_2 != a_3, ..., a_n != a_1 }
//! and p(G) is the smallest prime divisor of |G| (or infinity if G = Z).

use rayon::prelude::*;
use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

/// Fast pseudo-random number generator (Xorshift64*)
#[derive(Clone, Copy)]
pub struct FastRng {
    state: u64,
}

impl FastRng {
    #[inline(always)]
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x5a17_d00d_cafe_babe } else { seed },
        }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    #[inline(always)]
    pub fn gen_range(&mut self, min: usize, max: usize) -> usize {
        if max <= min {
            return min;
        }
        min + (self.next_u64() as usize % (max - min))
    }
}

/// Computes the conjectured lower bound:
/// Bound = min{ p(G), sum |A_i| - 2n + (-1)^n (1 + {n}_2) }
pub fn compute_conjectured_bound(sizes: &[usize], n: usize, p_g: usize) -> usize {
    let sum_sizes: usize = sizes.iter().sum();
    let term = if n % 2 == 0 { 1i32 } else { -2i32 };
    let linear_bound = (sum_sizes as i32) - 2 * (n as i32) + term;
    if linear_bound <= 0 {
        return 0;
    }
    (linear_bound as usize).min(p_g)
}

// -----------------------------------------------------------------------------
// Fast Bitmask Evaluators for Z_p
// -----------------------------------------------------------------------------

#[inline(always)]
pub fn eval_circular_zp_4(a1: &[u8], a2: &[u8], a3: &[u8], a4: &[u8], p: u8) -> u32 {
    let mut mask = 0u64;
    for &x1 in a1 {
        for &x2 in a2 {
            if x2 == x1 { continue; }
            for &x3 in a3 {
                if x3 == x2 { continue; }
                for &x4 in a4 {
                    if x4 == x3 || x4 == x1 { continue; }
                    let sum = (x1 + x2 + x3 + x4) % p;
                    mask |= 1u64 << sum;
                }
            }
        }
    }
    mask.count_ones()
}

#[inline(always)]
pub fn eval_circular_zp_5(a1: &[u8], a2: &[u8], a3: &[u8], a4: &[u8], a5: &[u8], p: u8) -> u32 {
    let mut mask = 0u64;
    for &x1 in a1 {
        for &x2 in a2 {
            if x2 == x1 { continue; }
            for &x3 in a3 {
                if x3 == x2 { continue; }
                for &x4 in a4 {
                    if x4 == x3 { continue; }
                    for &x5 in a5 {
                        if x5 == x4 || x5 == x1 { continue; }
                        let sum = (x1 + x2 + x3 + x4 + x5) % p;
                        mask |= 1u64 << sum;
                    }
                }
            }
        }
    }
    mask.count_ones()
}

#[inline(always)]
pub fn eval_circular_zp_6(a1: &[u8], a2: &[u8], a3: &[u8], a4: &[u8], a5: &[u8], a6: &[u8], p: u8) -> u32 {
    let mut mask = 0u64;
    for &x1 in a1 {
        for &x2 in a2 {
            if x2 == x1 { continue; }
            for &x3 in a3 {
                if x3 == x2 { continue; }
                for &x4 in a4 {
                    if x4 == x3 { continue; }
                    for &x5 in a5 {
                        if x5 == x4 { continue; }
                        for &x6 in a6 {
                            if x6 == x5 || x6 == x1 { continue; }
                            let sum = (x1 + x2 + x3 + x4 + x5 + x6) % p;
                            mask |= 1u64 << sum;
                        }
                    }
                }
            }
        }
    }
    mask.count_ones()
}

// -----------------------------------------------------------------------------
// Fast Bitmask Evaluators for Z (Integers)
// -----------------------------------------------------------------------------

#[inline(always)]
pub fn eval_circular_z_4(a1: &[u8], a2: &[u8], a3: &[u8], a4: &[u8]) -> u32 {
    let mut mask = 0u128;
    for &x1 in a1 {
        for &x2 in a2 {
            if x2 == x1 { continue; }
            for &x3 in a3 {
                if x3 == x2 { continue; }
                for &x4 in a4 {
                    if x4 == x3 || x4 == x1 { continue; }
                    let sum = x1 + x2 + x3 + x4;
                    mask |= 1u128 << sum;
                }
            }
        }
    }
    mask.count_ones()
}

#[inline(always)]
pub fn eval_circular_z_5(a1: &[u8], a2: &[u8], a3: &[u8], a4: &[u8], a5: &[u8]) -> u32 {
    let mut mask = 0u128;
    for &x1 in a1 {
        for &x2 in a2 {
            if x2 == x1 { continue; }
            for &x3 in a3 {
                if x3 == x2 { continue; }
                for &x4 in a4 {
                    if x4 == x3 { continue; }
                    for &x5 in a5 {
                        if x5 == x4 || x5 == x1 { continue; }
                        let sum = x1 + x2 + x3 + x4 + x5;
                        mask |= 1u128 << sum;
                    }
                }
            }
        }
    }
    mask.count_ones()
}

#[inline(always)]
pub fn eval_circular_z_6(a1: &[u8], a2: &[u8], a3: &[u8], a4: &[u8], a5: &[u8], a6: &[u8]) -> u32 {
    let mut mask = 0u128;
    for &x1 in a1 {
        for &x2 in a2 {
            if x2 == x1 { continue; }
            for &x3 in a3 {
                if x3 == x2 { continue; }
                for &x4 in a4 {
                    if x4 == x3 { continue; }
                    for &x5 in a5 {
                        if x5 == x4 { continue; }
                        for &x6 in a6 {
                            if x6 == x5 || x6 == x1 { continue; }
                            let sum = x1 + x2 + x3 + x4 + x5 + x6;
                            mask |= 1u128 << sum;
                        }
                    }
                }
            }
        }
    }
    mask.count_ones()
}

// -----------------------------------------------------------------------------
// Combinatorics Helpers
// -----------------------------------------------------------------------------

pub fn combinations(universe: &[u8], k: usize) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut current = Vec::new();

    fn comb(
        start: usize,
        k: usize,
        universe: &[u8],
        current: &mut Vec<u8>,
        result: &mut Vec<Vec<u8>>,
    ) {
        if current.len() == k {
            result.push(current.clone());
            return;
        }
        for i in start..universe.len() {
            current.push(universe[i]);
            comb(i + 1, k, universe, current, result);
            current.pop();
        }
    }

    comb(0, k, universe, &mut current, &mut result);
    result
}

// -----------------------------------------------------------------------------
// Exhaustive Parallel Search
// -----------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct SearchResult {
    pub group_desc: String,
    pub n: usize,
    pub sizes: Vec<usize>,
    pub bound: usize,
    pub tested_configs: u64,
    pub min_cardinality: usize,
    pub elapsed_secs: f64,
    pub counterexample_found: bool,
}

/// Exhaustive scan over all combinations of subset sizes in Z_p
pub fn scan_exhaustive_zp(p: u8, sizes: &[usize]) -> SearchResult {
    let n = sizes.len();
    let bound = compute_conjectured_bound(sizes, n, p as usize);
    let u: Vec<u8> = (0..p).collect();
    let c_subsets: Vec<Vec<Vec<u8>>> = sizes.iter().map(|&s| combinations(&u, s)).collect();

    let total_tested = AtomicU64::new(0);
    let min_card = AtomicUsize::new(usize::MAX);
    let ce_flag = AtomicBool::new(false);
    let ce_detail = Mutex::new(None);

    let start = Instant::now();

    match n {
        4 => {
            c_subsets[0].par_iter().for_each(|a1| {
                let mut local_min = usize::MAX;
                let mut local_count = 0u64;

                for a2 in &c_subsets[1] {
                    for a3 in &c_subsets[2] {
                        for a4 in &c_subsets[3] {
                            local_count += 1;
                            let c = eval_circular_zp_4(a1, a2, a3, a4, p) as usize;
                            if c < local_min {
                                local_min = c;
                            }
                            if c < bound {
                                ce_flag.store(true, Ordering::SeqCst);
                                let mut lk = ce_detail.lock().unwrap();
                                if lk.is_none() {
                                    *lk = Some((vec![a1.clone(), a2.clone(), a3.clone(), a4.clone()], c));
                                }
                            }
                        }
                    }
                }
                total_tested.fetch_add(local_count, Ordering::Relaxed);
                min_card.fetch_min(local_min, Ordering::Relaxed);
            });
        }
        5 => {
            c_subsets[0].par_iter().for_each(|a1| {
                let mut local_min = usize::MAX;
                let mut local_count = 0u64;

                for a2 in &c_subsets[1] {
                    for a3 in &c_subsets[2] {
                        for a4 in &c_subsets[3] {
                            for a5 in &c_subsets[4] {
                                local_count += 1;
                                let c = eval_circular_zp_5(a1, a2, a3, a4, a5, p) as usize;
                                if c < local_min {
                                    local_min = c;
                                }
                                if c < bound {
                                    ce_flag.store(true, Ordering::SeqCst);
                                    let mut lk = ce_detail.lock().unwrap();
                                    if lk.is_none() {
                                        *lk = Some((
                                            vec![a1.clone(), a2.clone(), a3.clone(), a4.clone(), a5.clone()],
                                            c,
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
                total_tested.fetch_add(local_count, Ordering::Relaxed);
                min_card.fetch_min(local_min, Ordering::Relaxed);
            });
        }
        6 => {
            c_subsets[0].par_iter().for_each(|a1| {
                let mut local_min = usize::MAX;
                let mut local_count = 0u64;

                for a2 in &c_subsets[1] {
                    for a3 in &c_subsets[2] {
                        for a4 in &c_subsets[3] {
                            for a5 in &c_subsets[4] {
                                for a6 in &c_subsets[5] {
                                    local_count += 1;
                                    let c = eval_circular_zp_6(a1, a2, a3, a4, a5, a6, p) as usize;
                                    if c < local_min {
                                        local_min = c;
                                    }
                                    if c < bound {
                                        ce_flag.store(true, Ordering::SeqCst);
                                        let mut lk = ce_detail.lock().unwrap();
                                        if lk.is_none() {
                                            *lk = Some((
                                                vec![
                                                    a1.clone(), a2.clone(), a3.clone(), a4.clone(),
                                                    a5.clone(), a6.clone(),
                                                ],
                                                c,
                                            ));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                total_tested.fetch_add(local_count, Ordering::Relaxed);
                min_card.fetch_min(local_min, Ordering::Relaxed);
            });
        }
        _ => panic!("Unsupported n = {}", n),
    }

    let elapsed = start.elapsed().as_secs_f64();
    let is_ce = ce_flag.load(Ordering::SeqCst);

    if is_ce {
        if let Some((subsets, card)) = ce_detail.lock().unwrap().as_ref() {
            println!(
                "\x1b[1;31m🚨🚨🚨 DISPROVED! COUNTEREXAMPLE FOUND IN Z_{}: {:?} with |C| = {} < bound {}! 🚨🚨🚨\x1b[0m",
                p, subsets, card, bound
            );
            let _ = OpenOptions::new()
                .create(true)
                .append(true)
                .open("discovered_counterexamples.txt")
                .map(|mut f| {
                    writeln!(
                        f,
                        "WANG & SUN CIRCULAR SUMSET COUNTEREXAMPLE in Z_{}: Subsets = {:?}, |C| = {}, Bound = {}",
                        p, subsets, card, bound
                    )
                });
        }
    }

    SearchResult {
        group_desc: format!("Z_{}", p),
        n,
        sizes: sizes.to_vec(),
        bound,
        tested_configs: total_tested.load(Ordering::Relaxed),
        min_cardinality: min_card.load(Ordering::Relaxed),
        elapsed_secs: elapsed,
        counterexample_found: is_ce,
    }
}

/// Exhaustive scan over all combinations of subset sizes in Z with universe {0..u_max-1}
pub fn scan_exhaustive_z(u_max: u8, sizes: &[usize]) -> SearchResult {
    let n = sizes.len();
    let bound = compute_conjectured_bound(sizes, n, usize::MAX);
    let u: Vec<u8> = (0..u_max).collect();
    let c_subsets: Vec<Vec<Vec<u8>>> = sizes.iter().map(|&s| combinations(&u, s)).collect();

    let total_tested = AtomicU64::new(0);
    let min_card = AtomicUsize::new(usize::MAX);
    let ce_flag = AtomicBool::new(false);
    let ce_detail = Mutex::new(None);

    let start = Instant::now();

    match n {
        4 => {
            c_subsets[0].par_iter().for_each(|a1| {
                let mut local_min = usize::MAX;
                let mut local_count = 0u64;

                for a2 in &c_subsets[1] {
                    for a3 in &c_subsets[2] {
                        for a4 in &c_subsets[3] {
                            local_count += 1;
                            let c = eval_circular_z_4(a1, a2, a3, a4) as usize;
                            if c < local_min {
                                local_min = c;
                            }
                            if c < bound {
                                ce_flag.store(true, Ordering::SeqCst);
                                let mut lk = ce_detail.lock().unwrap();
                                if lk.is_none() {
                                    *lk = Some((vec![a1.clone(), a2.clone(), a3.clone(), a4.clone()], c));
                                }
                            }
                        }
                    }
                }
                total_tested.fetch_add(local_count, Ordering::Relaxed);
                min_card.fetch_min(local_min, Ordering::Relaxed);
            });
        }
        5 => {
            c_subsets[0].par_iter().for_each(|a1| {
                let mut local_min = usize::MAX;
                let mut local_count = 0u64;

                for a2 in &c_subsets[1] {
                    for a3 in &c_subsets[2] {
                        for a4 in &c_subsets[3] {
                            for a5 in &c_subsets[4] {
                                local_count += 1;
                                let c = eval_circular_z_5(a1, a2, a3, a4, a5) as usize;
                                if c < local_min {
                                    local_min = c;
                                }
                                if c < bound {
                                    ce_flag.store(true, Ordering::SeqCst);
                                    let mut lk = ce_detail.lock().unwrap();
                                    if lk.is_none() {
                                        *lk = Some((
                                            vec![a1.clone(), a2.clone(), a3.clone(), a4.clone(), a5.clone()],
                                            c,
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
                total_tested.fetch_add(local_count, Ordering::Relaxed);
                min_card.fetch_min(local_min, Ordering::Relaxed);
            });
        }
        6 => {
            c_subsets[0].par_iter().for_each(|a1| {
                let mut local_min = usize::MAX;
                let mut local_count = 0u64;

                for a2 in &c_subsets[1] {
                    for a3 in &c_subsets[2] {
                        for a4 in &c_subsets[3] {
                            for a5 in &c_subsets[4] {
                                for a6 in &c_subsets[5] {
                                    local_count += 1;
                                    let c = eval_circular_z_6(a1, a2, a3, a4, a5, a6) as usize;
                                    if c < local_min {
                                        local_min = c;
                                    }
                                    if c < bound {
                                        ce_flag.store(true, Ordering::SeqCst);
                                        let mut lk = ce_detail.lock().unwrap();
                                        if lk.is_none() {
                                            *lk = Some((
                                                vec![
                                                    a1.clone(), a2.clone(), a3.clone(), a4.clone(),
                                                    a5.clone(), a6.clone(),
                                                ],
                                                c,
                                            ));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                total_tested.fetch_add(local_count, Ordering::Relaxed);
                min_card.fetch_min(local_min, Ordering::Relaxed);
            });
        }
        _ => panic!("Unsupported n = {}", n),
    }

    let elapsed = start.elapsed().as_secs_f64();
    let is_ce = ce_flag.load(Ordering::SeqCst);

    if is_ce {
        if let Some((subsets, card)) = ce_detail.lock().unwrap().as_ref() {
            println!(
                "\x1b[1;31m🚨🚨🚨 DISPROVED! COUNTEREXAMPLE FOUND IN Z: {:?} with |C| = {} < bound {}! 🚨🚨🚨\x1b[0m",
                subsets, card, bound
            );
            let _ = OpenOptions::new()
                .create(true)
                .append(true)
                .open("discovered_counterexamples.txt")
                .map(|mut f| {
                    writeln!(
                        f,
                        "WANG & SUN CIRCULAR SUMSET COUNTEREXAMPLE in Z: Subsets = {:?}, |C| = {}, Bound = {}",
                        subsets, card, bound
                    )
                });
        }
    }

    SearchResult {
        group_desc: format!("Z [0..{}]", u_max - 1),
        n,
        sizes: sizes.to_vec(),
        bound,
        tested_configs: total_tested.load(Ordering::Relaxed),
        min_cardinality: min_card.load(Ordering::Relaxed),
        elapsed_secs: elapsed,
        counterexample_found: is_ce,
    }
}

// -----------------------------------------------------------------------------
// Stochastic Tabu Local Search Engine for Deep Exploration
// -----------------------------------------------------------------------------

pub fn run_tabu_local_search_z(u_max: u8, sizes: &[usize], num_walkers: usize, steps_per_walker: usize) -> SearchResult {
    let n = sizes.len();
    let bound = compute_conjectured_bound(sizes, n, usize::MAX);
    let universe: Vec<u8> = (0..u_max).collect();

    let total_steps = AtomicU64::new(0);
    let global_min = AtomicUsize::new(usize::MAX);
    let ce_flag = AtomicBool::new(false);
    let start = Instant::now();

    (0..num_walkers).into_par_iter().for_each(|walker_id| {
        let mut rng = FastRng::new(walker_id as u64 * 999_983 + 12345);
        let mut a: Vec<Vec<u8>> = Vec::with_capacity(n);

        // Random initialization
        for &s in sizes {
            let mut set = Vec::with_capacity(s);
            let mut pool = universe.clone();
            for _ in 0..s {
                let idx = rng.gen_range(0, pool.len());
                set.push(pool.swap_remove(idx));
            }
            set.sort_unstable();
            a.push(set);
        }

        let eval_fn = |sets: &[Vec<u8>]| -> usize {
            match n {
                4 => eval_circular_z_4(&sets[0], &sets[1], &sets[2], &sets[3]) as usize,
                5 => eval_circular_z_5(&sets[0], &sets[1], &sets[2], &sets[3], &sets[4]) as usize,
                6 => eval_circular_z_6(&sets[0], &sets[1], &sets[2], &sets[3], &sets[4], &sets[5]) as usize,
                _ => 0,
            }
        };

        let mut curr_score = eval_fn(&a);
        let mut best_score = curr_score;

        for _ in 0..steps_per_walker {
            if ce_flag.load(Ordering::Relaxed) {
                break;
            }

            // Pick a set to mutate
            let set_idx = rng.gen_range(0, n);
            let elem_idx = rng.gen_range(0, a[set_idx].len());
            let old_val = a[set_idx][elem_idx];

            // Pick a new element not in a[set_idx]
            let available: Vec<u8> = universe.iter().cloned().filter(|x| !a[set_idx].contains(x)).collect();
            if available.is_empty() {
                continue;
            }
            let new_val = available[rng.gen_range(0, available.len())];

            a[set_idx][elem_idx] = new_val;
            a[set_idx].sort_unstable();

            let new_score = eval_fn(&a);

            if new_score < bound {
                ce_flag.store(true, Ordering::SeqCst);
                println!(
                    "\x1b[1;31m🚨 [TABU LOCAL SEARCH] DISPROVED IN Z: {:?} with |C| = {} < bound {}!\x1b[0m",
                    a, new_score, bound
                );
                break;
            }

            // Simulated annealing acceptance
            if new_score <= curr_score || (rng.gen_range(0, 100) < 5) {
                curr_score = new_score;
                if new_score < best_score {
                    best_score = new_score;
                    global_min.fetch_min(best_score, Ordering::Relaxed);
                }
            } else {
                // Revert
                a[set_idx][elem_idx] = old_val;
                a[set_idx].sort_unstable();
            }
        }

        total_steps.fetch_add(steps_per_walker as u64, Ordering::Relaxed);
        global_min.fetch_min(best_score, Ordering::Relaxed);
    });

    let elapsed = start.elapsed().as_secs_f64();
    SearchResult {
        group_desc: format!("Z [0..{}] (Tabu Local Search)", u_max - 1),
        n,
        sizes: sizes.to_vec(),
        bound,
        tested_configs: total_steps.load(Ordering::Relaxed),
        min_cardinality: global_min.load(Ordering::Relaxed),
        elapsed_secs: elapsed,
        counterexample_found: ce_flag.load(Ordering::SeqCst),
    }
}

// -----------------------------------------------------------------------------
// Main Hunter Orchestrator
// -----------------------------------------------------------------------------

pub fn hunt_wang_sun_circular() {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF 2: WANG & SUN (EJC 2026, #P2.64) CIRCULAR RESTRICTED SUMSET BREAKER\x1b[0m");
    println!("   Conjecture: |C(A_1, ..., A_n)| >= min{{ p(G), sum |A_i| - 2n + (-1)^n (1 + {{n}}_2) }}");
    println!("   Targeting Asymmetric Configurations across n = 4, 5, 6 in G = Z and G = Z_p");
    println!("   Engine: Rayon 8-Core Parallel Lock-Free Bitmask Vectorizer & Tabu Annealer");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let mut all_results: Vec<SearchResult> = Vec::new();

    // -------------------------------------------------------------------------
    // Phase 1: Finite Cyclic Groups Z_p (p = 5, 7, 11, 13)
    // -------------------------------------------------------------------------
    println!("\x1b[1;33m>>> PHASE 1: EXHAUSTIVE SCAN ON FINITE CYCLIC GROUPS Z_p (p in {{5, 7, 11, 13}}) <<<\x1b[0m\n");

    let zp_targets: Vec<(u8, Vec<usize>)> = vec![
        // n = 4
        (5, vec![3, 2, 3, 2]),
        (7, vec![3, 2, 3, 2]),
        (7, vec![4, 2, 4, 2]),
        (7, vec![4, 2, 3, 2]),
        (7, vec![4, 3, 4, 2]),
        (7, vec![3, 3, 3, 3]),
        (11, vec![3, 2, 3, 2]),
        (11, vec![4, 2, 4, 2]),
        // n = 5
        (5, vec![3, 2, 3, 2, 3]),
        (7, vec![3, 2, 3, 2, 3]),
        (7, vec![4, 2, 4, 2, 3]),
        (7, vec![3, 3, 3, 3, 3]),
        // n = 6
        (5, vec![3, 2, 3, 2, 3, 2]),
        (7, vec![3, 2, 3, 2, 3, 2]),
    ];

    for (p, sizes) in zp_targets {
        let res = scan_exhaustive_zp(p, &sizes);
        let speed = (res.tested_configs as f64) / res.elapsed_secs.max(1e-6) / 1e6;
        println!(
            "  [{}] n={} sizes={:?} | Bound={:>2} | MinObs={:>2} | Tested={:>10} | Speed={:>6.2} M/s | Time={:>6.3}s | {}",
            res.group_desc,
            res.n,
            res.sizes,
            res.bound,
            res.min_cardinality,
            res.tested_configs,
            speed,
            res.elapsed_secs,
            if res.counterexample_found { "\x1b[1;31m🚨 DISPROVED!\x1b[0m" } else { "\x1b[1;32m✓ SHARP/PASSED\x1b[0m" }
        );
        all_results.push(res);
    }

    // -------------------------------------------------------------------------
    // Phase 2: Infinite Torsion-Free Group Z (Integers)
    // -------------------------------------------------------------------------
    println!("\n\x1b[1;33m>>> PHASE 2: EXHAUSTIVE SCAN IN Z (INTEGERS) <<<\x1b[0m\n");

    let z_targets: Vec<(u8, Vec<usize>)> = vec![
        // n = 4
        (8, vec![3, 2, 3, 2]),
        (8, vec![4, 2, 4, 2]),
        (8, vec![4, 2, 3, 2]),
        (8, vec![4, 3, 4, 2]),
        (8, vec![3, 3, 3, 3]),
        (9, vec![4, 3, 4, 3]),
        // n = 5
        (7, vec![3, 2, 3, 2, 3]),
        (8, vec![3, 2, 3, 2, 3]),
        (7, vec![4, 2, 4, 2, 3]),
        (7, vec![3, 3, 3, 3, 3]),
        // n = 6
        (7, vec![3, 2, 3, 2, 3, 2]),
    ];

    for (u_max, sizes) in z_targets {
        let res = scan_exhaustive_z(u_max, &sizes);
        let speed = (res.tested_configs as f64) / res.elapsed_secs.max(1e-6) / 1e6;
        println!(
            "  [{}] n={} sizes={:?} | Bound={:>2} | MinObs={:>2} | Tested={:>10} | Speed={:>6.2} M/s | Time={:>6.3}s | {}",
            res.group_desc,
            res.n,
            res.sizes,
            res.bound,
            res.min_cardinality,
            res.tested_configs,
            speed,
            res.elapsed_secs,
            if res.counterexample_found { "\x1b[1;31m🚨 DISPROVED!\x1b[0m" } else { "\x1b[1;32m✓ SHARP/PASSED\x1b[0m" }
        );
        all_results.push(res);
    }

    // -------------------------------------------------------------------------
    // Phase 3: High-Depth Tabu Simulated Annealing Search (U = 12..14)
    // -------------------------------------------------------------------------
    println!("\n\x1b[1;33m>>> PHASE 3: STOCHASTIC TABU / SIMULATED ANNEALING DEEP EXPLORATION (U = 12..14) <<<\x1b[0m\n");

    let tabu_targets: Vec<(u8, Vec<usize>)> = vec![
        (13, vec![4, 3, 4, 3]),
        (13, vec![4, 3, 4, 2]),
        (13, vec![4, 2, 4, 2, 4]),
        (13, vec![4, 3, 4, 3, 4]),
        (14, vec![3, 2, 3, 2, 3, 2]),
        (14, vec![4, 2, 4, 2, 4, 2]),
        (14, vec![4, 3, 4, 3, 4, 3]),
    ];

    for (u_max, sizes) in tabu_targets {
        let res = run_tabu_local_search_z(u_max, &sizes, 64, 50_000);
        let speed = (res.tested_configs as f64) / res.elapsed_secs.max(1e-6) / 1e6;
        println!(
            "  [{}] n={} sizes={:?} | Bound={:>2} | MinObs={:>2} | Steps={:>11} | Speed={:>6.2} M/s | Time={:>6.3}s | {}",
            res.group_desc,
            res.n,
            res.sizes,
            res.bound,
            res.min_cardinality,
            res.tested_configs,
            speed,
            res.elapsed_secs,
            if res.counterexample_found { "\x1b[1;31m🚨 DISPROVED!\x1b[0m" } else { "\x1b[1;32m✓ BOUND RESPECTED\x1b[0m" }
        );
        all_results.push(res);
    }

    // -------------------------------------------------------------------------
    // Summary
    // -------------------------------------------------------------------------
    let total_configs: u64 = all_results.iter().map(|r| r.tested_configs).sum();
    let total_time: f64 = all_results.iter().map(|r| r.elapsed_secs).sum();
    let ce_count = all_results.iter().filter(|r| r.counterexample_found).count();

    println!("\n\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 WANG & SUN CIRCULAR SUMSET HUNT COMPLETED");
    println!("   Total Configurations / Steps Tested : {}", total_configs);
    println!("   Total Execution Time                : {:.2} seconds", total_time);
    println!("   Average Testing Throughput          : {:.2} M configs/s", total_configs as f64 / total_time.max(1e-6) / 1e6);
    println!("   Counterexamples Verified            : {}", ce_count);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bound_formula() {
        // n = 4, sizes (3, 2, 3, 2), sum = 10, linear = 10 - 8 + 1 = 3
        assert_eq!(compute_conjectured_bound(&[3, 2, 3, 2], 4, usize::MAX), 3);
        assert_eq!(compute_conjectured_bound(&[3, 2, 3, 2], 4, 5), 3);

        // n = 5, sizes (3, 2, 3, 2, 3), sum = 13, linear = 13 - 10 - 2 = 1
        assert_eq!(compute_conjectured_bound(&[3, 2, 3, 2, 3], 5, usize::MAX), 1);
        assert_eq!(compute_conjectured_bound(&[3, 2, 3, 2, 3], 5, 5), 1);

        // n = 6, sizes (3, 2, 3, 2, 3, 2), sum = 15, linear = 15 - 12 + 1 = 4
        assert_eq!(compute_conjectured_bound(&[3, 2, 3, 2, 3, 2], 6, usize::MAX), 4);
    }

    #[test]
    fn test_exhaustive_small_zp() {
        let res = scan_exhaustive_zp(5, &[3, 2, 3, 2]);
        assert_eq!(res.bound, 3);
        assert!(res.min_cardinality >= res.bound);
        assert!(!res.counterexample_found);
    }
}
