//! GÖREV 2: Low-Bound and Finite-Ceiling Prize Bounties Hunter
//!
//! Investigates 4 fragile, highly-targeted official prize conjectures:
//! 1. Hedef A: OEIS A247824 — Zhi-Wei Sun's $500 Prize Conjecture on (m+n) | (p_m + p_n) for n < m(m-1).
//! 2. Hedef B: OEIS A181901 — Zhi-Wei Sun's Prime Self-Recurrence / Modulo Squares Distinctness.
//! 3. Hedef C: Sun Conjecture 2.22(ii) — Inverted Square Prime Counting Sieve (pi(pi(k*n)) = x^2).
//! 4. Hedef D: Sun's Finite Abelian Group Cycle Conjecture (|G| >= 27) with distinct adjacent sums.

use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

/// Fast linear sieve for generating all primes up to `limit`.
pub fn sieve_primes(limit: usize) -> Vec<u32> {
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;
    let sqrt = (limit as f64).sqrt() as usize;
    for i in 2..=sqrt {
        if is_prime[i] {
            let mut j = i * i;
            while j <= limit {
                is_prime[j] = false;
                j += i;
            }
        }
    }
    is_prime
        .iter()
        .enumerate()
        .filter_map(|(idx, &p)| if p { Some(idx as u32) } else { None })
        .collect()
}

// ============================================================================
// HEDEF A: OEIS A247824 ($500 Sun Prize: (m+n) | (p_m + p_n))
// ============================================================================

/// Record of champion m values where n/m or n/(m(m-1)) reaches a new high.
#[derive(Debug, Clone)]
pub struct DivisibilityChampion {
    pub m: u64,
    pub n: u64,
    pub ratio_m: f64,
    pub ratio_bound: f64,
}

/// Hunt for counterexamples to Sun's $500 Prize Conjecture on OEIS A247824.
pub fn hunt_hedef_a_divisibility(max_m: u64, prime_limit: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF A: OEIS A247824 — SUN'S $500 PRIZE PRIME DIVISIBILITY CONJECTURE\x1b[0m");
    println!("   Conjecture: For each m > 2, exists 1 <= n < m(m-1) such that (m + n) | (p_m + p_n)");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Range      : m = 3 .. {}", max_m);
    println!("⚡ Sieve Prime Pool  : {} primes (up to {})", prime_limit, prime_limit);
    println!("🧵 Multi-threading   : Rayon 8-Core Work Stealing");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let sieve_start = Instant::now();
    let primes = sieve_primes(prime_limit);
    println!("✓ Sieve completed: {} primes in {:.2?}\n", primes.len(), sieve_start.elapsed());

    let num_primes = primes.len() as u64;
    let counterexamples_found = AtomicUsize::new(0);
    let total_m_tested = AtomicU64::new(0);

    let start_time = Instant::now();
    let mut last_telemetry = Instant::now();

    // Scan m in smaller chunks for responsive telemetry
    let chunk_size = 5_000u64;
    let mut current_m = 3u64;

    let mut top_champions: Vec<DivisibilityChampion> = Vec::new();
    let mut max_ratio_m = 0.0f64;
    let mut max_ratio_bound = 0.0f64;

    while current_m <= max_m {
        let chunk_end = (current_m + chunk_size).min(max_m + 1);
        let m_range: Vec<u64> = (current_m..chunk_end).collect();

        let chunk_results: Vec<(u64, Option<u64>)> = m_range
            .into_par_iter()
            .map(|m| {
                if m > num_primes {
                    return (m, None);
                }
                let p_m = primes[(m - 1) as usize] as u64;
                let limit = m.saturating_mul(m - 1);
                let max_n_check = limit.min(num_primes);

                let mut found_n = None;
                for n in 1..max_n_check {
                    let p_n = primes[(n - 1) as usize] as u64;
                    if (p_m + p_n) % (m + n) == 0 {
                        found_n = Some(n);
                        break;
                    }
                }
                (m, found_n)
            })
            .collect();

        for (m, res) in chunk_results {
            total_m_tested.fetch_add(1, Ordering::Relaxed);
            match res {
                Some(n) => {
                    let limit = m * (m - 1);
                    if n >= limit {
                        counterexamples_found.fetch_add(1, Ordering::SeqCst);
                        println!(
                            "\x1b[1;31m🚨🚨 $500 BOUNTY WON: COUNTEREXAMPLE AT m = {}! n = {} >= m(m-1) = {}! 🚨🚨\x1b[0m",
                            m, n, limit
                        );
                    } else {
                        let r_m = (n as f64) / (m as f64);
                        let r_b = (n as f64) / (limit as f64);
                        if r_m > max_ratio_m || r_b > max_ratio_bound {
                            if r_m > max_ratio_m {
                                max_ratio_m = r_m;
                            }
                            if r_b > max_ratio_bound {
                                max_ratio_bound = r_b;
                            }
                            top_champions.push(DivisibilityChampion {
                                m,
                                n,
                                ratio_m: r_m,
                                ratio_bound: r_b,
                            });
                        }
                    }
                }
                None => {
                    if m * (m - 1) <= num_primes {
                        // Truly exhausted entire domain m(m-1)!
                        counterexamples_found.fetch_add(1, Ordering::SeqCst);
                        println!(
                            "\x1b[1;31m🚨🚨 $500 BOUNTY WON: NO n < m(m-1) EXISTS FOR m = {}! 🚨🚨\x1b[0m",
                            m
                        );
                    }
                }
            }
        }

        let tested = total_m_tested.load(Ordering::Relaxed);
        if last_telemetry.elapsed().as_secs() >= 5 || chunk_end > max_m {
            let elapsed = start_time.elapsed().as_secs_f64();
            let progress = (tested as f64) / (max_m as f64) * 100.0;
            let rate = (tested as f64) / elapsed / 1000.0;
            println!(
                "📡 [TELEMETRY] [{:5.1}%] m = {:8} | {:6.1} km/s | Max n/m: {:7.2} | Max n/m(m-1): {:8.4} | Tested: {}",
                progress, chunk_end - 1, rate, max_ratio_m, max_ratio_bound, tested
            );
            last_telemetry = Instant::now();
        }

        current_m = chunk_end;
    }

    let elapsed = start_time.elapsed().as_secs_f64();
    let tested = total_m_tested.load(Ordering::Relaxed);

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 HEDEF A: SCAN COMPLETED");
    println!("   Total m Tested            : {}", tested);
    println!("   Time Elapsed              : {:.2}s ({:.1} km/s)", elapsed, (tested as f64) / elapsed / 1000.0);
    println!("   Counterexamples Found     : {}", counterexamples_found.load(Ordering::SeqCst));
    println!("\n🏆 TOP NOTABLE CHAMPIONS (Extreme n/m and n/m(m-1) Ratios):");
    println!("   {:>8} | {:>10} | {:>12} | {:>16}", "m", "n = a(m)", "n / m", "n / (m(m-1))");
    println!("   {:-<54}", "");
    for c in top_champions.iter().rev().take(10) {
        println!("   {:>8} | {:>10} | {:>12.2} | {:>16.6}", c.m, c.n, c.ratio_m, c.ratio_bound);
    }
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

// ============================================================================
// HEDEF B: OEIS A181901 (Sun's Prime Self-Recurrence / Modulo Squares)
// ============================================================================

/// Hunt for collisions in Sun's prime alternating sequence s_k and test A181901.
pub fn hunt_hedef_b_collisions(max_k: usize, prime_limit: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF B: OEIS A181901 — SUN'S PRIME SELF-RECURRENCE & COLLISION HUNTER\x1b[0m");
    println!("   Sequence: s_k = sum_{{j=1}}^k (-1)^{{k-j}} p_j = p_k - s_{{k-1}}");
    println!("   Conjecture: p_{{n+1}} is least modulus making {{2*s_1^2, ..., 2*s_n^2}} mod m distinct");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Depth      : k = 1 .. {}", max_k);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let sieve_start = Instant::now();
    let primes = sieve_primes(prime_limit);
    println!("✓ Sieved {} primes in {:.2?}\n", primes.len(), sieve_start.elapsed());

    let actual_k = max_k.min(primes.len());

    // 1. Compute s_k sequence
    let start_gen = Instant::now();
    let mut s = Vec::with_capacity(actual_k);
    let mut curr_s = 0i64;
    for k in 0..actual_k {
        let p_k = primes[k] as i64;
        curr_s = p_k - curr_s;
        s.push(curr_s);
    }
    println!("✓ Generated {} terms of s_k in {:.2?}\n", s.len(), start_gen.elapsed());

    // 2. Check for Exact Collisions (s_k = s_j)
    println!("⚡ Scanning for Exact Collisions (s_k == s_j)...");
    let max_val = *s.iter().max().unwrap_or(&0) as usize;
    let mut seen = vec![0u32; max_val + 1];
    let mut exact_collision: Option<(usize, usize, i64)> = None;

    for (k, &val) in s.iter().enumerate() {
        let v = val as usize;
        if seen[v] != 0 {
            exact_collision = Some((seen[v] as usize, k + 1, val));
            break;
        }
        seen[v] = (k + 1) as u32;
    }

    match exact_collision {
        Some((j, k, val)) => {
            println!(
                "\x1b[1;31m🚨🚨 EXACT COLLISION FOUND: s_{} == s_{} == {}! Conjecture REFUTED for all n >= {}! 🚨🚨\x1b[0m\n",
                j, k, val, k
            );
        }
        None => {
            println!("✓ Zero exact collisions in {} terms! (All s_k values are pairwise distinct)\n", actual_k);
        }
    }

    // 3. Check for Modulo p_{n+1} Sum Collisions (s_j + s_k = p_{n+1})
    let test_n_limit = actual_k.min(50_000);
    println!("⚡ Scanning for Modulo p_{{n+1}} Sum Collisions (s_j + s_k = p_{{n+1}}) up to n = {}...", test_n_limit);

    let sum_collisions = AtomicUsize::new(0);

    // Use two-pointer approach on sorted prefix of s
    let mut sorted_s = Vec::with_capacity(test_n_limit);
    for n in 1..test_n_limit {
        let new_elem = s[n - 1];
        let ins_pos = sorted_s.binary_search(&new_elem).unwrap_or_else(|e| e);
        sorted_s.insert(ins_pos, new_elem);

        let p_next = primes[n] as i64;

        // Two-pointer check if any two distinct elements sum to p_next
        let mut left = 0usize;
        let mut right = sorted_s.len().saturating_sub(1);
        let mut found = false;

        while left < right {
            let sum = sorted_s[left] + sorted_s[right];
            if sum == p_next {
                found = true;
                sum_collisions.fetch_add(1, Ordering::Relaxed);
                break;
            } else if sum < p_next {
                left += 1;
            } else {
                right -= 1;
            }
        }

        if found {
            println!(
                "\x1b[1;33m⚠️ Sum Collision at n = {}: s_j + s_k = p_{{n+1}} = {}\x1b[0m",
                n, p_next
            );
        }
    }

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 HEDEF B: SCAN COMPLETED");
    println!("   Total Terms s_k Evaluated : {}", actual_k);
    println!("   Exact Collisions          : {}", if exact_collision.is_some() { "1 (REFUTED)" } else { "0 (STABLE)" });
    println!("   Mod p_{{n+1}} Sum Collisions: {}", sum_collisions.load(Ordering::Relaxed));
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

// ============================================================================
// HEDEF C: Sun Conjecture 2.22(ii) (Inverted Square Prime-Counting Sieve)
// ============================================================================

/// Interval I_x = [L_x, R_x] where pi(pi(m)) = x^2.
#[derive(Debug, Clone, Copy)]
pub struct PiPiSquareInterval {
    pub x: u64,
    pub l: u64,
    pub r: u64,
}

/// Hunt for counterexamples to Sun's Conjecture 2.22(ii) using inverted interval hitting.
pub fn hunt_hedef_c_square_pipi(max_n: u64, prime_limit: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF C: SUN CONJECTURE 2.22(ii) — INVERTED SQUARE PRIME-COUNTING SIEVE\x1b[0m");
    println!("   Conjecture: For each n >= 1, exists 1 <= k <= n such that pi(pi(k*n)) = x^2");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Range      : n = 1 .. {}", max_n);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let sieve_start = Instant::now();
    let primes = sieve_primes(prime_limit);
    println!("✓ Sieved {} primes in {:.2?}", primes.len(), sieve_start.elapsed());

    // Build intervals [L_x, R_x] where pi(pi(m)) = x^2:
    // x = 0: pi(pi(1)) = 0 = 0^2, pi(pi(2)) = 0 = 0^2 => [1, 2]
    let mut intervals = vec![PiPiSquareInterval { x: 0, l: 1, r: 2 }];

    let num_primes = primes.len();
    let mut x = 1u64;

    while (x * x + 1) as usize <= num_primes {
        let idx1 = (x * x) as usize;
        let idx2 = (x * x + 1) as usize;

        let p_sq = primes[idx1 - 1] as usize;
        let p_sq_plus_1 = primes[idx2 - 1] as usize;

        if p_sq_plus_1 > num_primes {
            break;
        }

        let l_x = primes[p_sq - 1] as u64;
        let r_x = (primes[p_sq_plus_1 - 1] - 1) as u64;

        intervals.push(PiPiSquareInterval { x, l: l_x, r: r_x });
        x += 1;
    }

    let max_interval_r = intervals.last().map(|i| i.r).unwrap_or(0);
    println!("✓ Generated {} inverted intervals I_x up to R = {}\n", intervals.len(), max_interval_r);

    let start_time = Instant::now();
    let counterexamples_found = AtomicUsize::new(0);
    let total_tested = AtomicU64::new(0);

    // Parallel scan over n up to max_n
    let chunk_size = 50_000u64;
    let mut current_n = 1u64;

    while current_n <= max_n {
        let chunk_end = (current_n + chunk_size).min(max_n + 1);
        let n_range: Vec<u64> = (current_n..chunk_end).collect();

        n_range.into_par_iter().for_each(|n| {
            total_tested.fetch_add(1, Ordering::Relaxed);
            let n_sq = n.saturating_mul(n);

            let mut hit = false;
            for iv in &intervals {
                if iv.l > n_sq && iv.l > max_interval_r {
                    break;
                }
                // Check if any multiple k*n falls in [iv.l, iv.r] with 1 <= k <= n
                let q_r = iv.r / n;
                let q_l_minus_1 = (iv.l - 1) / n;
                if q_r > q_l_minus_1 {
                    let k_min = q_l_minus_1 + 1;
                    if k_min <= n {
                        hit = true;
                        break;
                    }
                }
            }

            if !hit && n_sq <= max_interval_r {
                counterexamples_found.fetch_add(1, Ordering::SeqCst);
                println!(
                    "\x1b[1;31m🚨🚨 OFFICIAL COUNTEREXAMPLE FOUND TO CONJECTURE 2.22(ii): n = {} has no valid k! 🚨🚨\x1b[0m",
                    n
                );
            }
        });

        current_n = chunk_end;
    }


    let elapsed = start_time.elapsed().as_secs_f64();
    let tested = total_tested.load(Ordering::Relaxed);

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 HEDEF C: SCAN COMPLETED");
    println!("   Total Integers n Tested   : {}", tested);
    println!("   Time Elapsed              : {:.2}s ({:.1} knums/s)", elapsed, (tested as f64) / elapsed / 1000.0);
    println!("   Counterexamples Found     : {}", counterexamples_found.load(Ordering::SeqCst));
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

// ============================================================================
// HEDEF D: Sun's Finite Abelian Group Cycle Conjecture (|G| >= 27)
// ============================================================================

/// Trait representing an abelian group of odd order.
pub trait OddAbelianGroup: Sync {
    fn order(&self) -> usize;
    fn add(&self, a: usize, b: usize) -> usize;
    fn name(&self) -> &str;
}

/// Elementary abelian group (Z_3)^3 of order 27.
pub struct Z3Cubed;
impl OddAbelianGroup for Z3Cubed {
    fn order(&self) -> usize { 27 }
    fn add(&self, a: usize, b: usize) -> usize {
        let (a0, a1, a2) = (a % 3, (a / 3) % 3, a / 9);
        let (b0, b1, b2) = (b % 3, (b / 3) % 3, b / 9);
        let c0 = (a0 + b0) % 3;
        let c1 = (a1 + b1) % 3;
        let c2 = (a2 + b2) % 3;
        c0 + c1 * 3 + c2 * 9
    }
    fn name(&self) -> &str { "(Z_3)^3 [Order 27]" }
}

/// Cyclic group Z_N of odd order N.
pub struct ZOdd {
    pub n: usize,
    pub label: String,
}
impl ZOdd {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            label: format!("Z_{}", n),
        }
    }
}
impl OddAbelianGroup for ZOdd {
    fn order(&self) -> usize { self.n }
    fn add(&self, a: usize, b: usize) -> usize { (a + b) % self.n }
    fn name(&self) -> &str { &self.label }
}

/// Find a Hamiltonian cycle in subset A of group G with pairwise distinct adjacent sums.
pub fn find_distinct_adjacent_sum_cycle(group: &dyn OddAbelianGroup, subset: &[usize]) -> bool {
    let k = subset.len();
    if k <= 2 {
        return true;
    }

    let mut path = vec![0usize; k];
    let used_elem_mask = 1u32; // bit 0 (subset[0]) is fixed as starting vertex
    path[0] = subset[0];

    // Track used sum values in group
    let mut used_sums = vec![false; group.order()];

    fn backtrack(
        step: usize,
        k: usize,
        subset: &[usize],
        path: &mut [usize],
        used_elem_mask: u32,
        used_sums: &mut [bool],
        group: &dyn OddAbelianGroup,
    ) -> bool {
        if step == k {
            // Check closing edge from path[k-1] to path[0]
            let closing_sum = group.add(path[k - 1], path[0]);
            return !used_sums[closing_sum];
        }

        let prev = path[step - 1];

        for i in 1..k {
            if (used_elem_mask & (1 << i)) == 0 {
                let candidate = subset[i];
                let adj_sum = group.add(prev, candidate);
                if !used_sums[adj_sum] {
                    used_sums[adj_sum] = true;
                    path[step] = candidate;
                    if backtrack(
                        step + 1,
                        k,
                        subset,
                        path,
                        used_elem_mask | (1 << i),
                        used_sums,
                        group,
                    ) {
                        return true;
                    }
                    used_sums[adj_sum] = false;
                }
            }
        }
        false
    }

    backtrack(1, k, subset, &mut path, used_elem_mask, &mut used_sums, group)
}

/// Hunt for counterexample subsets to Sun's Abelian Group Cycle Conjecture.
pub fn hunt_hedef_d_group_cycles(num_samples_per_size: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF D: SUN'S FINITE ABELIAN GROUP CYCLE CONJECTURE (|G| >= 27)\x1b[0m");
    println!("   Conjecture: Every subset A of odd-order abelian group G admits a cyclic permutation");
    println!("               with pairwise distinct adjacent sums a_i + a_{{i+1}}");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let groups: Vec<Box<dyn OddAbelianGroup>> = vec![
        Box::new(Z3Cubed),
        Box::new(ZOdd::new(31)),
        Box::new(ZOdd::new(33)),
    ];

    for g in groups {
        println!("⚡ Testing Group: \x1b[1;33m{}\x1b[0m", g.name());
        let ord = g.order();
        let mut failed_subsets = 0usize;
        let mut total_tested = 0usize;

        // Test subset sizes k from 4 up to min(14, ord)
        for k in 4..=14.min(ord) {
            // Generate deterministic combinations & pseudorandom subsets
            let mut rng_state = (k as u64) * 6364136223846793005 + 1;
            for _ in 0..num_samples_per_size {
                total_tested += 1;
                // Knuth shuffle to select k distinct elements
                let mut pool: Vec<usize> = (0..ord).collect();
                for i in 0..k {
                    rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let j = i + (rng_state as usize) % (ord - i);
                    pool.swap(i, j);
                }
                let subset = &pool[..k];

                if !find_distinct_adjacent_sum_cycle(g.as_ref(), subset) {
                    failed_subsets += 1;
                    println!(
                        "\x1b[1;31m🚨 COUNTEREXAMPLE SUBSET FOUND in {}: {:?} (size {}) has NO valid cycle!\x1b[0m",
                        g.name(), subset, k
                    );
                    break;
                }
            }
        }

        println!(
            "   ✓ Subsets Tested: {} | Failed (Counterexamples): {}\n",
            total_tested, failed_subsets
        );
    }

    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hedef_a_small() {
        let primes = sieve_primes(5000);
        for m in 3..=10 {
            let p_m = primes[m - 1] as u64;
            let m_u = m as u64;
            let max_n = m_u * (m_u - 1);
            let mut found = false;
            for n in 1..max_n {
                let p_n = primes[n as usize - 1] as u64;
                if (p_m + p_n) % (m_u + n) == 0 {
                    found = true;
                    break;
                }
            }
            assert!(found, "Failed for m = {}", m);
        }
    }

    #[test]
    fn test_hedef_b_s_sequence() {
        let primes = sieve_primes(100);
        let mut s = Vec::new();
        let mut curr = 0i64;
        for &p in primes.iter().take(10) {
            curr = p as i64 - curr;
            s.push(curr);
        }
        // p = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        // s1 = 2
        // s2 = 3 - 2 = 1
        // s3 = 5 - 1 = 4
        // s4 = 7 - 4 = 3
        // s5 = 11 - 3 = 8
        // s6 = 13 - 8 = 5
        // s7 = 17 - 5 = 12
        // s8 = 19 - 12 = 7
        // s9 = 23 - 7 = 16
        // s10 = 29 - 16 = 13
        assert_eq!(s, vec![2, 1, 4, 3, 8, 5, 12, 7, 16, 13]);
    }

    #[test]
    fn test_hedef_d_z3_cubed_cycle() {
        let group = Z3Cubed;
        let subset = vec![0, 1, 3, 9, 13];
        let has_cycle = find_distinct_adjacent_sum_cycle(&group, &subset);
        assert!(has_cycle);
    }
}
