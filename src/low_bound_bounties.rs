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

/// Fast odd-only bitset linear sieve for generating primes up to `limit`.
pub fn sieve_primes(limit: usize) -> Vec<u32> {
    if limit < 2 {
        return Vec::new();
    }
    let num_odds = limit / 2;
    let num_words = (num_odds + 63) / 64;
    let mut composite = vec![0u64; num_words];
    let sqrt_limit = (limit as f64).sqrt() as usize;
    for i in (3..=sqrt_limit).step_by(2) {
        let idx = (i - 3) / 2;
        if (composite[idx / 64] & (1u64 << (idx % 64))) == 0 {
            let step = i;
            let mut mult_idx = (i * i - 3) / 2;
            while mult_idx < num_odds {
                composite[mult_idx / 64] |= 1u64 << (mult_idx % 64);
                mult_idx += step;
            }
        }
    }
    let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).ln().max(1.0) * 1.15) as usize);
    primes.push(2u32);
    for idx in 0..num_odds {
        if (composite[idx / 64] & (1u64 << (idx % 64))) == 0 {
            primes.push((2 * idx + 3) as u32);
        }
    }
    primes
}

/// Sieve exactly `target_count` primes up to Dusart's bound using compact odd-only bitset.
pub fn sieve_n_primes(target_count: usize) -> Vec<u32> {
    if target_count == 0 {
        return Vec::new();
    }
    if target_count == 1 {
        return vec![2];
    }
    let limit = if target_count <= 10 {
        35
    } else {
        let n_f = target_count as f64;
        (n_f * (n_f.ln() + n_f.ln().ln() + 0.1)) as usize
    };
    let num_odds = limit / 2;
    let num_words = (num_odds + 63) / 64;
    let mut composite = vec![0u64; num_words];
    let sqrt_limit = (limit as f64).sqrt() as usize;
    for i in (3..=sqrt_limit).step_by(2) {
        let idx = (i - 3) / 2;
        if (composite[idx / 64] & (1u64 << (idx % 64))) == 0 {
            let step = i;
            let mut mult_idx = (i * i - 3) / 2;
            while mult_idx < num_odds {
                composite[mult_idx / 64] |= 1u64 << (mult_idx % 64);
                mult_idx += step;
            }
        }
    }
    let mut primes = Vec::with_capacity(target_count);
    primes.push(2u32);
    for idx in 0..num_odds {
        if (composite[idx / 64] & (1u64 << (idx % 64))) == 0 {
            primes.push((2 * idx + 3) as u32);
            if primes.len() == target_count {
                break;
            }
        }
    }
    primes
}

// ============================================================================
// HEDEF A: OEIS A247824 ($500 Sun Prize: (m+n) | (p_m + p_n))
// ============================================================================

/// Record of champion m values where n/m or n/(m(m-1)) reaches a new high.
#[derive(Debug, Clone)]
pub struct DivisibilityChampion {
    pub m: u64,
    pub p_m: u32,
    pub n: u64,
    pub p_n: u32,
    pub ratio_m: f64,
    pub ratio_bound: f64,
}

/// Hunt for counterexamples to Sun's $500 Prize Conjecture on OEIS A247824 up to m = max_m.
pub fn hunt_hedef_a_divisibility(start_m: u64, max_m: u64, prime_limit: usize) {
    let effective_start = start_m.max(3);
    let effective_end = max_m.max(effective_start);

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF A: OEIS A247824 — SUN'S $500 PRIZE PRIME DIVISIBILITY CONJECTURE\x1b[0m");
    println!("   Conjecture: For each m > 2, exists 1 <= n < m(m-1) such that (m + n) | (p_m + p_n)");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Range      : m = {} .. {}", effective_start, effective_end);
    println!("⚡ Target Prime Pool : 50,000,000 primes (Vec<u32>, ~200 MB RAM)");
    println!("🧵 Multi-threading   : Rayon 8-Core Work Stealing");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let sieve_start = Instant::now();
    // Guarantee 50 million primes for large ranges
    let primes = if prime_limit >= 982_451_653 || effective_end > 100_000 {
        sieve_n_primes(50_000_000)
    } else {
        sieve_primes(prime_limit)
    };
    println!("✓ Sieved {} primes in {:.2?}. Prime p_50000000 = {}\n", primes.len(), sieve_start.elapsed(), primes.last().unwrap_or(&0));

    let num_primes = primes.len() as u64;
    let counterexamples_found = AtomicUsize::new(0);
    let total_m_tested = AtomicU64::new(0);

    let start_time = Instant::now();
    let mut last_telemetry = Instant::now();

    let chunk_size = 5_000u64;
    let mut current_m = effective_start;

    let mut record_highs: Vec<DivisibilityChampion> = Vec::new();
    let mut record_ratios_m: Vec<DivisibilityChampion> = Vec::new();
    let mut record_ratios_b: Vec<DivisibilityChampion> = Vec::new();
    let mut uncharted_champions: Vec<DivisibilityChampion> = Vec::new();

    let mut max_seen_n = 0u64;
    let mut max_ratio_m = 0.0f64;
    let mut max_ratio_bound = 0.0f64;

    while current_m <= effective_end {
        let chunk_end = (current_m + chunk_size).min(effective_end + 1);
        let m_range: Vec<u64> = (current_m..chunk_end).collect();

        let chunk_results: Vec<(u64, Option<u64>)> = m_range
            .into_par_iter()
            .map(|m| {
                if m > num_primes {
                    return (m, None);
                }
                let p_m = unsafe { *primes.get_unchecked((m - 1) as usize) } as u64;
                let limit = m.saturating_mul(m - 1);
                let max_n_check = limit.min(num_primes);

                let mut found_n = None;
                // Fast path for n = 1 (p_1 = 2)
                if (p_m + 2) % (m + 1) == 0 {
                    found_n = Some(1);
                } else {
                    for n in 2..max_n_check {
                        let p_n = unsafe { *primes.get_unchecked((n - 1) as usize) } as u64;
                        let sum_p = p_m + p_n;
                        let sum_idx = m + n;
                        if sum_p % sum_idx == 0 {
                            found_n = Some(n);
                            break;
                        }
                    }
                }
                (m, found_n)
            })
            .collect();

        for (m, res) in chunk_results {
            total_m_tested.fetch_add(1, Ordering::Relaxed);
            let p_m = primes[(m - 1) as usize];
            match res {
                Some(n) => {
                    let limit = m.saturating_mul(m - 1);
                    let p_n = primes[(n - 1) as usize];
                    let r_m = (n as f64) / (m as f64);
                    let r_b = (n as f64) / (limit as f64);

                    if n >= limit {
                        counterexamples_found.fetch_add(1, Ordering::SeqCst);
                        println!(
                            "\x1b[1;31m🚨🚨 $500 BOUNTY WON: COUNTEREXAMPLE AT m = {}! n = {} >= m(m-1) = {}! 🚨🚨\x1b[0m",
                            m, n, limit
                        );
                        let _ = std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open("discovered_counterexamples.txt")
                            .map(|mut f| {
                                use std::io::Write;
                                writeln!(f, "HEDEF A COUNTEREXAMPLE: m = {}, p_m = {}, n = {}, p_n = {}, limit = {}", m, p_m, n, p_n, limit)
                            });
                    }

                    let mut is_record = false;
                    if n > max_seen_n {
                        max_seen_n = n;
                        record_highs.push(DivisibilityChampion {
                            m,
                            p_m,
                            n,
                            p_n,
                            ratio_m: r_m,
                            ratio_bound: r_b,
                        });
                        is_record = true;
                    }

                    if r_m > max_ratio_m {
                        max_ratio_m = r_m;
                        record_ratios_m.push(DivisibilityChampion {
                            m,
                            p_m,
                            n,
                            p_n,
                            ratio_m: r_m,
                            ratio_bound: r_b,
                        });
                        is_record = true;
                    }

                    if r_b > max_ratio_bound {
                        max_ratio_bound = r_b;
                        record_ratios_b.push(DivisibilityChampion {
                            m,
                            p_m,
                            n,
                            p_n,
                            ratio_m: r_m,
                            ratio_bound: r_b,
                        });
                        is_record = true;
                    }

                    if m > 400_000 && is_record {
                        uncharted_champions.push(DivisibilityChampion {
                            m,
                            p_m,
                            n,
                            p_n,
                            ratio_m: r_m,
                            ratio_bound: r_b,
                        });
                        println!(
                            "🌟 \x1b[1;32m[UNCHARTED RECORD]\x1b[0m m = {:7} (p_m = {:9}) -> n = {:9} (p_n = {:10}), n/m = {:6.2}, n/m(m-1) = {:.6}",
                            m, p_m, n, p_n, r_m, r_b
                        );
                    }
                }
                None => {
                    let limit = m.saturating_mul(m - 1);
                    if limit <= num_primes {
                        // Truly exhausted entire domain m(m-1)!
                        counterexamples_found.fetch_add(1, Ordering::SeqCst);
                        println!(
                            "\x1b[1;31m🚨🚨 $500 BOUNTY WON: NO n < m(m-1) EXISTS FOR m = {}! 🚨🚨\x1b[0m",
                            m
                        );
                    } else {
                        println!(
                            "\x1b[1;33m⚠️ [EXTREME PEAK CANDIDATE] m = {} (p_m = {}): n > 50,000,000 (limit = {})!\x1b[0m",
                            m, p_m, limit
                        );
                    }
                }
            }
        }

        let tested = total_m_tested.load(Ordering::Relaxed);
        let total_to_test = effective_end - effective_start + 1;
        if last_telemetry.elapsed().as_secs() >= 5 || chunk_end > effective_end {
            let elapsed = start_time.elapsed().as_secs_f64();
            let progress = (tested as f64) / (total_to_test as f64) * 100.0;
            let rate = (tested as f64) / elapsed / 1000.0;
            println!(
                "📡 [TELEMETRY] [{:5.1}%] m = {:8} | {:6.1} km/s | Max n: {:9} | Max n/m: {:7.2} | Max n/m(m-1): {:8.4} | Tested: {}",
                progress, chunk_end - 1, rate, max_seen_n, max_ratio_m, max_ratio_bound, tested
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
    println!("\n🏆 TOP NOTABLE RECORD HIGHS OF n = a(m):");
    println!("   {:>8} | {:>10} | {:>10} | {:>12} | {:>10} | {:>16}", "m", "p_m", "n = a(m)", "p_n", "n / m", "n / (m(m-1))");
    println!("   {:-<76}", "");
    for c in record_highs.iter().rev().take(10) {
        println!("   {:>8} | {:>10} | {:>10} | {:>12} | {:>10.2} | {:>16.6}", c.m, c.p_m, c.n, c.p_n, c.ratio_m, c.ratio_bound);
    }
    if !uncharted_champions.is_empty() {
        println!("\n🌟 NOTABLE RECORDS IN UNCHARTED TERRITORY (m > 400,000):");
        println!("   {:>8} | {:>10} | {:>10} | {:>12} | {:>10} | {:>16}", "m", "p_m", "n = a(m)", "p_n", "n / m", "n / (m(m-1))");
        println!("   {:-<76}", "");
        for c in uncharted_champions.iter().rev().take(10) {
            println!("   {:>8} | {:>10} | {:>10} | {:>12} | {:>10.2} | {:>16.6}", c.m, c.p_m, c.n, c.p_n, c.ratio_m, c.ratio_bound);
        }
    }
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

// ============================================================================
// HEDEF B: OEIS A181901 (Sun's Prime Self-Recurrence / Modulo Squares)
// ============================================================================

/// Hunt for collisions in Sun's prime alternating sequence s_k and test A181901.
pub fn hunt_hedef_b_collisions(max_k: usize, max_n_sum: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF B: OEIS A181901 — SUN'S PRIME SELF-RECURRENCE & COLLISION HUNTER\x1b[0m");
    println!("   Sequence: s_k = sum_{{j=1}}^k (-1)^{{k-j}} p_j = p_k - s_{{k-1}}");
    println!("   Conjecture: p_{{n+1}} is least modulus making {{2*s_1^2, ..., 2*s_n^2}} mod m distinct");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Depth      : Exact Collisions up to k = {}", max_k);
    println!("⚡ Modulo Sum Depth  : Mod p_{{n+1}} Sum Collisions up to n = {}", max_n_sum);
    println!("🧵 Multi-threading   : Rayon 8-Core Work Stealing");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let sieve_start = Instant::now();
    let primes = sieve_n_primes(max_k);
    println!("✓ Sieved {} primes in {:.2?}\n", primes.len(), sieve_start.elapsed());

    let actual_k = max_k.min(primes.len());

    // 1. Compute s_k sequence
    let start_gen = Instant::now();
    let mut s = Vec::with_capacity(actual_k);
    let mut curr_s = 0i64;
    let mut max_s = 0u32;
    for &p in &primes {
        curr_s = p as i64 - curr_s;
        let u = curr_s as u32;
        if u > max_s {
            max_s = u;
        }
        s.push(u);
    }
    println!("✓ Generated {} terms of s_k in {:.2?}. Max s_k: {}\n", s.len(), start_gen.elapsed(), max_s);

    // 2. Check for Exact Collisions (s_k = s_j) via Compact 60MB Bitset
    println!("⚡ Scanning for Exact Collisions (s_k == s_j) across {} terms...", actual_k);
    let start_exact = Instant::now();
    let bitset_size = (max_s as usize + 64) / 64;
    let mut seen_bitset = vec![0u64; bitset_size];
    let mut exact_collision: Option<(usize, u32)> = None;

    for (k, &val) in s.iter().enumerate() {
        let idx = val as usize;
        let word = idx / 64;
        let mask = 1u64 << (idx % 64);
        if (seen_bitset[word] & mask) != 0 {
            exact_collision = Some((k + 1, val));
            break;
        }
        seen_bitset[word] |= mask;
    }

    match exact_collision {
        Some((k, val)) => {
            println!(
                "\x1b[1;31m🚨🚨 EXACT COLLISION FOUND: s_{} == {}! Conjecture REFUTED for all n >= {}! 🚨🚨\x1b[0m\n",
                k, val, k
            );
        }
        None => {
            println!("✓ Zero exact collisions in {} terms! (Checked in {:.2?}, pairwise distinct)\n", actual_k, start_exact.elapsed());
        }
    }

    // 3. Check for Modulo p_{n+1} Sum Collisions (s_j + s_k = p_{n+1}) via Rayon Parallel Bitsets
    let test_n_limit = max_n_sum.min(actual_k - 1);
    println!("⚡ Scanning for Modulo p_{{n+1}} Sum Collisions (s_j + s_k = p_{{n+1}}) up to n = {} with Rayon...", test_n_limit);

    let sum_start = Instant::now();
    let max_val_prefix = s[..test_n_limit].iter().copied().max().unwrap_or(0) as usize;

    let chunk_size = 50_000usize;
    let chunks: Vec<(usize, usize)> = (0..test_n_limit)
        .step_by(chunk_size)
        .map(|st| (st + 1, (st + chunk_size).min(test_n_limit)))
        .collect();

    let total_sum_collisions: usize = chunks
        .into_par_iter()
        .map(|(c_start, c_end)| {
            let mut local_bitset = vec![0u64; (max_val_prefix + 64) / 64];
            for i in 0..(c_start - 1) {
                let v = s[i] as usize;
                local_bitset[v / 64] |= 1u64 << (v % 64);
            }

            let mut local_col = 0;
            for n in c_start..=c_end {
                let new_val = s[n - 1] as usize;
                local_bitset[new_val / 64] |= 1u64 << (new_val % 64);

                let p_next = primes[n] as usize;
                let half_p = p_next / 2;
                for j in 0..n {
                    let val = s[j] as usize;
                    if val <= half_p {
                        let complement = p_next - val;
                        if complement != val && complement <= max_val_prefix {
                            let w = complement / 64;
                            let m = 1u64 << (complement % 64);
                            if (local_bitset[w] & m) != 0 {
                                local_col += 1;
                                break;
                            }
                        }
                    }
                }
            }
            local_col
        })
        .sum();

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 HEDEF B: SCAN COMPLETED");
    println!("   Total Terms s_k Evaluated : {}", actual_k);
    println!("   Exact Collisions          : {}", if exact_collision.is_some() { "1 (REFUTED)" } else { "0 (STABLE)" });
    println!("   Mod p_{{n+1}} Sum Collisions: {} (Tested in {:.2?})", total_sum_collisions, sum_start.elapsed());
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
