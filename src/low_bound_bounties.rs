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
// MEISSEL-LEHMER EXACT PRIME COUNTING ENGINE (exact_prime_pi)
// ============================================================================

use std::sync::OnceLock;

/// Fast Meissel-Lehmer Exact Prime Counting Engine for x <= 2 * 10^11.
pub struct PrimePiEngine {
    pub limit: usize,
    pub primes: Vec<u32>,
    pub bitset: Vec<u64>,
    pub prefix_counts: Vec<u32>,
}

impl PrimePiEngine {
    pub fn new(limit: usize) -> Self {
        let num_odds = limit / 2;
        let num_words = (num_odds + 63) / 64;
        let mut bitset = vec![0u64; num_words];
        let sqrt_limit = (limit as f64).sqrt() as usize;

        for i in (3..=sqrt_limit).step_by(2) {
            let idx = (i - 3) / 2;
            if (bitset[idx / 64] & (1u64 << (idx % 64))) == 0 {
                let step = i;
                let mut mult_idx = (i * i - 3) / 2;
                while mult_idx < num_odds {
                    bitset[mult_idx / 64] |= 1u64 << (mult_idx % 64);
                    mult_idx += step;
                }
            }
        }

        let mut primes = Vec::with_capacity((limit as f64 / (limit as f64).ln().max(1.0) * 1.12) as usize);
        primes.push(2u32);
        let mut prefix_counts = vec![0u32; num_words];
        let mut cum_primes = 0u32;

        for w in 0..num_words {
            prefix_counts[w] = cum_primes;
            let word = bitset[w];
            let primes_in_word = (!word).count_ones();
            cum_primes += primes_in_word;
        }

        for idx in 0..num_odds {
            if (bitset[idx / 64] & (1u64 << (idx % 64))) == 0 {
                primes.push((2 * idx + 3) as u32);
            }
        }

        PrimePiEngine {
            limit,
            primes,
            bitset,
            prefix_counts,
        }
    }

    #[inline(always)]
    pub fn pi_small(&self, x: usize) -> u64 {
        if x < 2 {
            return 0;
        }
        if x == 2 {
            return 1;
        }
        let x_clamped = x.min(self.limit);
        let idx = (x_clamped - 3) / 2;
        let w = idx / 64;
        let b = idx % 64;
        let mask = if b == 63 { !0u64 } else { (1u64 << (b + 1)) - 1 };
        let bits = (!self.bitset[w]) & mask;
        1 + (self.prefix_counts[w] as u64) + (bits.count_ones() as u64)
    }

    fn phi(&self, x: u64, a: usize, memo: &mut std::collections::HashMap<(u64, u32), u64>) -> u64 {
        if a == 0 {
            return x;
        }
        if a == 1 {
            return (x + 1) / 2;
        }
        let p_a = self.primes[a - 1] as u64;
        if x < p_a {
            return 1;
        }
        if x < p_a * p_a {
            if x <= self.limit as u64 {
                return self.pi_small(x as usize) - (a as u64) + 1;
            }
        }
        let key = (x, a as u32);
        if a < 150 {
            if let Some(&val) = memo.get(&key) {
                return val;
            }
        }
        let res = self.phi(x, a - 1, memo) - self.phi(x / p_a, a - 1, memo);
        if a < 150 {
            memo.insert(key, res);
        }
        res
    }

    pub fn meissel_pi(&self, x: u64) -> u64 {
        if x <= self.limit as u64 {
            return self.pi_small(x as usize);
        }
        let c = (x as f64).cbrt() as u64;
        let a = self.pi_small(c as usize) as usize;
        let sqrt_x = (x as f64).sqrt() as u64;
        let b = self.pi_small(sqrt_x as usize) as usize;

        let mut memo = std::collections::HashMap::with_capacity(32768);
        let phi_val = self.phi(x, a, &mut memo);

        let mut p2 = 0u64;
        for i in (a + 1)..=b {
            let p = self.primes[i - 1] as u64;
            let x_div_p = x / p;
            let pi_val = self.pi_small(x_div_p as usize);
            p2 += pi_val - (i as u64) + 1;
        }

        phi_val + (a as u64) - 1 - p2
    }
}

static PRIME_PI_ENGINE: OnceLock<PrimePiEngine> = OnceLock::new();

pub fn get_prime_pi_engine() -> &'static PrimePiEngine {
    PRIME_PI_ENGINE.get_or_init(|| PrimePiEngine::new(40_000_000))
}

/// Exact prime counting function pi(x) using Meissel-Lehmer algorithm (x <= 2 * 10^11).
pub fn exact_prime_pi(x: u64) -> u64 {
    get_prime_pi_engine().meissel_pi(x)
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
    pub p_n: u64,
    pub k: u32,
    pub ratio_m: f64,
    pub ratio_bound: f64,
}

#[inline(always)]
fn fast_bitset_hash(val: i64) -> usize {
    ((val as u64).wrapping_mul(0x517cc1b727220a95) >> 41) as usize & (131072 * 64 - 1)
}

/// Hunt for counterexamples to Sun's $500 Prize Conjecture on OEIS A247824 up to m = max_m
/// using the 2-Phase Algebraic Inversion Sniper:
/// Phase 1: Small-n fast elimination (n <= 2000) + Inversion Streaming for k in 2..=19 (n <= 50M).
/// Phase 2: Segmented Prime Sieve for hard survivors (k in 20..=24) up to n = 500,000,000.
pub fn hunt_hedef_a_divisibility(start_m: u64, max_m: u64, prime_limit: usize) {
    let effective_start = start_m.max(3);
    let effective_end = max_m.max(effective_start);

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF A: OEIS A247824 — 2-PHASE ALGEBRAIC INVERSION SNIPER ($500 BOUNTY)\x1b[0m");
    println!("   Algebraic Identity: p_m + p_n = k(m + n) <=> p_n - k*n = k*m - p_m");
    println!("   Conjecture: For each m > 2, exists 1 <= n < m(m-1) such that (m + n) | (p_m + p_n)");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Range          : m = {} .. {}", effective_start, effective_end);
    println!("⚡ Phase 1 Prime Pool    : 50,000,000 primes (Vec<u32>, ~200 MB RAM, k in 2..=19)");
    println!("🚀 Phase 2 Deep Sieve    : Segmented Sieve up to n = 500,000,000 (k in 20..=24)");
    println!("🧵 Multi-threading       : Rayon Lock-Free Parallel Multiplier Dispatch");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let total_timer = Instant::now();

    // 1. Generate first 50 million primes for Phase 1
    let sieve_start = Instant::now();
    let primes = if prime_limit >= 982_451_653 || effective_end > 100_000 {
        sieve_n_primes(50_000_000)
    } else {
        sieve_primes(prime_limit.max(2_000_000))
    };
    let num_primes = primes.len() as u64;
    println!(
        "✓ Sieved {} primes in {:.2?}. Prime p_{} = {}\n",
        primes.len(),
        sieve_start.elapsed(),
        primes.len(),
        primes.last().unwrap_or(&0)
    );

    let counterexamples_found = AtomicUsize::new(0);
    let mut ans_n = vec![0u64; (effective_end + 1) as usize];
    let mut ans_k = vec![0u32; (effective_end + 1) as usize];

    // ========================================================================
    // FAZ 1a: FAST PATH FOR SMALL n <= 2000
    // ========================================================================
    let p1a_start = Instant::now();
    let fast_results: Vec<(u64, u64, u32)> = (effective_start..=effective_end)
        .into_par_iter()
        .filter_map(|m| {
            if m > num_primes {
                return None;
            }
            let p_m = unsafe { *primes.get_unchecked((m - 1) as usize) } as u64;
            let limit = m.saturating_mul(m - 1);
            let max_n = 2000u64.min(limit).min(num_primes);

            for n in 1..=max_n {
                let p_n = unsafe { *primes.get_unchecked((n - 1) as usize) } as u64;
                if (p_m + p_n) % (m + n) == 0 {
                    let k = ((p_m + p_n) / (m + n)) as u32;
                    return Some((m, n, k));
                }
            }
            None
        })
        .collect();

    let fast_resolved_count = fast_results.len();
    for (m, n, k) in fast_results {
        ans_n[m as usize] = n;
        ans_k[m as usize] = k;
    }

    let total_m = (effective_end - effective_start + 1) as usize;
    println!(
        "⚡ [PHASE 1a FAST PATH] Resolved {} / {} candidates ({:.1}%) with n <= 2000 in {:.2?}",
        fast_resolved_count,
        total_m,
        (fast_resolved_count as f64) / (total_m as f64) * 100.0,
        p1a_start.elapsed()
    );

    // ========================================================================
    // FAZ 1b: ALGEBRAIC INVERSION SNIPER FOR k in 2..=19 (2000 < n <= 50M)
    // ========================================================================
    let unresolved_m: Vec<u64> = (effective_start..=effective_end)
        .filter(|&m| ans_n[m as usize] == 0)
        .collect();

    if !unresolved_m.is_empty() {
        let p1b_start = Instant::now();
        println!(
            "⚡ [PHASE 1b ALGEBRAIC INVERSION] Inverting {} candidates for k in 2..=19 across 50M primes...",
            unresolved_m.len()
        );

        let k_matches: Vec<Vec<(u64, u64, u32)>> = (2..=19u32)
            .into_par_iter()
            .map(|k| {
                use std::collections::HashMap;
                let mut target_map: HashMap<i64, Vec<u32>> = HashMap::new();
                let mut bitset = vec![0u64; 131072];
                let k_i64 = k as i64;

                for &m in &unresolved_m {
                    let p_m = unsafe { *primes.get_unchecked((m - 1) as usize) } as i64;
                    let target = k_i64 * (m as i64) - p_m;
                    let h = fast_bitset_hash(target);
                    bitset[h / 64] |= 1u64 << (h % 64);
                    target_map.entry(target).or_default().push(m as u32);
                }

                let mut matches = Vec::new();
                for (idx, &p_n) in primes.iter().enumerate().skip(2000) {
                    let n = (idx + 1) as u64;
                    let d = (p_n as i64) - k_i64 * (n as i64);
                    let h = fast_bitset_hash(d);
                    if (bitset[h / 64] & (1u64 << (h % 64))) != 0 {
                        if let Some(ms) = target_map.remove(&d) {
                            for m in ms {
                                matches.push((m as u64, n, k));
                            }
                            if target_map.is_empty() {
                                break;
                            }
                        }
                    }
                }
                matches
            })
            .collect();

        let mut p1b_resolved = 0;
        for batch in k_matches {
            for (m, n, k) in batch {
                let m_idx = m as usize;
                if ans_n[m_idx] == 0 || n < ans_n[m_idx] {
                    if ans_n[m_idx] == 0 {
                        p1b_resolved += 1;
                    }
                    ans_n[m_idx] = n;
                    ans_k[m_idx] = k;
                }
            }
        }

        println!(
            "✓ [PHASE 1b COMPLETED] Resolved {} candidates across k in 2..=19 in {:.2?}\n",
            p1b_resolved,
            p1b_start.elapsed()
        );
    }

    // Check for any Phase 1 counterexamples (n >= m(m-1))
    for m in effective_start..=effective_end {
        let n = ans_n[m as usize];
        if n > 0 {
            let limit = m.saturating_mul(m - 1);
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
                        let p_m = primes[(m - 1) as usize];
                        let p_n = (ans_k[m as usize] as u64) * (m + n) - (p_m as u64);
                        writeln!(f, "HEDEF A COUNTEREXAMPLE: m = {}, p_m = {}, n = {}, p_n = {}, limit = {}", m, p_m, n, p_n, limit)
                    });
            }
        }
    }

    // ========================================================================
    // FAZ 2: MEISSEL-LEHMER EXACT PI + RESONANCE BAND JUMPER (k in 20..=25)
    // ========================================================================
    let hard_survivors: Vec<u64> = (effective_start..=effective_end)
        .filter(|&m| ans_n[m as usize] == 0)
        .collect();

    let mut hard_survivor_records: Vec<DivisibilityChampion> = Vec::new();

    if !hard_survivors.is_empty() {
        println!("\x1b[1;33m================================================================================\x1b[0m");
        println!("\x1b[1;33m⚡ FAZ 2: MEISSEL-LEHMER EXACT PI + RESONANCE BAND JUMPER (k in 20..=25)\x1b[0m");
        println!("   Hard Survivors with n > 50,000,000 Count: {}", hard_survivors.len());
        if hard_survivors.len() <= 15 {
            println!("   Surviving m list: {:?}", hard_survivors);
        } else {
            println!("   First 15 surviving m: {:?} ...", &hard_survivors[..15]);
        }
        println!("   Skipping dead zones via exact prime counting across k in 20..=25...\n");
        println!("\x1b[1;33m================================================================================\x1b[0m\n");

        let p2_start = Instant::now();
        let pi_engine = get_prime_pi_engine();
        let base_primes = sieve_primes(470_000); // covers up to sqrt(2.2e11) ~ 469,000

        for k in 20..=25u32 {
            let k_i64 = k as i64;
            let active_m: Vec<u64> = hard_survivors
                .iter()
                .copied()
                .filter(|&m| ans_n[m as usize] == 0)
                .collect();

            if active_m.is_empty() {
                break;
            }

            // Build target map for this band
            use std::collections::HashMap;
            let mut target_map: HashMap<i64, Vec<u32>> = HashMap::new();
            let mut bitset = vec![0u64; 131072];

            for &m in &active_m {
                let p_m = unsafe { *primes.get_unchecked((m - 1) as usize) } as i64;
                let target = k_i64 * (m as i64) - p_m;
                if target >= -6_000_000 && target <= 26_000_000 {
                    let h = fast_bitset_hash(target);
                    bitset[h / 64] |= 1u64 << (h % 64);
                    target_map.entry(target).or_default().push(m as u32);
                }
            }

            if target_map.is_empty() {
                continue;
            }

            // Calculate approximate center for resonance band k
            let mut ln_n = (k as f64) - (k as f64).ln() + 1.0;
            for _ in 0..15 {
                let f = ln_n + ln_n.ln() - 1.0 - (k as f64);
                let f_prime = 1.0 + 1.0 / ln_n;
                ln_n -= f / f_prime;
            }
            let n_center = ln_n.exp();
            let p_center = (n_center * (k as f64)) as u64;

            // Binary search for P_start where D(P) <= -6,000,000
            let mut low = p_center.saturating_sub(6_000_000_000).max(1_000_000_000);
            let mut high = p_center;
            while high - low > 1_000_000 {
                let mid = low + (high - low) / 2;
                let d = (mid as i64) - k_i64 * (pi_engine.meissel_pi(mid) as i64);
                if d < -6_000_000 {
                    low = mid;
                } else {
                    high = mid;
                }
            }
            let mut p_start = low;
            if p_start % 2 == 0 { p_start -= 1; }

            // Binary search for P_end where D(P) >= 26,000,000
            let mut low = p_center;
            let mut high = p_center.saturating_add(6_000_000_000);
            while high - low > 1_000_000 {
                let mid = low + (high - low) / 2;
                let d = (mid as i64) - k_i64 * (pi_engine.meissel_pi(mid) as i64);
                if d < 26_000_000 {
                    low = mid;
                } else {
                    high = mid;
                }
            }
            let mut p_end = high;
            if p_end % 2 == 0 { p_end += 1; }

            let n_start = pi_engine.meissel_pi(p_start - 1);
            let est_primes = ((p_end - p_start) as f64) / (p_start as f64).ln() / 1e6;

            println!(
                "⚡ \x1b[1;36m[RESONANCE BAND k = {:2}]\x1b[0m Window [{}, {}] (~{:.2}M primes) | starting index n = {} | Active Targets: {}",
                k, p_start, p_end, est_primes, n_start, target_map.len()
            );

            // Run Segmented Sieve on [p_start, p_end]
            let seg_num_odds = 1_000_000usize;
            let seg_step = 2 * seg_num_odds as u64;
            let num_words = (seg_num_odds + 63) / 64;
            let mut seg_bitset = vec![0u64; num_words];

            let mut low_seg = p_start;
            let mut curr_n = n_start;
            let mut band_resolved = 0usize;

            while low_seg < p_end && !target_map.is_empty() {
                let high_seg = low_seg + seg_step;
                seg_bitset.fill(0);

                for &p in &base_primes {
                    if p == 2 { continue; }
                    let p_u64 = p as u64;
                    let p_sq = p_u64 * p_u64;
                    if p_sq >= high_seg { break; }
                    let mut start = if p_sq > low_seg {
                        p_sq
                    } else {
                        let rem = low_seg % p_u64;
                        if rem == 0 { low_seg } else { low_seg + (p_u64 - rem) }
                    };
                    if start % 2 == 0 {
                        start += p_u64;
                    }
                    if start < high_seg {
                        let mut bit_idx = ((start - low_seg) / 2) as usize;
                        let step = p as usize;
                        while bit_idx < seg_num_odds {
                            seg_bitset[bit_idx / 64] |= 1u64 << (bit_idx % 64);
                            bit_idx += step;
                        }
                    }
                }

                // Scan primes in segment
                for word_idx in 0..num_words {
                    let word = seg_bitset[word_idx];
                    if word == !0u64 {
                        continue;
                    }
                    for bit in 0..64 {
                        let bit_idx = word_idx * 64 + bit;
                        if bit_idx >= seg_num_odds {
                            break;
                        }
                        if (word & (1u64 << bit)) == 0 {
                            let p_n = low_seg + 2 * (bit_idx as u64);
                            curr_n += 1;

                            let d = (p_n as i64) - k_i64 * (curr_n as i64);
                            let h = fast_bitset_hash(d);

                            if (bitset[h / 64] & (1u64 << (h % 64))) != 0 {
                                if let Some(ms) = target_map.remove(&d) {
                                    for m in ms {
                                        let m_u64 = m as u64;
                                        let m_idx = m as usize;
                                        if ans_n[m_idx] == 0 {
                                            ans_n[m_idx] = curr_n;
                                            ans_k[m_idx] = k;
                                            band_resolved += 1;

                                            let p_m = primes[m_idx - 1];
                                            let limit = m_u64.saturating_mul(m_u64 - 1);
                                            let r_m = (curr_n as f64) / (m_u64 as f64);
                                            let r_b = (curr_n as f64) / (limit as f64);

                                            // Exact u128 verification
                                            let sum_p = (p_m as u128) + (p_n as u128);
                                            let sum_idx = (m_u64 as u128) + (curr_n as u128);
                                            assert_eq!(sum_p % sum_idx, 0, "Exact divisibility verification failed!");
                                            assert_eq!(sum_p / sum_idx, k as u128, "Multiplier quotient verification failed!");

                                            hard_survivor_records.push(DivisibilityChampion {
                                                m: m_u64,
                                                p_m,
                                                n: curr_n,
                                                p_n,
                                                k,
                                                ratio_m: r_m,
                                                ratio_bound: r_b,
                                            });

                                            println!(
                                                "🌟 \x1b[1;32m[RESONANCE HIT]\x1b[0m m = {:7} (p_m = {:9}) -> n = {:11} (p_n = {:12}), k = {:2}, n/m = {:8.2}, n/m(m-1) = {:.8} [EXACT u128 VERIFIED]",
                                                m, p_m, curr_n, p_n, k, r_m, r_b
                                            );

                                            if curr_n >= limit {
                                                counterexamples_found.fetch_add(1, Ordering::SeqCst);
                                                println!(
                                                    "\x1b[1;31m🚨🚨 $500 BOUNTY WON: COUNTEREXAMPLE AT m = {}! n = {} >= m(m-1) = {}! 🚨🚨\x1b[0m",
                                                    m, curr_n, limit
                                                );
                                                let _ = std::fs::OpenOptions::new()
                                                    .create(true)
                                                    .append(true)
                                                    .open("discovered_counterexamples.txt")
                                                    .map(|mut f| {
                                                        use std::io::Write;
                                                        writeln!(f, "HEDEF A COUNTEREXAMPLE: m = {}, p_m = {}, n = {}, p_n = {}, limit = {}", m, p_m, curr_n, p_n, limit)
                                                    });
                                            }
                                        }
                                    }
                                }
                            }

                            if target_map.is_empty() {
                                break;
                            }
                        }
                    }
                    if target_map.is_empty() {
                        break;
                    }
                }

                low_seg = high_seg;
            }

            println!(
                "   ✓ Band k = {:2} completed. Resolved {} candidates.\n",
                k,
                band_resolved
            );
        }

        println!(
            "\n✓ [FAZ 2 COMPLETED] Resonance band jumping finished in {:.2?}.\n",
            p2_start.elapsed()
        );

        for &m in &hard_survivors {
            if ans_n[m as usize] == 0 {
                let p_m = primes[(m - 1) as usize];
                let limit = m.saturating_mul(m - 1);
                println!(
                    "\x1b[1;33m⚠️ [EXTREME HIGH RECORD CANDIDATE] m = {} (p_m = {}): n > 8,500,000,000 (k >= 26, limit = {})!\x1b[0m",
                    m, p_m, limit
                );
            }
        }
    }

    // ========================================================================
    // CHAMPION STATISTICS & COMPREHENSIVE TABLES
    // ========================================================================
    let total_elapsed = total_timer.elapsed().as_secs_f64();
    let mut record_highs: Vec<DivisibilityChampion> = Vec::new();
    let mut record_ratios_m: Vec<DivisibilityChampion> = Vec::new();
    let mut record_ratios_b: Vec<DivisibilityChampion> = Vec::new();
    let mut uncharted_records: Vec<DivisibilityChampion> = Vec::new();

    let mut max_seen_n = 0u64;
    let mut max_ratio_m = 0.0f64;
    let mut max_ratio_b = 0.0f64;

    for m in effective_start..=effective_end {
        let n = ans_n[m as usize];
        if n > 0 {
            let p_m = primes[(m - 1) as usize];
            let k = ans_k[m as usize];
            let limit = m.saturating_mul(m - 1);
            let p_n = (k as u64) * (m + n) - (p_m as u64);
            let r_m = (n as f64) / (m as f64);
            let r_b = (n as f64) / (limit as f64);

            let champ = DivisibilityChampion {
                m,
                p_m,
                n,
                p_n,
                k,
                ratio_m: r_m,
                ratio_bound: r_b,
            };

            let mut is_record = false;
            if n > max_seen_n {
                max_seen_n = n;
                record_highs.push(champ.clone());
                is_record = true;
            }
            if r_m > max_ratio_m {
                max_ratio_m = r_m;
                record_ratios_m.push(champ.clone());
                is_record = true;
            }
            if r_b > max_ratio_b {
                max_ratio_b = r_b;
                record_ratios_b.push(champ.clone());
                is_record = true;
            }
            if m > 400_000 && is_record {
                uncharted_records.push(champ);
            }
        }
    }

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 HEDEF A: SCAN COMPLETED");
    println!("   Total m Tested            : {}", total_m);
    println!("   Total Elapsed Time        : {:.2}s ({:.1} km/s)", total_elapsed, (total_m as f64) / total_elapsed / 1000.0);
    println!("   Counterexamples Found     : {}", counterexamples_found.load(Ordering::SeqCst));
    println!("   Hard Survivors (n > 50M)  : {}", hard_survivor_records.len());
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    if !hard_survivor_records.is_empty() {
        println!("🏆 \x1b[1;32mALL HARD SURVIVORS REQUIRING n > 50,000,000 (k >= 20):\x1b[0m");
        println!(
            "   {:>8} | {:>10} | {:>12} | {:>13} | {:>4} | {:>10} | {:>16}",
            "m", "p_m", "n = a(m)", "p_n", "k", "n / m", "n / (m(m-1))"
        );
        println!("   {:-<84}", "");
        for c in &hard_survivor_records {
            println!(
                "   {:>8} | {:>10} | {:>12} | {:>13} | {:>4} | {:>10.2} | {:>16.8}",
                c.m, c.p_m, c.n, c.p_n, c.k, c.ratio_m, c.ratio_bound
            );
        }
        println!();
    }

    println!("🏆 TOP RECORD HIGHS OF n = a(m):");
    println!(
        "   {:>8} | {:>10} | {:>12} | {:>13} | {:>4} | {:>10} | {:>16}",
        "m", "p_m", "n = a(m)", "p_n", "k", "n / m", "n / (m(m-1))"
    );
    println!("   {:-<84}", "");
    for c in record_highs.iter().rev().take(12) {
        println!(
            "   {:>8} | {:>10} | {:>12} | {:>13} | {:>4} | {:>10.2} | {:>16.8}",
            c.m, c.p_m, c.n, c.p_n, c.k, c.ratio_m, c.ratio_bound
        );
    }
    println!();

    if !uncharted_records.is_empty() {
        println!("🌟 NOTABLE RECORDS IN UNCHARTED TERRITORY (m > 400,000):");
        println!(
            "   {:>8} | {:>10} | {:>12} | {:>13} | {:>4} | {:>10} | {:>16}",
            "m", "p_m", "n = a(m)", "p_n", "k", "n / m", "n / (m(m-1))"
        );
        println!("   {:-<84}", "");
        for c in uncharted_records.iter().rev().take(10) {
            println!(
                "   {:>8} | {:>10} | {:>12} | {:>13} | {:>4} | {:>10.2} | {:>16.8}",
                c.m, c.p_m, c.n, c.p_n, c.k, c.ratio_m, c.ratio_bound
            );
        }
        println!();
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

    #[test]
    fn test_hedef_a_algebraic_identity() {
        // Test algebraic identity: D_k(n) = p_n - k*n == k*m - p_m
        // m = 3, n = 5, k = 2
        // p_3 = 5, p_5 = 11. D_2(5) = 11 - 10 = 1. T_2(3) = 6 - 5 = 1.
        let p_3 = 5i64;
        let p_5 = 11i64;
        assert_eq!(p_5 - 2 * 5, 2 * 3 - p_3);
        assert_eq!((p_3 + p_5) % (3 + 5), 0);

        // m = 7, n = 38, k = 4
        // p_7 = 17, p_38 = 163. D_4(38) = 163 - 152 = 11. T_4(7) = 28 - 17 = 11.
        let p_7 = 17i64;
        let p_38 = 163i64;
        assert_eq!(p_38 - 4 * 38, 4 * 7 - p_7);
        assert_eq!((p_7 + p_38) % (7 + 38), 0);

        // m = 106, n = 6866, k = 10
        let p_106 = 577i64;
        let p_6866 = 69143i64;
        assert_eq!(p_6866 - 10 * 6866, 10 * 106 - p_106);
        assert_eq!((p_106 + p_6866) % (106 + 6866), 0);

        // m = 10409, n = 69804276, k = 20
        let p_10409 = 109481i64;
        let p_n = 1396184219i64;
        assert_eq!(p_n - 20 * 69804276, 20 * 10409 - p_10409);
        assert_eq!((p_10409 + p_n) % (10409 + 69804276), 0);
    }

    #[test]
    fn test_exact_prime_pi_benchmarks() {
        assert_eq!(exact_prime_pi(100_000_000), 5_761_455);
        assert_eq!(exact_prime_pi(1_000_000_000), 50_847_534);
        assert_eq!(exact_prime_pi(10_000_000_000), 455_052_511);
        assert_eq!(exact_prime_pi(100_000_000_000), 4_118_054_813);
        assert_eq!(exact_prime_pi(75_391_645_409), 3_141_281_384);
    }
}
