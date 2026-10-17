//! HEDEF 1: PPL 043 ($1,000 Cash Prize — Zhi-Wei Sun)
//!
//! Investigates Sun's conjecture on representations of positive integers
//! by alternating sums of consecutive primes with ceiling p_n <= 2m + 2.2*sqrt(m).
//!
//! Identity: For consecutive primes p_k, ..., p_n (1 <= k < n):
//!   m = S(k, n) = sum_{j=k}^n (-1)^{n-j} p_j = s_n + (-1)^{n-k} s_{k-1}
//! where s_n = sum_{j=1}^n (-1)^{n-j} p_j (s_0 = 0).

use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

/// Fast odd-only bitset linear sieve for primes up to `limit`.
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

/// Champion record of extreme c(m) = (p_n - 2m) / sqrt(m)
#[derive(Debug, Clone)]
pub struct Ppl043Champion {
    pub m: u64,
    pub p_n: u64,
    pub c: f64,
    pub k: usize,
    pub n: usize,
}

/// Hunt for counterexamples to PPL 043 up to m = max_m.
pub fn hunt_ppl043_breaker(max_m: usize, max_j_depth: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 HEDEF 1: PPL 043 ($1,000 PRIZE) — ZHI-WEI SUN ALTERNATING SUMS BREAKER\x1b[0m");
    println!("   Conjecture: Every m >= 1 can be represented as m = sum_{{j=k}}^n (-1)^{{n-j}} p_j");
    println!("               with consecutive primes p_k, ..., p_n (1 <= k < n) and p_n <= 2m + 2.2*sqrt(m)");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("🔍 Search Range      : m = 1 .. {}", max_m);
    println!("⚡ Exploration Depth : j <= {} (k <= {})", max_j_depth, max_j_depth + 1);
    println!("🧵 Multi-threading   : Rayon 8-Core Work Stealing");
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let sieve_start = Instant::now();
    let limit_p = (2.0 * max_m as f64 + 2.2 * (max_m as f64).sqrt() + 1000.0) as usize;
    let primes = sieve_primes(limit_p);
    println!("✓ Sieved {} primes up to {} in {:.2?}", primes.len(), primes.last().unwrap_or(&0), sieve_start.elapsed());

    let num_primes = primes.len();
    let gen_start = Instant::now();
    let mut s = Vec::with_capacity(num_primes + 1);
    s.push(0i64);
    let mut curr = 0i64;
    for &p in &primes {
        curr = p as i64 - curr;
        s.push(curr);
    }
    println!("✓ Generated s_n sequence (length {}) in {:.2?}\n", s.len(), gen_start.elapsed());

    let words = (max_m + 64) / 64;
    let represented: Vec<AtomicU64> = (0..words).map(|_| AtomicU64::new(0)).collect();

    let scan_start = Instant::now();
    println!("⚡ Executing Parallel Rayon Bitmask Marking across all cores...");

    let total_marked = AtomicUsize::new(0);

    // Parallel marking over n
    (2..=num_primes).into_par_iter().for_each(|n| {
        let p_n = primes[n - 1] as f64;
        let s_n = s[n];

        // Solve 2m + 2.2*sqrt(m) >= p_n for m_min
        let u = (-2.2 + (4.84 + 8.0 * p_n).sqrt()) / 4.0;
        let m_min = (u * u).ceil() as i64;
        if m_min > max_m as i64 {
            return;
        }

        let max_j = (n - 2).min(max_j_depth);
        let mut local_count = 0usize;

        for j in 0..=max_j {
            let s_j = s[j];
            let m = if (n % 2) == (j % 2) {
                s_n - s_j
            } else {
                s_n + s_j
            };

            if m >= m_min && m <= max_m as i64 && m > 0 {
                let m_u = m as usize;
                let word = m_u / 64;
                let mask = 1u64 << (m_u % 64);
                let prev = represented[word].fetch_or(mask, Ordering::Relaxed);
                if (prev & mask) == 0 {
                    local_count += 1;
                }
            }
        }

        if local_count > 0 {
            total_marked.fetch_add(local_count, Ordering::Relaxed);
        }
    });

    let marking_elapsed = scan_start.elapsed();
    println!("✓ Parallel marking completed in {:.2?}\n", marking_elapsed);

    // Scan for unrepresented integers
    println!("🔍 Inspecting bitset for unrepresented integers (counterexample candidates)...");
    let mut unrepresented = Vec::new();
    for m in 1..=max_m {
        let word = m / 64;
        let mask = 1u64 << (m % 64);
        if (represented[word].load(Ordering::Relaxed) & mask) == 0 {
            unrepresented.push(m as u64);
        }
    }

    if !unrepresented.is_empty() {
        println!("\x1b[1;33m⚠️ {} candidate unrepresented integers found in fast pass. Running exhaustive search...\x1b[0m", unrepresented.len());
        let mut true_counterexamples = Vec::new();

        for &m in &unrepresented {
            let m_i = m as i64;
            let bound = 2.0 * m as f64 + 2.2 * (m as f64).sqrt();
            let mut found = false;

            for n in 2..num_primes {
                let p_n = primes[n - 1];
                if (p_n as f64) > bound {
                    break;
                }
                let s_n = s[n];
                for k in 1..n {
                    let s_prev = s[k - 1];
                    let val = if (n - k) % 2 == 0 {
                        s_n + s_prev
                    } else {
                        s_n - s_prev
                    };
                    if val == m_i {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }

            if !found {
                true_counterexamples.push(m);
                println!(
                    "\x1b[1;31m🚨🚨🚨 $1,000 BOUNTY WON! DEFINITIVE COUNTEREXAMPLE FOUND: m = {} HAS NO VALID REPRESENTATION! 🚨🚨🚨\x1b[0m",
                    m
                );
                let _ = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open("discovered_counterexamples.txt")
                    .map(|mut f| {
                        use std::io::Write;
                        writeln!(f, "PPL 043 COUNTEREXAMPLE: m = {}, ceiling = 2m + 2.2*sqrt(m)", m)
                    });
            }
        }

        if true_counterexamples.is_empty() {
            println!("✓ All candidate integers resolved under deeper search! Zero counterexamples.");
        }
    } else {
        println!("✓ ALL {} integers are REPRESENTED! Zero counterexamples detected!", max_m);
    }

    // Benchmark sample champions for c(m)
    println!("\n🏆 COMPUTING CHAMPION c(m) = (p_n - 2m) / sqrt(m) VALUES...");
    let sample_points = [
        10u64, 50, 100, 500, 968, 1_000, 5_000, 10_000, 50_000, 100_000,
        500_000, 1_000_000, 5_000_000, 10_000_000, 25_000_000, 50_000_000,
    ];

    let mut champions: Vec<Ppl043Champion> = Vec::new();
    for &m in &sample_points {
        if m as usize > max_m { continue; }
        let m_i = m as i64;
        let bound = 2.0 * m as f64 + 2.2 * (m as f64).sqrt();

        for n in 2..num_primes {
            let p_n = primes[n - 1] as u64;
            let s_n = s[n];
            let mut hit = false;
            let mut best_k = 0;

            for k in 1..n.min(2000) {
                let s_prev = s[k - 1];
                let val = if (n - k) % 2 == 0 {
                    s_n + s_prev
                } else {
                    s_n - s_prev
                };
                if val == m_i {
                    hit = true;
                    best_k = k;
                    break;
                }
            }

            if hit {
                let c = (p_n as f64 - 2.0 * m as f64) / (m as f64).sqrt();
                champions.push(Ppl043Champion {
                    m,
                    p_n,
                    c,
                    k: best_k,
                    n,
                });
                break;
            }
            if (p_n as f64) > bound {
                break;
            }
        }
    }

    println!("   {:>10} | {:>12} | {:>8} | {:>8} | {:>12} | {:>10}", "m", "p_n", "k", "n", "c(m)", "Ceiling");
    println!("   {:-<70}", "");
    for ch in &champions {
        println!(
            "   {:>10} | {:>12} | {:>8} | {:>8} | {:>12.4} | {:>10.2}",
            ch.m, ch.p_n, ch.k, ch.n, ch.c, 2.20
        );
    }

    let elapsed = scan_start.elapsed().as_secs_f64();
    println!("\n\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 PPL 043: SCAN COMPLETED");
    println!("   Total Integers Verified   : {}", max_m);
    println!("   Processing Speed          : {:.1} Mnums/s", (max_m as f64) / elapsed / 1e6);
    println!("   Max c(m) Observed         : {:.4} (Bound: 2.2000)", champions.iter().map(|c| c.c).fold(0.0f64, f64::max));
    println!("   Counterexamples Found     : {}", unrepresented.len());
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_m_representations() {
        let primes = sieve_primes(5000);
        let mut s = Vec::new();
        s.push(0i64);
        let mut curr = 0i64;
        for &p in &primes {
            curr = p as i64 - curr;
            s.push(curr);
        }

        for m in 1..=50 {
            let m_i = m as i64;
            let bound = 2.0 * m as f64 + 2.2 * (m as f64).sqrt();
            let mut found = false;
            for n in 2..primes.len() {
                let p_n = primes[n - 1] as f64;
                if p_n > bound { break; }
                let s_n = s[n];
                for k in 1..n {
                    let s_prev = s[k - 1];
                    let val = if (n - k) % 2 == 0 { s_n + s_prev } else { s_n - s_prev };
                    if val == m_i {
                        found = true;
                        break;
                    }
                }
                if found { break; }
            }
            assert!(found, "m = {} should be representable within bound!", m);
        }
    }
}
