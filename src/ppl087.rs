//! HEDEF 3: PPL 087 / OEIS A185636 — Zhi-Wei Sun's $100 Triple Simultaneous Prime Conjecture
//!
//! Official Problem Statement (PPL 087 / Sun 2011):
//! For every integer n > 1456, there exists 0 <= k < n such that:
//! n - k, n + k, and n + k^2 are ALL prime numbers.
//!
//! Known exceptions below 2000:
//! [16, 76, 166, 316, 341, 361, 411, 481, 556, 656, 766, 866, 1456].
//! The last known exception is exactly n = 1456.

use crate::verifier::is_prime_fast;
use rayon::prelude::*;
use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

/// Check if a given n has a valid k in [0, n) making (n - k), (n + k), (n + k^2) all prime.
/// Returns Some((k, checks_done)) if found, or None if no valid k exists (counterexample!).
#[inline(always)]
pub fn find_triple_prime_k(n: u64) -> Option<(u64, usize)> {
    let k_start = if n % 2 == 0 { 1 } else { 0 };
    let k_step = 2;
    let n_mod3 = n % 3;

    let mut checks = 0usize;
    let mut k = k_start;

    while k < n {
        // If 3 does not divide n, then (n-k) and (n+k) cannot both be primes unless k is a multiple of 3,
        // or one of them equals 3.
        if n_mod3 != 0 && k % 3 != 0 && (n - k != 3) && (n + k != 3) {
            k += k_step;
            continue;
        }

        checks += 1;
        let p1 = n - k;
        let p2 = n + k;
        let p3 = n + k * k;

        if is_prime_fast(p1) && is_prime_fast(p2) && is_prime_fast(p3) {
            return Some((k, checks));
        }

        k += k_step;
    }

    None
}

/// Configuration for PPL 087 hunt.
pub struct Ppl087Config {
    pub start_n: u64,
    pub end_n: u64,
    pub telemetry_interval_secs: u64,
    pub output_file: String,
}

/// Statistics for PPL 087 hunt.
pub struct Ppl087Stats {
    pub n_tested: AtomicU64,
    pub k_checks_total: AtomicU64,
    pub max_k_checked: AtomicU64,
    pub worst_n: AtomicU64,
    pub counterexamples_found: AtomicUsize,
    pub halt_flag: AtomicBool,
}

/// Execute multi-core parallel hunt for counterexamples to Sun's Triple Prime Conjecture (PPL 087).
pub fn hunt_ppl087(config: Ppl087Config) -> bool {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 PPL 087: ZHI-WEI SUN $100 TRIPLE SIMULTANEOUS PRIME HUNTER\x1b[0m");
    println!("   Conjecture: For every n > 1456, exists k in [0, n) with (n-k), (n+k), (n+k^2) all prime");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1m🔍 Search Space      :\x1b[0m {} .. {}", config.start_n, config.end_n);
    println!("\x1b[1m⚡ Known Last Bound  :\x1b[0m n = 1456 (All exceptions < 1457 verified)");
    println!("\x1b[1m🧵 Threading Engine  :\x1b[0m Rayon lock-free chunked work stealing");
    println!("\x1b[1m📁 Proof Log File    :\x1b[0m {}", config.output_file);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let stats = Ppl087Stats {
        n_tested: AtomicU64::new(0),
        k_checks_total: AtomicU64::new(0),
        max_k_checked: AtomicU64::new(0),
        worst_n: AtomicU64::new(0),
        counterexamples_found: AtomicUsize::new(0),
        halt_flag: AtomicBool::new(false),
    };

    let start_time = Instant::now();
    let total_n = config.end_n.saturating_sub(config.start_n);

    // Telemetry monitor thread
    let telemetry_interval = config.telemetry_interval_secs;
    let start_n = config.start_n;
    let end_n = config.end_n;

    // Process in chunks of 50,000 for work balancing and telemetry
    let chunk_size = 50_000u64;
    let num_chunks = (total_n + chunk_size - 1) / chunk_size;

    let mut last_telemetry = Instant::now();

    for chunk_idx in 0..num_chunks {
        let chunk_start = start_n + chunk_idx * chunk_size;
        let chunk_end = (chunk_start + chunk_size).min(end_n);

        let chunk_range: Vec<u64> = (chunk_start..chunk_end).collect();

        chunk_range.into_par_iter().for_each(|n| {
            if stats.halt_flag.load(Ordering::Relaxed) {
                return;
            }

            match find_triple_prime_k(n) {
                Some((k, checks)) => {
                    stats.n_tested.fetch_add(1, Ordering::Relaxed);
                    stats.k_checks_total.fetch_add(checks as u64, Ordering::Relaxed);

                    // Update worst n
                    let mut current_max = stats.max_k_checked.load(Ordering::Relaxed);
                    while (k as u64) > current_max {
                        match stats.max_k_checked.compare_exchange_weak(
                            current_max,
                            k as u64,
                            Ordering::Relaxed,
                            Ordering::Relaxed,
                        ) {
                            Ok(_) => {
                                stats.worst_n.store(n, Ordering::Relaxed);
                                break;
                            }
                            Err(actual) => current_max = actual,
                        }
                    }
                }
                None => {
                    // Counterexample discovered!
                    stats.counterexamples_found.fetch_add(1, Ordering::SeqCst);
                    println!(
                        "\x1b[1;31m🚨🚨 OFFICIAL $100 BOUNTY COUNTEREXAMPLE FOUND: n = {} has NO valid k! 🚨🚨\x1b[0m",
                        n
                    );
                    log_ppl087_counterexample(&config.output_file, n);
                    stats.halt_flag.store(true, Ordering::SeqCst);
                }
            }
        });

        if stats.halt_flag.load(Ordering::SeqCst) {
            break;
        }

        if last_telemetry.elapsed().as_secs() >= telemetry_interval || chunk_idx == num_chunks - 1 {
            let elapsed = start_time.elapsed().as_secs_f64();
            let tested = stats.n_tested.load(Ordering::Relaxed);
            let rate = (tested as f64) / elapsed / 1_000.0;
            let progress = (tested as f64) / (total_n as f64) * 100.0;
            let max_k = stats.max_k_checked.load(Ordering::Relaxed);
            let worst = stats.worst_n.load(Ordering::Relaxed);

            println!(
                "📡 [TELEMETRY] [{:5.1}%] n = {:10} | {:7.1} knums/s | Max k = {:6} (at n={:10}) | Total Tested: {}",
                progress, chunk_end, rate, max_k, worst, tested
            );
            last_telemetry = Instant::now();
        }
    }

    let total_elapsed = start_time.elapsed().as_secs_f64();
    let total_tested = stats.n_tested.load(Ordering::Relaxed);
    let total_checks = stats.k_checks_total.load(Ordering::Relaxed);
    let avg_checks = if total_tested > 0 { (total_checks as f64) / (total_tested as f64) } else { 0.0 };
    let max_k = stats.max_k_checked.load(Ordering::Relaxed);
    let worst = stats.worst_n.load(Ordering::Relaxed);
    let found = stats.counterexamples_found.load(Ordering::SeqCst) > 0;

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🏁 PPL 087 HUNT COMPLETED\x1b[0m");
    println!("   Total Integers n Tested   : {}", total_tested);
    println!("   Average Checks per n      : {:.2}", avg_checks);
    println!("   Maximum k Required        : {} (at n = {})", max_k, worst);
    println!("   Counterexamples Discovered: {}", if found { "\x1b[1;31m1 (VERIFIED)\x1b[0m" } else { "0 (Conjecture Holds)" });
    println!("   Total Time Elapsed        : {:.2} seconds", total_elapsed);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    found
}

fn log_ppl087_counterexample(path: &str, n: u64) {
    if let Ok(mut f) = OpenOptions::new().create(true).append(true).open(path) {
        let _ = writeln!(
            f,
            "=== OFFICIAL PPL 087 ($100 PRIZE) COUNTEREXAMPLE ===\nNumber n: {}\nSearch Domain: 0 <= k < {}\nCondition: No k exists such that (n-k), (n+k), (n+k^2) are all prime.\nTimestamp: {:?}\n====================================================\n",
            n, n, Instant::now()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_known_exceptions() {
        // Known exceptions below 1457:
        let exceptions = [16, 76, 166, 316, 341, 361, 411, 481, 556, 656, 766, 866, 1456];
        for &e in &exceptions {
            assert!(find_triple_prime_k(e).is_none(), "Expected {} to be an exception!", e);
        }

        // Test that 1457 has a solution:
        let sol = find_triple_prime_k(1457);
        assert!(sol.is_some(), "Expected 1457 to have a solution!");
        let (k, _) = sol.unwrap();
        assert!(is_prime_fast(1457 - k));
        assert!(is_prime_fast(1457 + k));
        assert!(is_prime_fast(1457 + k * k));
    }
}
