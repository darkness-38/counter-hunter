//! CRT Covering Sniper Engine
//! Multi-core, lock-free search using Chinese Remainder Theorem covering progressions.
//! Skips 80%-95% of shifts and leaps across 10^14 - 10^16 integer space at tera-scale speed.

use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use crate::conjectures::{PrimeKind, TargetConjecture};
use crate::crt::{is_composite_by_medium_primes, optimize_crt_covering};
use crate::verifier::{
    format_factorization, is_prime_fast, is_sophie_germain, is_twin_prime, smallest_prime_factor,
    CounterexampleReport, DifferenceAnalysis,
};

/// Global atomic counter for lock-free work-stealing across k-indices.
static GLOBAL_K_INDEX: AtomicU64 = AtomicU64::new(0);

/// Live telemetry stats.
#[derive(Default, Debug)]
pub struct TelemetryStats {
    pub k_steps_evaluated: AtomicU64,
    pub current_max_n: AtomicU64,
    pub candidates_tested: AtomicU64,
    pub counterexamples_found: AtomicU64,
    pub active_workers: AtomicU64,
}

/// Hunt engine configuration.
#[derive(Clone, Debug)]
pub struct HuntConfig {
    pub target: TargetConjecture,
    pub start: u64,
    pub end: u64,
    pub max_master_modulus: u64,
    pub num_threads: usize,
    pub telemetry_interval_secs: u64,
    pub output_file: String,
    pub max_counterexamples: usize,
    pub run_preflight: bool,
    pub verbose: bool,
}

impl Default for HuntConfig {
    fn default() -> Self {
        Self {
            target: TargetConjecture::Binom2aPrime,
            start: 100_000_000,
            end: 10_000_000_000_000, // 10^13
            max_master_modulus: 500_000_000_000, // 5 * 10^11
            num_threads: 0,
            telemetry_interval_secs: 10,
            output_file: "discovered_counterexamples.txt".to_string(),
            max_counterexamples: 1,
            run_preflight: true,
            verbose: false,
        }
    }
}

/// The CRT Sniper Hunter Engine.
pub struct Hunter {
    pub config: HuntConfig,
    pub stats: Arc<TelemetryStats>,
    pub counterexamples: Arc<Mutex<Vec<u64>>>,
}

impl Hunter {
    pub fn new(config: HuntConfig) -> Self {
        Self {
            config,
            stats: Arc::new(TelemetryStats::default()),
            counterexamples: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// Executes the full CRT Covering Sniper hunt.
    pub fn hunt(&self) -> bool {
        // Step 1: Pre-Flight Validation
        if self.config.run_preflight {
            println!("\x1b[1;36m>>> RUNNING MANDATORY PRE-FLIGHT VALIDATOR (10,000 .. 100,000) <<<\x1b[0m");
            let report = self.config.target.run_preflight_check(10_000, 100_000);
            println!("   Target                 : {}", self.config.target.name());
            println!("   Tested Candidates      : {}", report.tested_count);
            println!("   Small Exceptions Found : {}", report.exceptions_count);
            if !report.sample_exceptions.is_empty() {
                println!("   Sample Exceptions      : {:?}", report.sample_exceptions);
            }
            println!("   Validator Verdict      : \x1b[1m{}\x1b[0m\n", report.verdict);

            if !report.passed {
                eprintln!("\x1b[1;31m[ABORT] Target failed pre-flight validation due to insufficient density / trivial obstacles.\x1b[0m\n");
                return false;
            }
        }

        // Step 2: Generate All Shifts up to search limit
        let target = self.config.target;
        let prime_kind = target.prime_kind();
        let all_shifts = target.generate_shifts(self.config.end);
        let total_shifts_count = all_shifts.len();

        let (base_mod, base_residues) = target.base_modulus_and_residues();

        // Step 3: Build Optimal CRT Progressions for all base residues
        println!("\x1b[1;33m>>> OPTIMIZING CRT COVERING PROGRESSIONS (Multiplicative Covering) <<<\x1b[0m");
        let mut progressions = Vec::new();

        for &res in &base_residues {
            let crt = optimize_crt_covering(
                &all_shifts,
                res,
                base_mod,
                prime_kind,
                self.config.max_master_modulus,
            );
            println!(
                "   Base Class: n ≡ {} (mod {}) => Master Progression: n ≡ {} (mod {})",
                res, base_mod, crt.master_residue, crt.master_modulus
            );
            println!(
                "   Covering Primes Used : {:?}",
                crt.primes_used
            );
            println!(
                "   Total Shifts: {} | Survivors: {} | \x1b[1;32mKilled: {:.1}%\x1b[0m",
                total_shifts_count,
                crt.survivor_shifts.len(),
                crt.kill_percentage()
            );
            progressions.push(crt);
        }
        println!("\x1b[1;36m================================================================================\x1b[0m\n");

        let num_threads = if self.config.num_threads > 0 {
            self.config.num_threads
        } else {
            std::thread::available_parallelism()
                .map(|p| p.get())
                .unwrap_or(4)
        };

        let start_n = self.config.start;
        let end_n = self.config.end;
        let max_ce = self.config.max_counterexamples as u64;

        // Step 4: Scan Each Master Progression
        for (prog_idx, crt) in progressions.into_iter().enumerate() {
            if max_ce > 0 && self.stats.counterexamples_found.load(Ordering::Relaxed) >= max_ce {
                break;
            }

            println!(
                "\x1b[1;32m🎯 LAUNCHING CRT SNIPER STREAM {}/{} | Step Size M = {}\x1b[0m",
                prog_idx + 1,
                base_residues.len(),
                crt.master_modulus
            );

            let m = crt.master_modulus;
            let r = crt.master_residue;

            let k_start = if start_n > r {
                (start_n - r).div_ceil(m)
            } else {
                0
            };
            let k_end = if end_n >= r {
                (end_n - r) / m
            } else {
                0
            };

            if k_start > k_end {
                continue;
            }

            let total_k_steps = k_end - k_start + 1;
            println!(
                "   Search Range in n-Space : {} .. {}",
                format_number(r + k_start * m),
                format_number(r + k_end * m)
            );
            println!("   Total k Steps to Scan   : {}\n", total_k_steps);

            GLOBAL_K_INDEX.store(k_start, Ordering::SeqCst);
            self.stats.current_max_n.store(r + k_start * m, Ordering::Relaxed);

            let running = Arc::new(AtomicBool::new(true));

            // Spawn Telemetry Thread
            let stats_telemetry = Arc::clone(&self.stats);
            let running_telemetry = Arc::clone(&running);
            let interval = self.config.telemetry_interval_secs;
            let target_short = target.short_code();
            let survivor_count = crt.survivor_shifts.len();

            let telemetry_handle = thread::spawn(move || {
                let mut last_k = 0u64;
                let mut last_time = Instant::now();

                while running_telemetry.load(Ordering::Relaxed) {
                    thread::sleep(Duration::from_secs(interval));
                    if !running_telemetry.load(Ordering::Relaxed) {
                        break;
                    }

                    let now = Instant::now();
                    let current_k = stats_telemetry.k_steps_evaluated.load(Ordering::Relaxed);
                    let delta_k = current_k.saturating_sub(last_k);
                    let delta_t = now.duration_since(last_time).as_secs_f64();

                    let k_per_sec = if delta_t > 0.0 {
                        (delta_k as f64) / delta_t
                    } else {
                        0.0
                    };

                    let current_n = stats_telemetry.current_max_n.load(Ordering::Relaxed);
                    let active = stats_telemetry.active_workers.load(Ordering::Relaxed);
                    let candidates = stats_telemetry.candidates_tested.load(Ordering::Relaxed);
                    let found = stats_telemetry.counterexamples_found.load(Ordering::Relaxed);

                    // Speed in actual integer space scanned per second
                    let n_space_speed = k_per_sec * (m as f64);
                    let n_speed_str = if n_space_speed >= 1e12 {
                        format!("{:.2} Tnums/s", n_space_speed / 1e12)
                    } else if n_space_speed >= 1e9 {
                        format!("{:.2} Gnums/s", n_space_speed / 1e9)
                    } else {
                        format!("{:.2} Mnums/s", n_space_speed / 1e6)
                    };

                    let progress_pct = if total_k_steps > 0 {
                        ((current_k.min(total_k_steps) as f64) / (total_k_steps as f64)) * 100.0
                    } else {
                        100.0
                    };

                    println!(
                        "\x1b[1;36m[CRT SNIPER ({})]\x1b[0m n ≈ {:>14} ({:>5.1}%) | \
                         \x1b[1;32mSpeed:\x1b[0m {:>12} ({:>7.0} k/s) | \
                         \x1b[1mCores:\x1b[0m {:>2}/{} | \
                         \x1b[1;33mSurvivors:\x1b[0m {:>4} | \
                         \x1b[1mCandidates:\x1b[0m {} | \x1b[1;31mFound:\x1b[0m {}",
                        target_short,
                        format_number(current_n),
                        progress_pct,
                        n_speed_str,
                        k_per_sec,
                        active,
                        num_threads,
                        survivor_count,
                        candidates,
                        found
                    );

                    last_k = current_k;
                    last_time = now;
                }
            });

            // Multi-Core Rayon/Scoped Lock-Free Worker Threads
            let crt_arc = Arc::new(crt);
            let all_shifts_arc = Arc::new(all_shifts.clone());

            thread::scope(|s| {
                for _ in 0..num_threads {
                    let crt_ref = Arc::clone(&crt_arc);
                    let all_shifts_ref = Arc::clone(&all_shifts_arc);

                    s.spawn(move || {
                        self.stats.active_workers.fetch_add(1, Ordering::Relaxed);
                        const CHUNK_SIZE: u64 = 1024;

                        loop {
                            let k_chunk_start = GLOBAL_K_INDEX.fetch_add(CHUNK_SIZE, Ordering::Relaxed);
                            if k_chunk_start > k_end {
                                break;
                            }
                            let k_chunk_end = (k_chunk_start + CHUNK_SIZE).min(k_end + 1);

                            for k in k_chunk_start..k_chunk_end {
                                if max_ce > 0 && self.stats.counterexamples_found.load(Ordering::Relaxed) >= max_ce {
                                    break;
                                }

                                let n_k = r + k * m;
                                self.stats.current_max_n.fetch_max(n_k, Ordering::Relaxed);

                                // Fast Survivor-Only Early-Exit Scan
                                let mut represented = false;

                                for &s in &crt_ref.survivor_shifts {
                                    if s >= n_k {
                                        break;
                                    }
                                    let p = n_k - s;
                                    if p <= 1 {
                                        continue;
                                    }

                                    // Medium prime modular trial division filter (41..251)
                                    if is_composite_by_medium_primes(p, prime_kind) {
                                        continue;
                                    }

                                    // 7-base Montgomery Miller-Rabin Primality Test
                                    let is_valid = match prime_kind {
                                        PrimeKind::Prime => is_prime_fast(p),
                                        PrimeKind::TwinPrime => is_twin_prime(p),
                                        PrimeKind::SophieGermain => is_sophie_germain(p),
                                    };

                                    if is_valid {
                                        represented = true;
                                        break; // Nanosecond Early-Exit! Move to n_{k+1}
                                    }
                                }

                                if !represented {
                                    // Potential Counterexample Gap Detected!
                                    self.stats.candidates_tested.fetch_add(1, Ordering::Relaxed);

                                    // Perform exhaustive deterministic verification across ALL shifts
                                    if let Some(report) = self.verify_candidate_exhaustive(n_k, &all_shifts_ref, prime_kind) {
                                        let count = self.stats.counterexamples_found.fetch_add(1, Ordering::SeqCst);
                                        let cert = report.format_certificate(target.name());

                                        eprintln!("\n{}", "=".repeat(80));
                                        eprintln!("🚨 [COUNTEREXAMPLE FOUND VIA CRT SNIPER] (Occurrence #{}) 🚨", count + 1);
                                        eprintln!("Target: {}", target.name());
                                        eprintln!("Counterexample Integer: n = {}", n_k);
                                        eprintln!("Master Progression: n ≡ {} (mod {})", r, m);
                                        eprintln!("{}\n", cert);
                                        eprintln!("{}", "=".repeat(80));

                                        if let Ok(mut f) = OpenOptions::new().create(true).append(true).open(&self.config.output_file) {
                                            writeln!(f, "{}", cert).ok();
                                        }

                                        self.counterexamples.lock().unwrap().push(n_k);
                                    }
                                }
                            }

                            self.stats.k_steps_evaluated.fetch_add(k_chunk_end - k_chunk_start, Ordering::Relaxed);
                        }

                        self.stats.active_workers.fetch_sub(1, Ordering::Relaxed);
                    });
                }
            });

            running.store(false, Ordering::Relaxed);
            telemetry_handle.join().ok();
        }

        self.stats.counterexamples_found.load(Ordering::Relaxed) > 0
    }

    /// Full exhaustive mathematical verification of candidate `n` with SPF factorization.
    fn verify_candidate_exhaustive(
        &self,
        n: u64,
        all_shifts: &[u64],
        kind: PrimeKind,
    ) -> Option<CounterexampleReport> {
        // Quick primality pass over ALL shifts
        for &s in all_shifts {
            if s >= n {
                break;
            }
            let p = n - s;
            let valid = match kind {
                PrimeKind::Prime => is_prime_fast(p),
                PrimeKind::TwinPrime => is_twin_prime(p),
                PrimeKind::SophieGermain => is_sophie_germain(p),
            };
            if valid {
                return None;
            }
        }

        // Complete factorization of all differences for certificate
        let mut analyses = Vec::new();
        for &s in all_shifts {
            if s >= n {
                break;
            }
            let diff = n - s;
            let spf = smallest_prime_factor(diff);
            let fact = format_factorization(diff);

            analyses.push(DifferenceAnalysis {
                shift: s,
                difference: diff,
                smallest_prime_factor: spf,
                factorization: fact,
            });
        }

        Some(CounterexampleReport {
            candidate: n,
            total_shifts_tested: analyses.len(),
            analyses,
        })
    }
}

fn format_number(n: u64) -> String {
    if n >= 1_000_000_000_000_000 {
        format!("{:.3} P ({})", n as f64 / 1e15, n)
    } else if n >= 1_000_000_000_000 {
        format!("{:.3} T ({})", n as f64 / 1e12, n)
    } else if n >= 1_000_000_000 {
        format!("{:.3} B ({})", n as f64 / 1e9, n)
    } else if n >= 1_000_000 {
        format!("{:.3} M ({})", n as f64 / 1e6, n)
    } else {
        n.to_string()
    }
}
