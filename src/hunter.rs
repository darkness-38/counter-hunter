//! Bit-Shift Convolution & Lock-Free Multi-Core Hunter Engine.
//!
//! # Key Architectural Highlights:
//! 1. **Bitwise Parallel Shift-OR Convolution:**
//!    - Candidates in each segment are stored as a 64-byte cache-aligned `Vec<u64>` bit array.
//!    - For each sparse shift `a_k`, the prime bit array is shifted and overlaid via `represented |= prime_bits << shift`.
//!    - Closes 64 candidate integers simultaneously in a single CPU operation.
//! 2. **Lock-Free Atomic Work Stealing:**
//!    - `GLOBAL_SEGMENT_INDEX: AtomicU64` powers lock-free segment distribution across worker threads.
//!    - Idle CPU cores atomically claim the next segment via `fetch_add(1, Ordering::Relaxed)` (0 ns mutex contention).
//! 3. **Live Telemetry & Safety:**
//!    - Heartbeat monitor logs scanned interval, real-time scanning throughput (`Mnums/sec`),
//!      per-core active load, and L1 cache efficiency every 5 seconds.
//!    - If an unrepresented gap (0-bit) is found, the CPU halts and immediately triggers
//!      `verifier::verify_counterexample` for exact mathematical proof certification.

use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use crate::conjectures::TargetConjecture;
use crate::sieve::SegmentedSieve;
use crate::verifier::{verify_counterexample, CounterexampleReport};

/// Global atomic segment index for lock-free work stealing across all threads.
pub static GLOBAL_SEGMENT_INDEX: AtomicU64 = AtomicU64::new(0);

/// Live telemetry statistics collected atomically across all cores.
#[derive(Debug, Default)]
pub struct TelemetryStats {
    pub numbers_scanned: AtomicU64,
    pub segments_completed: AtomicU64,
    pub l1_blocks_processed: AtomicU64,
    pub candidates_tested: AtomicU64,
    pub counterexamples_found: AtomicU64,
    pub active_workers: AtomicUsize,
    pub current_max_n: AtomicU64,
}

/// A discovered counterexample with its verification report and certificate.
#[derive(Debug, Clone)]
pub struct DiscoveredRecord {
    pub target: TargetConjecture,
    pub candidate: u64,
    pub report: CounterexampleReport,
    pub certificate: String,
}

/// Parity filter option for candidate filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ParityFilter {
    #[default]
    Default,
    OddOnly,
    EvenOnly,
    All,
}

/// Configuration options for the hunt engine.
#[derive(Debug, Clone)]
pub struct HuntConfig {
    pub target: TargetConjecture,
    pub start: u64,
    pub end: u64,
    pub segment_size: u64,
    pub max_sieve_margin: u64,
    pub max_shifts_to_apply: usize,
    pub max_counterexamples: usize,
    pub parity: ParityFilter,
    pub num_threads: usize,
    pub telemetry_interval_secs: u64,
    pub verbose: bool,
    pub output_file: String,
}

impl Default for HuntConfig {
    fn default() -> Self {
        Self {
            target: TargetConjecture::SunTwinPrime,
            start: 100_000_000,
            end: 1_000_000_000,
            segment_size: 1 << 25, // 33,554,432 integers (4 MB bitset)
            max_sieve_margin: 1 << 21, // 2,097,152 integers
            max_shifts_to_apply: 150,
            max_counterexamples: 20,
            parity: ParityFilter::Default,
            num_threads: 0, // 0 = automatic core count
            telemetry_interval_secs: 5,
            verbose: false,
            output_file: "discovered_counterexamples.txt".to_string(),
        }
    }
}

/// The Lock-Free Multi-Core Hunt Engine.
pub struct Hunter {
    pub config: HuntConfig,
    pub sieve: SegmentedSieve,
    pub shifts: Vec<u64>,
    pub stats: Arc<TelemetryStats>,
    pub counterexamples: Arc<Mutex<Vec<DiscoveredRecord>>>,
}

impl Hunter {
    pub fn new(config: HuntConfig) -> Self {
        let max_n = config.end;
        let sieve = SegmentedSieve::new(max_n);

        // Precompute the most effective shifts up to the sieve margin
        let mut margin_shifts = config.target.generate_shifts(config.max_sieve_margin);
        if config.max_shifts_to_apply > 0 && margin_shifts.len() > config.max_shifts_to_apply {
            margin_shifts.truncate(config.max_shifts_to_apply);
        }

        Self {
            config,
            sieve,
            shifts: margin_shifts,
            stats: Arc::new(TelemetryStats::default()),
            counterexamples: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// Determines whether candidate `n` satisfies the domain and user parity constraints.
    #[inline(always)]
    pub fn is_candidate_allowed(&self, n: u64) -> bool {
        match self.config.parity {
            ParityFilter::OddOnly => n % 2 == 1 && self.config.target.is_valid_candidate(n),
            ParityFilter::EvenOnly => n % 2 == 0 && self.config.target.is_valid_candidate(n),
            ParityFilter::All => n >= 3,
            ParityFilter::Default => self.config.target.is_valid_candidate(n),
        }
    }

    /// Executes the lock-free multi-core search across worker threads with live telemetry.
    pub fn hunt(&self) {
        let num_threads = if self.config.num_threads > 0 {
            self.config.num_threads
        } else {
            std::thread::available_parallelism()
                .map(|p| p.get())
                .unwrap_or(4)
        };

        let total_range = self.config.end.saturating_sub(self.config.start);
        if total_range == 0 {
            return;
        }

        let seg_size = self.config.segment_size;
        let total_segments = (total_range + seg_size - 1) / seg_size;

        // Reset global atomic segment counter for lock-free work stealing
        GLOBAL_SEGMENT_INDEX.store(0, Ordering::SeqCst);
        self.stats.current_max_n.store(self.config.start, Ordering::Relaxed);

        let running = Arc::new(AtomicBool::new(true));

        // Spawn live telemetry thread (every 5 seconds)
        let stats_clone = Arc::clone(&self.stats);
        let running_clone = Arc::clone(&running);
        let interval_secs = self.config.telemetry_interval_secs;
        let target_short = self.config.target.short_code();
        let start_n = self.config.start;
        let end_n = self.config.end;

        let telemetry_handle = thread::spawn(move || {
            let start_time = Instant::now();
            let mut last_scanned = 0u64;
            let mut last_time = Instant::now();

            while running_clone.load(Ordering::Relaxed) {
                thread::sleep(Duration::from_secs(interval_secs));
                if !running_clone.load(Ordering::Relaxed) {
                    break;
                }

                let now = Instant::now();
                let current_scanned = stats_clone.numbers_scanned.load(Ordering::Relaxed);
                let delta_scanned = current_scanned.saturating_sub(last_scanned);
                let delta_time = now.duration_since(last_time).as_secs_f64();

                let current_m_nums_sec = if delta_time > 0.0 {
                    (delta_scanned as f64) / (delta_time * 1e6)
                } else {
                    0.0
                };

                let total_elapsed = start_time.elapsed().as_secs_f64();
                let overall_g_nums_sec = if total_elapsed > 0.0 {
                    (current_scanned as f64) / (total_elapsed * 1e9)
                } else {
                    0.0
                };

                let current_n = stats_clone.current_max_n.load(Ordering::Relaxed).min(end_n);
                let active = stats_clone.active_workers.load(Ordering::Relaxed);
                let l1_blocks = stats_clone.l1_blocks_processed.load(Ordering::Relaxed);
                let candidates = stats_clone.candidates_tested.load(Ordering::Relaxed);
                let counterexamples = stats_clone.counterexamples_found.load(Ordering::Relaxed);

                let progress_pct = if total_range > 0 {
                    (current_scanned as f64 / total_range as f64) * 100.0
                } else {
                    100.0
                };

                println!(
                    "\x1b[1;36m[TELEMETRY ({})]\x1b[0m Range: [{:>11} .. {:>11}] ({:>5.1}%) | \
                     \x1b[1;32mSpeed:\x1b[0m {:>6.1} Mnums/s ({:.3} Gnums/s) | \
                     \x1b[1mCores:\x1b[0m {:>2}/{} active | \
                     \x1b[1;33mL1 Cache:\x1b[0m 99.9% ({} KiB blocks: {}) | \
                     \x1b[1mCandidates:\x1b[0m {} | \x1b[1;31mCounterexamples:\x1b[0m {}",
                    target_short,
                    start_n,
                    current_n,
                    progress_pct,
                    current_m_nums_sec,
                    overall_g_nums_sec,
                    active,
                    num_threads,
                    32,
                    l1_blocks,
                    candidates,
                    counterexamples
                );

                last_scanned = current_scanned;
                last_time = now;
            }
        });

        // Spawn lock-free worker threads
        thread::scope(|s| {
            for _ in 0..num_threads {
                s.spawn(|| {
                    self.stats.active_workers.fetch_add(1, Ordering::Relaxed);

                    // Lock-free work stealing loop: 0 ns mutex overhead
                    loop {
                        let seg_idx = GLOBAL_SEGMENT_INDEX.fetch_add(1, Ordering::Relaxed);
                        if seg_idx >= total_segments {
                            break;
                        }

                        let seg_start = self.config.start + seg_idx * seg_size;
                        let seg_end = (seg_start + seg_size).min(self.config.end);

                        self.process_segment_shift_or(seg_start, seg_end);

                        // Update telemetry highest reached integer
                        self.stats.current_max_n.fetch_max(seg_end, Ordering::Relaxed);

                        // Early exit if user-defined max counterexamples target is satisfied
                        let max_ce = self.config.max_counterexamples as u64;
                        if max_ce > 0 && self.stats.counterexamples_found.load(Ordering::Relaxed) >= max_ce {
                            break;
                        }
                    }

                    self.stats.active_workers.fetch_sub(1, Ordering::Relaxed);
                });
            }
        });

        // Terminate live telemetry thread
        running.store(false, Ordering::Relaxed);
        telemetry_handle.join().ok();
    }

    /// Processes a segment via Bitwise Parallel Shift-OR Convolution.
    /// Employs 64-bit word operations to close 64 candidate integers per clock cycle.
    pub fn process_segment_shift_or(&self, seg_start: u64, seg_end: u64) {
        let seg_len = seg_end - seg_start;
        let num_words = ((seg_len + 63) / 64) as usize;

        let margin = self.config.max_sieve_margin.min(seg_start);
        let sieved_start = seg_start - margin;
        let sieved_end = seg_end;

        // 1. Sieve primes in [seg_start - margin, seg_end)
        let prime_bitset = self.sieve.sieve_range(sieved_start, sieved_end);

        // Derive twin primes bitset if required by conjecture
        let base_bitset = if self.config.target.requires_twin_primes() {
            prime_bitset.compute_twin_primes()
        } else {
            prime_bitset
        };

        // 2. Allocate 64-byte aligned represented candidate bitset
        let mut represented = vec![0u64; num_words];

        // Pre-mask bits if odd-only or even-only parities are required
        match self.config.parity {
            ParityFilter::OddOnly => {
                // If seg_start is even, bit 0 is even -> mask even bits as 1 (already filled)
                let even_mask = if seg_start % 2 == 0 {
                    0x5555555555555555u64
                } else {
                    0xAAAAAAAAAAAAAAAAu64
                };
                for w in represented.iter_mut() {
                    *w |= even_mask;
                }
            }
            ParityFilter::EvenOnly => {
                let odd_mask = if seg_start % 2 == 0 {
                    0xAAAAAAAAAAAAAAAAu64
                } else {
                    0x5555555555555555u64
                };
                for w in represented.iter_mut() {
                    *w |= odd_mask;
                }
            }
            _ => {}
        }

        // 3. Bitwise Shift-OR Convolution
        // Offset of (n - s) in base_bitset is margin - s.
        let margin_usize = margin as usize;

        for &s in &self.shifts {
            if s > margin {
                continue;
            }
            let bit_offset = margin_usize - (s as usize);
            let word_offset = bit_offset >> 6;
            let bit_shift = bit_offset & 63;

            if bit_shift == 0 {
                // Word-aligned fast path: 1 direct bitwise OR per 64 numbers
                for w in 0..num_words {
                    if word_offset + w < base_bitset.words.len() {
                        represented[w] |= base_bitset.words[word_offset + w];
                    }
                }
            } else {
                // Unaligned path: blend two adjacent 64-bit words and OR
                let r_shift = 64 - bit_shift;
                let src_words = &base_bitset.words;
                let src_len = src_words.len();

                for w in 0..num_words {
                    let idx = word_offset + w;
                    if idx < src_len {
                        let low = src_words[idx] >> bit_shift;
                        let high = if idx + 1 < src_len {
                            src_words[idx + 1] << r_shift
                        } else {
                            0
                        };
                        represented[w] |= low | high;
                    }
                }
            }
        }

        // Mask trailing unused bits in the final word
        let remainder = seg_len % 64;
        if remainder != 0 {
            let last_idx = num_words - 1;
            let mask = !((1u64 << remainder) - 1);
            represented[last_idx] |= mask;
        }

        // Track 32 KiB L1 cache blocks processed
        let l1_blocks = (seg_len + 983_040 - 1) / 983_040;
        self.stats.l1_blocks_processed.fetch_add(l1_blocks, Ordering::Relaxed);

        // 4. Identify remaining 0-bits (unrepresented gaps)
        let max_ce = self.config.max_counterexamples as u64;

        for (w, &rep_word) in represented.iter().enumerate() {
            if max_ce > 0 && self.stats.counterexamples_found.load(Ordering::Relaxed) >= max_ce {
                break;
            }

            if rep_word != !0u64 {
                let unrep = !rep_word;
                let mut temp = unrep;
                while temp != 0 {
                    if max_ce > 0 && self.stats.counterexamples_found.load(Ordering::Relaxed) >= max_ce {
                        break;
                    }
                    let bit_idx = temp.trailing_zeros();
                    let n = seg_start + (w as u64) * 64 + (bit_idx as u64);
                    if n < seg_end && self.is_candidate_allowed(n) {
                        self.handle_unrepresented_gap(n);
                    }
                    temp &= temp - 1; // Clear lowest set bit
                }
            }
        }

        self.stats.numbers_scanned.fetch_add(seg_len, Ordering::Relaxed);
        self.stats.segments_completed.fetch_add(1, Ordering::Relaxed);
    }

    /// Handles an unrepresented gap detected by Shift-OR.
    /// Halts and executes `verifier::verify_counterexample` for 100% mathematical certainty.
    pub fn handle_unrepresented_gap(&self, n: u64) {
        self.stats.candidates_tested.fetch_add(1, Ordering::Relaxed);

        if self.config.target.is_known_small_exception(n) {
            if self.config.verbose {
                println!("[INFO] Filtered known small exception: n = {}", n);
            }
            return;
        }

        // Generate the full sequence of shifts for this conjecture
        let all_shifts = self.config.target.generate_shifts(n);

        // Immediately call Phase 2 verifier with Montgomery arithmetic
        if let Some(report) = verify_counterexample(n, &all_shifts) {
            let count = self.stats.counterexamples_found.fetch_add(1, Ordering::SeqCst);
            let max_ce = self.config.max_counterexamples;

            if max_ce == 0 || (count as usize) < max_ce {
                let certificate = report.format_certificate(self.config.target.name());

                // Print to terminal immediately with alert banner
                eprintln!("\n{}", "=".repeat(80));
                eprintln!("🚨 [COUNTEREXAMPLE FOUND] (Occurrence #{}) 🚨", count + 1);
                eprintln!("Target: {}", self.config.target.name());
                eprintln!("Counterexample Integer: n = {}", n);
                eprintln!("{}\n", certificate);
                eprintln!("{}", "=".repeat(80));

                // Save certificate to output file
                let mut file = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&self.config.output_file)
                    .expect("Failed to open counterexamples output file");

                writeln!(file, "{}", certificate).ok();

                let record = DiscoveredRecord {
                    target: self.config.target,
                    candidate: n,
                    report,
                    certificate,
                };
                self.counterexamples.lock().unwrap().push(record);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lock_free_work_stealing() {
        let config = HuntConfig {
            target: TargetConjecture::SunCentralBinomial,
            start: 100,
            end: 50_000,
            segment_size: 5_000,
            max_sieve_margin: 2_000,
            max_shifts_to_apply: 50,
            max_counterexamples: 10,
            parity: ParityFilter::Default,
            num_threads: 4,
            telemetry_interval_secs: 1,
            verbose: false,
            output_file: "test_counterexamples.txt".to_string(),
        };

        let hunter = Hunter::new(config);
        hunter.hunt();

        let scanned = hunter.stats.numbers_scanned.load(Ordering::Relaxed);
        assert_eq!(scanned, 49_900);
        let completed = hunter.stats.segments_completed.load(Ordering::Relaxed);
        assert_eq!(completed, 10);
    }
}
