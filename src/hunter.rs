//! Bit-Parallel Shift-OR Convolution Hunter Engine.
//! Executes ultra-fast vectorized bitwise sieving across Rayon parallel segments.

use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use rayon::prelude::*;

use crate::conjectures::{TargetConjecture, VerificationResult};
use crate::sieve::SegmentedSieve;

/// Statistics gathered during the search.
#[derive(Debug, Default)]
pub struct SearchStats {
    pub numbers_scanned: AtomicU64,
    pub candidates_tested: AtomicU64,
    pub counterexamples_found: AtomicU64,
}

/// A discovered counterexample record.
#[derive(Debug, Clone)]
pub struct DiscoveredCounterexample {
    pub target: TargetConjecture,
    pub n: u64,
    pub proof: String,
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

/// Configuration options for the hunt.
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
    pub verbose: bool,
    pub output_file: String,
}

impl Default for HuntConfig {
    fn default() -> Self {
        Self {
            target: TargetConjecture::SunTwinPrime,
            start: 10_000_000,
            end: 100_000_000,
            segment_size: 1 << 25, // 33,554,432 integers (4 MB bitset)
            max_sieve_margin: 1 << 21, // 2,097,152 integers
            max_shifts_to_apply: 150, // Optimal saturation threshold
            max_counterexamples: 20, // Stop generating proofs after this limit
            parity: ParityFilter::Default,
            num_threads: 0, // 0 = automatic
            verbose: false,
            output_file: "discovered_counterexamples.txt".to_string(),
        }
    }
}

/// The main Hunt Engine.
pub struct Hunter {
    pub config: HuntConfig,
    pub sieve: SegmentedSieve,
    pub shifts: Vec<u64>,
    pub stats: Arc<SearchStats>,
    pub counterexamples: Arc<Mutex<Vec<DiscoveredCounterexample>>>,
}

impl Hunter {
    pub fn new(config: HuntConfig) -> Self {
        let max_n = config.end;
        let sieve = SegmentedSieve::new(max_n);

        // Precompute shifts up to the sieve margin, taking the most effective shifts
        let mut margin_shifts = config.target.generate_shifts(config.max_sieve_margin);
        if config.max_shifts_to_apply > 0 && margin_shifts.len() > config.max_shifts_to_apply {
            margin_shifts.truncate(config.max_shifts_to_apply);
        }

        Self {
            config,
            sieve,
            shifts: margin_shifts,
            stats: Arc::new(SearchStats::default()),
            counterexamples: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// Executes the parallel hunt over the range `[start, end)`.
    pub fn hunt<F>(&self, on_segment_completed: F)
    where
        F: Fn(u64, u64, f64, &SearchStats) + Send + Sync,
    {
        if self.config.num_threads > 0 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(self.config.num_threads)
                .build_global()
                .ok();
        }

        let total_range = self.config.end.saturating_sub(self.config.start);
        if total_range == 0 {
            return;
        }

        let seg_size = self.config.segment_size;
        let num_segments = (total_range + seg_size - 1) / seg_size;

        // Build segment boundaries
        let segments: Vec<(u64, u64)> = (0..num_segments)
            .map(|i| {
                let s = self.config.start + i * seg_size;
                let e = (s + seg_size).min(self.config.end);
                (s, e)
            })
            .collect();

        let start_time = Instant::now();

        // Process segments in parallel batches or chunks
        segments.into_par_iter().for_each(|(seg_start, seg_end)| {
            self.process_segment(seg_start, seg_end);

            let scanned = self.stats.numbers_scanned.load(Ordering::Relaxed);
            let elapsed_sec = start_time.elapsed().as_secs_f64();
            let gnums_per_sec = if elapsed_sec > 0.0 {
                (scanned as f64) / (elapsed_sec * 1e9)
            } else {
                0.0
            };

            on_segment_completed(seg_start, seg_end, gnums_per_sec, &self.stats);
        });
    }

    /// Processes a single segment `[seg_start, seg_end)` using Bitwise Shift-OR Convolution.
    pub fn process_segment(&self, seg_start: u64, seg_end: u64) {
        let seg_len = seg_end - seg_start;
        let num_words = ((seg_len + 63) / 64) as usize;

        let margin = self.config.max_sieve_margin.min(seg_start);
        let sieved_start = seg_start - margin;
        let sieved_end = seg_end;

        // 1. Sieve primes in [seg_start - margin, seg_end)
        let prime_bitset = self.sieve.sieve_range(sieved_start, sieved_end);

        // For twin prime conjectures, compute twin prime bitset
        let base_bitset = if self.config.target.requires_twin_primes() {
            prime_bitset.compute_twin_primes()
        } else {
            prime_bitset
        };

        // 2. Initialize represented bitset for [seg_start, seg_end)
        let mut represented = vec![0u64; num_words];

        // 3. Bitwise Shift-OR Convolution
        // Bit `i` in `represented` represents `n = seg_start + i`.
        // Condition: (n - s) is prime => offset in base_bitset is (n - s) - sieved_start
        // = (seg_start + i - s) - (seg_start - margin) = margin - s + i.
        let margin_usize = margin as usize;

        for &s in &self.shifts {
            if s > margin {
                continue;
            }
            let bit_offset = margin_usize - (s as usize);
            let word_offset = bit_offset >> 6;
            let bit_shift = bit_offset & 63;

            if bit_shift == 0 {
                // Word-aligned fast path: 1 direct OR per word
                for w in 0..num_words {
                    if word_offset + w < base_bitset.words.len() {
                        represented[w] |= base_bitset.words[word_offset + w];
                    }
                }
            } else {
                // Unaligned path: combine two words and OR
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

        // Mask trailing unused bits in the last word
        let remainder = seg_len % 64;
        if remainder != 0 {
            let last_idx = num_words - 1;
            let mask = !((1u64 << remainder) - 1);
            represented[last_idx] |= mask;
        }

        // 4. Identify remaining 0-bits (unrepresented candidates)
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
                        self.handle_unrepresented_candidate(n);
                    }
                    temp &= temp - 1; // Clear lowest set bit
                }
            }
        }

        self.stats.numbers_scanned.fetch_add(seg_len, Ordering::Relaxed);
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

    /// Handles an unrepresented candidate `n` found by the bit-parallel sieve.
    /// Runs full deterministic Miller-Rabin across the entire domain.
    pub fn handle_unrepresented_candidate(&self, n: u64) {
        self.stats.candidates_tested.fetch_add(1, Ordering::Relaxed);

        let result = self.config.target.verify_candidate_deterministic(n);

        match result {
            VerificationResult::Represented { .. } => {
                // Was resolved by a larger shift outside sieve margin
            }
            VerificationResult::KnownSmallException => {
                if self.config.verbose {
                    println!("[INFO] Filtered known small exception: n = {}", n);
                }
            }
            VerificationResult::Counterexample { .. } => {
                let count = self.stats.counterexamples_found.fetch_add(1, Ordering::SeqCst);

                if self.config.max_counterexamples == 0 || (count as usize) < self.config.max_counterexamples {
                    // 🚨 AUTHENTIC COUNTEREXAMPLE FOUND! 🚨
                    let proof = self.config.target.generate_proof(n);

                    // Print to terminal immediately
                    eprintln!("\n{}", "=".repeat(80));
                    eprintln!("🚨 [COUNTEREXAMPLE FOUND] (Occurrence #{}) 🚨", count + 1);
                    eprintln!("Target: {}", self.config.target.name());
                    eprintln!("Counterexample Integer: n = {}", n);
                    eprintln!("{}\n", proof);
                    eprintln!("{}", "=".repeat(80));

                    // Save to output file
                    let mut file = OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&self.config.output_file)
                        .expect("Failed to open counterexamples output file");

                    writeln!(file, "{}", proof).ok();

                    let record = DiscoveredCounterexample {
                        target: self.config.target,
                        n,
                        proof,
                    };
                    self.counterexamples.lock().unwrap().push(record);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hunter_small_range() {
        let config = HuntConfig {
            target: TargetConjecture::SunCentralBinomial,
            start: 10,
            end: 1000,
            segment_size: 256,
            max_sieve_margin: 128,
            max_shifts_to_apply: 50,
            max_counterexamples: 20,
            parity: ParityFilter::Default,
            num_threads: 1,
            verbose: false,
            output_file: "test_counterexamples.txt".to_string(),
        };

        let hunter = Hunter::new(config);
        hunter.hunt(|_start, _end, _speed, _stats| {});

        let scanned = hunter.stats.numbers_scanned.load(Ordering::Relaxed);
        assert_eq!(scanned, 990);
    }
}
