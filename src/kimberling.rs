//! HEDEF 1: PPL 150 & PPL 143 — Clark Kimberling #13 Difference Sequence Hunter
//!
//! Investigates Clark Kimberling's Problem #13 (OEIS A131388 / A131389)
//! and the associated Prize Problem Ledger conjectures:
//! - PPL 143 (Kimberling #13.3): Whenever d(k) > 0, at least one of d(k+1), d(k+2), d(k+3) is positive.
//!   (Counterexample: d(k) > 0 followed by 3 consecutive non-positive terms, i.e., run of >= 3 negatives).
//! - PPL 150 (Kimberling #13.4): Whenever d(k) < 0, at least one of d(k+1), d(k+2), d(k+3) is negative.
//!   (Counterexample: d(k) < 0 followed by 3 consecutive non-negative terms, i.e., run of >= 3 positives).

use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

/// Dynamic flat bitset optimized for high-speed sequential lookups and bit tests.
pub struct FlatBitSet {
    words: Vec<u64>,
}

impl FlatBitSet {
    pub fn new(capacity_bits: usize) -> Self {
        let num_words = (capacity_bits + 63) / 64 + 4;
        Self {
            words: vec![0u64; num_words],
        }
    }

    #[inline(always)]
    pub fn is_set(&self, idx: usize) -> bool {
        let w = idx >> 6;
        if w < self.words.len() {
            ((self.words[w] >> (idx & 63)) & 1) != 0
        } else {
            false
        }
    }

    #[inline(always)]
    pub fn set(&mut self, idx: usize) {
        let w = idx >> 6;
        if w >= self.words.len() {
            self.words.resize(w + 1024, 0);
        }
        self.words[w] |= 1u64 << (idx & 63);
    }
}

/// Statistics and findings from a Kimberling sequence search.
#[derive(Debug, Clone)]
pub struct KimberlingReport {
    pub terms_evaluated: usize,
    pub max_consecutive_pos: usize,
    pub max_consecutive_neg: usize,
    pub ppl143_counterexamples: Vec<(usize, usize)>, // (step k, run length)
    pub ppl150_counterexamples: Vec<(usize, usize)>, // (step k, run length)
    pub final_a: usize,
    pub min_unused_a: usize,
    pub elapsed_secs: f64,
}

/// Run Kimberling #13 sequence generator up to `max_terms` and check PPL 143 and PPL 150.
pub fn hunt_kimberling(max_terms: usize, telemetry_interval_secs: u64, output_file: &str) -> KimberlingReport {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 PPL 143 & PPL 150: KIMBERLING #13 SEQUENCE COUNTEREXAMPLE HUNTER\x1b[0m");
    println!("   Clark Kimberling's 2-Sequence Permutation Algorithm (OEIS A131388 / A131389)");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1m🔍 Target Length     :\x1b[0m {} terms", max_terms);
    println!("\x1b[1m🎯 PPL 143 Condition :\x1b[0m d(k) > 0 => at least one of next 3 terms is positive (no run >= 3 negatives)");
    println!("\x1b[1m🎯 PPL 150 Condition :\x1b[0m d(k) < 0 => at least one of next 3 terms is negative (no run >= 3 positives)");
    println!("\x1b[1m📁 Proof Log File    :\x1b[0m {}", output_file);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let capacity = max_terms * 3 + 100_000;
    let mut used_a = FlatBitSet::new(capacity);
    let mut used_d_pos = FlatBitSet::new(capacity);
    let mut used_d_neg = FlatBitSet::new(capacity);

    // Initial state: a(1) = 1, d(1) = 0
    used_a.set(1);
    let mut curr = 1usize;
    let mut min_unused_a = 2usize;
    let mut min_unused_d_pos = 1usize;

    let mut max_pos_run = 0usize;
    let mut max_neg_run = 0usize;
    let mut cur_pos_run = 0usize;
    let mut cur_neg_run = 0usize;

    let mut ppl143_violations = Vec::new();
    let mut ppl150_violations = Vec::new();

    let start_time = Instant::now();
    let mut last_telemetry = Instant::now();

    for k in 1..max_terms {
        // Step 1: Look for greatest negative h in [-(curr-1), -1]
        // That is, smallest |h| >= 1 such that |h| not in used_d_neg and curr - |h| not in used_a.
        // Since curr - |h| must be an unused positive integer, curr - |h| >= min_unused_a.
        // Therefore, |h| <= curr - min_unused_a.
        let mut found_neg = false;
        let mut chosen_h = 0usize;

        if curr > min_unused_a {
            let max_h = curr - min_unused_a;
            let mut h = 1usize;
            while h <= max_h {
                if !used_d_neg.is_set(h) && !used_a.is_set(curr - h) {
                    found_neg = true;
                    chosen_h = h;
                    break;
                }
                h += 1;
            }
        }

        let diff: i64;
        if found_neg {
            used_d_neg.set(chosen_h);
            curr -= chosen_h;
            used_a.set(curr);
            diff = -(chosen_h as i64);
        } else {
            // Step 2: Choose least positive integer h not in used_d_pos such that curr + h not in used_a
            let mut h_p = min_unused_d_pos;
            while used_d_pos.is_set(h_p) || used_a.is_set(curr + h_p) {
                h_p += 1;
            }
            used_d_pos.set(h_p);
            curr += h_p;
            used_a.set(curr);
            diff = h_p as i64;

            while used_d_pos.is_set(min_unused_d_pos) {
                min_unused_d_pos += 1;
            }
        }

        while used_a.is_set(min_unused_a) {
            min_unused_a += 1;
        }

        // Track difference runs
        if diff > 0 {
            cur_pos_run += 1;
            cur_neg_run = 0;
            if cur_pos_run > max_pos_run {
                max_pos_run = cur_pos_run;
            }
        } else {
            cur_neg_run += 1;
            cur_pos_run = 0;
            if cur_neg_run > max_neg_run {
                max_neg_run = cur_neg_run;
            }
        }

        // Check for counterexamples
        if cur_neg_run >= 3 {
            println!(
                "\x1b[1;31m🚨 PPL 143 COUNTEREXAMPLE DISCOVERED at k = {}! Run length = {} negatives!\x1b[0m",
                k, cur_neg_run
            );
            ppl143_violations.push((k, cur_neg_run));
            log_counterexample(output_file, "PPL 143", k, cur_neg_run, curr);
        }

        if cur_pos_run >= 3 {
            println!(
                "\x1b[1;31m🚨 PPL 150 COUNTEREXAMPLE DISCOVERED at k = {}! Run length = {} positives!\x1b[0m",
                k, cur_pos_run
            );
            ppl150_violations.push((k, cur_pos_run));
            log_counterexample(output_file, "PPL 150", k, cur_pos_run, curr);
        }

        // Telemetry pulse
        if last_telemetry.elapsed().as_secs() >= telemetry_interval_secs || k == max_terms - 1 {
            let elapsed = start_time.elapsed().as_secs_f64();
            let rate = (k as f64) / elapsed / 1_000.0;
            let progress = (k as f64) / (max_terms as f64) * 100.0;
            println!(
                "📡 [TELEMETRY] [{:5.1}%] k = {:10} | {:7.1} kterms/s | Max Pos Run: {} | Max Neg Run: {} | a(k) = {:10} | min_unused = {:10}",
                progress, k, rate, max_pos_run, max_neg_run, curr, min_unused_a
            );
            last_telemetry = Instant::now();
        }
    }

    let total_elapsed = start_time.elapsed().as_secs_f64();

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🏁 KIMBERLING HUNT FINISHED\x1b[0m");
    println!("   Total Steps Evaluated      : {}", max_terms);
    println!("   Max Consecutive Positive d : {} (Conjecture PPL 150 expects <= 2)", max_pos_run);
    println!("   Max Consecutive Negative d : {} (Conjecture PPL 143 expects <= 2)", max_neg_run);
    println!("   PPL 143 Violations Found   : {}", ppl143_violations.len());
    println!("   PPL 150 Violations Found   : {}", ppl150_violations.len());
    println!("   Final a(N) Value           : {}", curr);
    println!("   Total Time Elapsed         : {:.2} seconds", total_elapsed);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    KimberlingReport {
        terms_evaluated: max_terms,
        max_consecutive_pos: max_pos_run,
        max_consecutive_neg: max_neg_run,
        ppl143_counterexamples: ppl143_violations,
        ppl150_counterexamples: ppl150_violations,
        final_a: curr,
        min_unused_a,
        elapsed_secs: total_elapsed,
    }
}

fn log_counterexample(path: &str, problem: &str, step_k: usize, run_len: usize, curr_a: usize) {
    if let Ok(mut f) = OpenOptions::new().create(true).append(true).open(path) {
        let _ = writeln!(
            f,
            "=== OFFICIAL BOUNTY COUNTEREXAMPLE FOUND ===\nProblem: {}\nStep k: {}\nConsecutive Run: {}\nCurrent a(k): {}\nTimestamp: {:?}\n============================================\n",
            problem, step_k, run_len, curr_a, Instant::now()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kimberling_initial_terms() {
        // Verify OEIS A131388 first terms:
        // [1, 2, 4, 3, 6, 10, 8, 5, 11, 7, 12, 19, 14, 22, 16, 9, 18, 28, 20, 31]
        let report = hunt_kimberling(20, 100, "test_log.txt");
        assert_eq!(report.final_a, 31);
        assert_eq!(report.max_consecutive_pos, 2);
        assert_eq!(report.max_consecutive_neg, 2);
    }
}
