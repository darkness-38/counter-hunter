//! Chinese Remainder Theorem (CRT) Covering and Modular Sniper Optimizer
//! Computes optimal residue progressions R (mod M) that pre-kill 80%-95% of shifts,
//! enabling high-velocity sniper scans across 10^14 - 10^16 integer space.

use crate::conjectures::PrimeKind;

/// Extended Euclidean Algorithm returning (gcd, x, y) such that a*x + b*y = gcd.
pub fn ext_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = ext_gcd(b, a % b);
        let x = y1;
        let y = x1 - (a / b) * y1;
        (g, x, y)
    }
}

/// Modular inverse of `a` modulo `m` using extended Euclidean algorithm.
pub fn mod_inverse(a: i128, m: i128) -> Option<i128> {
    let (g, x, _) = ext_gcd(a, m);
    if g == 1 {
        Some((x % m + m) % m)
    } else {
        None
    }
}

/// Combines two modular congruences n = r1 (mod m1) and n = r2 (mod m2) via Chinese Remainder Theorem.
/// Returns (r, m) where n = r (mod m) and m = m1 * m2, assuming gcd(m1, m2) = 1.
pub fn crt_combine(r1: u64, m1: u64, r2: u64, m2: u64) -> Option<(u64, u64)> {
    let inv = mod_inverse(m1 as i128, m2 as i128)?;
    let diff = ((r2 as i128) - ((r1 % m2) as i128)).rem_euclid(m2 as i128);
    let k = (diff * inv).rem_euclid(m2 as i128);
    let r = (r1 as u128) + (m1 as u128) * (k as u128);
    let m = (m1 as u128) * (m2 as u128);
    if m > u64::MAX as u128 {
        return None;
    }
    Some((r as u64, m as u64))
}

/// Medium prime pool used for ultra-fast modular screening before Montgomery Miller-Rabin tests.
pub const MEDIUM_PRIMES: [u64; 42] = [
    41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
    211, 223, 227, 229, 233, 239, 241, 251,
];

/// Small prime set tested for greedy modular covering.
pub const COVERING_PRIMES: [u64; 10] = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Determines whether choosing n = r (mod q) mathematically forces difference (n - s) to be composite.
#[inline(always)]
pub fn is_shift_killed(s: u64, r: u64, q: u64, kind: PrimeKind) -> bool {
    let s_mod = s % q;
    match kind {
        PrimeKind::Prime => s_mod == r,
        PrimeKind::TwinPrime => s_mod == r || s_mod == (r + 2) % q,
        PrimeKind::SophieGermain => s_mod == r || s_mod == (r + q.div_ceil(2)) % q,
    }
}

/// A fully optimized master CRT arithmetic progression: n_k = master_residue + k * master_modulus.
#[derive(Debug, Clone)]
pub struct CrtProgression {
    pub base_residue: u64,
    pub base_modulus: u64,
    pub master_residue: u64,
    pub master_modulus: u64,
    pub primes_used: Vec<(u64, u64)>, // (q, r_q)
    pub total_shifts: usize,
    pub survivor_shifts: Vec<u64>,
    pub killed_shifts_count: usize,
}

impl CrtProgression {
    /// Percentage of total shifts killed by this CRT master residue.
    pub fn kill_percentage(&self) -> f64 {
        if self.total_shifts == 0 {
            0.0
        } else {
            (self.killed_shifts_count as f64) / (self.total_shifts as f64) * 100.0
        }
    }
}

/// Builds an optimal CRT covering progression by selecting the most lethal residue mod each prime q in Q.
pub fn optimize_crt_covering(
    all_shifts: &[u64],
    base_residue: u64,
    base_modulus: u64,
    prime_kind: PrimeKind,
    max_master_modulus: u64,
) -> CrtProgression {
    let total_shifts = all_shifts.len();
    let mut survivor_shifts = all_shifts.to_vec();
    let mut current_r = base_residue;
    let mut current_m = base_modulus;
    let mut primes_used = Vec::new();

    for &q in &COVERING_PRIMES {
        // Skip primes dividing the base modulus
        if base_modulus.is_multiple_of(q) {
            continue;
        }

        // Check if multiplying by q exceeds max_master_modulus
        if let Some(next_m) = current_m.checked_mul(q) {
            if next_m > max_master_modulus {
                continue;
            }
        } else {
            continue;
        }

        // Find residue r in 0..q-1 that kills the highest number of surviving shifts
        let mut best_r = 0u64;
        let mut best_kill = 0usize;

        for r in 0..q {
            let mut kill_cnt = 0usize;
            for &s in &survivor_shifts {
                if is_shift_killed(s, r, q, prime_kind) {
                    kill_cnt += 1;
                }
            }
            if kill_cnt > best_kill {
                best_kill = kill_cnt;
                best_r = r;
            }
        }

        if best_kill > 0
            && let Some((new_r, new_m)) = crt_combine(current_r, current_m, best_r, q) {
                current_r = new_r;
                current_m = new_m;
                primes_used.push((q, best_r));

                // Filter out killed shifts from survivor set
                survivor_shifts.retain(|&s| !is_shift_killed(s, best_r, q, prime_kind));
            }
    }

    let killed_shifts_count = total_shifts.saturating_sub(survivor_shifts.len());

    CrtProgression {
        base_residue,
        base_modulus,
        master_residue: current_r,
        master_modulus: current_m,
        primes_used,
        total_shifts,
        survivor_shifts,
        killed_shifts_count,
    }
}

/// Medium prime trial division test. Returns true if difference `p = n - s` is divisible by any q in MEDIUM_PRIMES.
#[inline(always)]
pub fn is_composite_by_medium_primes(p: u64, kind: PrimeKind) -> bool {
    for &q in &MEDIUM_PRIMES {
        if p <= q {
            return false;
        }
        match kind {
            PrimeKind::Prime => {
                if p.is_multiple_of(q) {
                    return true;
                }
            }
            PrimeKind::TwinPrime => {
                if p.is_multiple_of(q) || (p + 2).is_multiple_of(q) {
                    return true;
                }
            }
            PrimeKind::SophieGermain => {
                if p.is_multiple_of(q) || (2 * p + 1).is_multiple_of(q) {
                    return true;
                }
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crt_combine() {
        // n = 2 mod 3, n = 3 mod 5 => n = 8 mod 15
        let (r, m) = crt_combine(2, 3, 3, 5).unwrap();
        assert_eq!(m, 15);
        assert_eq!(r % 3, 2);
        assert_eq!(r % 5, 3);
        assert_eq!(r, 8);
    }

    #[test]
    fn test_optimize_crt_covering() {
        let shifts = vec![2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
        let crt = optimize_crt_covering(&shifts, 1, 2, PrimeKind::Prime, 1000);
        assert!(crt.master_modulus > 2);
        assert!(crt.killed_shifts_count > 0);
        assert!(crt.survivor_shifts.len() < shifts.len());
    }
}

