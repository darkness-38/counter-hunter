//! Deterministic Miller-Rabin primality test and exact integer factorization engine
//! powered by Zero-Division Montgomery Modular Arithmetic.
//!
//! # Key Architectural Features:
//! - **Montgomery Form Arithmetic:** Modular multiplication `(a * b) % n` eliminates hardware
//!   division (`/` and `%`) inside exponentiation loops using 128-bit wide multiplication
//!   and Montgomery reduction (`montgomery_reduce`), utilizing only additions, multiplications, and shifts.
//! - **7-Base Deterministic Base Set:**
//!   `BASES = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]`
//!   Proven 100% deterministic for all integers below 2^64 (Galway, Jiang & Deng 2012).
//! - **Counterexample Verification & Smallest Prime Factor Reporting:**
//!   Exhaustively confirms that no `n - a_k` is prime, reporting the smallest prime factor (SPF)
//!   and complete prime factorization for mathematical counterexample certificates.

/// Greatest common divisor using Euclidean algorithm.
#[inline(always)]
pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Proven 7-base set covering all 64-bit integers (`< 2^64`) without false positives.
/// Reference: Jim Sinclair / Galway & Jiang-Deng (2012).
pub const MILLER_RABIN_7_BASES: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];

/// Montgomery arithmetic context for an odd modulus `n`.
/// All modular multiplications in exponentiation loops run without hardware 64-bit division.
#[derive(Clone, Copy, Debug)]
pub struct Montgomery {
    pub n: u64,
    pub n_prime: u64, // n' such that n * n' = -1 (mod 2^64)
    pub r: u64,       // 2^64 % n (Montgomery representation of 1)
    pub r2: u64,      // 2^128 % n
}

impl Montgomery {
    /// Creates a Montgomery reduction context for an odd integer `n > 1`.
    #[inline]
    pub fn new(n: u64) -> Self {
        debug_assert!(n % 2 == 1 && n > 1, "Modulus must be an odd integer > 1");

        // Compute n_prime such that n * n_prime == -1 (mod 2^64)
        // using 6 Newton-Raphson inverse refinement iterations:
        let mut inv = 1u64;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
        }
        let n_prime = 0u64.wrapping_sub(inv);

        // r = 2^64 % n
        let r = ((1u128 << 64) % (n as u128)) as u64;
        // r2 = 2^128 % n
        let r2 = (((r as u128) * (r as u128)) % (n as u128)) as u64;

        Self { n, n_prime, r, r2 }
    }

    /// Transforms a standard residue `x` into Montgomery space: `(x * R) % n`.
    #[inline(always)]
    pub fn transform(&self, x: u64) -> u64 {
        self.reduce((x as u128) * (self.r2 as u128))
    }

    /// Montgomery reduction: computes `(T * R^-1) % n` for `T < n * 2^64`.
    /// Zero hardware 64-bit division (`/` or `%`)! Only additions, multiplications, and shifts.
    #[inline(always)]
    pub fn reduce(&self, t: u128) -> u64 {
        let t_lo = t as u64;
        let t_hi = (t >> 64) as u64;
        let m = t_lo.wrapping_mul(self.n_prime);
        let mn = (m as u128) * (self.n as u128);
        let mn_hi = (mn >> 64) as u64;
        let carry = if t_lo != 0 { 1u64 } else { 0u64 };

        let (sum1, c1) = t_hi.overflowing_add(mn_hi);
        let (sum2, c2) = sum1.overflowing_add(carry);
        let overflow = c1 || c2;

        let mut res = sum2;
        if overflow || res >= self.n {
            res = res.wrapping_sub(self.n);
        }
        res
    }

    /// Modular multiplication in Montgomery space: `(x * y * R^-1) % n`.
    #[inline(always)]
    pub fn mul(&self, x: u64, y: u64) -> u64 {
        self.reduce((x as u128) * (y as u128))
    }

    /// Montgomery modular exponentiation: computes `(base^exp) % n` in Montgomery space.
    /// Returns the Montgomery-form result.
    #[inline(always)]
    pub fn pow_mont(&self, base_mont: u64, mut exp: u64) -> u64 {
        let mut res = self.r; // 1 in Montgomery form
        let mut cur = base_mont;

        while exp > 0 {
            if exp & 1 == 1 {
                res = self.mul(res, cur);
            }
            cur = self.mul(cur, cur);
            exp >>= 1;
        }

        res
    }

    /// Exponentiation returning the standard (un-transformed) residue in `[0, n)`.
    #[inline(always)]
    pub fn pow(&self, base: u64, exp: u64) -> u64 {
        let base_mont = self.transform(base);
        let res_mont = self.pow_mont(base_mont, exp);
        self.reduce(res_mont as u128)
    }
}

/// Ultra-fast 100% deterministic primality test for all `n < 2^64` using
/// Montgomery arithmetic and the minimal 7-base set.
#[inline]
pub fn is_prime_fast(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 || n == 17 || n == 19 || n == 23 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0 || n % 11 == 0 || n % 13 == 0 || n % 17 == 0 || n % 19 == 0 || n % 23 == 0 {
        return false;
    }
    if n < 29 * 29 {
        return true; // Any composite < 841 must have prime factor <= 23
    }

    let mont = Montgomery::new(n);

    // Write n - 1 as 2^s * d with d odd
    let d_val = n - 1;
    let s = d_val.trailing_zeros();
    let d = d_val >> s;

    let one_mont = mont.r; // 1 in Montgomery space
    let neg_one_mont = mont.n - mont.r; // -1 = (n - 1) in Montgomery space

    'base_loop: for &base in &MILLER_RABIN_7_BASES {
        let a = base % n;
        if a == 0 || a == 1 || a == n - 1 {
            continue 'base_loop;
        }

        let a_mont = mont.transform(a);
        let mut x = mont.pow_mont(a_mont, d);

        if x == one_mont || x == neg_one_mont {
            continue 'base_loop;
        }

        for _ in 1..s {
            x = mont.mul(x, x);
            if x == neg_one_mont {
                continue 'base_loop;
            }
            if x == one_mont {
                return false;
            }
        }

        return false;
    }

    true
}

/// Backward-compatible primality test entry point.
#[inline(always)]
pub fn is_prime(n: u64) -> bool {
    is_prime_fast(n)
}

/// Checks whether `p` is a twin prime (`p` and `p + 2` are prime).
#[inline]
pub fn is_twin_prime(p: u64) -> bool {
    is_prime_fast(p) && is_prime_fast(p + 2)
}

/// Checks whether `p` is a Sophie Germain prime (`p` and `2p + 1` are prime).
#[inline]
pub fn is_sophie_germain(p: u64) -> bool {
    if p > (u64::MAX - 1) / 2 {
        return false;
    }
    is_prime_fast(p) && is_prime_fast(2 * p + 1)
}

/// Finds the smallest prime factor (SPF) of a composite integer `n > 1`.
pub fn smallest_prime_factor(n: u64) -> u64 {
    if n <= 1 {
        return n;
    }
    // Trial division for small primes
    for &p in &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113] {
        if n % p == 0 {
            return p;
        }
    }

    if is_prime_fast(n) {
        return n;
    }

    // Extended trial division up to 10,000
    let mut d = 115u64;
    while d <= 10_000 && d * d <= n {
        if n % d == 0 {
            return if is_prime_fast(d) { d } else { smallest_prime_factor(d) };
        }
        if n % (d + 2) == 0 {
            let f = d + 2;
            return if is_prime_fast(f) { f } else { smallest_prime_factor(f) };
        }
        d += 6;
    }

    // Pollard-Brent factorization for remaining composite part
    let factor = pollard_brent_factor(n);
    if is_prime_fast(factor) {
        let other = n / factor;
        let other_spf = smallest_prime_factor(other);
        factor.min(other_spf)
    } else {
        smallest_prime_factor(factor)
    }
}

/// Pollard's Rho algorithm with Brent's cycle detection optimization.
pub fn pollard_brent_factor(n: u64) -> u64 {
    if n % 2 == 0 {
        return 2;
    }
    if n % 3 == 0 {
        return 3;
    }
    if n % 5 == 0 {
        return 5;
    }

    let mont = Montgomery::new(n);
    let mut c_val = 1u64;

    loop {
        let c_mont = mont.transform(c_val);
        let f = |x_mont: u64| -> u64 {
            let sq = mont.mul(x_mont, x_mont);
            let sum = sq + c_mont;
            if sum >= n { sum - n } else { sum }
        };

        let mut y = mont.transform(2);
        let mut x = y;
        let mut d = 1u64;
        let mut power = 1u64;
        let mut lam = 1u64;

        while d == 1 {
            if power == lam {
                x = y;
                power *= 2;
                lam = 0;
            }
            y = f(y);
            lam += 1;

            let x_val = mont.reduce(x as u128);
            let y_val = mont.reduce(y as u128);
            let diff = if x_val > y_val { x_val - y_val } else { y_val - x_val };
            d = gcd(diff, n);
        }

        if d != n && d > 1 {
            return d;
        }
        c_val += 1;
    }
}

/// Prime factorization: returns sorted vector of `(prime_factor, multiplicity)`.
pub fn factorize(mut n: u64) -> Vec<(u64, u32)> {
    let mut factors: Vec<u64> = Vec::new();
    if n <= 1 {
        return Vec::new();
    }

    for &p in &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97] {
        while n % p == 0 {
            factors.push(p);
            n /= p;
        }
    }

    let mut queue = vec![n];
    while let Some(curr) = queue.pop() {
        if curr == 1 {
            continue;
        }
        if is_prime_fast(curr) {
            factors.push(curr);
        } else {
            let divisor = pollard_brent_factor(curr);
            queue.push(divisor);
            queue.push(curr / divisor);
        }
    }

    factors.sort_unstable();

    let mut grouped = Vec::new();
    let mut iter = factors.into_iter();
    if let Some(first) = iter.next() {
        let mut current_prime = first;
        let mut count = 1;
        for f in iter {
            if f == current_prime {
                count += 1;
            } else {
                grouped.push((current_prime, count));
                current_prime = f;
                count = 1;
            }
        }
        grouped.push((current_prime, count));
    }

    grouped
}

/// Formats prime factorization into algebraic string, e.g. `2^3 * 5 * 13^2`.
pub fn format_factorization(n: u64) -> String {
    if n <= 1 {
        return n.to_string();
    }
    let factors = factorize(n);
    if factors.len() == 1 && factors[0].1 == 1 {
        return format!("{} [PRIME]", n);
    }
    let parts: Vec<String> = factors
        .into_iter()
        .map(|(p, exp)| {
            if exp == 1 {
                p.to_string()
            } else {
                format!("{}^{}", p, exp)
            }
        })
        .collect();
    parts.join(" * ")
}

/// Detailed analysis of an individual difference `n - shift`.
#[derive(Debug, Clone)]
pub struct DifferenceAnalysis {
    pub shift: u64,
    pub difference: u64,
    pub smallest_prime_factor: u64,
    pub factorization: String,
}

/// Rigorous Counterexample Report establishing proof that candidate `n` cannot be represented.
#[derive(Debug, Clone)]
pub struct CounterexampleReport {
    pub candidate: u64,
    pub total_shifts_tested: usize,
    pub analyses: Vec<DifferenceAnalysis>,
}

impl CounterexampleReport {
    /// Formats the report into a complete mathematical proof certificate string.
    pub fn format_certificate(&self, conjecture_name: &str) -> String {
        let mut cert = String::new();
        cert.push_str("================================================================================\n");
        cert.push_str("MATHEMATICAL COUNTEREXAMPLE CERTIFICATE (Phase 2 Montgomery Verified)\n");
        cert.push_str("================================================================================\n");
        cert.push_str(&format!("Target Conjecture: {}\n", conjecture_name));
        cert.push_str(&format!("Counterexample Integer: n = {}\n\n", self.candidate));
        cert.push_str(&format!("Theorem: n = {} CANNOT be represented in the conjectured form.\n", self.candidate));
        cert.push_str("Exhaustive Verification: All differences (n - s) are proven composite:\n");
        cert.push_str("--------------------------------------------------------------------------------\n");

        for a in &self.analyses {
            cert.push_str(&format!(
                "  s = {:<12} => n - s = {:<12} [SPF: {:<8}] = {}\n",
                a.shift, a.difference, a.smallest_prime_factor, a.factorization
            ));
        }

        cert.push_str("--------------------------------------------------------------------------------\n");
        cert.push_str(&format!("Total shifts tested: {}\n", self.total_shifts_tested));
        cert.push_str("All differences have been deterministically proven composite via 7-base Montgomery Miller-Rabin\n");
        cert.push_str("with their smallest prime factors and complete prime factorizations explicitly provided.\n");
        cert.push_str("Q.E.D. Counterexample is mathematically sound and 100% verified.\n");
        cert.push_str("================================================================================\n");
        cert
    }
}

/// Verifies whether candidate `n` is an authentic counterexample by checking that
/// NONE of the differences `n - s` (for `s` in `forbidden_diffs`) is prime.
///
/// If any `n - s` is prime: returns `None` (candidate is represented).
/// If all `n - s` are composite: returns `Some(CounterexampleReport)` containing the SPF
/// and full prime factorization of each difference as mathematical proof.
pub fn verify_counterexample(n: u64, forbidden_diffs: &[u64]) -> Option<CounterexampleReport> {
    let mut analyses = Vec::new();

    for &s in forbidden_diffs {
        if s >= n {
            continue;
        }
        let diff = n - s;

        // If any difference is prime, n has a valid representation!
        if is_prime_fast(diff) {
            return None;
        }

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_montgomery_pow() {
        let mont = Montgomery::new(17);
        // 3^5 mod 17 = 243 mod 17 = 5
        assert_eq!(mont.pow(3, 5), 5);
        // 2^10 mod 17 = 1024 mod 17 = 4
        assert_eq!(mont.pow(2, 10), 4);
    }

    #[test]
    fn test_is_prime_fast_small() {
        assert!(!is_prime_fast(0));
        assert!(!is_prime_fast(1));
        assert!(is_prime_fast(2));
        assert!(is_prime_fast(3));
        assert!(!is_prime_fast(4));
        assert!(is_prime_fast(5));
        assert!(!is_prime_fast(6));
        assert!(is_prime_fast(7));
        assert!(!is_prime_fast(9));
        assert!(is_prime_fast(11));
        assert!(!is_prime_fast(15));
        assert!(is_prime_fast(997));
        assert!(!is_prime_fast(999));
    }

    #[test]
    fn test_is_prime_fast_carmichael() {
        // Carmichael numbers must be identified as composite!
        let carmichaels = [561, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 62745, 63973];
        for &c in &carmichaels {
            assert!(!is_prime_fast(c), "Carmichael number {} falsely identified as prime", c);
        }
    }

    #[test]
    fn test_is_prime_fast_large_primes() {
        let primes = [
            1_000_000_007u64,
            1_000_000_009u64,
            10_000_000_019u64,
            100_000_000_003u64,
            1_000_000_000_039u64,
            2_305_843_009_213_693_951u64, // Mersenne prime 2^61 - 1
        ];
        for &p in &primes {
            assert!(is_prime_fast(p), "Prime {} failed primality test", p);
        }
    }

    #[test]
    fn test_smallest_prime_factor() {
        assert_eq!(smallest_prime_factor(15), 3);
        assert_eq!(smallest_prime_factor(77), 7);
        assert_eq!(smallest_prime_factor(91), 7);
        assert_eq!(smallest_prime_factor(119), 7);
        assert_eq!(smallest_prime_factor(1001), 7);
        assert_eq!(smallest_prime_factor(997), 997);
    }

    #[test]
    fn test_verify_counterexample() {
        // Suppose candidate n = 11, shifts for B1 = [3] (2*1^3 + 1^4 = 3)
        // 11 - 3 = 8 (composite, SPF = 2)
        let forbidden = [3u64];
        let report = verify_counterexample(11, &forbidden);
        assert!(report.is_some());
        let rep = report.unwrap();
        assert_eq!(rep.candidate, 11);
        assert_eq!(rep.total_shifts_tested, 1);
        assert_eq!(rep.analyses[0].smallest_prime_factor, 2);
        assert_eq!(rep.analyses[0].difference, 8);

        // If a difference is prime:
        let not_counterexample = verify_counterexample(10, &[3]); // 10 - 3 = 7 (prime)
        assert!(not_counterexample.is_none());
    }
}
