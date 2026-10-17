use crate::verifier::{
    format_factorization, is_prime, is_sophie_germain, is_twin_prime, verify_counterexample,
};

/// Type of prime required by the conjecture target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimeKind {
    /// Ordinary prime p.
    Prime,
    /// Twin prime: both p and p+2 are prime.
    TwinPrime,
    /// Sophie Germain prime: both p and 2p+1 are prime.
    SophieGermain,
}

/// The supported open conjecture targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetConjecture {
    // === NEW MATHEMATICALLY RIGOROUS CRT SNIPER TARGETS ===
    /// 1. binom-2a-prime: n = p + binom(2k, k) + 2^a (k >= 1, a >= 1, p prime), n = 1 mod 2
    Binom2aPrime,
    /// 2. binom-2a-twin: n = p + binom(2k, k) + 2^a (k >= 1, a >= 1, p twin prime), n = 5 mod 6
    Binom2aTwin,
    /// 3. twin-2a-3b: n = p + 2^a + 3^b (a >= 1, b >= 0, p twin prime), n = 0, 4 mod 6
    Twin2a3b,
    /// 4. sophie-2a-3b: n = p + 2^a + 3^b (a >= 1, b >= 0, p Sophie Germain), n = 0, 4 mod 6
    Sophie2a3b,
    /// 5. erdos-2a-2b-twin: n = p + 2^a + 2^b (1 <= a <= b, p twin prime), n = 5 mod 6
    Erdos2a2bTwin,

    // === LEGACY BRUTE-FORCE TARGETS ===
    SunTwinPrime,
    SunSophieGermain,
    SunCubic2X3Y4,
    SunCubicX32Y3,
    SunCentralBinomial,
    SunFibonacciLucas,
    SunPrimeSquareRomanov,
}

/// Result of deterministic verification of an unrepresented candidate.
#[derive(Debug, Clone)]
pub enum VerificationResult {
    /// A valid representation was found: n = prime_part + shift.
    Represented { prime_part: u64, shift: u64, details: String },
    /// 100% Verified Counterexample: NO representation exists in the entire domain!
    Counterexample { failed_shifts: Vec<u64> },
    /// Candidate is an established OEIS small exception (< 10^8).
    KnownSmallException,
}

/// Pre-Flight Validator Report evaluating candidate density and small exception behavior.
#[derive(Debug, Clone)]
pub struct PreflightReport {
    pub target: TargetConjecture,
    pub start: u64,
    pub end: u64,
    pub tested_count: usize,
    pub exceptions_count: usize,
    pub sample_exceptions: Vec<u64>,
    pub passed: bool,
    pub verdict: &'static str,
}

impl TargetConjecture {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Binom2aPrime => "Target 1: binom-2a-prime (n = p + binom(2k,k) + 2^a, n ≡ 1 mod 2)",
            Self::Binom2aTwin => "Target 2: binom-2a-twin (n = p + binom(2k,k) + 2^a, p twin prime, n ≡ 5 mod 6)",
            Self::Twin2a3b => "Target 3: twin-2a-3b (n = p + 2^a + 3^b, p twin prime, n ≡ 0, 4 mod 6)",
            Self::Sophie2a3b => "Target 4: sophie-2a-3b (n = p + 2^a + 3^b, p Sophie Germain, n ≡ 0, 4 mod 6)",
            Self::Erdos2a2bTwin => "Target 5: erdos-2a-2b-twin (n = p + 2^a + 2^b, p twin prime, n ≡ 5 mod 6)",

            Self::SunTwinPrime => "Target A1: Sun's Twin-Prime Power Conjecture (n = p + 2^a + 3^b)",
            Self::SunSophieGermain => "Target A2: Sun's Sophie-Germain Power Conjecture (n = p + 2^a + 3^b)",
            Self::SunCubic2X3Y4 => "Target B1: Sun's Mixed Cubic-Prime Representation (n = p + 2x^3 + y^4)",
            Self::SunCubicX32Y3 => "Target B2: Sun's Mixed Cubic-Prime Representation (n = p + x^3 + 2y^3)",
            Self::SunCentralBinomial => "Target C: Sun's Central Binomial + Prime Conjecture (n = p + binom(2k,k) ± 2^m)",
            Self::SunFibonacciLucas => "Target D1: Sun's Fibonacci/Lucas + Prime Conjecture (n = p + F_k + L_m)",
            Self::SunPrimeSquareRomanov => "Target D2: Sun's Prime Square + Prime + Power of 2 (n = p^2 + q + 2^a)",
        }
    }

    pub fn short_code(&self) -> &'static str {
        match self {
            Self::Binom2aPrime => "binom-2a-prime",
            Self::Binom2aTwin => "binom-2a-twin",
            Self::Twin2a3b => "twin-2a-3b",
            Self::Sophie2a3b => "sophie-2a-3b",
            Self::Erdos2a2bTwin => "erdos-2a-2b-twin",

            Self::SunTwinPrime => "a-twin",
            Self::SunSophieGermain => "a-sophie",
            Self::SunCubic2X3Y4 => "b-cubic1",
            Self::SunCubicX32Y3 => "b-cubic2",
            Self::SunCentralBinomial => "c-binom",
            Self::SunFibonacciLucas => "d-fib",
            Self::SunPrimeSquareRomanov => "d-square",
        }
    }

    /// Prime kind required by the conjecture.
    pub fn prime_kind(&self) -> PrimeKind {
        match self {
            Self::Binom2aPrime
            | Self::SunCubic2X3Y4
            | Self::SunCubicX32Y3
            | Self::SunCentralBinomial
            | Self::SunFibonacciLucas
            | Self::SunPrimeSquareRomanov => PrimeKind::Prime,

            Self::Binom2aTwin
            | Self::Twin2a3b
            | Self::Erdos2a2bTwin
            | Self::SunTwinPrime => PrimeKind::TwinPrime,

            Self::Sophie2a3b
            | Self::SunSophieGermain => PrimeKind::SophieGermain,
        }
    }

    /// Base modulus and allowed residue classes.
    pub fn base_modulus_and_residues(&self) -> (u64, Vec<u64>) {
        match self {
            Self::Binom2aPrime => (2, vec![1]),
            Self::Binom2aTwin => (6, vec![5]),
            Self::Twin2a3b | Self::Sophie2a3b => (6, vec![0, 4]),
            Self::Erdos2a2bTwin => (6, vec![5]),

            Self::SunTwinPrime | Self::SunSophieGermain => (6, vec![0, 4]),
            Self::SunCubic2X3Y4 | Self::SunCubicX32Y3 => (1, vec![0]),
            Self::SunCentralBinomial => (2, vec![1]),
            Self::SunFibonacciLucas => (1, vec![0]),
            Self::SunPrimeSquareRomanov => (2, vec![0]),
        }
    }

    /// Whether the fast sieve should sieve twin primes instead of regular primes.
    pub fn requires_twin_primes(&self) -> bool {
        matches!(self.prime_kind(), PrimeKind::TwinPrime)
    }

    /// Checks if `n` satisfies the parity/domain constraints of the conjecture.
    pub fn is_valid_candidate(&self, n: u64) -> bool {
        match self {
            Self::Binom2aPrime => n > 4 && n % 2 == 1,
            Self::Binom2aTwin => n > 4 && n % 6 == 5,
            Self::Twin2a3b | Self::Sophie2a3b => n > 4 && n.is_multiple_of(2) && (n.is_multiple_of(6) || n % 6 == 4),
            Self::Erdos2a2bTwin => n > 4 && n % 6 == 5,

            Self::SunTwinPrime | Self::SunSophieGermain => n > 4 && n.is_multiple_of(2) && (n.is_multiple_of(6) || n % 6 == 4),
            Self::SunCubic2X3Y4 | Self::SunCubicX32Y3 => n >= 5,
            Self::SunCentralBinomial => n > 2 && n % 2 == 1,
            Self::SunFibonacciLucas => n >= 4,
            Self::SunPrimeSquareRomanov => n >= 10 && n.is_multiple_of(2),
        }
    }

    /// Generates the sequence of sparse shifts `a_k <= max_shift`.
    pub fn generate_shifts(&self, max_shift: u64) -> Vec<u64> {
        let mut shifts = Vec::new();

        match self {
            Self::Binom2aPrime | Self::Binom2aTwin => {
                // s = binom(2k, k) + 2^a (k >= 1, a >= 1)
                let mut k = 1u64;
                loop {
                    let c = central_binomial(k);
                    if c == 0 || c >= max_shift {
                        break;
                    }
                    let mut a = 1u32;
                    while let Some(p2) = 1u64.checked_shl(a) {
                        if let Some(sum) = c.checked_add(p2) {
                            if sum <= max_shift {
                                shifts.push(sum);
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                        a += 1;
                    }
                    k += 1;
                }
            }

            Self::Twin2a3b | Self::Sophie2a3b => {
                // s = 2^a + 3^b (a >= 1, b >= 0)
                let mut a = 1u32;
                while let Some(p2) = 1u64.checked_shl(a) {
                    if p2 >= max_shift {
                        break;
                    }
                    let mut p3 = 1u64;
                    loop {
                        if let Some(sum) = p2.checked_add(p3) {
                            if sum <= max_shift {
                                shifts.push(sum);
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                        match p3.checked_mul(3) {
                            Some(next_p3) => p3 = next_p3,
                            None => break,
                        }
                    }
                    a += 1;
                }
            }

            Self::Erdos2a2bTwin => {
                // s = 2^a + 2^b (1 <= a <= b)
                let mut a = 1u32;
                while let Some(p2a) = 1u64.checked_shl(a) {
                    if p2a >= max_shift {
                        break;
                    }
                    let mut b = a;
                    while let Some(p2b) = 1u64.checked_shl(b) {
                        if let Some(sum) = p2a.checked_add(p2b) {
                            if sum <= max_shift {
                                shifts.push(sum);
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                        b += 1;
                    }
                    a += 1;
                }
            }

            Self::SunTwinPrime | Self::SunSophieGermain => {
                // s = 2^a + 3^b (a >= 1, b >= 0 for even integers)
                let mut a = 1u32;
                while let Some(p2) = 1u64.checked_shl(a) {
                    if p2 >= max_shift {
                        break;
                    }
                    let mut p3 = 1u64;
                    loop {
                        if let Some(sum) = p2.checked_add(p3) {
                            if sum <= max_shift {
                                shifts.push(sum);
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                        match p3.checked_mul(3) {
                            Some(next_p3) => p3 = next_p3,
                            None => break,
                        }
                    }
                    a += 1;
                }
            }

            Self::SunCubic2X3Y4 => {
                // s = 2*x^3 + y^4 for x >= 1, y >= 1
                let max_x = (((max_shift as f64) / 2.0).cbrt().ceil() as u64).max(1);
                for x in 1..=max_x {
                    let x3_2 = 2 * x * x * x;
                    if x3_2 >= max_shift {
                        break;
                    }
                    let mut y = 1u64;
                    while let Some(val) = y.checked_mul(y).and_then(|y2| y2.checked_mul(y2)) {
                        let y4 = val;
                        let sum = x3_2 + y4;
                        if sum <= max_shift {
                            shifts.push(sum);
                        } else {
                            break;
                        }
                        y += 1;
                    }
                }
            }

            Self::SunCubicX32Y3 => {
                // s = x^3 + 2*y^3 for x >= 1, y >= 1
                let max_x = ((max_shift as f64).cbrt().ceil() as u64).max(1);
                for x in 1..=max_x {
                    let x3 = x * x * x;
                    if x3 >= max_shift {
                        break;
                    }
                    let mut y = 1u64;
                    while let Some(val) = y.checked_mul(y).and_then(|y2| y2.checked_mul(y)).and_then(|v| v.checked_mul(2)) {
                        let y3_2 = val;
                        let sum = x3 + y3_2;
                        if sum <= max_shift {
                            shifts.push(sum);
                        } else {
                            break;
                        }
                        y += 1;
                    }
                }
            }

            Self::SunCentralBinomial => {
                // s = binom(2k, k) +/- 2^m
                let max_k = 30u64;
                for k in 1..=max_k {
                    let c = central_binomial(k);
                    if c == 0 || c >= max_shift {
                        break;
                    }
                    let mut m = 0u32;
                    while let Some(p2) = 1u64.checked_shl(m) {
                        if let Some(sum) = c.checked_add(p2)
                            && sum <= max_shift {
                                shifts.push(sum);
                            }
                        if c > p2 {
                            let diff = c - p2;
                            if diff <= max_shift && diff > 0 {
                                shifts.push(diff);
                            }
                        }
                        m += 1;
                        if p2 >= max_shift {
                            break;
                        }
                    }
                }
            }

            Self::SunFibonacciLucas => {
                // s = F_k + L_m (k >= 1, m >= 1)
                let fibs = generate_fibonacci(max_shift);
                let lucas = generate_lucas(max_shift);

                for &f in &fibs {
                    for &l in &lucas {
                        let sum = f + l;
                        if sum <= max_shift {
                            shifts.push(sum);
                        } else {
                            break;
                        }
                    }
                }
            }

            Self::SunPrimeSquareRomanov => {
                // s = p^2 + 2^a (p prime, a >= 0)
                let max_p = ((max_shift as f64).sqrt().ceil() as u64).max(2);
                let primes = crate::sieve::generate_base_primes(max_p);

                for &p in &primes {
                    let p2 = p * p;
                    if p2 >= max_shift {
                        break;
                    }
                    let mut p2_pow = 1u64;
                    loop {
                        let sum = p2 + p2_pow;
                        if sum <= max_shift {
                            shifts.push(sum);
                        } else {
                            break;
                        }
                        match p2_pow.checked_mul(2) {
                            Some(next_pow) => p2_pow = next_pow,
                            None => break,
                        }
                    }
                }
            }
        }

        shifts.sort_unstable();
        shifts.dedup();
        shifts
    }

    /// Known OEIS small exceptions below 10^8 that are mathematically expected.
    pub fn is_known_small_exception(&self, n: u64) -> bool {
        if n >= 100_000_000 {
            return false;
        }
        match self {
            Self::SunCubic2X3Y4 => {
                matches!(
                    n,
                    1..=4 | 7 | 9 | 11..=13 | 15 | 17..=18 | 27 | 33 | 38
                        | 42 | 52..=53 | 67 | 80 | 95 | 109 | 165 | 193 | 235 | 239 | 454
                        | 552 | 1227
                )
            }
            Self::SunCubicX32Y3 => {
                matches!(n, 1..=4 | 7 | 9 | 11 | 18 | 25 | 38 | 49 | 59 | 101)
            }
            Self::SunCentralBinomial | Self::Binom2aPrime => {
                matches!(n, 1..=2)
            }
            Self::SunFibonacciLucas => {
                matches!(n, 1..=3)
            }
            Self::SunPrimeSquareRomanov => {
                matches!(n, 1..=6)
            }
            _ => false,
        }
    }

    /// Runs the Mandatory Pre-Flight Validator on a small range (e.g. 10,000 .. 100,000).
    /// Tests candidate density to guarantee that the target is a genuine open problem without trivial obstacles.
    pub fn run_preflight_check(&self, start: u64, end: u64) -> PreflightReport {
        let shifts = self.generate_shifts(end);
        let mut tested_count = 0;
        let mut exceptions = Vec::new();
        let kind = self.prime_kind();

        for n in start..=end {
            if !self.is_valid_candidate(n) {
                continue;
            }
            tested_count += 1;
            let mut represented = false;

            for &s in &shifts {
                if s >= n {
                    break;
                }
                let p = n - s;
                let is_valid_prime = match kind {
                    PrimeKind::Prime => crate::verifier::is_prime_fast(p),
                    PrimeKind::TwinPrime => crate::verifier::is_twin_prime(p),
                    PrimeKind::SophieGermain => crate::verifier::is_sophie_germain(p),
                };
                if is_valid_prime {
                    represented = true;
                    break;
                }
            }

            if !represented {
                exceptions.push(n);
            }
        }

        let exceptions_count = exceptions.len();
        // Rare exceptions (<= 100 out of thousands tested) -> Passed!
        // > 200 exceptions -> Trivial / Density Insufficient
        let passed = exceptions_count <= 100;
        let verdict = if passed {
            "✅ [PRE-FLIGHT PASSED: Nadir İstisnalı / Temiz Hedef - CRT Avına Uygun]"
        } else {
            "⚠️ [PRE-FLIGHT FAILED: Trivial / Yoğunluk Yetersiz]"
        };

        let sample_exceptions = exceptions.into_iter().take(20).collect();

        PreflightReport {
            target: *self,
            start,
            end,
            tested_count,
            exceptions_count,
            sample_exceptions,
            passed,
            verdict,
        }
    }

    /// Deterministically verifies whether candidate `n` can be represented in conjectured form.
    pub fn verify_candidate_deterministic(&self, n: u64) -> VerificationResult {
        if self.is_known_small_exception(n) {
            return VerificationResult::KnownSmallException;
        }

        let shifts = self.generate_shifts(n);
        let kind = self.prime_kind();

        for &s in &shifts {
            if s >= n {
                break;
            }
            let p = n - s;
            let is_valid_prime = match kind {
                PrimeKind::Prime => is_prime(p),
                PrimeKind::TwinPrime => is_twin_prime(p),
                PrimeKind::SophieGermain => is_sophie_germain(p),
            };

            if is_valid_prime {
                return VerificationResult::Represented {
                    prime_part: p,
                    shift: s,
                    details: format!("n = {} = {} + {}", n, p, s),
                };
            }
        }

        VerificationResult::Counterexample {
            failed_shifts: shifts,
        }
    }

    /// Formats verified counterexample proof certificate.
    pub fn format_counterexample_proof(&self, n: u64) -> String {
        let shifts = self.generate_shifts(n);
        if let Some(report) = verify_counterexample(n, &shifts) {
            return report.format_certificate(self.name());
        }

        // Fallback formatting
        let mut proof = String::new();
        proof.push_str("================================================================================\n");
        proof.push_str("MATHEMATICAL COUNTEREXAMPLE CERTIFICATE\n");
        proof.push_str("================================================================================\n");
        proof.push_str(&format!("Target Conjecture: {}\n", self.name()));
        proof.push_str(&format!("Counterexample Integer: n = {}\n\n", n));
        for s in shifts {
            if s >= n {
                continue;
            }
            let diff = n - s;
            let fact_str = format_factorization(diff);
            proof.push_str(&format!("  s = {:<12} => n - s = {:<12} = {}\n", s, diff, fact_str));
        }
        proof.push_str("================================================================================\n");
        proof
    }
}

/// Computes central binomial coefficient binom(2k, k) with overflow checking.
pub fn central_binomial(k: u64) -> u64 {
    if k == 0 {
        return 1;
    }
    if k > 30 {
        return 0; // overflows 64-bit unsigned integer (binom(62, 31) > 2^64)
    }
    let mut result: u128 = 1;
    for i in 1..=k {
        result = result * (k + i) as u128 / i as u128;
    }
    if result > u64::MAX as u128 {
        0
    } else {
        result as u64
    }
}

/// Generates Fibonacci numbers up to `limit`.
pub fn generate_fibonacci(limit: u64) -> Vec<u64> {
    let mut fibs = vec![1, 2];
    let mut a = 1u64;
    let mut b = 2u64;
    loop {
        match a.checked_add(b) {
            Some(next) if next <= limit => {
                fibs.push(next);
                a = b;
                b = next;
            }
            _ => break,
        }
    }
    fibs
}

/// Generates Lucas numbers up to `limit` (L_1=1, L_2=3, L_3=4, L_4=7, ...).
pub fn generate_lucas(limit: u64) -> Vec<u64> {
    let mut lucas = vec![1, 3];
    let mut a = 1u64;
    let mut b = 3u64;
    loop {
        match a.checked_add(b) {
            Some(next) if next <= limit => {
                lucas.push(next);
                a = b;
                b = next;
            }
            _ => break,
        }
    }
    lucas
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_central_binomial() {
        assert_eq!(central_binomial(1), 2);
        assert_eq!(central_binomial(2), 6);
        assert_eq!(central_binomial(3), 20);
        assert_eq!(central_binomial(4), 70);
        assert_eq!(central_binomial(5), 252);
    }

    #[test]
    fn test_fibonacci_lucas() {
        let fibs = generate_fibonacci(100);
        assert_eq!(&fibs[..6], &[1, 2, 3, 5, 8, 13]);
        let lucas = generate_lucas(100);
        assert_eq!(&lucas[..6], &[1, 3, 4, 7, 11, 18]);
    }

    #[test]
    fn test_preflight_checks() {
        let targets = [
            (TargetConjecture::Binom2aPrime, true),
            (TargetConjecture::Twin2a3b, true),
            (TargetConjecture::Sophie2a3b, true),
            (TargetConjecture::Binom2aTwin, false), // Correctly identified as Trivial / Density Insufficient
            (TargetConjecture::Erdos2a2bTwin, true),
        ];
        for (t, should_pass) in targets {
            let report = t.run_preflight_check(10_000, 20_000);
            println!("{}: passed={}, exceptions={}", t.short_code(), report.passed, report.exceptions_count);
            assert_eq!(report.passed, should_pass, "Target {} preflight status mismatch", t.short_code());
        }
    }
}
