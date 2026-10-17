//! Deterministic Miller-Rabin primality test and exact integer factorization engine.
//! Used for 100% rigorous mathematical verification and counterexample proof certificates.

/// Modular multiplication `(a * b) % m` avoiding 64-bit overflow using 128-bit arithmetic.
#[inline(always)]
pub fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Modular exponentiation `(base ^ exp) % m`.
#[inline(always)]
pub fn pow_mod(mut base: u64, mut exp: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let mut result: u64 = 1;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mul_mod(result, base, m);
        }
        base = mul_mod(base, base, m);
        exp >>= 1;
    }
    result
}

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

/// Proven deterministic Miller-Rabin bases covering all 64-bit integers (`< 2^64`).
/// Reference: Jim Sinclair / Galway & Feitsma (verified against entire 64-bit space).
const MILLER_RABIN_BASES_64: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// 100% deterministic primality test for all `n < 2^64`.
#[inline]
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 || n == 5 || n == 7 || n == 11 || n == 13 || n == 17 || n == 19 || n == 23 || n == 29 || n == 31 || n == 37 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 || n % 5 == 0 {
        return false;
    }

    // Write n - 1 as 2^s * d with d odd
    let d_val = n - 1;
    let s = d_val.trailing_zeros();
    let d = d_val >> s;

    'base_loop: for &a in &MILLER_RABIN_BASES_64 {
        if a >= n {
            break;
        }

        let mut x = pow_mod(a, d, n);
        if x == 1 || x == n - 1 {
            continue 'base_loop;
        }

        for _ in 1..s {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'base_loop;
            }
            if x == 1 {
                return false;
            }
        }

        return false;
    }

    true
}

/// Checks whether `p` is a twin prime (specifically: `p` is prime and `p + 2` is prime).
#[inline]
pub fn is_twin_prime(p: u64) -> bool {
    is_prime(p) && is_prime(p + 2)
}

/// Checks whether `p` is a Sophie Germain prime (`p` is prime and `2p + 1` is prime).
#[inline]
pub fn is_sophie_germain(p: u64) -> bool {
    if p > (u64::MAX - 1) / 2 {
        return false;
    }
    is_prime(p) && is_prime(2 * p + 1)
}

/// Pollard's Rho algorithm with Brent's cycle detection optimization.
/// Finds a non-trivial factor of a composite integer `n`.
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

    let mut c = 1u64;
    loop {
        let f = |x: u64| -> u64 {
            let sq = mul_mod(x, x, n);
            let res = sq.wrapping_add(c);
            if res >= n { res - n } else { res }
        };

        let mut y = 2u64;
        let mut x = 2u64;
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
            let diff = if x > y { x - y } else { y - x };
            d = gcd(diff, n);
        }

        if d != n {
            return d;
        }
        c += 1; // Try different polynomial constant
    }
}

/// Computes complete prime factorization of `n`.
/// Returns a sorted list of `(prime_factor, multiplicity)` pairs.
pub fn factorize(mut n: u64) -> Vec<(u64, u32)> {
    let mut factors: Vec<u64> = Vec::new();
    if n <= 1 {
        return Vec::new();
    }

    // Trial division for small primes
    for &p in &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97] {
        while n % p == 0 {
            factors.push(p);
            n /= p;
        }
    }

    // Sieve or Pollard for remaining part
    let mut queue = vec![n];
    while let Some(curr) = queue.pop() {
        if curr == 1 {
            continue;
        }
        if is_prime(curr) {
            factors.push(curr);
        } else {
            let divisor = pollard_brent_factor(curr);
            queue.push(divisor);
            queue.push(curr / divisor);
        }
    }

    factors.sort_unstable();

    // Group into (factor, count)
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

/// Formats prime factorization into human-readable algebraic proof, e.g. `2^3 * 5 * 13^2`.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_miller_rabin() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
        assert!(is_prime(997));
        assert!(!is_prime(999));
        assert!(is_prime(1_000_000_007));
        assert!(is_prime(10_000_000_019));
    }

    #[test]
    fn test_twin_prime() {
        assert!(is_twin_prime(3)); // (3, 5)
        assert!(is_twin_prime(5)); // (5, 7)
        assert!(is_twin_prime(11)); // (11, 13)
        assert!(!is_twin_prime(7)); // (7, 9 - not prime)
    }

    #[test]
    fn test_sophie_germain() {
        assert!(is_sophie_germain(2)); // 2*2+1 = 5
        assert!(is_sophie_germain(3)); // 2*3+1 = 7
        assert!(is_sophie_germain(5)); // 2*5+1 = 11
        assert!(is_sophie_germain(11)); // 2*11+1 = 23
        assert!(!is_sophie_germain(7)); // 2*7+1 = 15 (composite)
    }

    #[test]
    fn test_factorization() {
        assert_eq!(factorize(60), vec![(2, 2), (3, 1), (5, 1)]);
        assert_eq!(format_factorization(60), "2^2 * 3 * 5");
        assert_eq!(format_factorization(17), "17 [PRIME]");
    }
}
