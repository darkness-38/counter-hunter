//! High-performance Segmented Bit-Sieve of Eratosthenes.
//! Designed for cache-friendly parallel sieve operations adhering to strict memory limits.

/// Precomputed base primes up to `limit` using a fast odd-only bitset.
pub fn generate_base_primes(limit: u64) -> Vec<u64> {
    if limit < 2 {
        return Vec::new();
    }
    let mut primes = vec![2];
    if limit == 2 {
        return primes;
    }

    let num_odds = ((limit - 1) / 2) as usize;
    // composite_bits[i] = true if 2*i + 3 is composite
    let mut composite_bits = vec![false; num_odds];
    let sqrt_limit = (limit as f64).sqrt() as u64;

    for i in 0..num_odds {
        if !composite_bits[i] {
            let p = (2 * i + 3) as u64;
            if p <= sqrt_limit {
                let start_idx = ((p * p - 3) / 2) as usize;
                let step = p as usize;
                for j in (start_idx..num_odds).step_by(step) {
                    composite_bits[j] = true;
                }
            }
            primes.push(p);
        }
    }

    primes
}

/// Bitset representing primality of integers in an interval `[start, end)`.
/// Memory usage: 1 bit per integer -> 8 MB per 67.1M numbers (`2^26`).
#[derive(Clone)]
pub struct PrimeBitset {
    pub start: u64,
    pub end: u64,
    pub words: Vec<u64>,
}

impl PrimeBitset {
    /// Creates an empty bitset for `[start, end)`.
    pub fn new(start: u64, end: u64) -> Self {
        assert!(end >= start, "end must be >= start");
        let len = (end - start) as usize;
        let num_words = (len + 63) / 64;
        Self {
            start,
            end,
            words: vec![0u64; num_words],
        }
    }

    /// Number of integers represented in this bitset.
    #[inline(always)]
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// Checks if a number in `[start, end)` is prime.
    #[inline(always)]
    pub fn is_prime(&self, n: u64) -> bool {
        if n < self.start || n >= self.end {
            return false;
        }
        let offset = (n - self.start) as usize;
        let word_idx = offset >> 6;
        let bit_idx = offset & 63;
        (self.words[word_idx] & (1u64 << bit_idx)) != 0
    }

    /// Sets or clears the prime bit for `n`.
    #[inline(always)]
    pub fn set_prime(&mut self, n: u64, is_prime: bool) {
        if n < self.start || n >= self.end {
            return;
        }
        let offset = (n - self.start) as usize;
        let word_idx = offset >> 6;
        let bit_idx = offset & 63;
        if is_prime {
            self.words[word_idx] |= 1u64 << bit_idx;
        } else {
            self.words[word_idx] &= !(1u64 << bit_idx);
        }
    }

    /// Extracts a 64-bit word starting at bit `bit_offset`.
    /// Handles unaligned boundaries by combining two adjacent words.
    #[inline(always)]
    pub fn get_u64_at(&self, bit_offset: usize) -> u64 {
        let word_idx = bit_offset >> 6;
        let shift = bit_offset & 63;

        if word_idx >= self.words.len() {
            return 0;
        }

        let low = self.words[word_idx] >> shift;
        if shift == 0 || word_idx + 1 >= self.words.len() {
            low
        } else {
            low | (self.words[word_idx + 1] << (64 - shift))
        }
    }

    /// Computes twin prime bitset where bit `i` is 1 if both `start + i` and `start + i + 2` are prime.
    pub fn compute_twin_primes(&self) -> Self {
        let mut twin_bitset = Self::new(self.start, self.end);
        let n_words = self.words.len();
        if n_words == 0 {
            return twin_bitset;
        }

        for i in 0..n_words {
            let curr = self.words[i];
            let next = if i + 1 < n_words { self.words[i + 1] } else { 0 };
            // Shift curr and next by 2 bits to align with prime + 2
            let shifted = (curr >> 2) | (next << 62);
            twin_bitset.words[i] = curr & shifted;
        }

        // 2 is prime but 2+2=4 is not, twin bit for 2 would be 0 anyway.
        twin_bitset
    }
}

/// Fast segmented sieve engine holding precomputed base primes.
pub struct SegmentedSieve {
    pub base_primes: Vec<u64>,
}

impl SegmentedSieve {
    /// Initializes base primes sufficient to sieve up to `max_n`.
    pub fn new(max_n: u64) -> Self {
        let sqrt_limit = ((max_n as f64).sqrt().ceil() as u64).max(10_000);
        let base_primes = generate_base_primes(sqrt_limit);
        Self { base_primes }
    }

    /// Sieves primes in `[start, end)` into a `PrimeBitset`.
    pub fn sieve_range(&self, start: u64, end: u64) -> PrimeBitset {
        let mut bitset = PrimeBitset::new(start, end);
        let len = (end - start) as usize;
        if len == 0 {
            return bitset;
        }

        // Initialize pattern: odd numbers = 1, even numbers = 0
        let odd_pattern = if start % 2 == 0 {
            0xAAAAAAAAAAAAAAAAu64 // even bit = 0, odd bit = 1
        } else {
            0x5555555555555555u64 // even bit = 1, odd bit = 0
        };

        for w in bitset.words.iter_mut() {
            *w = odd_pattern;
        }

        // Fix boundaries: 0 and 1 are not prime; 2 is prime
        if start <= 0 && end > 0 {
            bitset.set_prime(0, false);
        }
        if start <= 1 && end > 1 {
            bitset.set_prime(1, false);
        }
        if start <= 2 && end > 2 {
            bitset.set_prime(2, true);
        }

        let sqrt_end = (end as f64).sqrt().ceil() as u64;

        // Eliminate multiples of each base odd prime
        for &p in &self.base_primes {
            if p == 2 {
                continue;
            }
            if p > sqrt_end {
                break;
            }

            // Smallest multiple >= max(p*p, start)
            let mut first_mult = ((start + p - 1) / p) * p;
            if first_mult < p * p {
                first_mult = p * p;
            }

            // Ensure first multiple is odd
            if first_mult % 2 == 0 {
                first_mult += p;
            }

            let step = 2 * p; // step by 2p to only touch odd numbers
            let mut m = first_mult;
            while m < end {
                let offset = (m - start) as usize;
                bitset.words[offset >> 6] &= !(1u64 << (offset & 63));
                m += step;
            }
        }

        // Mask trailing unused bits in the last word
        let remainder = len % 64;
        if remainder != 0 {
            let last_idx = bitset.words.len() - 1;
            let mask = (1u64 << remainder) - 1;
            bitset.words[last_idx] &= mask;
        }

        bitset
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verifier::is_prime;

    #[test]
    fn test_base_primes() {
        let primes = generate_base_primes(30);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn test_sieve_range() {
        let sieve = SegmentedSieve::new(100_000);
        let bitset = sieve.sieve_range(100, 200);

        for n in 100..200 {
            assert_eq!(
                bitset.is_prime(n),
                is_prime(n),
                "Mismatch at n = {}",
                n
            );
        }
    }

    #[test]
    fn test_twin_primes_bitset() {
        let sieve = SegmentedSieve::new(100_000);
        let primes = sieve.sieve_range(1, 100);
        let twins = primes.compute_twin_primes();

        assert!(twins.is_prime(3)); // (3, 5)
        assert!(twins.is_prime(5)); // (5, 7)
        assert!(twins.is_prime(11)); // (11, 13)
        assert!(twins.is_prime(17)); // (17, 19)
        assert!(twins.is_prime(29)); // (29, 31)
        assert!(!twins.is_prime(7)); // (7, 9) is not twin prime
        assert!(!twins.is_prime(13)); // (13, 15) is not twin prime
    }
}
