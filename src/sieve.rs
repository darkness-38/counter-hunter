//! High-performance Segmented Bit-Sieve of Eratosthenes with Wheel-30 & L1D Cache Alignment.
//!
//! # Architecture:
//! - **Wheel-30 Modular Design:** Eliminates multiples of 2, 3, and 5.
//!   Only 8 coprime residues (1, 7, 11, 13, 17, 19, 23, 29) are tracked in bits.
//!   1 byte represents exactly 30 integers (8 bits = 8 coprime candidates).
//! - **L1D Cache Alignment (32 KiB):** Each segment is exactly 32,768 bytes (32 KiB),
//!   representing 983,040 integers (32,768 * 30), fitting precisely into CPU L1 Data Cache.
//! - **Cache-Line Alignment:** Memory structures are aligned to 64-byte CPU cache lines with `#[repr(align(64))]`.
//! - **Pre-Sieve Offset Tracking:** Base primes below `sqrt(R)` are generated once. Each prime's multiple
//!   offsets and wheel residue steps are tracked across consecutive 32 KiB segments without recalculation.

pub const L1D_CACHE_SIZE: usize = 32_768; // 32 KiB
pub const WHEEL_MODULUS: usize = 30;
pub const NUMBERS_PER_SEGMENT: u64 = (L1D_CACHE_SIZE * WHEEL_MODULUS) as u64; // 983,040 integers

/// The 8 coprime residues modulo 30.
pub const WHEEL30_RESIDUES: [u8; 8] = [1, 7, 11, 13, 17, 19, 23, 29];

/// Steps between successive coprime residues modulo 30:
/// (7-1=6, 11-7=4, 13-11=2, 17-13=4, 19-17=2, 23-19=4, 29-23=6, 31-29=2).
pub const WHEEL30_STEPS: [u8; 8] = [6, 4, 2, 4, 2, 4, 6, 2];

/// Maps a residue modulo 30 to its bit index 0..7 (or 0xFF if not coprime to 30).
pub const RESIDUE_TO_BIT: [u8; 30] = {
    let mut table = [0xFFu8; 30];
    table[1] = 0;
    table[7] = 1;
    table[11] = 2;
    table[13] = 3;
    table[17] = 4;
    table[19] = 5;
    table[23] = 6;
    table[29] = 7;
    table
};

/// Maps a residue modulo 30 to its index 0..7 in `WHEEL30_RESIDUES`.
pub const RESIDUE_TO_INDEX: [u8; 30] = {
    let mut table = [0xFFu8; 30];
    table[1] = 0;
    table[7] = 1;
    table[11] = 2;
    table[13] = 3;
    table[17] = 4;
    table[19] = 5;
    table[23] = 6;
    table[29] = 7;
    table
};

/// 64-byte cache-line aligned byte buffer for exactly 32 KiB L1 Data Cache.
#[repr(align(64))]
#[derive(Clone)]
pub struct AlignedBytes(pub [u8; L1D_CACHE_SIZE]);

impl Default for AlignedBytes {
    fn default() -> Self {
        Self([0xFFu8; L1D_CACHE_SIZE])
    }
}

/// A 32 KiB Wheel-30 sieved segment representing an interval `[start, end)`.
/// Size: 32 KiB buffer + 16 bytes metadata = exactly aligned to 64-byte cache lines.
#[repr(align(64))]
#[derive(Clone)]
pub struct Wheel30Segment {
    /// 32 KiB bitset, strictly aligned to 64-byte CPU cache lines.
    pub bytes: AlignedBytes,
    /// Starting integer (always a multiple of 30).
    pub start: u64,
    /// Ending integer (start + 983,040).
    pub end: u64,
}

impl Wheel30Segment {
    /// Creates a new uninitialized/all-ones segment for `start`.
    pub fn new(start: u64) -> Self {
        assert_eq!(start % 30, 0, "Segment start must be a multiple of 30");
        Self {
            bytes: AlignedBytes::default(),
            start,
            end: start + NUMBERS_PER_SEGMENT,
        }
    }

    /// Checks if a number in `[start, end)` is prime.
    #[inline(always)]
    pub fn is_prime(&self, n: u64) -> bool {
        if n < self.start || n >= self.end {
            return false;
        }
        if n == 2 || n == 3 || n == 5 {
            return true;
        }
        if n < 2 || n.is_multiple_of(2) || n.is_multiple_of(3) || n.is_multiple_of(5) {
            return false;
        }

        let offset = (n - self.start) as usize;
        let byte_idx = offset / 30;
        let rem = offset % 30;
        let bit = RESIDUE_TO_BIT[rem];
        if bit == 0xFF {
            return false;
        }
        (self.bytes.0[byte_idx] & (1u8 << bit)) != 0
    }

    /// Iterates over all prime numbers in this segment.
    pub fn for_each_prime<F: FnMut(u64)>(&self, mut f: F) {
        if self.start <= 2 && self.end > 2 {
            f(2);
        }
        if self.start <= 3 && self.end > 3 {
            f(3);
        }
        if self.start <= 5 && self.end > 5 {
            f(5);
        }

        for (byte_idx, &b) in self.bytes.0.iter().enumerate() {
            if b == 0 {
                continue;
            }
            let base = self.start + (byte_idx as u64) * 30;
            for (bit_idx, &res) in WHEEL30_RESIDUES.iter().enumerate() {
                if (b & (1u8 << bit_idx)) != 0 {
                    let p = base + (res as u64);
                    if p >= self.start && p < self.end && p > 1 {
                        f(p);
                    }
                }
            }
        }
    }

    /// Counts total primes in this segment.
    pub fn count_primes(&self) -> usize {
        let mut count = 0;
        if self.start <= 2 && self.end > 2 {
            count += 1;
        }
        if self.start <= 3 && self.end > 3 {
            count += 1;
        }
        if self.start <= 5 && self.end > 5 {
            count += 1;
        }

        // Count ones in the bitset
        for &b in self.bytes.0.iter() {
            count += b.count_ones() as usize;
        }

        // 1 was cleared if start == 0, but if start > 0 and 1 was not in range,
        // all bits represent numbers >= start > 5 coprime to 30.
        count
    }
}

/// Tracks pre-sieve state of a base prime across consecutive 32 KiB segments.
#[derive(Clone, Copy, Debug)]
pub struct BasePrimeState {
    pub prime: u64,
    pub next_multiple: u64,
    pub q_step_idx: u8,
}

/// The Wheel-30 Pre-Sieve Engine.
/// Generates base primes up to `sqrt(max_n)` and drives 32 KiB L1D-cache-friendly segmented sieving.
pub struct Wheel30Sieve {
    pub base_primes: Vec<u64>,
    pub max_n: u64,
}

impl Wheel30Sieve {
    /// Initializes base primes sufficient to sieve up to `max_n`.
    pub fn new(max_n: u64) -> Self {
        let sqrt_limit = ((max_n as f64).sqrt().ceil() as u64).max(10_000);
        let base_primes = generate_base_primes(sqrt_limit);
        Self { base_primes, max_n }
    }

    /// Initializes pre-sieve states for all base primes >= 7 for an interval starting at `start`.
    pub fn init_base_states(&self, start: u64) -> Vec<BasePrimeState> {
        let mut states = Vec::with_capacity(self.base_primes.len());

        for &p in &self.base_primes {
            if p < 7 {
                continue; // 2, 3, 5 are handled by Wheel-30 structure
            }

            // Find smallest multiple m = p * q >= max(p*p, start) with gcd(q, 30) == 1
            let min_target = (p * p).max(start);
            let mut q = min_target.div_ceil(p);
            if q < p {
                q = p;
            }

            // Advance q to the nearest integer coprime to 30
            while q % 2 == 0 || q % 3 == 0 || q % 5 == 0 {
                q += 1;
            }

            let m = p * q;
            let rem = (q % 30) as usize;
            let q_step_idx = RESIDUE_TO_INDEX[rem];

            states.push(BasePrimeState {
                prime: p,
                next_multiple: m,
                q_step_idx,
            });
        }

        states
    }

    /// Returns a segment-based bitset iterator producing consecutive 32 KiB Wheel-30 segments.
    pub fn iter_segments(&self, start: u64, end: u64) -> Wheel30SegmentIter {
        let aligned_start = (start / 30) * 30;
        let states = self.init_base_states(aligned_start);
        Wheel30SegmentIter {
            states,
            current_start: aligned_start,
            end,
        }
    }

    /// Sieves a single standalone 32 KiB segment starting at `seg_start`.
    pub fn sieve_segment(&self, seg_start: u64) -> Box<Wheel30Segment> {
        assert_eq!(seg_start % 30, 0, "Segment start must be multiple of 30");
        let mut segment = Box::new(Wheel30Segment::new(seg_start));
        let seg_end = segment.end;

        // Clear residue 1 (number 1) if in segment 0
        if seg_start == 0 {
            segment.bytes.0[0] &= !(1u8 << RESIDUE_TO_BIT[1]);
        }

        let sqrt_end = (seg_end as f64).sqrt().ceil() as u64;

        for &p in &self.base_primes {
            if p < 7 {
                continue;
            }
            if p > sqrt_end {
                break;
            }

            let min_target = (p * p).max(seg_start);
            let mut q = min_target.div_ceil(p);
            if q < p {
                q = p;
            }
            while q % 2 == 0 || q % 3 == 0 || q % 5 == 0 {
                q += 1;
            }

            let mut m = p * q;
            let mut step_idx = RESIDUE_TO_INDEX[(q % 30) as usize] as usize;

            let p_steps = [
                6 * p, 4 * p, 2 * p, 4 * p, 2 * p, 4 * p, 6 * p, 2 * p,
            ];

            while m < seg_end {
                let offset = (m - seg_start) as usize;
                let byte_idx = offset / 30;
                let rem = offset % 30;
                let bit = RESIDUE_TO_BIT[rem];
                segment.bytes.0[byte_idx] &= !(1u8 << bit);

                m += p_steps[step_idx];
                step_idx = (step_idx + 1) & 7;
            }
        }

        segment
    }
}

/// Segment-based iterator streaming 32 KiB Wheel-30 sieved segments across L1D cache.
pub struct Wheel30SegmentIter {
    pub states: Vec<BasePrimeState>,
    pub current_start: u64,
    pub end: u64,
}

impl Iterator for Wheel30SegmentIter {
    type Item = Box<Wheel30Segment>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_start >= self.end {
            return None;
        }

        let seg_start = self.current_start;
        let mut segment = Box::new(Wheel30Segment::new(seg_start));
        let seg_end = segment.end;

        // Number 1 is not prime
        if seg_start == 0 {
            segment.bytes.0[0] &= !(1u8 << RESIDUE_TO_BIT[1]);
        }

        let sqrt_end = (seg_end as f64).sqrt().ceil() as u64;

        // Strike out multiples of base primes using precomputed state
        for state in self.states.iter_mut() {
            let p = state.prime;
            if p > sqrt_end && state.next_multiple >= seg_end {
                break;
            }

            if state.next_multiple >= seg_end {
                continue;
            }

            let mut m = state.next_multiple;
            let mut step_idx = state.q_step_idx as usize;

            let p_steps = [
                6 * p, 4 * p, 2 * p, 4 * p, 2 * p, 4 * p, 6 * p, 2 * p,
            ];

            while m < seg_end {
                let offset = (m - seg_start) as usize;
                let byte_idx = offset / 30;
                let rem = offset % 30;
                let bit = RESIDUE_TO_BIT[rem];
                segment.bytes.0[byte_idx] &= !(1u8 << bit);

                m += p_steps[step_idx];
                step_idx = (step_idx + 1) & 7;
            }

            state.next_multiple = m;
            state.q_step_idx = step_idx as u8;
        }

        self.current_start = seg_end;
        Some(segment)
    }
}

// ============================================================================
// Legacy Segmented Bit-Sieve (Dense 1-bit per integer representation)
// Preserved for dense Shift-OR operations and backward compatibility
// ============================================================================

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
    pub fn new(start: u64, end: u64) -> Self {
        assert!(end >= start, "end must be >= start");
        let len = (end - start) as usize;
        let num_words = len.div_ceil(64);
        Self {
            start,
            end,
            words: vec![0u64; num_words],
        }
    }

    #[inline(always)]
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

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

    pub fn compute_twin_primes(&self) -> Self {
        let mut twin_bitset = Self::new(self.start, self.end);
        let n_words = self.words.len();
        if n_words == 0 {
            return twin_bitset;
        }

        for i in 0..n_words {
            let curr = self.words[i];
            let next = if i + 1 < n_words { self.words[i + 1] } else { 0 };
            let shifted = (curr >> 2) | (next << 62);
            twin_bitset.words[i] = curr & shifted;
        }

        twin_bitset
    }
}

pub struct SegmentedSieve {
    pub base_primes: Vec<u64>,
}

impl SegmentedSieve {
    pub fn new(max_n: u64) -> Self {
        let sqrt_limit = ((max_n as f64).sqrt().ceil() as u64).max(10_000);
        let base_primes = generate_base_primes(sqrt_limit);
        Self { base_primes }
    }

    pub fn sieve_range(&self, start: u64, end: u64) -> PrimeBitset {
        let mut bitset = PrimeBitset::new(start, end);
        let len = (end - start) as usize;
        if len == 0 {
            return bitset;
        }

        let odd_pattern = if start.is_multiple_of(2) {
            0xAAAAAAAAAAAAAAAAu64
        } else {
            0x5555555555555555u64
        };

        for w in bitset.words.iter_mut() {
            *w = odd_pattern;
        }

        if start == 0 && end > 0 {
            bitset.set_prime(0, false);
        }
        if start <= 1 && end > 1 {
            bitset.set_prime(1, false);
        }
        if start <= 2 && end > 2 {
            bitset.set_prime(2, true);
        }

        let sqrt_end = (end as f64).sqrt().ceil() as u64;

        for &p in &self.base_primes {
            if p == 2 {
                continue;
            }
            if p > sqrt_end {
                break;
            }

            let mut first_mult = start.div_ceil(p) * p;
            if first_mult < p * p {
                first_mult = p * p;
            }

            if first_mult.is_multiple_of(2) {
                first_mult += p;
            }

            let step = 2 * p;
            let mut m = first_mult;
            while m < end {
                let offset = (m - start) as usize;
                bitset.words[offset >> 6] &= !(1u64 << (offset & 63));
                m += step;
            }
        }

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
    fn test_wheel30_alignment_and_size() {
        assert_eq!(std::mem::size_of::<AlignedBytes>(), 32_768);
        assert_eq!(std::mem::align_of::<AlignedBytes>(), 64);
        assert_eq!(std::mem::align_of::<Wheel30Segment>(), 64);

        let seg = Box::new(Wheel30Segment::new(0));
        let bytes_addr = seg.bytes.0.as_ptr() as usize;
        assert_eq!(bytes_addr % 64, 0, "Segment bytes must be 64-byte cache-line aligned");
        assert_eq!(seg.end - seg.start, 983_040);
    }

    #[test]
    fn test_wheel30_single_segment() {
        let sieve = Wheel30Sieve::new(1_000_000);
        let seg = sieve.sieve_segment(0);

        // Verify known primes and composites below 1000
        for n in 0..1000 {
            assert_eq!(
                seg.is_prime(n),
                is_prime(n),
                "Primality mismatch at n = {}",
                n
            );
        }

        // Verify prime count up to 100,000 matches pi(100,000) = 9,592
        let mut count_100k = 0;
        seg.for_each_prime(|p| {
            if p <= 100_000 {
                count_100k += 1;
            }
        });
        assert_eq!(count_100k, 9_592, "pi(100,000) must equal 9,592");
    }

    #[test]
    fn test_wheel30_segment_iterator() {
        let sieve = Wheel30Sieve::new(2_000_000);
        let mut total_primes = 0;

        for seg in sieve.iter_segments(0, 1_000_000) {
            seg.for_each_prime(|p| {
                if p <= 1_000_000 {
                    assert!(is_prime(p), "Reported prime {} is not prime", p);
                    total_primes += 1;
                }
            });
        }

        assert_eq!(total_primes, 78_498, "pi(1,000,000) must equal 78,498");
    }

    #[test]
    fn test_wheel30_arbitrary_start() {
        let sieve = Wheel30Sieve::new(200_000_000);
        // Test a high range [100_000_000, 101_000_000]
        let start = 100_000_020;
        let mut seg_iter = sieve.iter_segments(start, start + 50_000);
        let seg = seg_iter.next().unwrap();

        for n in (start..start + 5_000).step_by(1) {
            assert_eq!(
                seg.is_prime(n),
                is_prime(n),
                "Primality mismatch at high range n = {}",
                n
            );
        }
    }
}
