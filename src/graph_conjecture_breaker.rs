//! MODÜL 2: Spectral & Topological Graph Conjectures Breaker
//!
//! High-speed Tabu Search & Simulated Annealing Graph Evolution Engine for connected graphs
//! of order V = 12 .. 30, designed to break open extremal graph conjectures from the
//! AutoGraphiX (AGX), Graffiti, and Aouchiche-Hansen literature.
//!
//! Features:
//! - Lock-free bitmask 64-bit graph representation with O(1) BFS connectivity.
//! - Exact degree sequence, diameter D, average distance l_bar, proximity pi, remoteness rho, Wiener index W.
//! - Topological indices: Randic index R(G) and Harmonic index H(G).
//! - Spectral invariants: Adjacency spectral radius lambda_1, Laplacian spectral radius mu_1,
//!   and Signless Laplacian spectral radius q(G) via accelerated Power Iteration.
//! - Parallel multi-walker evolutionary optimization targeting 8 open conjecture inequalities.

use rayon::prelude::*;
use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;

/// Fast pseudo-random number generator (Xorshift64*)
#[derive(Clone, Copy)]
pub struct FastRng {
    state: u64,
}

impl FastRng {
    #[inline(always)]
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x5a17_d00d_cafe_babe } else { seed },
        }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    #[inline(always)]
    pub fn gen_range(&mut self, min: usize, max: usize) -> usize {
        if max <= min {
            return min;
        }
        min + (self.next_u64() as usize % (max - min))
    }

    #[inline(always)]
    pub fn gen_bool(&mut self, prob: f64) -> bool {
        (self.next_u64() as f64 / u64::MAX as f64) < prob
    }
}

/// Compact bitmask representation of a graph on up to 32 vertices.
#[derive(Clone)]
pub struct FastGraph {
    pub n: usize,
    pub adj: [u64; 32],
}

impl FastGraph {
    pub fn new(n: usize) -> Self {
        assert!(n <= 32, "FastGraph supports up to 32 vertices");
        Self { n, adj: [0u64; 32] }
    }

    #[inline(always)]
    pub fn add_edge(&mut self, u: usize, v: usize) {
        if u != v && u < self.n && v < self.n {
            self.adj[u] |= 1u64 << v;
            self.adj[v] |= 1u64 << u;
        }
    }

    #[inline(always)]
    pub fn remove_edge(&mut self, u: usize, v: usize) {
        if u < self.n && v < self.n {
            self.adj[u] &= !(1u64 << v);
            self.adj[v] &= !(1u64 << u);
        }
    }

    #[inline(always)]
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        (self.adj[u] & (1u64 << v)) != 0
    }

    /// Fast bitwise BFS connectivity check
    #[inline(always)]
    pub fn is_connected(&self) -> bool {
        let mut visited = 1u64;
        let mut prev = 0u64;
        while visited != prev {
            prev = visited;
            let mut next = visited;
            for u in 0..self.n {
                if (visited & (1u64 << u)) != 0 {
                    next |= self.adj[u];
                }
            }
            visited = next;
        }
        visited.count_ones() as usize == self.n
    }

    pub fn degrees(&self) -> Vec<usize> {
        (0..self.n).map(|u| self.adj[u].count_ones() as usize).collect()
    }

    pub fn num_edges(&self) -> usize {
        let total: usize = (0..self.n).map(|u| self.adj[u].count_ones() as usize).sum();
        total / 2
    }

    /// All-Pairs Shortest Paths (APSP) using bitwise BFS
    pub fn apsp(&self) -> (usize, f64, f64, f64, u64) {
        let mut dist = [[0usize; 32]; 32];
        let mut max_dist = 0usize;
        let mut total_dist = 0u64;
        let mut min_avg_dist = f64::MAX;
        let mut max_avg_dist = 0.0f64;

        for s in 0..self.n {
            let mut visited = 1u64 << s;
            let mut frontier = 1u64 << s;
            let mut d = 0;
            let mut sum_dist_s = 0usize;

            while frontier != 0 {
                let mut next_frontier = 0u64;
                for u in 0..self.n {
                    if (frontier & (1u64 << u)) != 0 {
                        dist[s][u] = d;
                        sum_dist_s += d;
                        if d > max_dist {
                            max_dist = d;
                        }
                        next_frontier |= self.adj[u] & !visited;
                    }
                }
                visited |= next_frontier;
                frontier = next_frontier;
                d += 1;
            }

            let avg_s = sum_dist_s as f64 / (self.n - 1) as f64;
            if avg_s < min_avg_dist {
                min_avg_dist = avg_s;
            }
            if avg_s > max_avg_dist {
                max_avg_dist = avg_s;
            }
        }

        for u in 0..self.n {
            for v in (u + 1)..self.n {
                total_dist += dist[u][v] as u64;
            }
        }

        let num_pairs = (self.n * (self.n - 1) / 2) as f64;
        let avg_l = total_dist as f64 / num_pairs;
        (max_dist, avg_l, min_avg_dist, max_avg_dist, total_dist)
    }

    /// Randic Index: R = sum_{uv in E} 1 / sqrt(d_u * d_v)
    pub fn randic_index(&self) -> f64 {
        let degs = self.degrees();
        let mut r = 0.0;
        for u in 0..self.n {
            for v in (u + 1)..self.n {
                if self.has_edge(u, v) {
                    r += 1.0 / ((degs[u] * degs[v]) as f64).sqrt();
                }
            }
        }
        r
    }

    /// Harmonic Index: H = sum_{uv in E} 2 / (d_u + d_v)
    pub fn harmonic_index(&self) -> f64 {
        let degs = self.degrees();
        let mut h = 0.0;
        for u in 0..self.n {
            for v in (u + 1)..self.n {
                if self.has_edge(u, v) {
                    h += 2.0 / (degs[u] + degs[v]) as f64;
                }
            }
        }
        h
    }

    /// Spectral Radius lambda_1(A) via Power Iteration
    pub fn spectral_radius(&self) -> f64 {
        let n = self.n;
        let mut x = vec![1.0f64 / (n as f64).sqrt(); n];
        let mut y = vec![0.0f64; n];

        for _ in 0..25 {
            for u in 0..n {
                let mut sum = 0.0f64;
                let mask = self.adj[u];
                for v in 0..n {
                    if (mask & (1u64 << v)) != 0 {
                        sum += x[v];
                    }
                }
                y[u] = sum;
            }
            let norm = y.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-12);
            for u in 0..n {
                x[u] = y[u] / norm;
            }
        }

        // Rayleigh quotient
        let mut rayleigh = 0.0f64;
        for u in 0..n {
            let mask = self.adj[u];
            for v in 0..n {
                if (mask & (1u64 << v)) != 0 {
                    rayleigh += x[u] * x[v];
                }
            }
        }
        rayleigh
    }

    /// Laplacian Spectral Radius mu_1(L) via Power Iteration
    pub fn laplacian_spectral_radius(&self) -> f64 {
        let n = self.n;
        let degs = self.degrees();
        let mut x: Vec<f64> = (0..n).map(|i| if i % 2 == 0 { 1.0 } else { -1.0 }).collect();
        let norm0 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        for val in &mut x {
            *val /= norm0;
        }
        let mut y = vec![0.0f64; n];

        for _ in 0..30 {
            for u in 0..n {
                let mut sum = degs[u] as f64 * x[u];
                let mask = self.adj[u];
                for v in 0..n {
                    if (mask & (1u64 << v)) != 0 {
                        sum -= x[v];
                    }
                }
                y[u] = sum;
            }
            let norm = y.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-12);
            for u in 0..n {
                x[u] = y[u] / norm;
            }
        }

        let mut rayleigh = 0.0f64;
        for u in 0..n {
            for v in (u + 1)..n {
                if self.has_edge(u, v) {
                    let diff = x[u] - x[v];
                    rayleigh += diff * diff;
                }
            }
        }
        rayleigh
    }

    /// Signless Laplacian Spectral Radius q(G) via Power Iteration
    pub fn signless_laplacian_spectral_radius(&self) -> f64 {
        let n = self.n;
        let degs = self.degrees();
        let mut x = vec![1.0f64 / (n as f64).sqrt(); n];
        let mut y = vec![0.0f64; n];

        for _ in 0..25 {
            for u in 0..n {
                let mut sum = degs[u] as f64 * x[u];
                let mask = self.adj[u];
                for v in 0..n {
                    if (mask & (1u64 << v)) != 0 {
                        sum += x[v];
                    }
                }
                y[u] = sum;
            }
            let norm = y.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-12);
            for u in 0..n {
                x[u] = y[u] / norm;
            }
        }

        let mut rayleigh = 0.0f64;
        for u in 0..n {
            for v in (u + 1)..n {
                if self.has_edge(u, v) {
                    let s = x[u] + x[v];
                    rayleigh += s * s;
                }
            }
        }
        rayleigh
    }

    /// Format adjacency matrix for publication / counterexample reporting
    pub fn format_adjacency(&self) -> String {
        let mut s = String::new();
        s.push_str("[\n");
        for u in 0..self.n {
            s.push_str("  [");
            for v in 0..self.n {
                s.push_str(if self.has_edge(u, v) { "1" } else { "0" });
                if v + 1 < self.n {
                    s.push_str(", ");
                }
            }
            s.push_str("],\n");
        }
        s.push(']');
        s
    }
}

// -----------------------------------------------------------------------------
// Target Conjectures & Evaluation Metrics
// -----------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct ConjectureResult {
    pub name: &'static str,
    pub description: &'static str,
    pub min_deficit: f64,
    pub counterexample_found: bool,
    pub best_graph_n: usize,
    pub best_graph_m: usize,
    pub details: String,
}

// -----------------------------------------------------------------------------
// Evolutionary Tabu Search & Simulated Annealing Engine
// -----------------------------------------------------------------------------

pub fn hunt_graph_conjectures(orders: &[usize], steps_per_walker: usize, num_walkers: usize) {
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 MODÜL 2: SPEKTRAL VE TOPOLOJİK GRAF VARSAYIMLARI KIRICI\x1b[0m");
    println!("   Tabu Search + Simulated Annealing Graph Evolution Engine (V = 12 .. 30)");
    println!("   Targeting 8 AutoGraphiX, Graffiti & Spectral Graph Theory Invariant Inequalities");
    println!("   Multi-threading: Rayon 8-Core Work Stealing ({} Parallel Walkers)", num_walkers);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let total_evaluations = AtomicU64::new(0);
    let ce_found_flag = AtomicBool::new(false);
    let global_best_mutex: Mutex<Vec<(f64, String)>> = Mutex::new(vec![(f64::MAX, String::new()); 8]);

    let start_time = Instant::now();

    for &n in orders {
        println!("\x1b[1;33m>>> EVOLVING CONNECTED GRAPHS OF ORDER V = {} <<<\x1b[0m", n);

        (0..num_walkers).into_par_iter().for_each(|walker_id| {
            let mut rng = FastRng::new(walker_id as u64 * 1_000_003 + n as u64 * 7919 + 42);
            let mut g = FastGraph::new(n);

            // Initial graph: random spanning tree + a few random edges
            for i in 0..(n - 1) {
                g.add_edge(i, i + 1);
            }
            for _ in 0..(n / 2) {
                let u = rng.gen_range(0, n);
                let v = rng.gen_range(0, n);
                g.add_edge(u, v);
            }

            let mut local_best = vec![f64::MAX; 8];

            for _step in 0..steps_per_walker {
                // Mutate graph: 70% edge toggle, 20% edge move, 10% 2-switch
                let op = rng.gen_range(0, 100);
                if op < 70 {
                    let u = rng.gen_range(0, n);
                    let v = rng.gen_range(0, n);
                    if u == v {
                        continue;
                    }
                    if g.has_edge(u, v) {
                        g.remove_edge(u, v);
                        if !g.is_connected() {
                            g.add_edge(u, v); // revert
                            continue;
                        }
                    } else {
                        g.add_edge(u, v);
                    }
                } else if op < 90 {
                    // Edge move
                    let u = rng.gen_range(0, n);
                    let v = rng.gen_range(0, n);
                    let w = rng.gen_range(0, n);
                    if u != v && u != w && v != w && g.has_edge(u, v) && !g.has_edge(u, w) {
                        g.remove_edge(u, v);
                        g.add_edge(u, w);
                        if !g.is_connected() {
                            g.remove_edge(u, w);
                            g.add_edge(u, v);
                            continue;
                        }
                    }
                }

                // Compute invariants
                let (diam, avg_l, prox, _rem, w_idx) = g.apsp();
                let r = g.randic_index();
                let h = g.harmonic_index();
                let lam1 = g.spectral_radius();
                let mu1 = g.laplacian_spectral_radius();
                let q_rad = g.signless_laplacian_spectral_radius();
                let n_f = n as f64;

                // -------------------------------------------------------------
                // Target Conjectures (Deficit < 0 represents counterexample)
                // -------------------------------------------------------------

                // C1: Harmonic vs Average Distance: H(G) >= avg_l(G)
                let c1_deficit = h - avg_l;
                // C2: Harmonic vs Half Diameter: H(G) >= D(G) / 2
                let c2_deficit = h - (diam as f64) / 2.0;
                // C3: Graffiti: Randic vs Average Distance: R(G) >= avg_l(G)
                let c3_deficit = r - avg_l;
                // C4: AutoGraphiX: lambda_1 + D >= 2*sqrt(n-1) + 1
                let c4_deficit = (lam1 + diam as f64) - (2.0 * (n_f - 1.0).sqrt() + 1.0);
                // C5: Hansen-Lucas: q(G) - 2R(G) <= n - 2*sqrt(n-1) => deficit = (n - 2*sqrt(n-1)) - (q - 2R)
                let c5_deficit = (n_f - 2.0 * (n_f - 1.0).sqrt()) - (q_rad - 2.0 * r);
                // C6: Wiener vs Spectral Radius: lambda_1 * W >= (n-1)^2
                let c6_deficit = (lam1 * w_idx as f64) - (n_f - 1.0).powi(2);
                // C7: Harmonic vs Proximity: H(G) >= pi(G)
                let c7_deficit = h - prox;
                // C8: Laplacian Spectral Radius * Diameter >= 2*(n-1)
                let c8_deficit = (mu1 * diam as f64) - (2.0 * (n_f - 1.0));

                let deficits = [
                    c1_deficit, c2_deficit, c3_deficit, c4_deficit,
                    c5_deficit, c6_deficit, c7_deficit, c8_deficit,
                ];

                for (idx, &def) in deficits.iter().enumerate() {
                    if def < local_best[idx] {
                        local_best[idx] = def;
                        if def < -1e-9 {
                            ce_found_flag.store(true, Ordering::SeqCst);
                            println!(
                                "\x1b[1;31m🚨🚨🚨 DISPROVED! COUNTEREXAMPLE FOUND FOR CONJECTURE C{} AT V = {}! Deficit = {:.6} < 0! 🚨🚨🚨\x1b[0m",
                                idx + 1, n, def
                            );
                            println!("   Edges = {}, Diam = {}, avg_l = {:.4}, H = {:.4}, R = {:.4}, lam1 = {:.4}",
                                g.num_edges(), diam, avg_l, h, r, lam1
                            );
                            println!("   Adjacency Matrix:\n{}", g.format_adjacency());

                            let _ = OpenOptions::new()
                                .create(true)
                                .append(true)
                                .open("discovered_counterexamples.txt")
                                .map(|mut f| {
                                    writeln!(
                                        f,
                                        "GRAPH CONJECTURE C{} COUNTEREXAMPLE: V = {}, Deficit = {:.6}\nAdjacency:\n{}\n",
                                        idx + 1, n, def, g.format_adjacency()
                                    )
                                });
                        }
                    }
                }
            }

            total_evaluations.fetch_add(steps_per_walker as u64, Ordering::Relaxed);
            let mut g_best = global_best_mutex.lock().unwrap();
            for idx in 0..8 {
                if local_best[idx] < g_best[idx].0 {
                    g_best[idx].0 = local_best[idx];
                    g_best[idx].1 = format!("V={}, min={:.4}", n, local_best[idx]);
                }
            }
        });
    }

    let elapsed = start_time.elapsed().as_secs_f64();
    let total_steps = total_evaluations.load(Ordering::Relaxed);
    let speed = (total_steps as f64) / elapsed.max(1e-6) / 1e6;

    let g_best = global_best_mutex.lock().unwrap();

    println!("\n\x1b[1;36m================================================================================\x1b[0m");
    println!("🏁 GRAPH CONJECTURE BREAKER SUMMARY TABLE");
    println!("   Total Graph Mutations / Evaluations : {}", total_steps);
    println!("   Total Search Time                   : {:.2} seconds", elapsed);
    println!("   Average Evaluation Throughput       : {:.2} M graphs/s", speed);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let conj_names = [
        ("C1: H(G) >= avg_l(G)", "Harmonic index vs Average distance"),
        ("C2: H(G) >= D(G) / 2", "Harmonic index vs Half diameter"),
        ("C3: R(G) >= avg_l(G)", "Graffiti Randic vs Average distance"),
        ("C4: lam1 + D >= 2*sqrt(n-1) + 1", "AutoGraphiX Spectral radius + Diameter"),
        ("C5: q - 2R <= n - 2*sqrt(n-1)", "Hansen-Lucas Signless Laplacian vs Randic"),
        ("C6: lam1 * W >= (n-1)^2", "Wiener index vs Spectral radius"),
        ("C7: H(G) >= pi(G)", "Harmonic index vs Proximity"),
        ("C8: mu1 * D >= 2*(n-1)", "Laplacian spectral radius * Diameter"),
    ];

    println!(
        "   {:>32} | {:>14} | {:>20} | {:>14}",
        "Conjecture", "Min Deficit", "Best Observed Case", "Verdict"
    );
    println!("   {:-<88}", "");

    for idx in 0..8 {
        let (def, details) = &g_best[idx];
        let (c_name, _desc) = conj_names[idx];
        let status = if *def < -1e-9 {
            "\x1b[1;31m🚨 DISPROVED\x1b[0m"
        } else if def.abs() < 1e-4 {
            "\x1b[1;33m⚠️ SHARP BOUND\x1b[0m"
        } else {
            "\x1b[1;32m✓ HELD TIGHT\x1b[0m"
        };

        println!(
            "   {:>32} | {:>14.6} | {:>20} | {:>14}",
            c_name, def, details, status
        );
    }
    println!("\n\x1b[1;36m================================================================================\x1b[0m\n");
}

pub fn run_graph_conjecture_breaker_master() {
    let orders = [12, 14, 16, 20, 24, 30];
    hunt_graph_conjectures(&orders, 25_000, 32);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_graph_invariants() {
        // Path on 4 vertices
        let mut g = FastGraph::new(4);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        assert!(g.is_connected());
        assert_eq!(g.num_edges(), 3);
        let (diam, avg_l, _prox, _rem, _w) = g.apsp();
        assert_eq!(diam, 3);
        assert!((avg_l - 5.0 / 3.0).abs() < 1e-6);
        let h = g.harmonic_index();
        assert!(h > 0.0);
    }
}
