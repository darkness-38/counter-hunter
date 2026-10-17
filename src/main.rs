//! Counter-Hunter CLI Main Entry Point
//! Ultra-fast CRT Covering Sniper Engine for Zhi-Wei Sun & Erdős Conjectures.

use std::sync::atomic::Ordering;
use std::time::Instant;

use clap::Parser;
use counter_hunter::conjectures::TargetConjecture;
use counter_hunter::hunter::{HuntConfig, Hunter};

#[derive(Parser, Debug)]
#[command(
    name = "counter-hunter",
    author = "Advanced Agentic Coding Team",
    version = "0.2.0",
    about = "CRT Covering Sniper Engine for Zhi-Wei Sun & Erdős Conjectures"
)]
struct Cli {
    /// Conjecture target to hunt:
    /// 'power-2a-3b' (default), 'b-cubic1' (cubic-2x3-y4), 'binom-2a-prime', etc.
    #[arg(short = 't', long = "target", default_value = "power-2a-3b")]
    target: String,

    /// Start of search range (default: 100000000 = 10^8)
    #[arg(short = 's', long = "start", default_value = "100000000")]
    start: u64,

    /// End of search range (default: 10000000000 = 10^10, or up to 10^16)
    #[arg(short = 'e', long = "end", default_value = "10000000000")]
    end: u64,

    /// Maximum master CRT modulus M_master (default: 500000000000 = 5*10^11)
    #[arg(long = "max-modulus", default_value = "500000000000")]
    max_modulus: u64,

    /// Maximum counterexamples to fully factorize before halting (default: 1)
    #[arg(long = "max-counterexamples", default_value = "1")]
    max_counterexamples: usize,

    /// Live telemetry reporting interval in seconds (default: 10)
    #[arg(long = "telemetry-interval", default_value = "10")]
    telemetry_interval: u64,

    /// Number of worker threads (default: 0 = Rayon auto-detect all cores)
    #[arg(short = 'j', long = "threads", default_value = "0")]
    threads: usize,

    /// Output file for confirmed counterexamples with mathematical proof certificates
    #[arg(short = 'o', long = "output", default_value = "discovered_counterexamples.txt")]
    output: String,

    /// Skip mandatory Pre-Flight validator
    #[arg(long = "no-preflight")]
    no_preflight: bool,

    /// Enable verbose candidate inspection logs
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    /// Run Zhi-Wei Sun's core conjectures sequentially
    #[arg(long = "run-all")]
    run_all: bool,
}

fn parse_target(s: &str) -> TargetConjecture {
    match s.to_lowercase().replace('_', "-").as_str() {
        // Zhi-Wei Sun Primary Conjectures (Ordinary Primes)
        "power" | "power-2a-3b" | "sun-power" | "a-power" => TargetConjecture::SunPower2a3b,
        "b" | "b1" | "b-cubic1" | "cubic1" | "cubic-2x3-y4" | "cubic" => TargetConjecture::SunCubic2X3Y4,

        // Additional CRT Sniper Targets
        "1" | "binom-2a-prime" | "binom-prime" => TargetConjecture::Binom2aPrime,
        "2" | "binom-2a-twin" | "binom-twin" => TargetConjecture::Binom2aTwin,
        "3" | "twin-2a-3b" | "twin-3b" => TargetConjecture::Twin2a3b,
        "4" | "sophie-2a-3b" | "sophie-3b" => TargetConjecture::Sophie2a3b,
        "5" | "erdos-2a-2b-twin" | "erdos-twin" | "erdos" => TargetConjecture::Erdos2a2bTwin,

        // Legacy / Backwards Compatibility Aliases
        "a" | "a-twin" | "twin" => TargetConjecture::SunTwinPrime,
        "a-sophie" | "sophie" => TargetConjecture::SunSophieGermain,
        "b2" | "b-cubic2" | "cubic2" => TargetConjecture::SunCubicX32Y3,
        "c" | "c-binom" => TargetConjecture::SunCentralBinomial,
        "d" | "d1" | "d-fib" | "fib" => TargetConjecture::SunFibonacciLucas,
        "d2" | "d-square" | "romanov" => TargetConjecture::SunPrimeSquareRomanov,

        other => {
            eprintln!("Unknown target '{}', defaulting to power-2a-3b", other);
            TargetConjecture::SunPower2a3b
        }
    }
}

fn format_number(n: u64) -> String {
    if n >= 1_000_000_000_000_000 {
        format!("{:.3} P ({})", n as f64 / 1e15, n)
    } else if n >= 1_000_000_000_000 {
        format!("{:.3} T ({})", n as f64 / 1e12, n)
    } else if n >= 1_000_000_000 {
        format!("{:.3} B ({})", n as f64 / 1e9, n)
    } else if n >= 1_000_000 {
        format!("{:.3} M ({})", n as f64 / 1e6, n)
    } else {
        n.to_string()
    }
}

fn run_hunt_for_target(target: TargetConjecture, cli: &Cli) -> bool {
    let config = HuntConfig {
        target,
        start: cli.start,
        end: cli.end,
        max_master_modulus: cli.max_modulus,
        max_counterexamples: cli.max_counterexamples,
        num_threads: cli.threads,
        telemetry_interval_secs: cli.telemetry_interval,
        output_file: cli.output.clone(),
        run_preflight: !cli.no_preflight,
        verbose: cli.verbose,
    };

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🎯 CRT COVERING SNIPER ENGINE — Number Theory Counterexample Hunter\x1b[0m");
    println!("   Multiplicative Order Covering & Chinese Remainder Progression Leaper");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1m🎯 Target Conjecture :\x1b[0m \x1b[1;33m{}\x1b[0m", target.name());
    println!("\x1b[1m🔍 Search Space     :\x1b[0m {} .. {}", format_number(cli.start), format_number(cli.end));
    println!("\x1b[1m⚡ Modulus Cap      :\x1b[0m M_master ≤ {}", format_number(cli.max_modulus));
    println!("\x1b[1m🧵 Threading Engine :\x1b[0m {} workers (Lock-free work stealing)", if cli.threads == 0 { "All available cores".to_string() } else { cli.threads.to_string() });
    println!("\x1b[1m📡 Telemetry Pulse  :\x1b[0m Every {}s", cli.telemetry_interval);
    println!("\x1b[1m📁 Proof Log File   :\x1b[0m {}", cli.output);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    let start_time = Instant::now();
    let hunter = Hunter::new(config);
    let found_ce = hunter.hunt();

    let total_elapsed = start_time.elapsed().as_secs_f64();
    let total_k = hunter.stats.k_steps_evaluated.load(Ordering::Relaxed);
    let total_candidates = hunter.stats.candidates_tested.load(Ordering::Relaxed);
    let total_counterexamples = hunter.stats.counterexamples_found.load(Ordering::Relaxed);

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🏁 HUNT COMPLETED for {}\x1b[0m", target.short_code());
    println!("   Total k Steps Evaluated  : {}", total_k);
    println!("   Candidates Tested        : {}", total_candidates);
    println!("   Counterexamples Verified : \x1b[1;31m{}\x1b[0m", total_counterexamples);
    println!("   Total Time Elapsed       : {:.2} seconds", total_elapsed);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");

    found_ce
}

fn main() {
    let cli = Cli::parse();

    if cli.run_all {
        println!("\x1b[1;33m>>> EXECUTING SEQUENTIAL SCAN ACROSS ZHI-WEI SUN CORE CONJECTURES <<<\x1b[0m\n");
        let targets = [
            TargetConjecture::SunPower2a3b,
            TargetConjecture::SunCubic2X3Y4,
        ];
        for t in targets {
            let found = run_hunt_for_target(t, &cli);
            if found {
                println!("\x1b[1;31m🚨 HALTING RUN-ALL: Counterexample detected for {}!\x1b[0m\n", t.short_code());
                break;
            }
        }
    } else {
        let target = parse_target(&cli.target);
        run_hunt_for_target(target, &cli);
    }
}
