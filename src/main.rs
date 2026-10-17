//! Counter-Hunter CLI Main Entry Point
//! Ultra-fast bit-parallel counterexample hunter for Zhi-Wei Sun & Erdős conjectures.

use std::sync::atomic::Ordering;
use std::time::Instant;

use clap::Parser;
use counter_hunter::conjectures::TargetConjecture;
use counter_hunter::hunter::{HuntConfig, Hunter};

#[derive(Parser, Debug)]
#[command(
    name = "counter-hunter",
    author = "Advanced Agentic Coding Team",
    version = "0.1.0",
    about = "Bit-Parallel Counterexample Hunter for Zhi-Wei Sun & Erdős Conjectures"
)]
struct Cli {
    /// Conjecture target to hunt:
    /// 'a-twin', 'a-sophie', 'b-cubic1', 'b-cubic2', 'c-binom', 'd-fib', 'd-square'
    #[arg(short = 't', long = "target", default_value = "b-cubic1")]
    target: String,

    /// Start of search range (e.g. 100000000)
    #[arg(short = 's', long = "start", default_value = "100000000")]
    start: u64,

    /// End of search range (e.g. 1000000000)
    #[arg(short = 'e', long = "end", default_value = "1000000000")]
    end: u64,

    /// Segment size in number of integers (default: 33554432 = 2^25 integers, 4 MB bitset)
    #[arg(long = "segment-size", default_value = "33554432")]
    segment_size: u64,

    /// Sieve margin for fast Shift-OR (default: 2097152 = 2^21 integers)
    #[arg(long = "margin", default_value = "2097152")]
    margin: u64,

    /// Maximum number of shifts to apply in bitwise sieve phase (default: 150)
    #[arg(long = "max-shifts", default_value = "150")]
    max_shifts: usize,

    /// Maximum counterexamples to fully factorize and save to proof file (default: 20)
    #[arg(long = "max-counterexamples", default_value = "20")]
    max_counterexamples: usize,

    /// Parity filter: 'default', 'odd', 'even', 'all'
    #[arg(long = "parity", default_value = "default")]
    parity: String,

    /// Live telemetry reporting interval in seconds (default: 5)
    #[arg(long = "telemetry-interval", default_value = "5")]
    telemetry_interval: u64,

    /// Number of worker threads (default: 0 = Rayon auto-detect all cores)
    #[arg(short = 'j', long = "threads", default_value = "0")]
    threads: usize,

    /// Output file for confirmed counterexamples with mathematical proof certificates
    #[arg(short = 'o', long = "output", default_value = "discovered_counterexamples.txt")]
    output: String,

    /// Enable verbose candidate inspection logs
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    /// Benchmark all 4 conjecture targets sequentially
    #[arg(long = "benchmark-all")]
    benchmark_all: bool,
}

fn parse_target(s: &str) -> TargetConjecture {
    match s.to_lowercase().as_str() {
        "a" | "a-twin" | "twin" => TargetConjecture::SunTwinPrime,
        "a-sophie" | "sophie" => TargetConjecture::SunSophieGermain,
        "b" | "b1" | "b-cubic1" | "cubic1" => TargetConjecture::SunCubic2X3Y4,
        "b2" | "b-cubic2" | "cubic2" => TargetConjecture::SunCubicX32Y3,
        "c" | "c-binom" | "binom" => TargetConjecture::SunCentralBinomial,
        "d" | "d1" | "d-fib" | "fib" => TargetConjecture::SunFibonacciLucas,
        "d2" | "d-square" | "romanov" => TargetConjecture::SunPrimeSquareRomanov,
        other => {
            eprintln!("Unknown target '{}', defaulting to Target B1 (n = p + 2x^3 + y^4)", other);
            TargetConjecture::SunCubic2X3Y4
        }
    }
}

fn format_number(n: u64) -> String {
    if n >= 1_000_000_000_000 {
        format!("{:.3} T ({})", n as f64 / 1e12, n)
    } else if n >= 1_000_000_000 {
        format!("{:.3} B ({})", n as f64 / 1e9, n)
    } else if n >= 1_000_000 {
        format!("{:.3} M ({})", n as f64 / 1e6, n)
    } else {
        n.to_string()
    }
}

fn parse_parity(s: &str) -> counter_hunter::hunter::ParityFilter {
    match s.to_lowercase().as_str() {
        "odd" => counter_hunter::hunter::ParityFilter::OddOnly,
        "even" => counter_hunter::hunter::ParityFilter::EvenOnly,
        "all" => counter_hunter::hunter::ParityFilter::All,
        _ => counter_hunter::hunter::ParityFilter::Default,
    }
}

fn run_hunt_for_target(target: TargetConjecture, cli: &Cli) {
    let parity = parse_parity(&cli.parity);
    let config = HuntConfig {
        target,
        start: cli.start,
        end: cli.end,
        segment_size: cli.segment_size,
        max_sieve_margin: cli.margin,
        max_shifts_to_apply: cli.max_shifts,
        max_counterexamples: cli.max_counterexamples,
        parity,
        num_threads: cli.threads,
        telemetry_interval_secs: cli.telemetry_interval,
        verbose: cli.verbose,
        output_file: cli.output.clone(),
    };

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🚀 COUNTER-HUNTER — High-Throughput Bit-Parallel Number Theory Engine\x1b[0m");
    println!("   Zhi-Wei Sun & Paul Erdős Open Conjectures Hunter (Lock-Free Multi-Core)");
    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1m🎯 Target Conjecture :\x1b[0m \x1b[1;33m{}\x1b[0m", target.name());
    println!("\x1b[1m🔍 Search Space     :\x1b[0m {} .. {}", format_number(cli.start), format_number(cli.end));
    println!("\x1b[1m⚡ Segment Size     :\x1b[0m {} numbers ({:.1} MB bitset)", cli.segment_size, cli.segment_size as f64 / 8.0 / 1024.0 / 1024.0);
    println!("\x1b[1m🛡️ Sieve Margin     :\x1b[0m {} numbers ({:.1} MB)", cli.margin, cli.margin as f64 / 8.0 / 1024.0 / 1024.0);
    println!("\x1b[1m🧵 Threading Engine :\x1b[0m {} workers (Lock-free atomic work stealing)", if cli.threads == 0 { "All available cores".to_string() } else { cli.threads.to_string() });
    println!("\x1b[1m📡 Telemetry Pulse  :\x1b[0m Every {}s (heartbeat monitor)", cli.telemetry_interval);
    println!("\x1b[1m📁 Proof Log File   :\x1b[0m {}", cli.output);
    println!("\x1b[1;36m================================================================================\x1b[0m");

    let start_time = Instant::now();
    let hunter = Hunter::new(config);
    hunter.hunt();

    let total_elapsed = start_time.elapsed().as_secs_f64();
    let total_scanned = hunter.stats.numbers_scanned.load(Ordering::Relaxed);
    let total_candidates = hunter.stats.candidates_tested.load(Ordering::Relaxed);
    let total_counterexamples = hunter.stats.counterexamples_found.load(Ordering::Relaxed);
    let total_l1_blocks = hunter.stats.l1_blocks_processed.load(Ordering::Relaxed);
    let overall_speed = if total_elapsed > 0.0 {
        (total_scanned as f64) / (total_elapsed * 1e9)
    } else {
        0.0
    };

    println!("\x1b[1;36m================================================================================\x1b[0m");
    println!("\x1b[1;32m🏁 HUNT COMPLETED for {}\x1b[0m", target.short_code());
    println!("   Total Integers Scanned   : {}", total_scanned);
    println!("   Average Throughput       : \x1b[1;32m{:.3} Giga-integers / second\x1b[0m ({:.1} Mnums/s)", overall_speed, overall_speed * 1000.0);
    println!("   Total Time Elapsed       : {:.2} seconds", total_elapsed);
    println!("   L1 Cache 32 KiB Blocks   : {}", total_l1_blocks);
    println!("   Candidates Filtered      : {}", total_candidates);
    println!("   Counterexamples Verified : \x1b[1;31m{}\x1b[0m", total_counterexamples);
    println!("\x1b[1;36m================================================================================\x1b[0m\n");
}

fn main() {
    let cli = Cli::parse();

    if cli.benchmark_all {
        println!("\x1b[1;33m>>> EXECUTING SEQUENTIAL SCAN ACROSS ALL 4 CONJECTURE TARGETS <<<\x1b[0m\n");
        let targets = [
            TargetConjecture::SunTwinPrime,
            TargetConjecture::SunSophieGermain,
            TargetConjecture::SunCubic2X3Y4,
            TargetConjecture::SunCubicX32Y3,
            TargetConjecture::SunCentralBinomial,
            TargetConjecture::SunFibonacciLucas,
            TargetConjecture::SunPrimeSquareRomanov,
        ];
        for t in targets {
            run_hunt_for_target(t, &cli);
        }
    } else {
        let target = parse_target(&cli.target);
        run_hunt_for_target(target, &cli);
    }
}
