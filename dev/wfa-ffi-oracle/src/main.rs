//! Oracle for WFA2-lib ends-free + nonzero-match-score divergence (issue #102).
//!
//! Asks whether the *C++* WFA2-lib returns the OPTIMAL ends-free alignment with
//! match score = +5 (DADA2 scoring), or reproduces the suboptimality we observe
//! in the pure-Rust `wfa2lib-rs`. See README.md. Nothing here ships — it is a
//! standalone FFI oracle used only to decide whether the upstream fix is real
//! and therefore worth porting into the pure-Rust `termination.rs`.
//!
//! Status: STUB. `wfa_endsfree_score` is not yet wired to the C++ library, so
//! the harness prints the known values and the expected verdict shape. Wire the
//! FFI binding (README "Wiring the FFI binding") and replace the stub.

/// DADA2 ends-free scoring used throughout dada2-rs (`src/nwalign.rs`).
const MATCH: i32 = 5;
const MISMATCH: i32 = -4;
const GAP: i32 = -8;

/// A reproducer: pattern/text, plus the optimal ends-free score that
/// `align_endsfree` (the true positional optimum) achieves.
struct Case {
    name: &'static str,
    pattern: &'static str,
    text: &'static str,
    optimal: i32,
    note: &'static str,
}

/// The two committed dada2-rs reproducers (see `nwalign.rs::wfa_endsfree_known_divergence`).
fn cases() -> Vec<Case> {
    // Mode 1: text is the pattern minus its leading base; the optimum credits a
    // free leading end-gap so all 51 shared columns score as matches (51*5=255).
    // The pure-Rust WFA mis-places it as a penalized internal gap (255-8=247).
    let p1 = "AACAGCGCAAACCAACTCGCTAGCTAGCAAAATCTTGTGTTTCTGCCTAGCG"; // 52 nt
    let t1 = "ACAGCGCAAACCAACTCGCTAGCTAGCAAAATCTTGTGTTTCTGCCTAGCG"; //  51 nt (p1[1..])
    vec![
        Case {
            name: "leading-gap (Mode 1)",
            pattern: p1,
            text: t1,
            optimal: 255,
            note: "optimum = free leading end-gap; pure-Rust WFA = 247 (internal gap)",
        },
        // The #102 minimal reproducer from the upstream issue thread.
        // Optimal ends-free: GCGG / -CG-  (the CG aligns to the interior GG-run
        // boundary with free ends). WFA prefers -GCGG / CG-- (suboptimal).
        Case {
            name: "#102 minimal (GCGG/CG)",
            pattern: "GCGG",
            text: "CG",
            // 2 matches (C,G) under the optimal placement, no mismatch; free
            // end-gaps unpenalized → 2*5 = 10 in DADA2 scoring.
            optimal: 2 * MATCH,
            note: "upstream #102 thread's 4-base reproducer",
        },
    ]
}

/// Ends-free alignment score of `pattern` vs `text` under (MATCH, MISMATCH, GAP)
/// as computed by the **C++ WFA2-lib** via FFI.
///
/// STUB — returns None until the FFI binding is wired (see README). When wired,
/// configure WFA2-lib with: gap-affine, match = +5 (or the Eizenga-equivalent),
/// mismatch = 4, gap-opening = 0, gap-extension = 8, and free ends on both ends
/// of both sequences; then convert its CIGAR to the DADA2 score and return it.
fn wfa_endsfree_score(_pattern: &str, _text: &str) -> Option<i32> {
    None
}

fn main() {
    println!(
        "WFA2-lib ends-free oracle (match={MATCH}, mismatch={MISMATCH}, gap={GAP})\n\
         purpose: does the C++ lib give the OPTIMAL ends-free score with match != 0? (#102)\n"
    );

    let mut wired = true;
    for c in cases() {
        print!(
            "[{}]  pattern={} ({} nt)  text={} ({} nt)  optimal={}",
            c.name,
            c.pattern,
            c.pattern.len(),
            c.text,
            c.text.len(),
            c.optimal
        );
        match wfa_endsfree_score(c.pattern, c.text) {
            Some(s) if s == c.optimal => println!("  ->  C++ WFA={s}  PASS (fix is real)"),
            Some(s) => println!("  ->  C++ WFA={s}  FAIL (still suboptimal upstream)"),
            None => {
                wired = false;
                println!("  ->  C++ WFA=<stub>  ({})", c.note);
            }
        }
    }

    if !wired {
        eprintln!(
            "\nSTUB: wfa_endsfree_score is not wired to the C++ WFA2-lib yet.\n\
             Add an FFI binding (README) and record the WFA2-lib commit tested.\n\
             PASS on every case => port the match-aware ends-free termination into\n\
             the pure-Rust wfa2lib-rs::termination.rs. Any FAIL => #102 unfixed upstream."
        );
        std::process::exit(2);
    }
}
