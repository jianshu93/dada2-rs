//! k-mer-distance vs true-divergence calibration (hidden `kdist-calibrate` subcommand).
//!
//! DADA2's k-mer screen skips alignment for pairs with k-mer distance above
//! `KDIST_CUTOFF = 0.42` (nominally ~10% nucleotide divergence, calibrated on
//! Illumina 16S). The k-mer distance traces to ESPRIT (Sun et al. 2009, whose
//! reference implementation is gone); the 0.42/10% calibration to DADA2
//! (Callahan et al. 2016). This re-derives the relationship empirically on real
//! data: for sampled unique-sequence pairs it emits the k-mer distance
//! (`kmer_dist8`, our port of the ESPRIT metric) alongside the true UNBANDED
//! `align_endsfree` divergence, so the constant can be checked per dataset /
//! platform / k / pooling regime.
//!
//! Alignment is unbanded by default (`--band -1`): the curve must measure the
//! true divergence of distant pairs, which a band would truncate. It is the
//! expensive part, so it is parallelised across `--threads`.
//!
//! POOLING: with multiple derep JSONs, `--per-sample` computes pairs WITHIN each
//! sample (the per-sample / independent regime); the default pools all uniques
//! into one set and computes pairs across the union (the full-pool regime).
//! Pseudo's screen population is per-sample (priors change the partition, not
//! which pairs are screened), so model it with `--per-sample`.

use crate::kmers::{assign_kmer8, kmer_dist8};
use crate::misc::WithPath as _;
use crate::nwalign::{AlignBuffers, align_endsfree_with_buf};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const GAP: u8 = b'-';

/// Parameters for [`run`] (mirrors the CLI flags).
pub struct Params {
    pub k: usize,
    pub cutoff: f64,
    pub leak_pct: f64,
    pub band: i32,
    pub max_pairs: usize,
    pub max_uniques: usize,
    pub per_sample: bool,
    pub threads: usize,
    pub seed: u64,
    pub output: Option<PathBuf>,
    pub verbose: bool,
}

fn encode(seq: &str) -> Vec<u8> {
    seq.bytes()
        .map(|b| match b {
            b'A' | b'a' => 1,
            b'C' | b'c' => 2,
            b'G' | b'g' => 3,
            b'T' | b't' => 4,
            _ => 5, // N etc. — never a valid k-mer, never matches
        })
        .collect()
}

/// Internal edit divergence of an ends-free alignment: trim terminal gap
/// overhang (length difference, not divergence), then count substitution and
/// indel columns in the aligned core. Returns (edits, core_len).
fn aln_divergence(a: &[u8], b: &[u8]) -> (usize, usize) {
    let n = a.len();
    let mut lo = 0;
    while lo < n && (a[lo] == GAP || b[lo] == GAP) {
        lo += 1;
    }
    let mut hi = n;
    while hi > lo && (a[hi - 1] == GAP || b[hi - 1] == GAP) {
        hi -= 1;
    }
    let mut edits = 0;
    for k in lo..hi {
        if a[k] == GAP || b[k] == GAP || a[k] != b[k] {
            edits += 1;
        }
    }
    (edits, hi - lo)
}

/// One sample's uniques: encoded sequences, abundances, and k-mer screens.
struct Sample {
    name: String,
    enc: Vec<Vec<u8>>,
    counts: Vec<u64>,
    kmers: Vec<Vec<u8>>,
}

/// Load a derep JSON (`uniques[].sequence` + `count`), gzip-transparent. Only
/// derep JSONs are accepted (the screen operates on per-sample uniques).
fn load_derep(path: &Path, k: usize, max_uniques: usize, seed: u64) -> io::Result<Sample> {
    let f = File::open(path).with_path(path)?;
    let mut txt = String::new();
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        MultiGzDecoder::new(f).read_to_string(&mut txt)?;
    } else {
        io::BufReader::new(f).read_to_string(&mut txt)?;
    }
    let v: serde_json::Value =
        serde_json::from_str(&txt).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let uniques = v.get("uniques").and_then(|u| u.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: not a derep JSON (no `uniques`)", path.display()),
        )
    })?;
    let name = v
        .get("sample")
        .and_then(|s| s.as_str())
        .map(str::to_string)
        .unwrap_or_else(|| {
            path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("?")
                .to_string()
        });
    let mut seqs: Vec<(Vec<u8>, u64)> = uniques
        .iter()
        .filter_map(|e| {
            let s = e.get("sequence")?.as_str()?;
            let c = e.get("count").and_then(|c| c.as_u64()).unwrap_or(1);
            Some((encode(s), c))
        })
        .collect();
    // Optional per-sample subsample of uniques to bound the O(n^2) pair count.
    // Random (not abundance-top) keeps the divergence distribution unbiased.
    if max_uniques > 0 && seqs.len() > max_uniques {
        let mut st = seed ^ (seqs.len() as u64).wrapping_mul(0x9E37_79B9);
        // partial Fisher–Yates: move `max_uniques` random picks to the front.
        for i in 0..max_uniques {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            let j = i + ((st >> 33) as usize) % (seqs.len() - i);
            seqs.swap(i, j);
        }
        seqs.truncate(max_uniques);
    }
    let (enc, counts): (Vec<_>, Vec<_>) = seqs.into_iter().unzip();
    let kmers = enc.iter().map(|e| assign_kmer8(e, k)).collect();
    Ok(Sample {
        name,
        enc,
        counts,
        kmers,
    })
}

/// Build the (i, j) pair list for a population of `n` uniques: enumerate all if
/// `n*(n-1)/2 <= max_pairs`, else draw `max_pairs` random pairs (with possible
/// repeats — fine for a calibration scatter).
fn pairs_for(n: usize, max_pairs: usize, seed: u64) -> Vec<(usize, usize)> {
    let total = n.saturating_mul(n.saturating_sub(1)) / 2;
    if n < 2 {
        return Vec::new();
    }
    if total <= max_pairs {
        let mut v = Vec::with_capacity(total);
        for i in 0..n {
            for j in (i + 1)..n {
                v.push((i, j));
            }
        }
        return v;
    }
    let mut st = seed;
    let mut rnd = |m: usize| {
        st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((st >> 33) as usize) % m
    };
    (0..max_pairs)
        .map(|_| {
            let i = rnd(n);
            let mut j = rnd(n);
            if i == j {
                j = (j + 1) % n;
            }
            (i.min(j), i.max(j))
        })
        .collect()
}

pub fn run(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    if inputs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no input derep JSON(s) given",
        ));
    }
    let loaded: Vec<Sample> = inputs
        .iter()
        .map(|path| load_derep(path, p.k, p.max_uniques, p.seed))
        .collect::<io::Result<_>>()?;

    // Form populations: one per sample (per-sample) or a single merged pool.
    let pops: Vec<Sample> = if p.per_sample {
        loaded
    } else {
        let mut pool = Sample {
            name: "pool".into(),
            enc: Vec::new(),
            counts: Vec::new(),
            kmers: Vec::new(),
        };
        for s in loaded {
            pool.enc.extend(s.enc);
            pool.counts.extend(s.counts);
            pool.kmers.extend(s.kmers);
        }
        vec![pool]
    };

    let mut w: Box<dyn Write> = match &p.output {
        Some(path) => Box::new(BufWriter::new(File::create(path).with_path(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    writeln!(
        w,
        "sample,kdist,edits,core_len,pct_div,screened_in,ab_i,ab_j"
    )?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads)
        .build()
        .map_err(io::Error::other)?;

    let (mut tot, mut scr, mut leak) = (0u64, 0u64, 0u64);
    for (pi, s) in pops.iter().enumerate() {
        let n = s.enc.len();
        let seed = p.seed.wrapping_add(pi as u64).wrapping_mul(0x100000001B3);
        let pairs = pairs_for(n, p.max_pairs, seed);
        if p.verbose {
            let total = n.saturating_mul(n.saturating_sub(1)) / 2;
            eprintln!(
                "[kdist] {} : {n} uniques, {total} pairs -> {} computed (k={}, band={}, {} threads)",
                s.name,
                pairs.len(),
                p.k,
                p.band,
                p.threads,
            );
        }
        // Parallel per-pair: kdist (cheap) + unbanded NW (the cost). Per-thread
        // AlignBuffers reuse via map_init.
        let rows: Vec<(f64, usize, usize, f64, bool)> = pool.install(|| {
            pairs
                .par_iter()
                .map_init(AlignBuffers::new, |buf, &(i, j)| {
                    let kd = kmer_dist8(
                        &s.kmers[i],
                        s.enc[i].len(),
                        &s.kmers[j],
                        s.enc[j].len(),
                        p.k,
                    );
                    align_endsfree_with_buf(&s.enc[i], &s.enc[j], 5, -4, -8, p.band, buf);
                    let (edits, core) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (kd, edits, core, pct, kd < p.cutoff)
                })
                .collect()
        });
        for (idx, &(kd, edits, core, pct, sin)) in rows.iter().enumerate() {
            let (i, j) = pairs[idx];
            tot += 1;
            if sin {
                scr += 1;
                if pct > p.leak_pct {
                    leak += 1;
                }
            }
            writeln!(
                w,
                "{},{kd:.4},{edits},{core},{pct:.3},{},{},{}",
                s.name, sin as u8, s.counts[i], s.counts[j]
            )?;
        }
    }
    w.flush()?;
    if p.verbose && tot > 0 {
        eprintln!(
            "[kdist] {tot} pairs: screened-in (kdist<{}) {scr} ({:.1}%); of those {leak} are \
             >{}% divergent (leakage)",
            p.cutoff,
            100.0 * scr as f64 / tot as f64,
            p.leak_pct,
        );
    }
    Ok(())
}
