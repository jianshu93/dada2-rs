#!/usr/bin/env bash
# run_kmer_sweep.sh
# ---------------------------------------------------------------------------
# Sweep the k-mer pre-alignment screen size (--kmer-size) and measure its
# effect on (a) the inferred ASVs, (b) how aggressively the screen rejects
# pairs ("shrouding"), and (c) runtime / iterations. Built for issue #15.
#
# The pre-alignment screen rejects candidate pairs whose k-mer distance
# exceeds --kdist-cutoff, so only similar sequences reach the expensive
# Needleman-Wunsch alignment. On long reads (e.g. PacBio HiFi ~1.4 kb) the
# screen is effectively a no-op at the k=5 default (it shrouds ~nothing, so
# every pair is fully aligned); larger k makes it engage. The open question
# (issue #15) is whether a larger k changes the final ASVs vs. k=5 / k=6 —
# the hypothesis being that a larger k better isolates highly-similar
# sequences and may be an acceptable memory/speed tradeoff. This script lets
# you test that on any dataset.
#
# Usage:
#   bash run_kmer_sweep.sh <input.fastq[.gz] | derep.json> [outdir] [errfun] [band] [k-list]
#
# Examples:
#   # PacBio HiFi sample, default k list (5 6 8):
#   bash run_kmer_sweep.sh reads.trim.fastq.gz ./ksweep_out pacbio 32
#
#   # Illumina, custom k list:
#   bash run_kmer_sweep.sh F3D0_R1.fastq.gz ./ksweep_out loess 16 "5 6 7 8"
#
# Defaults: outdir=./kmer_sweep_out  errfun=pacbio  band=32  k-list="5 6 8"
#
# IMPORTANT for PacBio raw reads: reads must already be primer-trimmed AND
# uniformly oriented before this script. A sequence and its reverse complement
# share almost no k-mers, so mixed-orientation input defeats the screen at ANY
# k. If your reads are not oriented, run `remove-primers --orient` first.
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DADA2RS="${SCRIPT_DIR}/../target/release/dada2-rs"

INPUT="${1:?Usage: run_kmer_sweep.sh <input.fastq|derep.json> [outdir] [errfun] [band] [k-list]}"
OUTDIR="${2:-${SCRIPT_DIR}/kmer_sweep_out}"
ERRFUN="${3:-pacbio}"
BAND="${4:-32}"
KLIST="${5:-5 6 8}"

# Tunables (override via environment): held CONSTANT across all k so the only
# variable is the k-mer screen size.
MAX_CONSIST="${MAX_CONSIST:-10}"
KDIST_CUTOFF="${KDIST_CUTOFF:-0.42}"
THREADS="${THREADS:-1}"   # keep at 1 for clean, comparable wall times

if [[ ! -x "$DADA2RS" ]]; then
    echo "ERROR: binary not found/executable at $DADA2RS" >&2
    echo "       Run 'cargo build --release' first." >&2
    exit 1
fi

mkdir -p "$OUTDIR"
SAMPLE="$(basename "$INPUT")"
SAMPLE="${SAMPLE%.gz}"; SAMPLE="${SAMPLE%.fastq}"; SAMPLE="${SAMPLE%.json}"

# ---------------------------------------------------------------------------
# Step 0: get a derep JSON (skip if input is already one).
# ---------------------------------------------------------------------------
case "$INPUT" in
    *.json|*.json.gz)
        DEREP="$INPUT"
        echo "==> Using existing derep JSON: $DEREP"
        ;;
    *)
        DEREP="${OUTDIR}/${SAMPLE}_derep.json"
        echo "==> Dereplicating $INPUT"
        "$DADA2RS" derep "$INPUT" -o "$DEREP"
        echo "    -> $DEREP"
        ;;
esac

# ---------------------------------------------------------------------------
# Step 1: per-k learn-errors + dada. Capture iterations, shroud counts, ASVs.
# ---------------------------------------------------------------------------
# /usr/bin/time -l on macOS prints "<bytes> maximum resident set size"; GNU
# /usr/bin/time -v prints "Maximum resident set size (kbytes)". We grab peak
# RSS in a portable-ish way below.
TIME_BIN=/usr/bin/time

SUMMARY_CSV="${OUTDIR}/summary.csv"
echo "k,learn_iters,dada_aligns,dada_shrouded,shroud_pct,n_asv,wall_s_dada,maxrss_kb_dada" > "$SUMMARY_CSV"

for k in $KLIST; do
    echo ""
    echo "================  k = $k  ================"
    ERR_JSON="${OUTDIR}/${SAMPLE}_errors_k${k}.json"
    DADA_JSON="${OUTDIR}/${SAMPLE}_dada_k${k}.json"
    LEARN_LOG="${OUTDIR}/${SAMPLE}_learn_k${k}.log"
    DADA_LOG="${OUTDIR}/${SAMPLE}_dada_k${k}.log"
    DADA_TIME="${OUTDIR}/${SAMPLE}_dada_k${k}.time"

    echo "==> learn-errors (k=$k)"
    "$DADA2RS" errors-from-sample "$DEREP" \
        --errfun "$ERRFUN" --band "$BAND" \
        --kmer-size "$k" --kdist-cutoff "$KDIST_CUTOFF" \
        --max-consist "$MAX_CONSIST" --threads "$THREADS" \
        --verbose -o "$ERR_JSON" 2> "$LEARN_LOG"

    echo "==> dada (k=$k)"
    # Time the dada step; tolerate either BSD or GNU /usr/bin/time output.
    "$TIME_BIN" -l "$DADA2RS" dada "$DEREP" \
        --error-model "$ERR_JSON" \
        --band "$BAND" --kmer-size "$k" --kdist-cutoff "$KDIST_CUTOFF" \
        --threads "$THREADS" --verbose -o "$DADA_JSON" \
        > /dev/null 2> "$DADA_TIME" || {
            # Some systems lack `/usr/bin/time -l`; retry without timing.
            "$DADA2RS" dada "$DEREP" --error-model "$ERR_JSON" \
                --band "$BAND" --kmer-size "$k" --kdist-cutoff "$KDIST_CUTOFF" \
                --threads "$THREADS" --verbose -o "$DADA_JSON" 2> "$DADA_TIME"
        }
    # The dada --verbose progress (incl. the "ALIGN:" line) goes to stderr,
    # which /usr/bin/time also writes to; keep a copy for grepping.
    cp "$DADA_TIME" "$DADA_LOG"

    echo "    errors -> $ERR_JSON"
    echo "    dada   -> $DADA_JSON"
done

# ---------------------------------------------------------------------------
# Step 2: analyze — fill the CSV and compare ASV sets across k.
# All parsing is in python for robustness; ASV extraction is schema-tolerant.
# ---------------------------------------------------------------------------
python3 - "$OUTDIR" "$SAMPLE" "$SUMMARY_CSV" "$KLIST" <<'PY'
import json, os, re, sys

outdir, sample, csv_path, klist = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split()

def extract_asvs(path):
    """Return the set of inferred ASV sequences from a dada output JSON.
    Schema-tolerant: finds the largest list whose elements are dicts carrying
    a sequence-like string field."""
    try:
        d = json.load(open(path))
    except Exception as e:
        return None
    # Canonical dada output schema: {"asvs": [{"sequence","abundance","reads"}]}
    if isinstance(d, dict) and isinstance(d.get("asvs"), list):
        seqs = [el["sequence"].upper() for el in d["asvs"]
                if isinstance(el, dict) and isinstance(el.get("sequence"), str)]
        if seqs:
            return set(seqs)
    # Fallback heuristic for other/older schemas.
    best = None
    def walk(node):
        nonlocal best
        if isinstance(node, list):
            seqs = []
            for el in node:
                if isinstance(el, dict):
                    for key in ("sequence", "seq", "denoised", "asv"):
                        v = el.get(key)
                        if isinstance(v, str) and set(v.upper()) <= set("ACGTN-"):
                            seqs.append(v.upper())
                            break
            if seqs and (best is None or len(seqs) > len(best)):
                best = seqs
            for el in node:
                walk(el)
        elif isinstance(node, dict):
            # top-level "sequence":count maps are also common
            strkeys = [kk for kk in node
                       if isinstance(kk, str) and kk
                       and set(kk.upper()) <= set("ACGTN-") and len(kk) > 20]
            if len(strkeys) > 1 and (best is None or len(strkeys) > len(best)):
                best = [s.upper() for s in strkeys]
            for v in node.values():
                walk(v)
    walk(d)
    return set(best) if best else set()

def parse_iters(log):
    """Count self-consistency rounds. The learn-errors --verbose log emits one
    line per round like '[learn_errors] iter=N sample=1: K cluster(s)' and a
    final '[learn_errors] converged after N iteration(s)'. Prefer the explicit
    'converged after N' count; otherwise fall back to the max iter= seen
    (de-duped). Also tolerate older 'Selfconsist round N' / 'Iteration N' logs."""
    if not os.path.exists(log): return ""
    converged = None
    rounds = set()
    for ln in open(log, errors="replace"):
        m = re.search(r"converged after\s+(\d+)\s+iteration", ln)
        if m:
            converged = int(m.group(1))
            continue
        m = re.search(r"\biter=(\d+)", ln)
        if m:
            rounds.add(int(m.group(1)))
            continue
        m = re.search(r"Selfconsist round\s+(\d+)", ln) or re.search(r"Iteration\s+(\d+)", ln)
        if m:
            rounds.add(int(m.group(1)))
    if converged is not None:
        return converged
    return max(rounds) if rounds else 0

def parse_align(log):
    """Return (aligns, shrouded) from the dada 'ALIGN:' verbose line."""
    if not os.path.exists(log): return ("", "")
    a = s = ""
    for ln in open(log, errors="replace"):
        m = re.search(r"ALIGN:\s*(\d+)\s+aligns,\s*(\d+)\s+shrouded", ln)
        if m:
            a, s = m.group(1), m.group(2)
    return (a, s)

def parse_time(tf):
    """Return (wall_s, maxrss_kb) from BSD or GNU /usr/bin/time output."""
    if not os.path.exists(tf): return ("", "")
    wall = rss = ""
    txt = open(tf, errors="replace").read()
    # BSD: "        2.45 real ..."  ;  GNU: "Elapsed (wall clock) time ... m:ss"
    m = re.search(r"^\s*([\d.]+)\s+real", txt, re.M)
    if m: wall = m.group(1)
    else:
        m = re.search(r"wall clock.*?(\d+):([\d.]+)", txt)
        if m: wall = str(int(m.group(1))*60 + float(m.group(2)))
    # BSD: "<bytes> maximum resident set size" ; GNU: "Maximum resident set size (kbytes): N"
    m = re.search(r"(\d+)\s+maximum resident set size", txt)
    if m: rss = str(round(int(m.group(1))/1024))      # bytes -> kB
    else:
        m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)
        if m: rss = m.group(1)
    return (wall, rss)

asv_sets = {}
rows = []
for k in klist:
    learn_log = f"{outdir}/{sample}_learn_k{k}.log"
    dada_log  = f"{outdir}/{sample}_dada_k{k}.log"
    dada_time = f"{outdir}/{sample}_dada_k{k}.time"
    dada_json = f"{outdir}/{sample}_dada_k{k}.json"
    iters = parse_iters(learn_log)
    aligns, shrouded = parse_align(dada_log)
    wall, rss = parse_time(dada_time)
    asvs = extract_asvs(dada_json)
    asv_sets[k] = asvs
    n_asv = len(asvs) if asvs is not None else ""
    pct = ""
    try:
        tot = int(aligns) + int(shrouded)
        if tot: pct = f"{100*int(shrouded)/tot:.1f}"
    except ValueError:
        pass
    rows.append((k, iters, aligns, shrouded, pct, n_asv, wall, rss))

with open(csv_path, "w") as fh:
    fh.write("k,learn_iters,dada_aligns,dada_shrouded,shroud_pct,n_asv,wall_s_dada,maxrss_kb_dada\n")
    for r in rows:
        fh.write(",".join(str(x) for x in r) + "\n")

# Console summary
print("\n================  SUMMARY  ================")
hdr = f"{'k':>3} {'iters':>5} {'aligns':>9} {'shroud':>9} {'shroud%':>7} {'#ASV':>6} {'wall_s':>7} {'RSS_MB':>7}"
print(hdr); print("-"*len(hdr))
for k, iters, aligns, shrouded, pct, n_asv, wall, rss in rows:
    rss_mb = f"{int(rss)/1024:.0f}" if rss else ""
    print(f"{k:>3} {iters:>5} {aligns:>9} {shrouded:>9} {pct:>7} {n_asv:>6} {wall:>7} {rss_mb:>7}")

# ASV set comparison vs the smallest k in the list (usually k=5, the baseline)
ks = list(klist)
base = ks[0]
b = asv_sets.get(base) or set()
print(f"\nASV set comparison (baseline k={base}, {len(b)} ASVs):")
for k in ks[1:]:
    s = asv_sets.get(k) or set()
    only_base = len(b - s)
    only_k    = len(s - b)
    shared    = len(b & s)
    status = "IDENTICAL" if (only_base == 0 and only_k == 0) else "DIFFERS"
    print(f"  k={k}: shared={shared}  only_in_k{base}={only_base}  only_in_k{k}={only_k}   [{status}]")

print(f"\nWrote {csv_path}")
print("Per-k outputs: errors_k*.json, dada_k*.json, *_learn_k*.log, *_dada_k*.log")
PY

echo ""
echo "Done. CSV: ${SUMMARY_CSV}"
