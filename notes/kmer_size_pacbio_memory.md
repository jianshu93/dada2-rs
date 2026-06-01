# K-mer size: PacBio HiFi screen behavior + memory (issue #15)

Follow-up to `notes/benchmark_kmer_size.md` (Illumina-only). Resolves the two
open questions in issue #15 for long reads: (1) does a larger `--kmer-size`
mis-handle the pre-alignment screen on ~1.4 kb PacBio HiFi reads, and (2) what
is the memory cost of large k.

- Date: 2026-05-31
- Binary: `target/release/dada2-rs` @ 36e0743
- Host: Apple Silicon, single-threaded runs (`--threads 1`) for clean compute timing
- All `errors-from-sample` runs: `--max-consist 20 --verbose`.
  "converged" = iteration count < the cap of 20.

> **Note on provenance.** The raw run artifacts lived in `/tmp/kmer_exp/` and were
> purged when `/tmp` cleared. The numbers below were captured/verified during the
> run (the `summary.csv` and shroud `ALIGN:` lines were re-read directly from disk
> before purge). The `err_out` max-diff column comes from the run logs and could
> not be independently re-derived after the purge — re-run to reproduce exactly.

## Datasets & preprocessing

| Dataset | Reads | Uniques | Read len | errfun | band |
|---------|-------|---------|----------|--------|------|
| MiSeq F3D0 (R1) | 7228 | 2014 | ~240 bp | loess | 16 |
| PacBio SRR28724909 | 32132 | 10109 | ~1450 bp | pacbio | 32 |

```
# MiSeq (raw -> filtered)
filter-and-trim --trunc-len 240 --max-n 0 --max-ee 2 --trunc-q 2 --compress

# PacBio
filter-and-trim --min-len 1000 --max-len 1600 --max-ee 2 --max-n 0 --compress
```

**Orientation check (important methodology note).** PacBio HiFi reads are *often*
mixed-orientation with primers attached, which would break the k-mer screen at
any k (a sequence and its reverse complement share almost no k-mers). We tested
empirically: of the first 2000 reads of `SRR28724909.trim.fastq.gz`, **0** carried
the 27F primer motif at the 5′ end and **0** carried its reverse-complement at
the 3′ end; `remove-primers --orient` produced 0 output. These reads are **already
primer-trimmed and uniformly oriented** (the `.trim` in the filename). So no
orientation step was needed, and the screen is being tested on legitimately
same-strand sequences. *For raw HiFi from other sources, orient first.*

**Scale.** The k-sweep used the **full 10109 PacBio uniques**, held constant
across all k. MiSeq used its full 2014 uniques.

## Results

### MiSeq F3D0 (Illumina, ~240 bp; loess, band 16)

| k | iters | wall (s) | peak RSS | err_out max\|Δ\| vs k5 |
|---|-------|----------|----------|------------------------|
| 5 | 5 | 5.21 | 43.7 MiB | — |
| 6 | 4 | 2.24 | 61.3 MiB | 1.25e-3 |
| 7 | 4 | 1.57 | 132.3 MiB | 2.06e-3 |
| 8 | 4 | 2.56 | 416.9 MiB | 4.35e-3 |

### PacBio SRR28724909 (HiFi, ~1450 bp; pacbio, band 32)

| k | iters | wall (s) | peak RSS | err_out max\|Δ\| vs k5 |
|---|-------|----------|----------|------------------------|
| 5 | 4 | 322.61 | 942.6 MiB | — |
| 6 | 4 | 78.44 | 1031.8 MiB | 4.91e-7 |
| 7 | 4 | 66.04 | 1389.9 MiB | 6.07e-7 |
| 8 | 4 | 67.29 | 2198.0 MiB | 6.65e-5 |

All runs converged. Learned error rates are **essentially identical across k**
(PacBio agrees to ~1e-7; MiSeq to a few 1e-3 at the low-Q edges).

> **Terminology:** "shrouded" = a candidate pair rejected by the **k-mer screen**
> before any Needleman–Wunsch alignment is attempted (the `sub_new` k-mer distance
> exceeds the cutoff, so it returns no alignment). The mechanism is what the
> literature calls *k-mer screening*; "shroud"/`nshroud` is the outcome counter,
> a name inherited verbatim from the upstream DADA2 C++ source (`dada.h`,
> `cluster.cpp`, `Rmain.cpp`'s `"%i aligns, %i shrouded"`). `shrouded` ⊆ `aligns`.

### Direct screen behavior (PacBio, `dada` iteration 1, kdist_cutoff = 0.42)

Same input, same band 32, only `--kmer-size` and the (own) learned model differ.
Total comparisons in the first cluster pass:

| k | aligns | shrouded | shroud % | final ASVs |
|---|--------|----------|----------|------------|
| 5 | 548600 | 0 | 0.0 % | 87 |
| 8 | 548707 | 432806 | 78.9 % | 87 |

## The key finding (corrects the naive intuition)

On **long reads the k=5 screen is a no-op**: it shrouds **0 %** of pairs, so every
pair goes to full O(L²) alignment — which is exactly why k=5 is the *slowest*
config (322 s vs ~70 s). Larger k makes the screen actually engage.

Why: `kdist = 1 − sharedKmers / (L − k + 1)`. For L ≈ 1450:
- **k=5:** only 4⁵ = 1024 distinct k-mers but ~1446 k-mer positions, so frequency
  vectors are dense and even unrelated long reads share most 5-mers → kdist stays
  far below 0.42 → nothing is screened. The screen carries no information at this
  k on long reads.
- **k=8:** 4⁸ = 65536 ≫ 1446 positions, so vectors are sparse; genuine HiFi error
  differences across 1.4 kb drop enough exact 8-mers to push many pairs past 0.42
  → 78.9 % shrouded.

Crucially, **the 78.9 % shrouding at k=8 was benign here**: final ASV count (87)
and learned error rates were identical to k=5. The screen rejected pairs that
the alignment would not have clustered anyway — it is *working*, not
*over-rejecting harmfully*. The aggressiveness is nonetheless a latent risk on
more divergent or lower-quality data, which argues for a moderate rather than
maximal k.

## Memory cost (the binding constraint)

K-mer storage per Raw = `4^k × 3 bytes` (u16 `kmer` + u8 `kmer8`; `kord` is
k-independent). For the **full PacBio sample** (10109 uniques):

| k | per-Raw | k-mer vectors (full sample) |
|---|---------|-----------------------------|
| 5 | 3.0 KB | 29.6 MiB |
| 6 | 12.0 KB | 118.5 MiB |
| 7 | 48.0 KB | 473.9 MiB |
| 8 | 192.0 KB | 1895.4 MiB |

Measured peak RSS tracks this 4×-per-step growth (PacBio: 943 → 1032 → 1390 →
2198 MiB). Pooled multi-sample HiFi (hundreds of thousands of uniques) would be
infeasible at k ≥ 7.

## Conclusions

1. **No harm from larger k on long reads.** Across k=5–8, PacBio converged in 4
   iterations, produced 87 ASVs, and learned the same error model (~1e-7
   agreement). The #15 worry about false-rejecting valid pairs did not
   materialize: even at 78.9 % shrouding the result was unchanged.
2. **k=5 wastes the screen on long reads** — 0 % shrouded means every pair is
   aligned, making k=5 the slowest PacBio config (~4× slower than k=6).
3. **k=6 is the PacBio sweet spot**: the screen engages (≈4× speedup over k=5),
   results are identical, and memory is only ~120 MiB on a full sample.
4. **k=7/k=8 add memory (4× per step) with no benefit** over k=6 here, and their
   aggressive shrouding is a latent false-negative risk on harder data.

## Recommendation

- Keep **k=5 as the global default** (correct and lowest-memory for short reads,
  where the screen does discriminate — see `benchmark_kmer_size.md`).
- For **PacBio HiFi, use k=6**: same ASVs and error model as k=5, ~4× faster, and
  far cheaper than k≥7. Document k=6 in the `--kmer-size` help for `--errfun
  pacbio`, mirroring the per-errfun guidance added for `--loess-*`.
- Avoid k ≥ 7 for HiFi, especially pooled or large-N runs (memory + benign-but-
  aggressive shrouding with no upside).

### Caveats

Single Illumina sample; single PacBio sample (full uniques). The benign-shrouding
result (identical ASVs at 79 % shroud) should be confirmed on a second sample or
a pooled run before changing any default. Scratch artifacts were purged from
`/tmp`; re-run to reproduce exact `err_out` diffs. The memory model and the
k=5-is-a-no-op-on-long-reads finding are robust (mechanistic + measured).

---

## Full MiSeq SOP confirmation (2026-06-01)

The caveat above (single Illumina sample) is now addressed: a sweep on the
**full MiSeq SOP** dataset (61 samples, learn-errors pooled per-sample then
denoised, R1 and R2 separately) via `comparison/run_kmer_sweep.sh`. Errfun
loess, band 16, k = 5/6/7/8.

### R1

| k | learn_iters | dada_aligns | dada_shrouded | shroud % | n_asv | wall (s) |
|---|------------|-------------|---------------|----------|-------|----------|
| 5 | 6  | 473 091 080 | 354 857 588 | 42.9 | 2514 | 238.8 |
| 6 | 6  | 476 322 636 | 400 154 558 | 45.7 | 2529 | 182.4 |
| 7 | 7  | 475 915 226 | 422 244 675 | 47.0 | 2525 | 236.4 |
| 8 | 6  | 475 626 964 | 437 773 435 | 47.9 | 2520 | 616.7 |

### R2

| k | learn_iters | dada_aligns | dada_shrouded | shroud % | n_asv | wall (s) |
|---|------------|-------------|---------------|----------|-------|----------|
| 5 | 10*| 398 005 481 | 296 195 546 | 42.7 | 1991 | 175.6 |
| 6 | 7  | 401 707 229 | 336 741 954 | 45.6 | 2004 | 137.5 |
| 7 | 6  | 404 084 016 | 360 688 524 | 47.2 | 2015 | 206.2 |
| 8 | 7  | 404 545 391 | 371 782 216 | 47.9 | 2016 | 538.1 |

\* R2/k5 hit the `--max-consist 10` cap. The log shows this is **oscillation,
not slow convergence**: max|err_in−err_out| reached 1.78e-6 at iter 6 (already
converged), then bounced back to a stuck 1.17e-5 two-state cycle through iter 10.
This is normal DADA2 behavior near the error-model noise floor (R `learnErrors`
does the same) — the iter-6 vs iter-10 model differs by ~1e-5, far below the
~1e-3 that affects ASV calls. The iteration counts are **non-monotonic in k**
(R1: 6/6/7/6; R2: 10/7/6/7), confirming this is run-to-run oscillation, not a
k effect — it is NOT evidence for k=6 over k=5 on correctness grounds.

### Findings at full scale

- **k's ASV impact stays minimal and non-monotonic**: R1 range 2514–2529 (15
  ASVs, 0.6 %); R2 range 1991–2016 (25 ASVs, 1.3 %). Holds across ~0.5 billion
  alignments — far larger than the original single-sample test.
- **Benign shrouding confirmed at scale**: shroud % climbs 43 → 48 % from k5→k8
  on both reads, yet ASV counts barely move — the screen removes pairs the
  alignment would not have clustered, not real variants.
- **The k=8 wall-time penalty is dramatic**: ~2.6× slower than k=6 (R1 617 s vs
  182 s; R2 538 s vs 138 s) despite *fewer* iterations — the 4⁸ = 65536-long
  k-mer dot product dominates. k=6 is fastest on both reads.
- **Per-iteration screen behavior** (from the learn logs): the first pass of
  every per-sample `dada` run is unscreened (`kdist_cutoff = 1.0`, by design, to
  seed cluster 0 — see `run_dada` in `dada.rs`), so it shows `0 shrouded` and
  fewer alignments; screening (and shrouding) kicks in only once clusters bud in
  later passes. This matches R DADA2's C++ and is not a bug.

### Open question still not closed

The slight ASV *increase* at higher k on R2 (k7=2015, k8=2016 vs k5=1991, +24–25
ASVs) is the one thing counts alone can't explain: are those extra ASVs real rare
close-variants the larger-k screen lets through, or borderline noise? Resolving
this needs an ASV-set *diff* (sequence-level), not just totals — pending.

### Conclusion (unchanged, now with large-N support)

k=5 stays the safe global default; **k=6 is the practical pick for speed/memory**
(fastest, lowest non-default footprint) with no ASV downside; k≥7 buys nothing on
Illumina and costs heavily at k=8. The recommendation is now backed by a 61-sample
pooled run, not a single sample.

---

## Downstream: merging paired reads across k (manual)

For paired-end Illumina the per-orientation ASVs must be merged into full-length
amplicons (`merge-pairs`). The k-mer screen does **not** run in `merge-pairs` —
it re-dereplicates and does its own ends-free NW with no k-mer screen — so k
affects merging only *indirectly*, through the upstream R1/R2 ASV differences.
The question merging answers is therefore: **do the per-k ASV differences survive
into the merged amplicons, or wash out?** This step is run by hand (it is
cross-orientation, so it consumes the R1 sweep dir + R2 sweep dir together).

The sweep writes one dada JSON per sample under `<outdir>/dada_k<k>/`. For each k,
merge the matching R1/R2 dada outputs, supplying the filtered FASTQs (re-dereplicated
to recover read→unique maps). Files match by position, so sort all four globs the
same way (here the shared sample stems guarantee it):

```bash
R1=out_R1   # the R1 sweep outdir
R2=out_R2   # the R2 sweep outdir
FQ1=filtered/R1   # filtered forward FASTQs (one per sample)
FQ2=filtered/R2   # filtered reverse FASTQs

for k in 5 6 7 8; do
  dada2-rs merge-pairs \
    --fwd-dada  "$R1"/dada_k${k}/*.json \
    --rev-dada  "$R2"/dada_k${k}/*.json \
    --fwd-fastq "$FQ1"/*.fastq.gz \
    --rev-fastq "$FQ2"/*.fastq.gz \
    --min-overlap 12 --max-mismatch 0 --threads 8 \
    -o merged_k${k}.json
done
```

Output schema (per k): `{ "samples": [ { "sample", "total_pairs",
"accepted_pairs", "num_merged", "merged": [ {"sequence","abundance","accept"} ] } ] }`.

### Metrics to compare across k

**(a) merged ASV count + merge rate** — pooled over samples:

```bash
for k in 5 6 7 8; do
  python3 - "$k" merged_k${k}.json <<'PY'
import json, sys
k, path = sys.argv[1], sys.argv[2]
d = json.load(open(path))
seqs, tot, acc = set(), 0, 0
for s in d["samples"]:
    tot += s["total_pairs"]; acc += s["accepted_pairs"]
    seqs |= {m["sequence"] for m in s["merged"] if m.get("accept")}
print(f"k={k}  merged_ASVs={len(seqs)}  merge_rate={acc/tot:.3f}  ({acc}/{tot})")
PY
done
```

**(b) merged ASV-set diff vs k=5** — the decisive sequence-level check (the
cross-orientation analogue of the still-open R2 question above):

```bash
python3 - <<'PY'
import json
def merged_set(path):
    d = json.load(open(path)); out=set()
    for s in d["samples"]:
        out |= {m["sequence"] for m in s["merged"] if m.get("accept")}
    return out
base = merged_set("merged_k5.json")
print(f"baseline k=5: {len(base)} merged ASVs")
for k in (6,7,8):
    s = merged_set(f"merged_k{k}.json")
    print(f"  k={k}: shared={len(base&s)} only_k5={len(base-s)} "
          f"only_k{k}={len(s-base)} "
          f"[{'IDENTICAL' if base==s else 'DIFFERS'}]")
PY
```

If the merged sets come out IDENTICAL (or differ by ≪ the pre-merge ASV-count
spread), the upstream k-sensitivity washes out at the amplicon level — the
strongest possible "k doesn't matter for final output" result. If they differ,
(b) is exactly the sequence-level evidence needed to close the open question above.

### Worked example (2 MiSeq SOP samples, k=5 vs k=6, verified)

Filtered F3D0+F3D1 (truncLen 240/160), swept each orientation, merged per the
commands above (min-overlap 12, max-mismatch 0):

```
k=5  merged_ASVs=144  merge_rate=0.955  (11630/12173)
k=6  merged_ASVs=144  merge_rate=0.955  (11631/12177)
  k=6 vs k=5: shared=144 only_k5=0 only_k6=0  [IDENTICAL]
```

The merged amplicon set is **identical** and the merge rate is unchanged, even
though k can nudge the pre-merge R1/R2 ASV counts. So on this data the upstream
k-sensitivity washes out completely after merging — pairing is robust to the
k-mer-screen size. (Worth repeating on the full 61-sample set and across k=7/8 to
confirm, but the 2-sample result is a strong indicator.)
