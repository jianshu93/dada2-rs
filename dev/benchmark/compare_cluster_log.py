#!/usr/bin/env python3
"""Compare dada/learn-errors verbose cluster logs across alignment backends.

Quantifies E-M convergence behaviour from the `--verbose` "New Cluster" output
so backend A/B differences (e.g. NW vs experimental WFA) are read off the logs
instead of eyeballed. Two signals:

  births     — number of "New Cluster" lines = cluster-birth steps summed over
               all dada calls / self-consistency rounds. A large delta between
               backends is a convergence-TRAJECTORY difference (the slightly
               different alignments push learn-errors down a different number of
               rounds), not a per-step instability.
  shuffle-S  — each 'S' on a "New Cluster" line is one b_shuffle2 round (reads
               reassigned among clusters until the partition is stable). The
               textual tokens on those lines have no capital S, so the count is
               exact. mean S/birth is the per-step convergence cost; if it
               matches across backends, the backend is NOT destabilising E-M
               even if the raw 'S' runs look long.

Usage:
  compare_cluster_log.py LABEL=path/to.log [LABEL=path/to.log ...]
  # e.g. compare across arms and read directions:
  compare_cluster_log.py \\
      nw-fwd=nw/learn_fwd.log  e50-fwd=e50/learn_fwd.log \\
      nw-rev=nw/learn_rev.log  e50-rev=e50/learn_rev.log

Background: dada2-rs issue #51 (WFA edit-budget cap) scale validation. The MiSeq
finding was per-birth shuffle cost identical NW vs WFA (~2.3 fwd / ~2.5 rev);
the only real divergence was total births on the lower-quality REVERSE reads.
"""
import sys


def summarize(path):
    births = 0
    shuffle_s = 0
    max_run = 0
    long_runs = 0  # births that took >=4 shuffle rounds
    try:
        fh = open(path)
    except OSError as e:
        return None, str(e)
    with fh:
        for line in fh:
            if "New Cluster" in line:
                births += 1
                c = line.count("S")  # all S on this line are shuffle markers
                shuffle_s += c
                max_run = max(max_run, c)
                if c >= 4:
                    long_runs += 1
    mean = shuffle_s / births if births else 0.0
    return (births, shuffle_s, mean, max_run, long_runs), None


def main(argv):
    if len(argv) < 2 or any("=" not in a for a in argv[1:]):
        sys.exit(f"usage: {argv[0]} LABEL=path.log [LABEL=path.log ...]")
    rows = []
    width = max(len(a.split("=", 1)[0]) for a in argv[1:])
    for spec in argv[1:]:
        label, path = spec.split("=", 1)
        stats, err = summarize(path)
        if err:
            print(f"{label:>{width}}: <skip: {err}>")
            continue
        births, s, mean, mx, long = stats
        rows.append((label, births, s, mean, mx, long))

    if not rows:
        return
    hdr = f"{'label':>{width}}  {'births':>8} {'shuffleS':>9} {'mean S/birth':>13} {'max_run':>8} {'>=4S':>7}"
    print(hdr)
    print("-" * len(hdr))
    for label, births, s, mean, mx, long in rows:
        print(
            f"{label:>{width}}  {births:>8d} {s:>9d} {mean:>13.2f} {mx:>8d} {long:>7d}"
        )
    print(
        "\nread: matching 'mean S/birth' across backends = E-M equally stable "
        "per step\n      (long S runs are normal, not a backend artifact); a "
        "large 'births'\n      delta = convergence-trajectory difference, "
        "watch it on lower-quality\n      (reverse / noisier) reads where "
        "backends diverge most."
    )


if __name__ == "__main__":
    main(sys.argv)
