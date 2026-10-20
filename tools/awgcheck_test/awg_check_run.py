"""Offline tests of the AWG RF power check (lib/nacs-seq/awg_check.cpp).

    python awg_check_run.py <exe> <corpus dir> [<corpus dir> ...] --out <dir>

<exe>         build-awg/output-bin/nacs-seq-awg-check.exe
<corpus dir>  a coil_libnacs_parity* folder: <scan>/coil.json lists the points, each
              with its serialized sequence (s<id>.bin) and globals (p<point>.txt)

Every run uses `use_dummy_device: true` (the driver refuses anything else). Per point:
  default   the proposed limits: verdict, quick bound / swept sum, check time
  sweep     warn level 0, so the time sweep runs on every output; its peak sum(amp^2)
            is compared with an independent recount from get_nominal_output (the .seq
            dump path) of the same shot
  tight     max_sum_amp2 just below the swept peak: must refuse, naming that output
  freq      one output's window moved just above a frequency in use while its tone is
            on: must refuse
  off       check disabled: pre_run time without it
Rearrangement: the same branch with random loading patterns (move targets), default
limits, for verdicts and the check time on the branch.
"""

import argparse
import json
import os
import random
import subprocess
import sys

BASE = """tick_per_sec: 1000000000000
use_dummy_device: true
"""

WINDOWS = {"AWG1/OUT0": [100e6, 185e6], "AWG1/OUT1": [102e6, 171e6]}


def cfg_text(enabled=True, max_sum=0.125, warn=None, windows=WINDOWS, pts=None):
    s = BASE + "awg_power:\n  enabled: %s\n  devices: [AWG1]\n" % ("true" if enabled else "false")
    s += "  max_sum_amp2: %.17g\n" % max_sum
    if warn is not None:
        s += "  warn_sum_amp2: %.17g\n" % warn
    if pts is not None:
        s += "  pts_per_ramp: %d\n" % pts
    if windows:
        s += "  freq_windows:\n"
        for k, (lo, hi) in windows.items():
            s += "    %s: [%.17g, %.17g]\n" % (k, lo, hi)
    return s


def run(exe, binf, cfg, globs, outdir, tag, shots=20, trace=None):
    cpath = os.path.join(outdir, "cfg_%s.yml" % tag)
    with open(cpath, "w", newline="\n") as f:
        f.write(cfg)
    args = [exe, binf, cpath, globs or "-", str(shots)]
    if trace:
        args.append(trace)
    p = subprocess.run(args, capture_output=True, text=True, timeout=600)
    last = [l for l in p.stdout.splitlines() if l.startswith("{")]
    if p.returncode != 0 or not last:
        raise RuntimeError("%s failed (%d): %s" % (tag, p.returncode, p.stderr[-2000:]))
    return json.loads(last[-1])


def recount(trace_path):
    """Independent sum(amp^2) peak per bseq and output from get_nominal_output samples.
    Between two samples at different times the larger |amp| is used (the value ramps
    or holds between them)."""
    res = {}
    for b in json.load(open(trace_path)):
        outs = {}
        for name, pts in b["chns"].items():
            parts = name.split("/")
            if len(parts) != 4 or parts[3] != "AMP":
                continue
            outs.setdefault("/".join(parts[:2]), []).append(pts)
        per = {}
        for out, tones in outs.items():
            grid = sorted({t for pts in tones for t, _ in pts})
            best = (0.0, 0)
            idx = [0] * len(tones)
            for t in grid:
                s = 0.0
                for k, pts in enumerate(tones):
                    while idx[k] + 1 < len(pts) and pts[idx[k] + 1][0] <= t:
                        idx[k] += 1
                    a = abs(pts[idx[k]][1])
                    if idx[k] + 1 < len(pts) and pts[idx[k] + 1][0] > pts[idx[k]][0]:
                        a = max(a, abs(pts[idx[k] + 1][1]))
                    s += a * a
                if s > best[0]:
                    best = (s, t)
            per[out] = {"sum_amp2": best[0], "t_ms": best[1] / b["tick_per_sec"] * 1e3}
        res[b["bseq_idx"]] = per
    return res


def on_freq(trace_path):
    """(bseq_idx, out, freq) of a tone that is on (amp > 0) during the whole bseq."""
    for b in json.load(open(trace_path)):
        chns = b["chns"]
        for name, pts in sorted(chns.items()):
            parts = name.split("/")
            if len(parts) != 4 or parts[3] != "AMP":
                continue
            if min(abs(v) for _, v in pts) <= 0:
                continue
            fpts = chns.get("/".join(parts[:3]) + "/FREQ")
            if fpts:
                return b["bseq_idx"], "/".join(parts[:2]), max(v for _, v in fpts)
    return None


def corpus(dirs):
    for d in dirs:
        for scan in sorted(os.listdir(d)):
            cj = os.path.join(d, scan, "coil.json")
            if not os.path.exists(cj):
                continue
            for p in json.load(open(cj))["points"]:
                binf = os.path.join(d, scan, p["bin"])
                globs = os.path.join(d, scan, "p%d.txt" % p["point"])
                yield scan, p["point"], binf, globs if os.path.exists(globs) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exe")
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--patterns", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    bad = []
    for scan, point, binf, globs in corpus(a.dirs):
        tag = "%s_p%d" % (scan, point)
        wd = os.path.join(a.out, tag)
        os.makedirs(wd, exist_ok=True)
        row = {"scan": scan, "point": point}
        d = run(a.exe, binf, cfg_text(), globs, wd, "default")
        row["default"] = {k: d[k] for k in ("verdict", "refusal", "outs", "bseqs", "timing")}
        if d["verdict"] != "PASS":
            bad.append("%s default: %s %s" % (tag, d["verdict"], d["refusal"]))
        tr = os.path.join(wd, "trace.json")
        sw = run(a.exe, binf, cfg_text(warn=0), globs, wd, "sweep", shots=1, trace=tr)
        rc = recount(tr)
        worst = (0.0, None, None)
        cmp = []
        for bi, b in enumerate(sw["bseqs"]):
            for o in b["outs"]:
                ref = rc.get(bi, {}).get(o["name"])
                if ref is None:
                    bad.append("%s sweep: no recount for bseq %d %s" % (tag, bi, o["name"]))
                    continue
                diff = abs(ref["sum_amp2"] - o["sum_amp2"])
                cmp.append({"bseq": bi, "out": o["name"], "cpp": o["sum_amp2"],
                            "recount": ref["sum_amp2"], "bound": o["bound"]})
                if diff > 1e-12:
                    bad.append("%s sweep: bseq %d %s C++ %.17g recount %.17g" %
                               (tag, bi, o["name"], o["sum_amp2"], ref["sum_amp2"]))
                if o["sum_amp2"] > worst[0]:
                    worst = (o["sum_amp2"], bi, o["name"])
        row["sweep"] = cmp
        if worst[1] is not None:
            lim = worst[0] * (1 - 1e-6)
            t = run(a.exe, binf, cfg_text(max_sum=lim), globs, wd, "tight", shots=1)
            row["tight"] = {"limit": lim, "verdict": t["verdict"], "refusal": t["refusal"],
                            "refused_bseq_idx": t["refused_bseq_idx"]}
            if t["verdict"] != "FAIL" or worst[2] not in t["refusal"] or \
               t["refused_bseq_idx"] > worst[1]:
                bad.append("%s tight: expected a refusal of %s by bseq %d, got %s %s" %
                           (tag, worst[2], worst[1], t["verdict"], t["refusal"]))
        f = on_freq(tr)
        if f is not None:
            bi, out, freq = f
            w = dict(WINDOWS)
            w[out] = [freq + 1.0, 1e12]
            fr = run(a.exe, binf, cfg_text(windows=w), globs, wd, "freq", shots=1)
            row["freq"] = {"out": out, "freq": freq, "verdict": fr["verdict"],
                           "refusal": fr["refusal"], "refused_bseq_idx": fr["refused_bseq_idx"]}
            if fr["verdict"] != "FAIL" or "frequency" not in fr["refusal"] or \
               out not in fr["refusal"] or fr["refused_bseq_idx"] > bi:
                bad.append("%s freq: expected a refusal of %s, got %s %s" %
                           (tag, out, fr["verdict"], fr["refusal"]))
        off = run(a.exe, binf, cfg_text(enabled=False), globs, wd, "off")
        row["off_timing"] = off["timing"]
        rows.append(row)
        print("%-45s %s  bound/sum %s  check_us %s" % (
            tag, d["verdict"],
            ["%s %.4f/%.4f" % (c["out"][-4:], c["bound"], c["cpp"]) for c in cmp][:4],
            [round(x["check_us_median"], 1) for x in d["timing"]]), flush=True)

    # Rearrangement branch with random loading patterns.
    rearr = []
    for scan, point, binf, globs in corpus(a.dirs):
        if not scan.startswith("rearrScan") or globs is None:
            continue
        base = [l.split() for l in open(globs).read().split("\n") if l.strip()]
        ids = [int(i) for i, _ in base]
        rng = random.Random(1)
        wd = os.path.join(a.out, "rearr_patterns")
        os.makedirs(wd, exist_ok=True)
        for k in range(a.patterns):
            # Same layout as the export: two blocks of 16 sites (Na, Cs), value = target
            # index 1..8 for a loaded site that moves, 0 = dropped.
            vals = {}
            for blk in (ids[:16], ids[16:32]):
                loaded = sorted(rng.sample(range(16), rng.randint(8, 16)))
                chosen = loaded[:8]
                for j, gid in enumerate(blk):
                    vals[gid] = chosen.index(j) + 1 if j in chosen else 0
            gp = os.path.join(wd, "g%d.txt" % k)
            with open(gp, "w", newline="\n") as f:
                for gid in ids:
                    f.write("%d %s\n" % (gid, vals.get(gid, 0)))
            r = run(a.exe, binf, cfg_text(), gp, wd, "p%d" % k)
            rearr.append({"pattern": k, "verdict": r["verdict"], "refusal": r["refusal"],
                          "timing": r["timing"],
                          "sweep": [[o["name"], o["bound"], o["sum_amp2"]]
                                    for b in r["bseqs"] for o in b["outs"]]})
            if r["verdict"] != "PASS":
                bad.append("rearr pattern %d: %s %s" % (k, r["verdict"], r["refusal"]))
        br = [x["timing"][-1] for x in rearr]
        print("rearr patterns: %d, branch check_us median of medians %.1f, worst max %.1f" % (
            len(rearr), sorted(b["check_us_median"] for b in br)[len(br) // 2],
            max(b["check_us_max"] for b in br)))

    res = {"points": rows, "rearr_patterns": rearr, "problems": bad}
    with open(os.path.join(a.out, "results.json"), "w", newline="\n") as f:
        json.dump(res, f, indent=1)
    print("problems: %d" % len(bad))
    for b in bad:
        print("  " + b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
