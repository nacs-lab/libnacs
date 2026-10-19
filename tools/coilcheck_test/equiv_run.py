"""Byte-for-byte device output: live DLL (backup copy) vs the coil-check build.

For every exported point (C:\\ec-coil-libnacs\\debug\\coil_libnacs_parity*\\<case>\\coil.json),
run nacs-seq-devdump three times and compare every output file:
  old      live libnacs-utils.dll (backup), config without the coil block
  new_off  coil-check build, same config (check disabled)
  new_on   coil-check build, config + `feshbach_coil: enabled` (check active)
"""
import filecmp
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "nacs-seq-devdump.exe")
DLLS = {
    "old": r"C:\msys64\mingw64\bin\libnacs-backup-20261019-live92ebcee\libnacs-utils.dll",
    "new_off": r"C:\projects\libnacs-coilcheck\build\output-bin\libnacs-utils.dll",
    "new_on": r"C:\projects\libnacs-coilcheck\build\output-bin\libnacs-utils.dll",
}
CFG = {"old": "cfg_equiv_off.yml", "new_off": "cfg_equiv_off.yml", "new_on": "cfg_equiv_on.yml"}
OUT = os.path.join(HERE, "equiv")

total = same = 0
bad = []
for coil in sorted(glob.glob(r"C:\ec-coil-libnacs\debug\coil_libnacs_parity*\*\coil.json")):
    case_dir = os.path.dirname(coil)
    case = os.path.basename(case_dir)
    res = json.load(open(coil))
    pts = res.get("points") or []
    pts = pts if isinstance(pts, list) else [pts]
    for pt in pts:
        tag = "%s_p%s" % (case, pt["point"])
        dirs = {}
        for v in DLLS:
            d = os.path.join(OUT, tag, v)
            os.makedirs(d, exist_ok=True)
            r = subprocess.run([EXE, DLLS[v], os.path.join(case_dir, pt["bin"]),
                                os.path.join(HERE, CFG[v]),
                                os.path.join(case_dir, "p%s.txt" % pt["point"]), d],
                               capture_output=True, text=True, timeout=600)
            if r.returncode != 0 or "ok:" not in r.stdout:
                bad.append("%s %s: exit %d %s" % (tag, v, r.returncode, (r.stderr or r.stdout)[-300:]))
            dirs[v] = d
        files = sorted(os.listdir(dirs["old"]))
        nb = len([f for f in files if f.endswith("_zynq.bin")])
        nz = sum(os.path.getsize(os.path.join(dirs["old"], f)) for f in files if f.endswith("_zynq.bin"))
        nn = sum(os.path.getsize(os.path.join(dirs["old"], f)) for f in files if f.endswith("_nidaq.bin"))
        diffs = []
        for v in ("new_off", "new_on"):
            other = sorted(os.listdir(dirs[v]))
            if other != files:
                diffs.append("%s file list differs" % v)
                continue
            for f in files:
                total += 1
                a, b = os.path.join(dirs["old"], f), os.path.join(dirs[v], f)
                if f.endswith(".txt") and not f.startswith("bseq"):
                    # IR text dumps: the live DLL prints the per-bseq "Assumptions" set in
                    # memory order, which changes from run to run of the SAME DLL; compare
                    # the lines as a multiset.
                    ok = sorted(open(a, "rb").read().split(b"\n")) == sorted(open(b, "rb").read().split(b"\n"))
                else:
                    ok = filecmp.cmp(a, b, shallow=False)
                if ok:
                    same += 1
                else:
                    diffs.append("%s %s" % (v, f))
        print("%-44s %d bseqs, FPGA bytecode %7d B, NI-DAQ %8d B, files %2d: %s" % (
            tag[:44], nb, nz, nn, len(files), "IDENTICAL" if not diffs else "DIFF " + ", ".join(diffs[:6])),
            flush=True)
        if diffs:
            bad.append("%s: %s" % (tag, diffs))
print("\ncompared %d files, %d identical; problems: %d" % (total, same, len(bad)))
for b in bad:
    print("  " + b)
sys.exit(1 if bad else 0)
