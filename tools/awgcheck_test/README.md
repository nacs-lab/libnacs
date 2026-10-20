# AWG RF power check in `pre_run`: offline tests

Branch `awg-power-check` (on top of `coil-check` 87008f2). Code: `lib/nacs-seq/awg_check.cpp`,
hooks in `manager.cpp` / `manager.h`, `Error::BasicSeq::AwgLimit` in `error.h`. Inert unless
`config.yml` has an `awg_power:` block. Built and tested offline on 2026-10-19 only; nothing
was installed. The live DLL is still the coil-check build (md5 `cab68153...` / `e5f56282...`,
identical to `build/output-bin`).

## What it checks

Per basic sequence, after `host_seq.pre_run()` (so with this shot's globals, including the
rearrangement move targets) and after the coil check, before any device gets the bseq:

- **power**: for each output of each listed AWG (`AWG1/OUT0`, `AWG1/OUT1`), the sum over its
  tones of amp^2 stays <= `max_sum_amp2` (0.125) at every time. This is `awgInterlock.m`, but it now
  applies to every AMP write, including `seqFuns.rampAWG_amp` / `rampAWG_amp_sqrt` (interlock
  commented out there) and raw `seq.add('AWG1/OUT*/CHN*/AMP', ...)`.
- **frequency** (optional, per output): while a tone's amp is above `amp_on` (0), its
  frequency stays in the window (`awgFreqInterlock.m`: Na OUT0 100-185 MHz, Cs OUT1 102-171 MHz).
  A tone parked at 0 Hz with amp 0 is fine. The rearrangement relabel does exactly that.

A violation throws `Error(BasicSeq, AwgLimit)`. MATLAB sees
`AWG power check refused this sequence: AWG1/OUT0 sum(amp^2) 0.1400 > 0.1250 at 6.300 ms;`.
That is the same path as the coil refusal: no device is sent anything for that bseq.

How: a quick bound first, with each tone at its own peak, all at once. Only when it is above the warn
level does a time sweep run over the tones' change points (sorted runs merged, running sum). Ramps
are sampled at `pts_per_ramp` (32) points, and between two samples of one ramp the larger end
counts, so a monotonic ramp (`rampTo`, `rampToSqrt`) is never under-counted.

## Proposed config block (`matlab_new/config.yml`)

```yaml
# AWG RF power check in libnacs pre_run (awg-power-check build). Same rules as
# awgInterlock.m / awgFreqInterlock.m. Delete this block to turn it off.
awg_power:
  enabled: true
  devices: [AWG1]               # AWG2 is unused by Experiments/
  max_sum_amp2: 0.125           # per output, sum over its tones of amp^2
  freq_windows:                 # Hz, applied while a tone's amp > 0
    AWG1/OUT0: [100e6, 185e6]   # Na
    AWG1/OUT1: [102e6, 171e6]   # Cs
```

## Results (2026-10-19)

Corpus: the 20 exported points of the coil-check work (`C:\ec-coil-libnacs\debug\coil_libnacs_parity*`,
9 scans: Feshbach x3, Align, cooling, loading, merge, rearr, RSC rspec), plus 20 random loading
patterns on the sRearr branch (8 of 8-16 loaded sites moving per species).

| test | result |
|---|---|
| default limits | 20/20 points PASS, 20/20 rearr patterns PASS; real peak sum(amp^2) 0.075 on every output (60 % of the limit) |
| time sweep vs independent recount from `get_nominal_output` (the .seq dump path) | equal to 1e-12 on every bseq and output |
| limit just below the real peak | refused on every point, naming the right output, at or before the peak's bseq |
| frequency window moved above a used tone's frequency | refused on every point |
| device output (`equiv_run_awg.py`): installed coil-check build vs this build, check off and on | 590/590 files identical (FPGA bytecode, NI-DAQ samples, IR dumps) |
| coil check, this build vs the installed build (`nacs-seq-coil-check`, 86 configs) | identical results |

Cost. On the loading/imaging bseq the quick bound passes: 1-4 us. On the rearrangement
branch the bound does not pass, so the sweep runs. Tones hand over between sites (the
relabel steps some tones down while others step up), so each tone at its own peak gives
0.15, but the real sum is 0.075.

| sRearr branch (bseq 2), 500 shots, 4 loading patterns | check |
|---|---|
| median | 21-26 us |
| p99 | 41-52 us |
| worst after the first shot | 44-55 us |
| first shot after a compile (finds the AWG pulses, sizes the buffers) | 62-81 us |

Real-backend `pre_run` of the same branch (FPGA bytecode + NI-DAQ samples; dummy AWG, dead FPGA
url, no `start()`; `nacs-seq-devdump ... 300`): median ~105 us without the check, ~123 us with it.
The AWG upload round trip and the MATLAB image analysis before `pre_run` come on top of that and
are not part of this measurement.

How the cost was brought down: an earlier version (full grid x tones walk, frequency walk across
all tones every shot, per-shot allocation) cost ~60 us median on the branch. The frequency walk
ran every shot because the relabel parks unused tones at 0 Hz.

## Not covered

- Only what this sequence sends. Tones the sequence never touches, and whatever the AWG server
  plays between shots, are not seen. Each shot starts from the sequence's default values.
- A non-monotonic AMP ramp that peaks between two of its 32 samples could be under-counted
  slightly. All current AMP ramps are monotonic.
- AWG analog (`A<k>`) and PHASE channels are ignored; AWG2 is not listed (unused).

## Re-run

```bash
bash build_isolated_awg.sh nacs-seq-awg-check-tool nacs-seq-coil-check-tool   # msys2 MINGW64; own folder build-awg/
python tools/awgcheck_test/awg_check_run.py C:/projects/libnacs-coilcheck/build-awg/output-bin/nacs-seq-awg-check.exe C:/ec-coil-libnacs/debug/coil_libnacs_parity C:/ec-coil-libnacs/debug/coil_libnacs_parity2 --out <dir>
gcc -O2 -o tools/awgcheck_test/nacs-seq-devdump.exe tools/nacs_seq_devdump.c
python tools/awgcheck_test/equiv_run_awg.py <dir>
```

`nacs-seq-awg-check` refuses any config without `use_dummy_device: true`. The devdump configs here
(`cfg_real_*`, `cfg_equiv_*`) use dummy AWGs and a dead FPGA url, and devdump never calls `start()`.

## Install (operator, later)

Same as the coil check (docs/proposals/coil-check-in-libnacs.md in experiment-control):

1. Back up the installed `libnacs-seq.dll` / `libnacs-utils.dll` (the coil-check build).
2. Close every MATLAB holding the DLL (`tasklist /m libnacs-utils.dll`; generation via
   `restart_generation.ps1 -NoLaunch`, which runs AndorDisconnect), copy the two DLLs from
   `build-awg/output-bin`, relaunch. With no `awg_power` block, behaviour is unchanged (the coil
   check stays on).
3. Add the `awg_power` block and run one normal scan with rearrangement. Then do a deliberate refusal
   (`max_sum_amp2: 0.05`, below the real 0.075): the first shot must be refused at `pre_run` with
   0 shots taken. Then restore 0.125.

Rollback: copy the backup DLLs back, or delete the block.
