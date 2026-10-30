#!/usr/bin/env python3
"""Compare GaiaVario option variants on several receivers by what matters for a
displacement waveform: the INTEGRATED series against the table accelerometer,
the noise floor before the event, and the serial correlation of the increments.

    python3 fix_table.py ACCEL.mseed LAT LON T0 T1 OUTDIR ID [ID ...]

OUTDIR/<ID>_va.txt (PPPHVA) and OUTDIR/<ID>_<variant>.csv for every variant.
Prints one table per receiver and an average over receivers.
"""
import sys, os, math, glob, datetime
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mseed_float64 import read_mseed

accel_path, lat, lon, t0u, t1u, outdir = sys.argv[1:7]; ids = sys.argv[7:]
lat = float(lat); lon = float(lon); t0u = float(t0u); t1u = float(t1u)
LEAP = 18.0; fs = 10.0; MARGIN = 15.0; HP, LP = 0.1, 2.0
sp, cp, sl, cl = math.sin(math.radians(lat)), math.cos(math.radians(lat)), math.sin(math.radians(lon)), math.cos(math.radians(lon))
R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])

def read_va(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or not line.strip(): continue
        p = line.split(); d = np.array([float(p[0 + 2]), float(p[3]), float(p[4])])
        t.append(float(p[0]) - LEAP); v.append(R @ d)
    return np.array(t), np.array(v)

def read_vario(path):
    t, v, ns = [], [], []
    for line in open(path):
        if line.startswith("#") or line.startswith("sod,") or not line.strip(): continue
        p = line.rstrip("\n").split(","); t.append(float(p[0]) - LEAP); v.append([float(p[2]), float(p[3]), float(p[4])]); ns.append(int(p[6]))
    return np.array(t), np.array(v), np.array(ns)

g0, g1 = t0u - MARGIN, t1u + MARGIN
grid = np.arange(math.ceil(g0 * fs), math.floor(g1 * fs) + 1) / fs
def to_grid(t, v):
    out = np.full((len(grid), 3), np.nan)
    idx = np.round((t - grid[0]) * fs).astype(int); ok = (idx >= 0) & (idx < len(grid)); out[idx[ok]] = v[ok]; return out
def _bp(x, fs_, hp, lp, order=0):
    n = len(x); N = (1 << (n - 1).bit_length()) * 2
    X = np.fft.rfft(x - np.mean(x[:int(fs_ * 5)]), N); f = np.fft.rfftfreq(N, 1 / fs_)
    H = np.ones_like(f)
    if hp > 0: H *= (f / hp) ** 4 / (1 + (f / hp) ** 4)
    if lp > 0: H *= 1 / (1 + (f / lp) ** 4)
    if order:
        with np.errstate(divide="ignore", invalid="ignore"):
            W = (1j * 2 * np.pi * f) ** (-order); W[0] = 0
        H = H * W
    return np.fft.irfft(X * H, N)[:n]
tr = read_mseed(accel_path); chan = {c: k for k in tr for c in ("E", "N", "Z") if k.endswith(c)}
base = (datetime.datetime(2025, 10, 21) - datetime.datetime(1970, 1, 1)).total_seconds()
def accel_disp(name):
    t0, rate, x, _ = tr[name]; ts = t0 - base + np.arange(len(x)) / rate; m = (ts >= g0 - 60) & (ts <= g1 + 60)
    return np.interp(grid, ts[m], _bp(x[m], rate, HP, LP, 2))
truth = np.column_stack([accel_disp(chan["N"]), accel_disp(chan["E"]), accel_disp(chan["Z"])])
ROT = float(os.environ.get("ROT_DEG", "350"))
c_, s_ = math.cos(math.radians(ROT)), math.sin(math.radians(ROT))
truth = np.column_stack([c_ * truth[:, 0] - s_ * truth[:, 1], s_ * truth[:, 0] + c_ * truth[:, 1], truth[:, 2]])
inw = (grid >= t0u) & (grid <= t1u)
env = np.sqrt(truth[:, 0] ** 2 + truth[:, 1] ** 2); k2 = int(fs * 10)
sm = np.convolve(env ** 2, np.ones(2 * k2 + 1) / (2 * k2 + 1), mode="same"); tpk = grid[np.argmax(np.where(inw, sm, 0))]
strong = (grid >= tpk - 25) & (grid <= tpk + 35) & inw
quiet = (grid >= t0u + 3) & (grid <= tpk - 30)
tail = (grid >= tpk + 90) & (grid <= t1u)               # quiet tail for the increment statistics

def ac1(x):
    x = x - x.mean(); return float(np.dot(x[:-1], x[1:]) / np.dot(x, x))

def analyse(t, v):
    x = to_grid(t, v); xs = np.where(np.isfinite(x), x, 0.0); cum = np.cumsum(xs, axis=0)
    d = np.column_stack([_bp(cum[:, c], fs, HP, LP) for c in range(3)])
    corr = [float(np.corrcoef(d[strong, c], truth[strong, c])[0, 1]) for c in range(3)]
    rmse = [1000 * math.sqrt(np.mean((d[strong, c] - truth[strong, c]) ** 2)) for c in range(3)]
    noise = [1000 * math.sqrt(np.nanmean(d[quiet, c] ** 2)) for c in range(3)]
    inc = x[tail]; inc = inc[np.isfinite(inc[:, 0])]
    a1 = [ac1(inc[:, c]) for c in range(3)]
    sig = [1000 * 1.4826 * np.median(np.abs(inc[:, c] - np.median(inc[:, c]))) for c in range(3)]
    return corr, rmse, noise, a1, sig, float(np.isfinite(x[inw, 0]).mean())

order = os.environ.get("ORDER", "va,base,k6,huber,hyst5,taper5,admit5,admit10,comboH,comboN,mean,meanc,l1,l1c").split(",")
agg = {}
for rid in ids:
    print(f"\n### {rid}\n")
    print("| variant | corr N/E/U (strong) | RMSE N/E/U mm (strong) | noise floor N/E/U mm | per-epoch sigma N/E/U mm | lag-1 autocorr | epochs solved | sats used |")
    print("|---|---|---|---|---|---|---|---|")
    for name in order:
        path = f"{outdir}/{rid}_va.txt" if name == "va" else f"{outdir}/{rid}_{name}.csv"
        if not os.path.exists(path): continue
        if name == "va":
            t, v = read_va(path); ns = None
        else:
            t, v, ns = read_vario(path)
        corr, rmse, noise, a1, sig, solved = analyse(t, v)
        used = "" if ns is None else f"{np.median(ns):.0f}"
        print(f"| {name} | " + " / ".join(f"{x:.2f}" for x in corr) + " | " + " / ".join(f"{x:.1f}" for x in rmse) + " | " +
              " / ".join(f"{x:.1f}" for x in noise) + " | " + " / ".join(f"{x:.2f}" for x in sig) + " | " + " / ".join(f"{x:+.2f}" for x in a1) +
              f" | {100*solved:.0f} % | {used} |")
        agg.setdefault(name, []).append((np.mean(corr), np.mean(rmse), np.mean(noise), np.mean(a1), np.mean(sig)))
print("\n### average over receivers (mean of N/E/U)\n")
print("| variant | corr | RMSE mm | noise floor mm | per-epoch sigma mm | lag-1 autocorr |")
print("|---|---|---|---|---|---|")
for name in order:
    if name in agg:
        a = np.mean(agg[name], axis=0)
        print(f"| {name} | {a[0]:.3f} | {a[1]:.2f} | {a[2]:.2f} | {a[4]:.2f} | {a[3]:+.3f} |")
