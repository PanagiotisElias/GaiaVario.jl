#!/usr/bin/env python3
"""PPPHVA in dual-frequency (DFV) and single-frequency (SFV) mode against GaiaVario L1 and the mean of two phases,
on the shaking-table replay: correlation with the table accelerometer, RMSE, noise floor, per-epoch sigma, sigma of
10-s sums (the telescoping effect), lag-1 autocorrelation and PGD (horizontal and 3-D, strong motion).

    python3 sfv_table.py ACCEL.mseed LAT LON T0 T1 DFV_DIR SFV_DIR ID [ID ...]

DFV_DIR/<ID>_va.txt (PPPHVA DFV), DFV_DIR/<ID>_l1.csv, <ID>_mean.csv; SFV_DIR/<ID>_va.txt (PPPHVA SFV).
"""
import sys, os, math, datetime
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mseed_float64 import read_mseed

accel_path, lat, lon, t0u, t1u, dfv, sfv = sys.argv[1:8]; ids = sys.argv[8:]
lat = float(lat); lon = float(lon); t0u = float(t0u); t1u = float(t1u)
LEAP = 18.0; fs = 10.0; MARGIN = 15.0; HP, LP = 0.1, 2.0
sp, cp, sl, cl = math.sin(math.radians(lat)), math.cos(math.radians(lat)), math.sin(math.radians(lon)), math.cos(math.radians(lon))
R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])
g0, g1 = t0u - MARGIN, t1u + MARGIN
grid = np.arange(math.ceil(g0 * fs), math.floor(g1 * fs) + 1) / fs

def read_va(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or not line.strip(): continue
        p = line.split(); t.append(float(p[0]) - LEAP); v.append(R @ np.array([float(p[2]), float(p[3]), float(p[4])]))
    return np.array(t), np.array(v)
def read_vario(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or line.startswith("sod,") or not line.strip(): continue
        p = line.rstrip("\n").split(","); t.append(float(p[0]) - LEAP); v.append([float(p[2]), float(p[3]), float(p[4])])
    return np.array(t), np.array(v)
def to_grid(t, v):
    out = np.full((len(grid), 3), np.nan); idx = np.round((t - grid[0]) * fs).astype(int); ok = (idx >= 0) & (idx < len(grid)); out[idx[ok]] = v[ok]; return out
def _bp(x, fs_, hp, lp, order=0):
    n = len(x); N = (1 << (n - 1).bit_length()) * 2
    X = np.fft.rfft(x - np.mean(x[:int(fs_ * 5)]), N); f = np.fft.rfftfreq(N, 1 / fs_)
    H = (f / hp) ** 4 / (1 + (f / hp) ** 4) / (1 + (f / lp) ** 4)
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
c_, s_ = math.cos(math.radians(350)), math.sin(math.radians(350))
truth = np.column_stack([c_ * truth[:, 0] - s_ * truth[:, 1], s_ * truth[:, 0] + c_ * truth[:, 1], truth[:, 2]])
inw = (grid >= t0u) & (grid <= t1u)
env = np.hypot(truth[:, 0], truth[:, 1]); k2 = int(fs * 10)
sm = np.convolve(env ** 2, np.ones(2 * k2 + 1) / (2 * k2 + 1), mode="same"); tpk = grid[np.argmax(np.where(inw, sm, 0))]
strong = (grid >= tpk - 25) & (grid <= tpk + 35) & inw
quiet = (grid >= t0u + 3) & (grid <= tpk - 30); tail = (grid >= tpk + 90) & (grid <= t1u)
pgd_t = (np.max(np.hypot(truth[strong, 0], truth[strong, 1])), np.max(np.sqrt((truth[strong] ** 2).sum(axis=1))))

def rs(x):
    x = x[np.isfinite(x)]; return 1.4826 * np.median(np.abs(x - np.median(x)))
def analyse(t, v):
    x = to_grid(t, v); xs = np.where(np.isfinite(x), x, 0.0); cum = np.cumsum(xs, axis=0)
    d = np.column_stack([_bp(cum[:, c], fs, HP, LP) for c in range(3)])
    corr = np.mean([np.corrcoef(d[strong, c], truth[strong, c])[0, 1] for c in range(3)])
    rmse = np.mean([1000 * math.sqrt(np.mean((d[strong, c] - truth[strong, c]) ** 2)) for c in range(3)])
    noise = np.mean([1000 * math.sqrt(np.mean(d[quiet, c] ** 2)) for c in range(3)])
    inc = x[tail]; inc = inc[np.isfinite(inc[:, 0])]
    sig = np.mean([1000 * rs(inc[:, c]) for c in range(3)])
    ac = np.mean([np.dot(inc[:-1, c] - inc[:, c].mean(), inc[1:, c] - inc[:, c].mean()) / np.dot(inc[:, c] - inc[:, c].mean(), inc[:, c] - inc[:, c].mean()) for c in range(3)])
    s10 = []
    for c in range(3):
        n = (len(inc) // 100) * 100; sums = inc[:n, c].reshape(-1, 100).sum(axis=1) if n else np.array([0.0]); s10.append(1000 * rs(sums) if len(sums) > 3 else float("nan"))
    ph = np.max(np.hypot(d[strong, 0], d[strong, 1])); p3 = np.max(np.sqrt((d[strong] ** 2).sum(axis=1)))
    return corr, rmse, noise, sig, float(ac), np.mean(s10), 100 * ph / pgd_t[0] - 100, 100 * p3 / pgd_t[1] - 100

names = ["PPPHVA DFV", "PPPHVA SFV", "GaiaVario L1", "GaiaVario mean"]
agg = {n: [] for n in names}
for rid in ids:
    print(f"\n### {rid}\n")
    print("| series | corr | RMSE mm | noise floor mm | per-epoch sigma mm | sigma of 10-s sums mm | lag-1 autocorr | PGD H vs table | PGD 3-D vs table |")
    print("|---|---|---|---|---|---|---|---|---|")
    srcs = {"PPPHVA DFV": read_va(f"{dfv}/{rid}_va.txt"), "PPPHVA SFV": read_va(f"{sfv}/{rid}_va.txt"),
            "GaiaVario L1": read_vario(f"{dfv}/{rid}_l1.csv"), "GaiaVario mean": read_vario(f"{dfv}/{rid}_mean.csv")}
    for n in names:
        r = analyse(*srcs[n]); agg[n].append(r)
        print(f"| {n} | {r[0]:.3f} | {r[1]:.2f} | {r[2]:.2f} | {r[3]:.2f} | {r[5]:.2f} | {r[4]:+.2f} | {r[6]:+.0f} % | {r[7]:+.0f} % |")
print("\n### average over receivers\n")
print("| series | corr | RMSE mm | noise floor mm | per-epoch sigma mm | sigma of 10-s sums mm | lag-1 autocorr | PGD H vs table | PGD 3-D vs table |")
print("|---|---|---|---|---|---|---|---|---|")
for n in names:
    a = np.mean(agg[n], axis=0)
    print(f"| {n} | {a[0]:.3f} | {a[1]:.2f} | {a[2]:.2f} | {a[3]:.2f} | {a[5]:.2f} | {a[4]:+.2f} | {a[6]:+.0f} % | {a[7]:+.0f} % |")
