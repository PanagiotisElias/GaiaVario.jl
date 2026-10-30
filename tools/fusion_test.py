#!/usr/bin/env python3
"""Complementary fusion of two GaiaVario displacement series: the low frequencies of
the ionosphere-free series plus the high frequencies of a quieter but
ionosphere-sensitive series (L1 or the mean of two), with complementary zero-phase
gains  LP = 1/(1+(f/fc)^4),  HP = (f/fc)^4/(1+(f/fc)^4)  (LP + HP = 1).

    python3 fusion_test.py ACCEL.mseed LAT LON T0 T1 OUTDIR IF_NAME ID [ID ...]

OUTDIR/<ID>_<IF_NAME>.csv, _mean.csv, _l1.csv and _va.txt. Evaluated against the table
accelerometer in several bands (HP corner of the evaluation band: 0.1, 0.05, 0.02 Hz).
"""
import sys, os, math, datetime
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mseed_float64 import read_mseed

accel_path, lat, lon, t0u, t1u, outdir, ifname = sys.argv[1:8]; ids = sys.argv[8:]
lat = float(lat); lon = float(lon); t0u = float(t0u); t1u = float(t1u)
LEAP = 18.0; fs = 10.0; MARGIN = 15.0; LP = 2.0
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
def cum_series(t, v):
    x = np.full((len(grid), 3), 0.0); idx = np.round((t - grid[0]) * fs).astype(int); ok = (idx >= 0) & (idx < len(grid))
    x[idx[ok]] = np.where(np.isfinite(v[ok]), v[ok], 0.0)
    return np.cumsum(x, axis=0)

def spec(x):
    n = len(x); N = (1 << (n - 1).bit_length()) * 2
    return np.fft.rfft(x - np.mean(x[:int(fs * 5)]), N), np.fft.rfftfreq(N, 1 / fs), n
def hp_gain(f, fc): return (f / fc) ** 4 / (1 + (f / fc) ** 4) if fc > 0 else np.ones_like(f)
def lp_gain(f, fc): return 1 / (1 + (f / fc) ** 4)
def filt(x, hp, lp, order=0, fs_=fs):
    n = len(x); N = (1 << (n - 1).bit_length()) * 2
    X = np.fft.rfft(x - np.mean(x[:int(fs_ * 5)]), N); f = np.fft.rfftfreq(N, 1 / fs_)
    H = hp_gain(f, hp) * lp_gain(f, lp)
    if order:
        with np.errstate(divide="ignore", invalid="ignore"):
            W = (1j * 2 * np.pi * f) ** (-order); W[0] = 0
        H = H * W
    return np.fft.irfft(X * H, N)[:n]
def fuse(low, high, fc):
    Xl, f, n = spec(low); Xh, _, _ = spec(high)
    return np.fft.irfft(Xl * lp_gain(f, fc) + Xh * hp_gain(f, fc), len(f) * 2 - 2)[:n]

tr = read_mseed(accel_path); chan = {c: k for k in tr for c in ("E", "N", "Z") if k.endswith(c)}
base = (datetime.datetime(2025, 10, 21) - datetime.datetime(1970, 1, 1)).total_seconds()
def accel_raw(name):
    t0, rate, x, _ = tr[name]; ts = t0 - base + np.arange(len(x)) / rate; m = (ts >= g0 - 120) & (ts <= g1 + 120)
    return ts[m], x[m], rate
acc = {c: accel_raw(chan[c]) for c in "NEZ"}
def truth_band(hp):
    out = []
    for c in "NEZ":
        ts, x, rate = acc[c]; out.append(np.interp(grid, ts, filt(x, hp, LP, 2, rate)))
    d = np.column_stack(out)
    ROT = 350.0; c_, s_ = math.cos(math.radians(ROT)), math.sin(math.radians(ROT))
    return np.column_stack([c_ * d[:, 0] - s_ * d[:, 1], s_ * d[:, 0] + c_ * d[:, 1], d[:, 2]])

inw = (grid >= t0u) & (grid <= t1u)
t_ref = truth_band(0.1); env = np.hypot(t_ref[:, 0], t_ref[:, 1]); k2 = int(fs * 10)
sm = np.convolve(env ** 2, np.ones(2 * k2 + 1) / (2 * k2 + 1), mode="same"); tpk = grid[np.argmax(np.where(inw, sm, 0))]
strong = (grid >= tpk - 25) & (grid <= tpk + 35) & inw
quiet = (grid >= t0u + 3) & (grid <= tpk - 30)

HPS = [0.1, 0.05, 0.02]; FCS = [0.03, 0.05, 0.1]
variants = ["va", ifname, "l1", "mean"] + [f"{h}+{ifname}@{fc}" for h in ("l1", "mean") for fc in FCS]
EXTRA = [v for v in os.environ.get("EXTRA", "").split(",") if v]       # additional plain series (csv names) to evaluate
if os.environ.get("ONLY_EXTRA"):
    variants = ["va", ifname, "l1", "mean"] + EXTRA
else:
    variants = variants + EXTRA
res = {hp: {v: [] for v in variants} for hp in HPS}
for rid in ids:
    series = {}
    series["va"] = cum_series(*read_va(f"{outdir}/{rid}_va.txt"))
    for nm in (ifname, "l1", "mean") + tuple(EXTRA):
        series[nm] = cum_series(*read_vario(f"{outdir}/{rid}_{nm}.csv"))
    if not os.environ.get("ONLY_EXTRA"):
        for h in ("l1", "mean"):
            for fc in FCS:
                series[f"{h}+{ifname}@{fc}"] = np.column_stack([fuse(series[ifname][:, c], series[h][:, c], fc) for c in range(3)])
    for hp in HPS:
        tr_b = truth_band(hp)
        for v in variants:
            d = np.column_stack([filt(series[v][:, c], hp, LP) for c in range(3)])
            corr = np.mean([np.corrcoef(d[strong, c], tr_b[strong, c])[0, 1] for c in range(3)])
            rmse = np.mean([1000 * math.sqrt(np.mean((d[strong, c] - tr_b[strong, c]) ** 2)) for c in range(3)])
            noise = np.mean([1000 * math.sqrt(np.mean(d[quiet, c] ** 2)) for c in range(3)])
            res[hp][v].append((corr, rmse, noise))
for hp in HPS:
    print(f"\n### evaluation band {hp}-{LP} Hz (average over {len(ids)} receivers; corr and RMSE over the 60 s of strong motion, noise floor before the event)\n")
    print("| series | corr | RMSE mm | noise floor mm |"); print("|---|---|---|---|")
    for v in variants:
        a = np.mean(res[hp][v], axis=0); print(f"| {v} | {a[0]:.3f} | {a[1]:.2f} | {a[2]:.2f} |")
