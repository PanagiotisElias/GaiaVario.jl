#!/usr/bin/env python3
"""Low-frequency behaviour of per-epoch solutions: integrate dN, dE, dU over
the run and report the end value, the largest excursion of the integrated
series after removing a linear trend, and the sigma of 60-s mean velocities.
Compares several GaiaVario CSVs and, optionally, a PPPHVA log (needs the
station from the first CSV).

    python3 drift_table.py [--ppp ppphva.txt] a.csv b.csv ...
"""
import sys, math
import numpy as np

args = sys.argv[1:]
ppp_path = None
if args and args[0] == "--ppp":
    ppp_path = args[1]; args = args[2:]

def load_gv(path):
    station = None; sod = []; neu = []
    with open(path) as f:
        for line in f:
            if line.startswith("# station ECEF"):
                station = tuple(float(v) for v in line.split()[3:6]); continue
            if line.startswith("#") or line.startswith("sod,"): continue
            p = line.rstrip("\n").split(",")
            sod.append(float(p[0])); neu.append((float(p[2]), float(p[3]), float(p[4])))
    return station, np.array(sod), np.array(neu)

def load_ppp(path, station):
    x, y, z = station
    lon = math.atan2(y, x); r = math.hypot(x, y); e2 = 0.00669437999014
    lat = math.atan2(z, r * (1 - e2))
    for _ in range(5):
        Nn = 6378137.0 / math.sqrt(1 - e2 * math.sin(lat) ** 2); h = r / math.cos(lat) - Nn
        lat = math.atan2(z, r * (1 - e2 * Nn / (Nn + h)))
    sp, cp, sl, cl = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])
    sod = []; neu = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.split()
            sod.append(float(p[0])); neu.append(R @ np.array([float(p[2]), float(p[3]), float(p[4])]))
    return np.array(sod), np.array(neu)

def analyse(sod, neu, t0, t1):
    m = (sod >= t0) & (sod <= t1)
    sod = sod[m]; neu = neu[m]
    out = []
    for c in range(3):
        v = neu[:, c].copy()
        v[~np.isfinite(v)] = 0.0                     # a blank epoch adds no displacement
        x = np.cumsum(v)
        t = sod - sod[0]
        A = np.vstack([t, np.ones_like(t)]).T
        coef, *_ = np.linalg.lstsq(A, x, rcond=None)
        resid = x - A @ coef
        # 60-s mean velocity (mm/s) from blocks of 600 epochs
        nb = len(v) // 600
        blocks = v[:nb * 600].reshape(nb, 600).sum(axis=1) / 60.0 if nb > 1 else np.array([np.nan])
        out.append((x[-1], coef[0] * 60, np.max(np.abs(resid)), np.std(blocks)))
    return out, len(sod)

runs = []
station = None
for path in args:
    st, sod, neu = load_gv(path)
    station = station or st
    runs.append((path.split("/")[-1].rsplit(".", 1)[0].split("_", 1)[-1], sod, neu))
if ppp_path:
    sod, neu = load_ppp(ppp_path, station)
    runs.insert(0, ("PPPHVA", sod, neu))
t0 = max(r[1][0] for r in runs); t1 = min(r[1][-1] for r in runs)
print(f"window SOD {t0:.1f} .. {t1:.1f} s ({t1 - t0:.0f} s)")
print()
print("| series | epochs | integrated N / E / U at the end [mm] | trend N / E / U [mm/min] | max detrended excursion N / E / U [mm] | σ of 60-s mean velocity N / E / U [mm/s] |")
print("|---|---|---|---|---|---|")
for name, sod, neu in runs:
    res, n = analyse(sod, neu, t0, t1)
    f = lambda i, k=1000: " / ".join(f"{k*r[i]:+.1f}" if i < 3 else f"{k*r[i]:.3f}" for r in res)
    print(f"| {name} | {n} | {f(0)} | {f(1)} | " + " / ".join(f"{1000*r[2]:.1f}" for r in res) + " | " + " / ".join(f"{1000*r[3]:.3f}" for r in res) + " |")
