#!/usr/bin/env python3
"""Summarise replay variants against one PPPHVA log: robust sigma per
component, ratio GaiaVario/PPPHVA, sigma of the difference, satellites used,
rejection and slip rates. Markdown table on stdout.

    python3 variant_table.py ppphva.txt prefix_base.csv prefix_sin.csv ...
"""
import sys, os, math
import numpy as np

ppp_path = sys.argv[1]; files = sys.argv[2:]

def load_gv(path):
    station = None; d = {}
    with open(path) as f:
        for line in f:
            if line.startswith("# station ECEF"):
                station = tuple(float(v) for v in line.split()[3:6]); continue
            if line.startswith("#") or line.startswith("sod,"): continue
            p = line.rstrip("\n").split(",")
            d[round(float(p[0]), 1)] = (float(p[2]), float(p[3]), float(p[4]), int(p[6]), p[15], p[16])
    return station, d

def load_ppp(path, station):
    x, y, z = station
    lon = math.atan2(y, x); r = math.hypot(x, y); e2 = 0.00669437999014
    lat = math.atan2(z, r * (1 - e2))
    for _ in range(5):
        Nn = 6378137.0 / math.sqrt(1 - e2 * math.sin(lat) ** 2); h = r / math.cos(lat) - Nn
        lat = math.atan2(z, r * (1 - e2 * Nn / (Nn + h)))
    sp, cp, sl, cl = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])
    d = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.split()
            d[round(float(p[0]), 1)] = tuple(R @ np.array([float(p[2]), float(p[3]), float(p[4])])) + (int(p[1]),)
    return d

def rs(v):
    v = np.asarray(v); v = v[np.isfinite(v)]
    if len(v) < 10: return float("nan")
    m = np.median(v); return 1.4826 * np.median(np.abs(v - m))

station, first = load_gv(files[0])
ppp = load_ppp(ppp_path, station)
rows = []
for path in files:
    _, gv = load_gv(path)
    common = [s for s in sorted(set(gv) & set(ppp)) if np.isfinite(gv[s][0]) and np.isfinite(ppp[s][0])]
    name = os.path.basename(path).rsplit(".", 1)[0].split("_", 1)[-1]
    sg = [rs([gv[s][c] for s in common]) for c in range(3)]
    spp = [rs([ppp[s][c] for s in common]) for c in range(3)]
    sd = [rs([gv[s][c] - ppp[s][c] for s in common]) for c in range(3)]
    nused = np.median([gv[s][3] for s in common])
    rej = 100 * sum(1 for s in common if gv[s][4]) / len(common)
    slp = 100 * sum(1 for s in common if gv[s][5]) / len(common)
    rows.append((name, len(common), sg, spp, sd, nused, rej, slp))

print("PPPHVA σ N/E/U over the common epochs: " + " / ".join(f"{1000*v:.2f}" for v in rows[0][3]) + " mm")
print()
print("| variant | epochs | GaiaVario σ N / E / U [mm] | ratio to PPPHVA N / E / U | σ(diff) N / E / U [mm] | sats used | rejected % | slips % |")
print("|---|---|---|---|---|---|---|---|")
for name, n, sg, spp, sd, nused, rej, slp in rows:
    print(f"| {name} | {n} | " + " / ".join(f"{1000*v:.2f}" for v in sg) + " | " +
          " / ".join(f"{a/b:.2f}" for a, b in zip(sg, spp)) + " | " +
          " / ".join(f"{1000*v:.2f}" for v in sd) + f" | {nused:.0f} | {rej:.1f} | {slp:.1f} |")
