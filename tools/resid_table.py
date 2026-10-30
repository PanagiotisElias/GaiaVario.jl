#!/usr/bin/env python3
"""Post-fit residuals of several bridge runs on their common epochs. Each
argument is label=bridge_log.txt (the per-epoch log of GaiaSeis, columns
SOD nsat dX dY dZ DTG DTR DTE nused rms_w rms_u since 2026-10-30).

    python3 resid_table.py "PPPHVA DFV=va_dfv.txt" "GaiaVario IF=gv_if.txt" ...
"""
import sys, numpy as np
runs = []
for arg in sys.argv[1:]:
    lab, path = arg.split("=", 1)
    d = {}
    for line in open(path):
        if line.startswith("#") or not line.strip(): continue
        p = line.split()
        if len(p) < 11 or p[2] == "NaN": continue
        d[round(float(p[0]), 1)] = (int(p[8]), float(p[9]), float(p[10]))
    runs.append((lab, d))
common = sorted(set.intersection(*[set(d) for _, d in runs]))
print(f"common epochs with a solution in every run: {len(common)}")
print()
print("| run | satellites kept (median) | weighted post-fit RMS [mm], median (p90) | unweighted post-fit RMS [mm], median (p90) |")
print("|---|---|---|---|")
for lab, d in runs:
    A = np.array([d[s] for s in common])
    f = lambda c: f"{1e3*np.nanmedian(A[:,c]):.2f} ({1e3*np.nanpercentile(A[:,c],90):.2f})"
    print(f"| {lab} | {np.median(A[:,0]):.0f} | {f(1)} | {f(2)} |")
