#!/usr/bin/env python3
"""Align a GaiaSeis (PPPHVA) per-epoch log with a GaiaVario comparison CSV and
print the statistics of both engines plus their difference; draw a PNG.

    python3 compare_report.py logs/UPC10Z_cmp_<stamp>.txt logs/gaiavario_cmp_<stamp>.csv [out.png]

Needs numpy; the figure is drawn with Pillow (no matplotlib on this box).
"""
import sys, math, re
import numpy as np

ppp_path, gv_path = sys.argv[1], sys.argv[2]
png = sys.argv[3] if len(sys.argv) > 3 else None

# ---- GaiaVario CSV
station = None
gv = {}
# second of day unwrapped across a GPS midnight: keys keep increasing through a run
def unwrapper():
    st = {"prev": None, "off": 0.0}
    def f(s):
        if st["prev"] is not None and s < st["prev"] - 43200: st["off"] += 86400.0
        st["prev"] = s
        return round(s + st["off"], 1)
    return f
uw = unwrapper()
with open(gv_path) as f:
    for line in f:
        if line.startswith("# station ECEF"):
            station = tuple(float(v) for v in line.split()[3:6])
            continue
        if line.startswith("#") or line.startswith("sod,"):
            continue
        p = line.rstrip("\n").split(",")
        sod = uw(float(p[0]))
        gv[sod] = dict(N=float(p[2]), E=float(p[3]), U=float(p[4]), nin=int(p[5]), nused=int(p[6]),
                       nsys=int(p[7]), pdop=float(p[8]), rms=float(p[9]), sig=(float(p[10]), float(p[11]), float(p[12])),
                       ms=float(p[13]), status=p[14], rej=p[15], slips=p[16])
if station is None:
    sys.exit("no station position in the GaiaVario CSV")

# ---- PPPHVA log: SOD nsat dX dY dZ ... (ECEF) -> NEU
x, y, z = station
lon = math.atan2(y, x); r = math.hypot(x, y); lat = math.atan2(z, r * (1 - 0.00669437999014))
for _ in range(5):
    Nn = 6378137.0 / math.sqrt(1 - 0.00669437999014 * math.sin(lat) ** 2)
    h = r / math.cos(lat) - Nn
    lat = math.atan2(z, r * (1 - 0.00669437999014 * Nn / (Nn + h)))
sφ, cφ, sλ, cλ = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
R = np.array([[-sφ * cλ, -sφ * sλ, cφ], [-sλ, cλ, 0.0], [cφ * cλ, cφ * sλ, sφ]])
ppp = {}
uw = unwrapper()
with open(ppp_path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        sod = uw(float(p[0]))
        d = np.array([float(p[2]), float(p[3]), float(p[4])])
        neu = R @ d
        ppp[sod] = dict(N=neu[0], E=neu[1], U=neu[2], nsat=int(p[1]))

def robust_sigma(v):
    v = np.asarray(v); v = v[np.isfinite(v)]
    if len(v) < 10:
        return float("nan")
    m = np.median(v)
    return 1.4826 * np.median(np.abs(v - m))

def stats(name, series, key_n="nsat"):
    sods = sorted(series)
    out = {}
    for c in "NEU":
        v = np.array([series[s][c] for s in sods])
        fin = v[np.isfinite(v)]
        out[c] = dict(n=len(fin), sigma_robust=robust_sigma(fin), std=float(np.std(fin)) if len(fin) else float("nan"),
                      maxabs=float(np.max(np.abs(fin))) if len(fin) else float("nan"),
                      n_gt_1cm=int(np.sum(np.abs(fin) > 0.01)), n_gt_5cm=int(np.sum(np.abs(fin) > 0.05)))
    out["epochs"] = len(sods)
    out["blank"] = sum(1 for s in sods if not np.isfinite(series[s]["N"]))
    return out

sp = stats("PPPHVA", ppp)
sg = stats("GaiaVario", gv)
common = sorted(set(ppp) & set(gv))
both = [s for s in common if np.isfinite(ppp[s]["N"]) and np.isfinite(gv[s]["N"])]

def fmt_mm(v): return "nan" if not np.isfinite(v) else f"{1000*v:.2f}"

print(f"PPPHVA log: {ppp_path}")
print(f"GaiaVario : {gv_path}")
print(f"station ECEF {station}")
t0, t1 = (min(common), max(common)) if common else (float('nan'), float('nan'))
print(f"common epochs {len(common)} (SOD {t0:.1f} .. {t1:.1f}), both solved {len(both)}")
print(f"PPPHVA   epochs {sp['epochs']} blank {sp['blank']}")
print(f"GaiaVario epochs {sg['epochs']} blank {sg['blank']}")
print()
print("| component | engine | n | robust σ [mm] | std [mm] | max |d| [mm] | >1 cm | >5 cm |")
print("|---|---|---|---|---|---|---|---|")
for c in "NEU":
    for name, st in (("PPPHVA", sp), ("GaiaVario", sg)):
        d = st[c]
        print(f"| {c} | {name} | {d['n']} | {fmt_mm(d['sigma_robust'])} | {fmt_mm(d['std'])} | {fmt_mm(d['maxabs'])} | {d['n_gt_1cm']} | {d['n_gt_5cm']} |")
print()
print("| component | difference GaiaVario − PPPHVA: robust σ [mm] | std [mm] | mean [mm] | max |Δ| [mm] | correlation |")
print("|---|---|---|---|---|---|")
for c in "NEU":
    a = np.array([ppp[s][c] for s in both]); b = np.array([gv[s][c] for s in both])
    d = b - a
    corr = float(np.corrcoef(a, b)[0, 1]) if len(both) > 2 else float("nan")
    print(f"| {c} | {fmt_mm(robust_sigma(d))} | {fmt_mm(float(np.std(d)))} | {fmt_mm(float(np.mean(d)))} | {fmt_mm(float(np.max(np.abs(d))))} | {corr:.3f} |")
print()
nin = np.array([gv[s]["nin"] for s in gv]); nused = np.array([gv[s]["nused"] for s in gv])
nsys = np.array([gv[s]["nsys"] for s in gv]); pdop = np.array([gv[s]["pdop"] for s in gv])
ms = np.array([gv[s]["ms"] for s in gv]); rms = np.array([gv[s]["rms"] for s in gv])
nsat_p = np.array([ppp[s]["nsat"] for s in ppp])
print(f"GaiaVario: satellites in/used median {np.median(nin):.0f}/{np.median(nused):.0f}, systems {np.median(nsys):.0f}, "
      f"PDOP median {np.nanmedian(pdop):.2f}, fit RMS median {1000*np.nanmedian(rms):.1f} mm, step! median {np.median(ms):.2f} ms p95 {np.percentile(ms,95):.2f} ms")
print(f"PPPHVA   : satellites used median {np.median(nsat_p):.0f}")
rej = sum(1 for s in gv if gv[s]["rej"]); slp = sum(1 for s in gv if gv[s]["slips"])
print(f"GaiaVario: epochs with a rejected satellite {rej}, with a slip-flagged satellite {slp}")
sig = np.array([gv[s]["sig"] for s in gv if np.isfinite(gv[s]["N"])])
if len(sig):
    print(f"GaiaVario formal σ median N/E/U {1000*np.median(sig[:,0]):.2f}/{1000*np.median(sig[:,1]):.2f}/{1000*np.median(sig[:,2]):.2f} mm")
# largest differences
dd = sorted(((abs(gv[s]["N"] - ppp[s]["N"]) + abs(gv[s]["E"] - ppp[s]["E"]) + abs(gv[s]["U"] - ppp[s]["U"]), s) for s in both), reverse=True)[:8]
print("\nlargest |ΔN|+|ΔE|+|ΔU| epochs (SOD, PPPHVA NEU mm, GaiaVario NEU mm, GaiaVario rejected/slips):")
for _, s in dd:
    print(f"  {s:9.1f}  P {1000*ppp[s]['N']:+7.1f} {1000*ppp[s]['E']:+7.1f} {1000*ppp[s]['U']:+7.1f}   G {1000*gv[s]['N']:+7.1f} {1000*gv[s]['E']:+7.1f} {1000*gv[s]['U']:+7.1f}   {gv[s]['rej']} / {gv[s]['slips']}")

# ---- figure (Pillow)
if png and both:
    from PIL import Image, ImageDraw, ImageFont
    W, H = 1400, 1180
    img = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16)
        small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 13)
    except Exception:
        font = small = ImageFont.load_default()
    sods = np.array(common)
    tt = sods - sods[0]
    panels = [("N", 0), ("E", 1), ("U", 2), ("ΔN", 3), ("ΔE", 4), ("ΔU", 5)]
    left, right, top, ph, gap = 90, W - 30, 50, 160, 22
    dr.text((left, 12), f"UPC10Z live, {len(common)} epochs at 10 Hz — PPPHVA (red) and GaiaVario (teal), per-epoch displacement [mm]; Δ = GaiaVario − PPPHVA (purple)", fill="black", font=font)
    for name, i in panels:
        y0 = top + i * (ph + gap)
        dr.rectangle([left, y0, right, y0 + ph], outline="#888888")
        if i < 3:
            a = np.array([ppp[s][name] if np.isfinite(ppp[s][name]) else np.nan for s in common]) * 1000
            b = np.array([gv[s][name] if np.isfinite(gv[s][name]) else np.nan for s in common]) * 1000
            lim = max(5.0, float(np.nanpercentile(np.abs(np.concatenate([a, b])), 99.5)) * 1.3)
            series = [(a, "#D85A30"), (b, "#1D9E75")]
        else:
            c = name[1]
            d = np.array([(gv[s][c] - ppp[s][c]) if (np.isfinite(gv[s][c]) and np.isfinite(ppp[s][c])) else np.nan for s in common]) * 1000
            series = [(d, "#534AB7")]
            lim = max(2.0, float(np.nanpercentile(np.abs(d), 99.5)) * 1.3)
        ym = y0 + ph / 2
        dr.line([left, ym, right, ym], fill="#cccccc")
        for v, col in series:
            pts = []
            for k in range(len(tt)):
                if np.isfinite(v[k]):
                    px = left + (right - left) * tt[k] / max(tt[-1], 1)
                    py = ym - (ph / 2) * max(-1, min(1, v[k] / lim))
                    pts.append((px, py))
            if len(pts) > 1:
                dr.line(pts, fill=col, width=1)
        dr.text((8, y0 + 2), name, fill="black", font=font)
        dr.text((8, y0 + ph - 18), f"±{lim:.0f} mm", fill="#555555", font=small)
    dr.text((left, H - 25), f"SOD {sods[0]:.1f} to {sods[-1]:.1f} s; blanks are gaps in the lines", fill="#555555", font=small)
    img.save(png)
    print("figure:", png)
