#!/usr/bin/env python3
"""Waveform comparison of PPPHVA and GaiaVario against the table accelerometers.

    python3 zenodo_waveforms.py PREFIX ACCEL.mseed LOCATION LAT LON T0_UTC_SOD T1_UTC_SOD OUT_PREFIX

PREFIX_va.txt          PPPHVA bridge log (SOD GPS scale, nsat, dX dY dZ ...)
PREFIX_<variant>.csv   GaiaVario solutions (per-epoch NEU increments, m): vario, single, ...
ACCEL.mseed            float64 miniSEED, three components (E, N, Z), m/s^2

Per-epoch increments are summed to displacement. All series are band-passed with
the same zero-phase filter (HP_HZ, LP_HZ; default 0.1-2 Hz); the accelerometer is
integrated twice with the same filter. Prints Markdown tables; writes a PNG and a
JSON. Environment: HP_HZ, LP_HZ, SHOW (comma list of series to draw).
"""
import sys, os, math, json, datetime, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mseed_float64 import read_mseed

prefix, accel_path, loc, lat, lon, t0u, t1u, outp = sys.argv[1:9]
lat = float(lat); lon = float(lon); t0u = float(t0u); t1u = float(t1u)
LEAP = 18.0
HP = float(os.environ.get("HP_HZ", "0.1"))
LP = float(os.environ.get("LP_HZ", "2.0"))
MARGIN = 15.0                                        # processed beyond the window to keep filter edges out of it

sp, cp, sl, cl = math.sin(math.radians(lat)), math.cos(math.radians(lat)), math.sin(math.radians(lon)), math.cos(math.radians(lon))
R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])      # ECEF -> NEU

def read_va(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or not line.strip(): continue
        p = line.split()
        d = np.array([float(p[2]), float(p[3]), float(p[4])])
        t.append(float(p[0]) - LEAP)                  # UTC second of day
        v.append(R @ d if np.all(np.isfinite(d)) else np.full(3, np.nan))
    return np.array(t), np.array(v)

def read_vario(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or line.startswith("sod,") or not line.strip(): continue
        p = line.rstrip("\n").split(",")
        t.append(float(p[0]) - LEAP)
        v.append([float(p[2]), float(p[3]), float(p[4])])
    return np.array(t), np.array(v)

LABEL = {"vario": "GaiaVario", "single": "GaiaVario L1", "dualW": "GaiaVario L2 P-code",
         "dualL5": "GaiaVario GPS L5", "dualGAL7": "GaiaVario Galileo E5b"}
series = {}
series["PPPHVA"] = read_va(prefix + "_va.txt")
for fn in sorted(glob.glob(prefix + "_*.csv")):
    key = os.path.basename(fn)[len(os.path.basename(prefix)) + 1:-4]
    series[LABEL.get(key, "GaiaVario " + key)] = read_vario(fn)
SHOW = os.environ.get("SHOW", "PPPHVA,GaiaVario,GaiaVario L1").split(",")

fs = 10.0
g0, g1 = t0u - MARGIN, t1u + MARGIN
grid = np.arange(math.ceil(g0 * fs), math.floor(g1 * fs) + 1) / fs
def to_grid(t, v):
    out = np.full((len(grid), 3), np.nan)
    idx = np.round((t - grid[0]) * fs).astype(int)
    ok = (idx >= 0) & (idx < len(grid))
    out[idx[ok]] = v[ok]
    return out
inc = {k: to_grid(*s) for k, s in series.items()}

def _bp(x, fs, hp, lp, order=0):
    n = len(x); N = 1 << (n - 1).bit_length(); N *= 2
    X = np.fft.rfft(x - np.mean(x[:int(fs * 5)]), N); f = np.fft.rfftfreq(N, 1 / fs)
    H = np.ones_like(f)
    if hp > 0: H *= (f / hp) ** 4 / (1 + (f / hp) ** 4)
    if lp > 0: H *= 1 / (1 + (f / lp) ** 4)
    if order:
        with np.errstate(divide="ignore", invalid="ignore"):
            W = (1j * 2 * np.pi * f) ** (-order); W[0] = 0
        H = H * W
    return np.fft.irfft(X * H, N)[:n]

disp = {}
inw = (grid >= t0u) & (grid <= t1u)
for k, x in inc.items():
    ok = np.isfinite(x[:, 0]); xs = np.where(np.isfinite(x), x, 0.0)
    cum = np.cumsum(xs, axis=0)
    disp[k] = np.column_stack([_bp(cum[:, c], fs, HP, LP) for c in range(3)])
    print(f"{k}: epochs in window {int(ok[inw].sum())} of {int(inw.sum())} ({100*ok[inw].mean():.1f} %)")

# ---- accelerometer: E, N, Z in m/s^2 -> displacement on the same grid
tr = read_mseed(accel_path)
chan = {c: k for k in tr for c in ("E", "N", "Z") if k.endswith(c)}
base = (datetime.datetime(2025, 10, 21) - datetime.datetime(1970, 1, 1)).total_seconds()
def accel_disp(name):
    t0, rate, x, _ = tr[name]
    ts = t0 - base + np.arange(len(x)) / rate
    m = (ts >= g0 - 60) & (ts <= g1 + 60)
    return np.interp(grid, ts[m], _bp(x[m], rate, HP, LP, 2))
d_neu = np.column_stack([accel_disp(chan["N"]), accel_disp(chan["E"]), accel_disp(chan["Z"])])

# strong-motion window: 60 s around the peak of the smoothed horizontal truth envelope
env = np.sqrt(d_neu[:, 0] ** 2 + d_neu[:, 1] ** 2)
k2 = int(fs * 10); sm = np.convolve(env ** 2, np.ones(2 * k2 + 1) / (2 * k2 + 1), mode="same")
tpk = grid[np.argmax(np.where(inw, sm, 0))]
strong = (grid >= tpk - 25) & (grid <= tpk + 35) & inw

# ---- unknown table orientation / sign: horizontal rotation maximising the total correlation of ALL series
def best_rotation(ests, tru):
    best = None
    for ang in np.arange(0, 360, 1.0):
        c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        rot = np.column_stack([c * tru[:, 0] - s * tru[:, 1], s * tru[:, 0] + c * tru[:, 1]])
        sc = sum(sum(np.dot(e[:, i], rot[:, i]) for i in range(2)) / math.sqrt(np.sum(e[:, :2] ** 2) * np.sum(rot ** 2) + 1e-30) for e in ests)
        if best is None or sc > best[0]: best = (sc / len(ests), ang)
    return best
sc, ang = best_rotation([disp[k][inw] for k in series], d_neu[inw])
c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
d_neu[:, 0], d_neu[:, 1] = c * d_neu[:, 0] - s * d_neu[:, 1], s * d_neu[:, 0] + c * d_neu[:, 1]
print(f"table axes rotated by {ang:.0f} deg relative to the accelerometer E/N (mean correlation over the series {sc:.3f}); strong motion around UTC {tpk:.0f} s")

def lag_of(x, y, maxlag=3.0):
    best = (-2, 0)
    for kk in range(-int(maxlag * fs), int(maxlag * fs) + 1):
        xs, ys = (x[kk:], y[:len(y) - kk]) if kk >= 0 else (x[:kk], y[-kk:])
        cc = np.corrcoef(xs[:, :2].ravel(), ys[:, :2].ravel())[0, 1]
        if cc > best[0]: best = (cc, kk)
    return best[1] / fs, best[0]
lags = {k: lag_of(disp[k][inw], d_neu[inw]) for k in series}
print("best lag GNSS vs accelerometer [s] (corr): " + ", ".join(f"{k} {v[0]:+.1f} ({v[1]:.3f})" for k, v in lags.items()))

# ---- statistics (displacement) over the whole event window and the 60-s strong-motion window
def stat_rows(mask):
    rows = {}
    for k in series:
        for ci, cn in enumerate("NEU"):
            e = disp[k][mask, ci] * 1000; t = d_neu[mask, ci] * 1000
            ok = np.isfinite(e) & np.isfinite(t); e, t = e[ok], t[ok]
            rows[(k, cn)] = dict(corr=float(np.corrcoef(e, t)[0, 1]), rmse=float(math.sqrt(np.mean((e - t) ** 2))),
                                 rms_truth=float(math.sqrt(np.mean(t ** 2))), peak_est=float(np.max(np.abs(e))), peak_truth=float(np.max(np.abs(t))))
    return rows
stats_all = stat_rows(inw); stats_strong = stat_rows(strong)
print()
print("| series | comp | corr (event) | RMSE (event) mm | corr (strong 60 s) | RMSE (strong) mm | RMS truth (strong) mm | peak est / truth (strong) mm |")
print("|---|---|---|---|---|---|---|---|")
for k in series:
    for cn in "NEU":
        a, b = stats_all[(k, cn)], stats_strong[(k, cn)]
        print(f"| {k} | {cn} | {a['corr']:.3f} | {a['rmse']:.2f} | {b['corr']:.3f} | {b['rmse']:.2f} | {b['rms_truth']:.2f} | {b['peak_est']:.1f} / {b['peak_truth']:.1f} |")
dd = (disp["GaiaVario"] - disp["PPPHVA"])[inw] * 1000
print("\nGaiaVario - PPPHVA displacement rms N/E/U = " + " / ".join(f"{np.sqrt(np.nanmean(dd[:, i] ** 2)):.2f}" for i in range(3)) + " mm (event window)")
quiet = (grid >= t0u + 3) & (grid <= tpk - 30)
if quiet.sum() > 100:
    print("noise-only window before the strong motion, displacement rms N/E/U [mm]: " + "; ".join(
        f"{k}: " + " / ".join(f"{1000*np.sqrt(np.nanmean(disp[k][quiet, i]**2)):.2f}" for i in range(3)) for k in series) +
        f";  truth: " + " / ".join(f"{1000*np.sqrt(np.nanmean(d_neu[quiet, i]**2)):.2f}" for i in range(3)))

json.dump(dict(grid=grid.tolist(), truth=(d_neu * 1000).round(3).tolist(), strong=[float(tpk - 25), float(tpk + 35)],
               series={k: (disp[k] * 1000).round(3).tolist() for k in series}, rot_deg=float(ang), hp=HP, lp=LP,
               stats_all={f"{k}|{c}": v for (k, c), v in stats_all.items()}, stats_strong={f"{k}|{c}": v for (k, c), v in stats_strong.items()}),
          open(outp + "_waveforms.json", "w"))

# ---- figure
from PIL import Image, ImageDraw, ImageFont
W, H = 1500, 1330
img = Image.new("RGB", (W, H), "white"); dr = ImageDraw.Draw(img)
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 15); small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
except Exception:
    font = small = ImageFont.load_default()
cols = {"truth": "#222222", "PPPHVA": "#D85A30", "GaiaVario": "#1D9E75", "GaiaVario L1": "#534AB7",
        "GaiaVario L2 P-code": "#BA7517", "GaiaVario GPS L5": "#D4537E", "GaiaVario Galileo E5b": "#378ADD"}
dr.text((90, 8), f"Samos 2020 on the Timisoara shaking table (station {loc}), displacement: table accelerometer (black) vs " + ", ".join(SHOW) + f"; band {HP}-{LP} Hz", fill="black", font=font)
left, right, top, ph, gap = 90, W - 30, 40, 190, 16
panels = [(f"{'NEU'[ci]} displacement [mm], whole event", ci, inw) for ci in range(3)] + \
         [(f"{'NEU'[ci]} displacement [mm], 60 s of strong motion", ci, strong) for ci in range(3)]
for pi, (title, ci, mask) in enumerate(panels):
    y0 = top + pi * (ph + gap)
    dr.rectangle([left, y0, right, y0 + ph], outline="#888888")
    tru = d_neu[:, ci] * 1000
    lim = max(1.0, float(np.percentile(np.abs(tru[mask]), 99.5)) * 1.35)
    ym = y0 + ph / 2
    dr.line([left, ym, right, ym], fill="#dddddd")
    tt0, tt1 = grid[mask][0], grid[mask][-1]
    def draw(y, col, w=1):
        pts = [(left + (right - left) * (tt - tt0) / (tt1 - tt0), ym - (ph / 2) * max(-1.05, min(1.05, vv / lim)))
               for tt, vv in zip(grid[mask], y[mask]) if np.isfinite(vv)]
        if len(pts) > 1: dr.line(pts, fill=col, width=w)
    for k in SHOW:
        draw(disp[k][:, ci] * 1000, cols[k], 1)
    draw(tru, cols["truth"], 2)
    dr.text((8, y0 + 3), title, fill="black", font=small); dr.text((8, y0 + ph - 16), f"±{lim:.0f}", fill="#555555", font=small)
    st = stats_all if mask is inw else stats_strong
    for j, k in enumerate(SHOW):
        r = st[(k, 'NEU'[ci])]
        dr.text((right - 470 + 160 * j, y0 + 3), f"{'VA' if k == 'PPPHVA' else k.replace('GaiaVario', 'Vario')} r={r['corr']:.2f}", fill=cols[k], font=small)
dr.text((left, H - 22), f"UTC seconds of day {grid[inw][0]:.0f} to {grid[inw][-1]:.0f} (strong motion {tpk-25:.0f} to {tpk+35:.0f}); table axes rotated {ang:.0f} deg; GNSS lag {lags['PPPHVA'][0]:+.1f} s; r = correlation with the accelerometer", fill="#555555", font=small)
img.save(outp + "_waveforms.png")
print("figure:", outp + "_waveforms.png")
