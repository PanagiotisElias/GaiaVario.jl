#!/usr/bin/env python3
"""Interactive web viewer (one self-contained HTML file) of the shaking-table replay:
several receivers x engines against the table accelerometer, with zoom, pan,
series selection, band selection and statistics of the visible window.

    python3 zenodo_viewer.py OUT.html ACCEL.mseed LAT LON T0_UTC_SOD T1_UTC_SOD  ID:LABEL:PREFIX [ID:LABEL:PREFIX ...]

For every PREFIX the files PREFIX_va.txt (PPPHVA bridge log), PREFIX_vario.csv and
PREFIX_single.csv (GaiaVario dual-frequency and L1) are read. Also prints a summary
table (correlation with the accelerometer in the strong-motion window).
"""
import sys, os, math, json, datetime
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mseed_float64 import read_mseed

out_html, accel_path, lat, lon, t0u, t1u = sys.argv[1:7]
lat = float(lat); lon = float(lon); t0u = float(t0u); t1u = float(t1u)
recv = [a.split(":", 2) for a in sys.argv[7:]]
LEAP = 18.0; fs = 10.0; MARGIN = 15.0
BANDS = {"0.1-2 Hz": (0.1, 2.0), "0.05-1 Hz": (0.05, 1.0), "0.1-4 Hz": (0.1, 4.0)}
ENGINES = [("va", "PPPHVA"), ("single", "GaiaVario L1 (default)"), ("gf30", "GaiaVario L1 + GF correction (option)"), ("vario", "GaiaVario dual-frequency (ionosphere-free)")]

sp, cp, sl, cl = math.sin(math.radians(lat)), math.cos(math.radians(lat)), math.sin(math.radians(lon)), math.cos(math.radians(lon))
R = np.array([[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]])

def read_va(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or not line.strip(): continue
        p = line.split(); d = np.array([float(p[2]), float(p[3]), float(p[4])])
        t.append(float(p[0]) - LEAP); v.append(R @ d if np.all(np.isfinite(d)) else np.full(3, np.nan))
    return np.array(t), np.array(v)

def read_vario(path):
    t, v = [], []
    for line in open(path):
        if line.startswith("#") or line.startswith("sod,") or not line.strip(): continue
        p = line.rstrip("\n").split(","); t.append(float(p[0]) - LEAP); v.append([float(p[2]), float(p[3]), float(p[4])])
    return np.array(t), np.array(v)

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

# ---- truth per band
tr = read_mseed(accel_path)
chan = {c: k for k in tr for c in ("E", "N", "Z") if k.endswith(c)}
base = (datetime.datetime(2025, 10, 21) - datetime.datetime(1970, 1, 1)).total_seconds()
def accel_disp(name, hp, lp):
    t0, rate, x, _ = tr[name]; ts = t0 - base + np.arange(len(x)) / rate
    m = (ts >= g0 - 60) & (ts <= g1 + 60)
    return np.interp(grid, ts[m], _bp(x[m], rate, hp, lp, 2))
truth = {b: np.column_stack([accel_disp(chan["N"], *bw), accel_disp(chan["E"], *bw), accel_disp(chan["Z"], *bw)]) for b, bw in BANDS.items()}

# ---- GNSS series per receiver, engine, band
inw = (grid >= t0u) & (grid <= t1u)
data = {}      # receiver -> engine -> band -> [N,E,U] lists (mm)
summary = []
raw_disp = {}
for rid, label, prefix in recv:
    data[rid] = {}
    for key, ename in ENGINES:
        path = prefix + ("_va.txt" if key == "va" else f"_{key}.csv")
        t, v = read_va(path) if key == "va" else read_vario(path)
        x = to_grid(t, v); xs = np.where(np.isfinite(x), x, 0.0); cum = np.cumsum(xs, axis=0)
        data[rid][key] = {}
        for b, (hp, lp) in BANDS.items():
            d = np.column_stack([_bp(cum[:, c], fs, hp, lp) for c in range(3)])
            data[rid][key][b] = d
# ---- table orientation: one rotation for all (fit on the first receiver's three engines in the default band)
b0 = "0.1-2 Hz"
def best_rotation(ests, tru):
    best = None
    for ang in np.arange(0, 360, 1.0):
        c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        rot = np.column_stack([c * tru[:, 0] - s * tru[:, 1], s * tru[:, 0] + c * tru[:, 1]])
        sc = sum(sum(np.dot(e[:, i], rot[:, i]) for i in range(2)) / math.sqrt(np.sum(e[:, :2] ** 2) * np.sum(rot ** 2) + 1e-30) for e in ests)
        if best is None or sc > best[0]: best = (sc, ang)
    return best[1]
# the orientation is a property of the set-up, so fit it on the series with the best signal-to-noise (L1) of all receivers
ests = [data[rid]['gf30'][b0][inw] for rid, _, _ in recv]
ang = best_rotation(ests, truth[b0][inw])
def obj(a_):
    c_, s_ = math.cos(math.radians(a_)), math.sin(math.radians(a_)); tt = truth[b0][inw]
    rot_ = np.column_stack([c_ * tt[:, 0] - s_ * tt[:, 1], s_ * tt[:, 0] + c_ * tt[:, 1]])
    return sum(sum(np.dot(e[:, i], rot_[:, i]) for i in range(2)) / math.sqrt(np.sum(e[:, :2] ** 2) * np.sum(rot_ ** 2)) for e in ests) / len(ests)
print("rotation objective at fit-10, fit, fit+10 deg: " + " / ".join(f"{obj((ang + d) % 360):.3f}" for d in (-10, 0, 10)))
c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
for b in truth:
    t_ = truth[b]; truth[b] = np.column_stack([c * t_[:, 0] - s * t_[:, 1], s * t_[:, 0] + c * t_[:, 1], t_[:, 2]])
env = np.sqrt(truth[b0][:, 0] ** 2 + truth[b0][:, 1] ** 2)
k2 = int(fs * 10); sm = np.convolve(env ** 2, np.ones(2 * k2 + 1) / (2 * k2 + 1), mode="same")
tpk = float(grid[np.argmax(np.where(inw, sm, 0))])
strong = (grid >= tpk - 25) & (grid <= tpk + 35) & inw

print(f"table rotation {ang:.0f} deg; strong motion {tpk-25:.0f}-{tpk+35:.0f} UTC s")
print("\n| receiver | engine | corr N | corr E | corr U | RMSE N / E / U [mm] (strong 60 s, band 0.1-2 Hz) |")
print("|---|---|---|---|---|---|")
for rid, label, prefix in recv:
    for key, ename in ENGINES:
        d = data[rid][key][b0]; cc = []; rm = []
        for ci in range(3):
            e = d[strong, ci] * 1000; t_ = truth[b0][strong, ci] * 1000
            cc.append(np.corrcoef(e, t_)[0, 1]); rm.append(math.sqrt(np.mean((e - t_) ** 2)))
        print(f"| {label} | {ename} | {cc[0]:.2f} | {cc[1]:.2f} | {cc[2]:.2f} | {rm[0]:.1f} / {rm[1]:.1f} / {rm[2]:.1f} |")

# ---- export (mm, one decimal) for the page
def pack(a): return [[round(float(x) * 1000, 1) for x in a[:, ci]] for ci in range(3)]
payload = dict(t0=float(grid[0]), fs=fs, event=[t0u, t1u], strong=[tpk - 25, tpk + 35], rot=float(ang),
               bands=list(BANDS), engines=[dict(key=k, name=n) for k, n in ENGINES],
               receivers=[dict(id=rid, label=label) for rid, label, _ in recv],
               truth={b: pack(truth[b]) for b in BANDS},
               series={rid: {k: {b: pack(data[rid][k][b]) for b in BANDS} for k, _ in ENGINES} for rid, _, _ in recv})

html = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zenodo_viewer_template.html")).read()
open(out_html, "w").write(html.replace("__DATA__", json.dumps(payload, separators=(",", ":"))))
print("wrote", out_html, os.path.getsize(out_html) // 1024, "kB")
