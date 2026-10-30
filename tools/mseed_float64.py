#!/usr/bin/env python3
"""Minimal miniSEED reader for the Zenodo accelerometer files (encoding 5,
IEEE float64, big-endian, 512-byte records). Stdlib + numpy.

    from mseed_float64 import read_mseed
    traces = read_mseed(path)   # {channel: (t0_epoch_s, rate_hz, samples ndarray)}
"""
import struct, datetime
import numpy as np

def read_mseed(path, rec_len=512):
    out = {}
    with open(path, "rb") as f:
        data = f.read()
    for off in range(0, len(data) - rec_len + 1, rec_len):
        d = data[off:off + rec_len]
        if d[6:7] not in (b"D", b"R", b"Q", b"M"):
            continue
        cha = d[15:18].decode().strip(); loc = d[13:15].decode().strip(); sta = d[8:13].decode().strip()
        yr, doy, hh, mm, ss, _, ms = struct.unpack(">HHBBBBH", d[20:30])
        ns, sf, sm = struct.unpack(">Hhh", d[30:36])
        dataoff, blkoff = struct.unpack(">HH", d[44:48])
        micro = 0
        b = blkoff; enc = None
        while b:
            t, nxt = struct.unpack(">HH", d[b:b + 4])
            if t == 1000: enc = d[b + 4]
            if t == 1001: micro = struct.unpack(">b", d[b + 5:b + 6])[0]
            b = nxt
        if enc != 5:
            raise ValueError(f"{path}: encoding {enc} not supported")
        rate = sf if sf > 0 else 1.0 / -sf
        if sm > 0: rate *= sm
        elif sm < 0: rate /= -sm
        t0 = (datetime.datetime(yr, 1, 1) + datetime.timedelta(days=doy - 1, hours=hh, minutes=mm,
              seconds=ss, microseconds=ms * 100 + micro) - datetime.datetime(1970, 1, 1)).total_seconds()
        x = np.frombuffer(d[dataoff:dataoff + 8 * ns], dtype=">f8").astype(float)
        key = f"{sta}.{loc}.{cha}"
        out.setdefault(key, []).append((t0, rate, x))
    res = {}
    for k, recs in out.items():
        recs.sort(key=lambda r: r[0])
        t0 = recs[0][0]; rate = recs[0][1]
        res[k] = (t0, rate, np.concatenate([r[2] for r in recs]), [r[0] for r in recs])
    return res
