#!/usr/bin/env python3
"""Fetch single members out of the zip archives of a Zenodo record with HTTP range requests,
without downloading the archives (stdlib only).

Data set used by the GaiaVario examples:
    Kudlacik, I., Kaplon, J., Nastase, E. I., Tiganescu, A., Kaczmarek, A., Elias, P.,
    Karakonstantis, A., Zuliani, D., Nistor, S. (2026). High-Rate GNSS and Multi-Sensor Benchmark
    Dataset from Earthquake Simulation Experiments (v1) [Data set]. Zenodo.
    https://doi.org/10.5281/zenodo.18847292   Licence: CC BY 4.0
    (Timisoara shaking-table experiments, October 2025; the data are not part of this repository:
    cite the record when you use them.)

Usage:
    python3 fetch_members.py [--out DIR] --list ARCHIVE
    python3 fetch_members.py [--out DIR] ARCHIVE MEMBER [MEMBER ...]

ARCHIVE is a file name of the record (for example GNSS_observations.zip), MEMBER a path inside it as
printed by --list. Files are written to DIR (default ./zenodo_18847292), without the directory part.

The members used by the Samos 2020 replay on the NET-G5 10 Hz log (about 116 MB in total):
    python3 fetch_members.py GNSS_observations.zip GNSS_observations/TM0/TM0J13ROU_R_20252940834_04H_10Z_MO.crx.gz
    python3 fetch_members.py RealTime_RTCM_otherFiles.zip RealTime_RTCM_otherFiles/EPH_NTRIP/TM0/BRDC00WRD_S_20252940000_01D_MN.rnx
    python3 fetch_members.py Accelerometer_records.zip Accelerometer_records/Corrected/TM0/RO_TIME_2025_294_TM0_cor.mseed
"""
import sys, io, os, zipfile, urllib.request

RECORD = "18847292"
BASE = f"https://zenodo.org/records/{RECORD}/files/%s?download=1"


class HttpFile(io.RawIOBase):
    """Read-only seekable file on top of HTTP range requests."""
    def __init__(self, url):
        self.url = url; self.pos = 0
        self.size = int(self._get(0, 0)[1].split("/")[1])

    def _get(self, a, b):
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={a}-{b}"})
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read(), r.headers.get("Content-Range")

    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos

    def seek(self, o, w=0):
        self.pos = o if w == 0 else self.pos + o if w == 1 else self.size + o
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        d, _ = self._get(self.pos, self.pos + n - 1)
        b[:len(d)] = d
        self.pos += len(d)
        return len(d)


def open_archive(name, buffer_size=8 << 20):
    return zipfile.ZipFile(io.BufferedReader(HttpFile(BASE % name), buffer_size=buffer_size))


def main(argv):
    out = "zenodo_18847292"
    if argv[:1] == ["--out"]:
        out = argv[1]; argv = argv[2:]
    if not argv:
        sys.exit(__doc__)
    if argv[0] == "--list":
        z = open_archive(argv[1], 1 << 20)
        for i in z.infolist():
            if not i.is_dir():
                print(f"{i.file_size:>12} {i.compress_size:>12} {i.filename}")
        return
    archive, members = argv[0], argv[1:]
    os.makedirs(out, exist_ok=True)
    z = open_archive(archive)
    for m in members:
        info = z.getinfo(m)
        dest = os.path.join(out, os.path.basename(m))
        with z.open(info) as src, open(dest, "wb") as dst:
            n = 0
            while True:
                chunk = src.read(4 << 20)
                if not chunk:
                    break
                dst.write(chunk); n += len(chunk)
        print(f"{dest}  {n} bytes (expected {info.file_size})", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
