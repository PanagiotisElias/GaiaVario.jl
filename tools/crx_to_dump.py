#!/usr/bin/env python3
"""Decode a Compact RINEX 3 (Hatanaka) observation file and write an observation
dump for obs_dump.jl, only for a time window. Stdlib only.

    python3 crx_to_dump.py FILE.crx[.gz] OUT.dump GPS_SOW_START GPS_SOW_END [SYSTEMS]

Compact RINEX 3 (format description: Hatanaka 2008, "A compression format and
tools for GNSS observation data", Bull. Geospat. Inf. Auth. Japan 55): after the
header, every epoch is an epoch line (text-differenced against the previous one),
a receiver-clock line, and one line per satellite with the observables as
space-separated integers (units 0.001) coded as order-3 differences ("3&value"
starts an arc), followed by the LLI/SSI characters (text-differenced).
"""
import sys, gzip, datetime

def open_text(path):
    return gzip.open(path, "rt", encoding="latin1") if path.endswith(".gz") else open(path, encoding="latin1")

class Arc:
    """Order-m differential arc of one observable."""
    __slots__ = ("order", "y", "n")
    def __init__(self, order, value):
        self.order = order; self.y = [value, 0, 0, 0]; self.n = 1
    def update(self, d):
        m = min(self.n, self.order)
        if m == 0:
            self.y[0] = d
        else:
            self.y[m] = d
            for j in range(m - 1, -1, -1):
                self.y[j] += self.y[j + 1]
        self.n += 1
        return self.y[0]

def overlay(prev, new):
    """Text difference: ' ' keeps the previous character, '&' makes a blank; positions beyond the new line keep the previous."""
    out = list(prev)
    if len(new) > len(out):
        out.extend(" " * (len(new) - len(out)))
    for i, c in enumerate(new):
        if c == " ":
            continue
        out[i] = " " if c == "&" else c
    return "".join(out)

def main():
    path, outp, s0, s1 = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    keep = sys.argv[5] if len(sys.argv) > 5 else "GREC"
    f = open_text(path)
    # ---- header
    header = []
    for line in f:
        header.append(line.rstrip("\n"))
        if "END OF HEADER" in line:
            break
    if "CRINEX VERS" not in header[0]:
        sys.exit("not a Compact RINEX file")
    obs_types = {}
    cur = None
    for l in header:
        if "SYS / # / OBS TYPES" in l:
            if l[0] != " ":
                cur = l[0]; n = int(l[3:6]); obs_types[cur] = l[7:60].split()
            else:
                obs_types[cur] += l[7:60].split()
    first_date = None
    for l in header:
        if "TIME OF FIRST OBS" in l:
            p = l.split(); first_date = (int(p[0]), int(p[1]), int(p[2]))
    if first_date is None:
        sys.exit("no TIME OF FIRST OBS")
    d0 = datetime.date(*first_date); gps0 = datetime.date(1980, 1, 6)

    prev_epoch = ""
    state = {}                 # (sat, obsindex) -> Arc
    flags_prev = {}            # sat -> flag string
    out = open(outp, "w")
    nep = nwr = 0
    done = False
    while True:
        line = f.readline()
        if not line:
            break
        line = line.rstrip("\n")
        if not line.strip():
            continue
        # ---- epoch line
        if line[0] == "&":                           # re-initialisation of all arcs
            state.clear(); flags_prev.clear()
            ep = ">" + line[1:]
        else:
            ep = overlay(prev_epoch, line)
        try:
            flag = int(ep[31]); nrec = int(ep[32:35])
        except (ValueError, IndexError):
            sys.exit(f"bad epoch line: {ep!r}")
        if flag > 1:
            # Event record (RINEX epoch flag 2-5 start of moving antenna, new site, header information,
            # external event; 6 is a cycle-slip record with satellite lines and is handled below):
            # the NREC lines that follow are text, not observations. Skip them; the differencing state of
            # the epoch lines and of the observation arcs is left untouched. (Assumption, not verified
            # against a file that contains such records: the epoch line of an event is not differenced
            # against and does not become the reference of the next epoch line.)
            if flag != 6:
                for _ in range(nrec):
                    f.readline()
                continue
        prev_epoch = ep
        clk = f.readline()                           # receiver clock line (ignored, but consumed)
        try:
            y, mo, d, h, mi = int(ep[2:6]), int(ep[7:9]), int(ep[10:12]), int(ep[13:15]), int(ep[16:18])
            sec = float(ep[18:29])
            nsat = nrec
        except ValueError:
            sys.exit(f"bad epoch line: {ep!r}")
        sats = [ep[41 + 3 * i: 44 + 3 * i] for i in range(nsat)]
        days = (datetime.date(y, mo, d) - gps0).days
        sow = (days % 7) * 86400 + h * 3600 + mi * 60 + sec
        nep += 1
        rec = {}
        for sat in sats:
            sl = f.readline().rstrip("\n")
            sysc = sat[0]
            types = obs_types.get(sysc)
            if types is None:
                continue
            n = len(types)
            pos = 0; toks = []
            for i in range(n):
                j = sl.find(" ", pos)
                if j == -1:
                    toks.append(sl[pos:]); pos = len(sl) + 1
                else:
                    toks.append(sl[pos:j]); pos = j + 1
            toks += [""] * (n - len(toks))
            ftxt = sl[pos:] if pos <= len(sl) else ""
            fl = overlay(flags_prev.get(sat, ""), ftxt)
            fl = fl.ljust(2 * n)
            flags_prev[sat] = fl
            vals = {}
            for i, tk in enumerate(toks):
                key = (sat, i)
                if tk == "":
                    state.pop(key, None); continue
                if len(tk) >= 2 and tk[1] == "&":
                    arc = Arc(int(tk[0]), int(tk[2:])); state[key] = arc
                    v = arc.y[0]
                else:
                    arc = state.get(key)
                    if arc is None:
                        continue                      # difference without an arc: cannot be decoded
                    v = arc.update(int(tk))
                lli = fl[2 * i] if 2 * i < len(fl) else " "
                ssi = fl[2 * i + 1] if 2 * i + 1 < len(fl) else " "
                vals[types[i]] = (v / 1000.0, int(lli) if lli.isdigit() else 0)
            rec[sat] = vals
        if flag != 0:
            continue
        if sow > s1:
            break
        if sow >= s0:
            out.write(f"E {sow:.3f}\n")
            for sat, vals in rec.items():
                if sat[0] not in keep: continue
                out.write("S " + sat + " " + " ".join(
                    f"{code}={v:.3f}" + (f":{lli}" if code[0] == "L" else "") for code, (v, lli) in vals.items()) + "\n")
            nwr += 1
    out.close()
    print(f"decoded {nep} epochs, wrote {nwr} between GPS sow {s0} and {s1}; last epoch sow {sow}")

main()
