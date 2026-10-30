#!/usr/bin/env python3
"""Self-test of crx_to_dump.py on a hand-made Compact RINEX 3 file: three epochs of one satellite,
the order-3 arcs, a differenced epoch line, and an event record (flag 4) between the second and the
third epoch that must be skipped without disturbing the differencing.

    python3 test_crx_to_dump.py
"""
import os, subprocess, sys, tempfile

HEADER = """3.0                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE
test                                                        CRINEX PROG / DATE
     3.05           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G    2 C1C L1C                                              SYS / # / OBS TYPES
  2025    10    21     0     0    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
"""
BODY = (
    "> 2025 10 21  0  0  0.0000000  0  1      G01\n"            # epoch 1: full line
    "\n"                                                          # receiver clock line
    "3&20000000123 3&105000000456 \n"                             # arcs start: C1C 20000000.123, L1C 105000000.456
    "> 2025 10 21  0  0  1.0000000  0  1      G01\n"            # epoch 2
    "\n"
    "1000 2000 \n"                                                # first differences: +1.000 m, +2.000 cycles
    "> 2025 10 21  0  0  1.5000000  4  2\n"                      # event record, two text lines follow
    "THIS IS EVENT TEXT ONE                                      COMMENT\n"
    "THIS IS EVENT TEXT TWO                                      COMMENT\n"
    "                    2.0000000\n"                          # epoch 3, differenced against epoch 2 (the event line is not the reference)
    "\n"
    "0 0 \n"                                                      # second differences: no change of rate
)

def main():
    with tempfile.TemporaryDirectory() as d:
        crx = os.path.join(d, "t.crx"); out = os.path.join(d, "t.dump")
        open(crx, "w").write(HEADER + BODY)
        here = os.path.dirname(os.path.abspath(__file__))
        r = subprocess.run([sys.executable, os.path.join(here, "crx_to_dump.py"), crx, out, "0", "700000"], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr + r.stdout
        lines = [l.rstrip("\n") for l in open(out)]
    sows = [float(l.split()[1]) for l in lines if l.startswith("E ")]
    sat = [l for l in lines if l.startswith("S ")]
    assert sows == [172800.0, 172801.0, 172802.0], sows           # 2025-10-21 is a Tuesday: 2 days into the GPS week
    vals = [dict(t.split("=") for t in l.split()[2:]) for l in sat]
    want = [(20000000.123, 105000000.456), (20000001.123, 105000002.456), (20000002.123, 105000004.456)]
    for v, (c, p) in zip(vals, want):
        assert abs(float(v["C1C"]) - c) < 5e-4, (v, c)
        assert abs(float(v["L1C"].split(":")[0]) - p) < 5e-4, (v, p)
    print("crx_to_dump.py self-test passed: 3 epochs decoded, the event record skipped")

main()
