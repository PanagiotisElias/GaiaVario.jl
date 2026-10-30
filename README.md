# GaiaVario.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23810265.svg)](https://doi.org/10.5281/zenodo.23810265)

A variometric GNSS displacement engine for real-time seismogeodesy, written
as the default engine plug-in of the
[GaiaSeis](https://github.com/PanagiotisElias/GaiaSeis) platform. Version
0.1.0 (2026-10-30), licence GNU AGPL-3.0-only.

Deforming Gaia Lab, Department of Geology, University of Patras. Panagiotis
Elias (<pelias@upatras.gr>).

## What it does

For every satellite seen at two consecutive epochs with continuous lock, the
observation is the change of the carrier phase minus the change of everything
modelled: geometric range to the a-priori station with signal transit time and
Earth rotation, satellite clock with relativity, troposphere, solid Earth
tide, phase wind-up. The remainder is the station displacement over the pair
projected on the line of sight, plus the receiver clock change of that
constellation. Weighted least squares (elevation weights) with iterative
rejection on standardised residuals gives dX, dY, dZ and one clock change per
constellation, with per-satellite residuals.

Design rules:

- **One time scale.** GPS time is kept as (week, seconds of week); differences
  are exact, there is no seconds-of-day arithmetic, so a day boundary cannot
  wrap. Broadcast references are compared with the ICD half-week rule.
- **One record per pair.** Both epochs of a pair are modelled with the
  broadcast record that was current when the previous epoch was solved; a
  record arriving between the two epochs is used from the next pair on. A
  changeover therefore never enters a difference (test
  `ephemeris changeover produces no step`).
- **Flag, do not blank.** A satellite with a lock loss, a geometry-free jump or
  a large standardised residual is left out of that pair only; the solution
  reports why (`SatResult.reason`).
- **Windowed noise.** The robust σ̂ is the median of the last `sigma_window`
  weighted post-fit RMS values.

## Observables

| `Options().mode` | observable | when |
|---|---|---|
| `:l1gf` (default) | L1 phase corrected, per satellite, with a smoothed rate of the geometry-free phase (L1 − L2): the ionospheric change of L1 over an epoch is estimated from the geometry-free rate (exponential average, `gf_tau_s` = 30 s, kept through gaps up to `gf_gap_s` = 10 s) and added back. In a quiet ionosphere it is plain L1; in a disturbed one it keeps the ionospheric drift small. Satellites without a second frequency are used as plain L1. | general use |
| `:single` | L1 (E1, B1) alone; on every receiver tested 20–35 % quieter per epoch than the ionosphere-free combination | quiet ionosphere, 10 Hz |
| `:dual` | ionosphere-free pair (L1/L2, E1/E5a, B1I/B3I) | disturbed ionosphere, 1 Hz |
| `:mean` | average of the two phases | experiments |

Weights: `sin2` (sin² of elevation, default), `sin`, `uniform`, `sin2floor`,
`cnr`, `sin2cnr`. Rejection: `reject` (hard, `reject_k` = 4 robust sigmas on
the signed weighted residuals), `none`, `huber`, `hyst`. Civil codes are
preferred inside a frequency slot (GPS L2C `2X` before `2W`); BeiDou uses
B1I + B3I by default.

## Models and their sources

| term | source |
|---|---|
| GPS, Galileo, BeiDou MEO/IGSO orbits and clocks | IS-GPS-200 table 20-IV; Galileo OS SIS ICD; BeiDou SIS ICD (B1I) |
| BeiDou GEO | BeiDou SIS ICD 5.2.4.12 (−5° rotation) |
| GLONASS orbits | GLONASS ICD 5.1 appendix 3 (RK4 of the PZ-90 equations) |
| Signal transit time, Earth rotation | Sagnac iteration |
| Troposphere | Saastamoinen zenith delays from a standard atmosphere; Niell (1996) mapping functions |
| Solid Earth tide | IERS Conventions 2010 eq. 7.5, degree 2, nominal h₂ l₂; Sun/Moon from Montenbruck & Gill low-precision series |
| Phase wind-up | Wu et al. (1993), nominal yaw attitude |
| Antenna phase centres | not modelled: their projection on the line of sight changes by micrometres between epochs |

Because the engine differences consecutive epochs, every model term enters
only as its change over 0.1–1 s. This is why a standard atmosphere and Niell
mapping functions suffice (the pair difference of the troposphere is below
1 mm at 10° elevation and the mapping-function error is 1 % of that), and why
no external data files (ANTEX, GPT3 grid, VMF coefficients) are needed.

The code was written from the sources above. No code of other variometric
packages was consulted while writing it; the comparison with other engines
is done by outputs only (synthetic truth, side-by-side runs on the same
stream; report `docs/GaiaVario_vs_PPPHVA_2026-10-29.md` in the GaiaSeis
repository: 25–27 % lower per-epoch noise than the ionosphere-free solution
on a 10 Hz station, no step at the GPS day boundary).

## Install and use

GaiaVario depends on the interface package `GaiaSeisEngineAPI`, which lives
in the GaiaSeis repository (sub-directory `GaiaSeisEngineAPI/`).

As the GaiaSeis engine (the normal case): in the GaiaSeis project,

```julia
julia --project=. -e 'using Pkg; Pkg.add(url="https://github.com/PanagiotisElias/GaiaVario.jl")'
```

and in the profile `engine = "gaiavario"`, `engine_mode = "l1gf"` (or clone
this repository beside `GaiaSeis.jl` and set `engine_path = "../GaiaVario.jl"`).
The GaiaSeis GUI shows the engine's parameters on the *Engine* page of the
profile editor; headless profiles set them in `engine_options`
(`elev_mask_deg`, `weight_model`, `reject_k`, `robust`, `min_sats`).

Stand-alone (the Manifest points at the interface package in a GaiaSeis
clone beside this one):

```bash
git clone https://github.com/PanagiotisElias/GaiaSeis.git
git clone https://github.com/PanagiotisElias/GaiaVario.jl.git
cd GaiaVario.jl && julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

```julia
using GaiaVario
nav = NavStore()                        # push!(nav, eph_from_row(sys, prn, row, t))
eng = Engine((x, y, z); nav = nav)      # Options(mode = :single) for plain L1, :dual for ionosphere-free
sol = step!(eng, Epoch(t_gps, sats))    # sats::Vector{SatObs}
sol.dneu, sol.clk, sol.rms_m, sol.sats
```

Tests: `julia --project=. test/runtests.jl` (unit tests and synthetic-truth
tests: motion, noise, slips, ephemeris changeovers, a real broadcast set).

## Layout

```
src/constants.jl   physical constants, frequencies
src/time.jl        GPS time, wraps, day of year, sidereal time
src/types.jl       SatObs, Epoch, KeplerEph, GlonassEph, NavStore, Options, Solution
src/orbits.jl      satellite position and clock, Sagnac range
src/models.jl      geodesy, troposphere, Sun/Moon, tide, wind-up
src/adapters.jl    typed records from the GaiaSeis decoder's 39-column rows
src/engine.jl      Engine, step!, least squares with rejection
src/plugin.jl      the GaiaSeis plug-in (GaiaSeisEngineAPI registration, parameters)
src/synth.jl       synthetic-truth generator (motion, noise, slips, changeovers)
test/runtests.jl   unit and synthetic-truth tests
tools/             replay and comparison tools (recorded streams, RINEX/Zenodo records, reports)
```

## Citing

See `CITATION.cff`. Zenodo DOI [10.5281/zenodo.23810265](https://doi.org/10.5281/zenodo.23810265) (all versions; this release: [10.5281/zenodo.23810266](https://doi.org/10.5281/zenodo.23810266)).

## Licence

GNU Affero General Public License v3.0 only (`LICENSE`).
Copyright (C) 2026 Panagiotis Elias, Deforming Gaia Lab, Department of Geology, University of Patras.
Contributions: `CONTRIBUTING.md`.

## Funding

This work is part of the [TRANSFORM²](https://www.transform2-project.eu/) project which aims to improve physical and digital infrastructure across Near-Fault Observatories (NFOs) in Europe. TRANSFORM² is funded by the European Union under project number 101188365 within the HORIZON-INFRA-2024-DEV-01-01 call.

<div align="center"><img src="https://www.transform2-project.eu/wp-content/uploads/2022/08/Logo_TRANSFORM2-round-logo-100x100-1.png" alt="TRANSFORM² logo"></div>

Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Executive Agency (REA). Neither the European Union nor the granting authority can be held responsible for them.
