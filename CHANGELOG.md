# GaiaVario changelog

## 0.1.0 (2026-10-30)

First public release, as the default engine plug-in of GaiaSeis 0.5.0.

- Variometric engine: GPS, GLONASS, Galileo, BeiDou from broadcast
  ephemerides; GPS time kept as (week, seconds of week); one broadcast record
  per epoch pair; models from the published sources listed in the README.
- Observables `l1gf` (default), `single`, `dual`, `mean`; weights `sin2`
  (default) and alternatives; per-satellite rejection on signed weighted
  residuals (k = 4), `huber` and `hyst` variants; civil-code preference inside
  a frequency slot; BeiDou B1I + B3I.
- GaiaSeis plug-in (`src/plugin.jl`): registration as `gaiavario`, parameters
  declared through `GaiaSeisEngineAPI.ParamSpec` for the platform's Engine page.
- Tests: unit tests and synthetic-truth tests (motion, noise, slips,
  ephemeris changeovers, a recorded broadcast set).
- Tools: replay of recorded RTCM streams and of RINEX/Zenodo records,
  comparison and report scripts.
