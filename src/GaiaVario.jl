"""
GaiaVario — a variometric GNSS displacement engine for real-time seismic
monitoring, written for GaiaSeis.

Epoch-to-epoch carrier-phase differences from four constellations are turned
into station displacement with weighted least squares. The engine is built
around continuous GPS time (no day wraps), one broadcast record per satellite
for both epochs of a pair, per-satellite residual rejection and a windowed
robust noise estimate.

Models are implemented from the public sources named in each file: the GNSS
interface control documents, IERS Conventions 2010, Niell (1996) and
Wu et al. (1993). No code of other variometric packages was consulted.
"""
module GaiaVario

using LinearAlgebra

include("constants.jl")
include("time.jl")
include("types.jl")
include("orbits.jl")
include("models.jl")
include("adapters.jl")
include("engine.jl")
include("synth.jl")
include("plugin.jl")

export Engine, Epoch, SatObs, NavStore, Options, Solution, SatResult, GpsTime,
       KeplerEph, GlonassEph, step!, reset!, set_position!, eph_from_row,
       gps_time, gps_week_tow, to_seconds, sat_state, sagnac_range, ecef2geodetic,
       synthesize, SynthOptions, synthetic_constellation,
       motion_static, motion_step, motion_sine, motion_pulse

end # module
