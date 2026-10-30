# GaiaVario as a GaiaSeis engine plug-in (interface: GaiaSeisEngineAPI).
# The platform selects it with `engine = "gaiavario"` in the profile; the
# engine options come from the profile as a flat Dict:
#   mode            "single" (default) | "dual" | "mean" | "l1gf"
#   elev_mask_deg   10.0
#   systems         "GREC" (subset)
#   weight_model    "sin2" | "sin" | "uniform" | ...
#   reject_k        4.0

import GaiaSeisEngineAPI as API

mutable struct GaiaVarioEngine <: API.AbstractEngine
    eng::Engine
end

function _options_from_dict(opts::Dict{String,Any})
    o = Options()
    haskey(opts, "mode") && (o.mode = Symbol(opts["mode"]))
    haskey(opts, "elev_mask_deg") && (o.elev_mask_deg = Float64(opts["elev_mask_deg"]))
    haskey(opts, "systems") && (o.systems = collect(Char, String(opts["systems"])))
    haskey(opts, "weight_model") && (o.weight_model = Symbol(opts["weight_model"]))
    haskey(opts, "reject_k") && (o.reject_k = Float64(opts["reject_k"]))
    haskey(opts, "min_sats") && (o.min_sats = Int(opts["min_sats"]))
    haskey(opts, "robust") && (o.robust = Symbol(opts["robust"]))
    return o
end

const PARAMS = [
    API.ParamSpec("mode", "Observable", :choice; default = "l1gf", choices = ["l1gf", "single", "dual", "mean"],
                  help = "single = L1 alone (quietest at 10 Hz); dual = ionosphere-free pair; mean = average of two phases; l1gf = L1 corrected with the smoothed geometry-free rate"),
    API.ParamSpec("elev_mask_deg", "Elevation mask (°)", :float; default = 10.0, min = 0.0, max = 60.0, step = 1.0,
                  help = "satellites below this elevation are left out of the fit"),
    API.ParamSpec("weight_model", "Weights", :choice; default = "sin2", choices = ["sin2", "sin", "uniform", "sin2floor", "cnr", "sin2cnr"],
                  help = "sin2 = sin²(elevation), the 1/sin noise model (best on the receivers tested)"),
    API.ParamSpec("reject_k", "Rejection threshold (σ)", :float; default = 4.0, min = 2.0, max = 20.0, step = 0.5,
                  help = "a satellite whose weighted residual exceeds this many robust sigmas is left out of the pair"),
    API.ParamSpec("robust", "Rejection rule", :choice; default = "reject", choices = ["reject", "none", "huber", "hyst"],
                  help = "reject = hard rejection; none = keep all; huber = down-weight; hyst = rejected satellites stay out for a few epochs"),
    API.ParamSpec("min_sats", "Minimum satellites", :int; default = 5, min = 4, max = 20, step = 1,
                  help = "fewer satellites than this: no solution for the epoch"),
]
const DESCRIPTION = "GaiaVario: in-house variometric engine (continuous GPS time, one broadcast record per epoch pair, per-satellite rejection). Default observable l1gf (L1 with the geometry-free ionospheric correction)."

function _create(station::API.Station, opts::Dict{String,Any})
    GaiaVarioEngine(Engine(station.xyz; opt = _options_from_dict(opts)))
end

API.engine_info(e::GaiaVarioEngine) = "GaiaVario mode=$(e.eng.opt.mode) mask=$(e.eng.opt.elev_mask_deg)° weights=$(e.eng.opt.weight_model)"

function API.engine_nav!(e::GaiaVarioEngine, nav::API.NavRow, t::API.EpochTime)
    push!(e.eng.nav, eph_from_row(nav.sys, nav.prn, nav.row, GpsTime(t.week, t.tow)))
    return nothing
end

_to_satobs(s::API.SatObservation) = SatObs(s.sys, s.prn, s.phase, s.code, s.lock, s.half, s.cnr)

function API.engine_step!(e::GaiaVarioEngine, ep::API.ObsEpoch)
    sol = step!(e.eng, Epoch(GpsTime(ep.t.week, ep.t.tow), [_to_satobs(s) for s in ep.sats]))
    lines = [API.SatLine(r.sys, r.prn, r.elev_deg, r.azim_deg, r.resid_m, r.weight, r.used, r.reason) for r in sol.sats]
    return API.EngineResult(sol.dxyz, sol.dneu, sol.clk, sol.nsat_in, sol.nsat_used, sol.nsys,
                            sol.pdop, sol.rms_m, sol.sig_neu, lines, sol.status)
end

API.engine_reset!(e::GaiaVarioEngine) = (reset!(e.eng); nothing)

function __init__()
    API.register_engine!(:gaiavario, _create; params = PARAMS, description = DESCRIPTION)
end
