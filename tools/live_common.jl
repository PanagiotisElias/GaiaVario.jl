# Shared by compare_live.jl and replay_gaiavario.jl: MSM epochs from the
# GaiaSeis decoder -> GaiaVario epochs, with epoch grouping across systems.
# Expects `GaiaSeis` and `GaiaVario` (as GV) to be loaded by the includer.

# MSM epoch encoding -> GPS time (all systems on the GPS scale) for GPS week `week`
function epoch_time(sys::Char, epoch_ms::Float64, week::Int)
    if sys == 'G' || sys == 'E'
        return GV.GpsTime(week, epoch_ms / 1000)
    elseif sys == 'C'
        return GV.GpsTime(week, epoch_ms / 1000 + 14.0)
    elseif sys == 'R'
        dow = floor(epoch_ms / 2^27); tod = mod(epoch_ms, 2^27) / 1000
        return GV.GpsTime(week, dow * 86400 + tod - 3 * 3600 + 18.0)
    end
    return GV.GpsTime(week, NaN)
end

# Signal preference inside a frequency slot: open civil codes before the
# semi-codeless / P codes (GPS 2X/2L before 2W). Lower rank wins.
# pref = :civil (default) or :pcode (semi-codeless W and P codes first).
function _code_rank(code::String, pref::Symbol = :civil)
    pref == :pcode && return (endswith(code, "W") || endswith(code, "P")) ? 1 : 2
    return endswith(code, "W") ? 3 : endswith(code, "P") ? 2 : 1
end

function to_satobs(ep::GaiaSeis.MsmEpoch, pref::Symbol = :civil)
    out = Dict{Int,GV.SatObs}()
    rank = Dict{Tuple{Int,Int},Int}()
    for o in ep.obs
        s = get!(out, o.prn) do; GV.SatObs(ep.sys, o.prn); end
        sl = o.freq_no
        sl <= 0 && continue
        isfinite(o.carrier_m) || continue
        r = _code_rank(o.rinex_code, pref)
        haskey(rank, (o.prn, sl)) && rank[(o.prn, sl)] <= r && continue
        rank[(o.prn, sl)] = r
        s.phase[sl] = o.carrier_m; s.code[sl] = o.pseudorange_m
        s.lock[sl] = o.lock_indicator; s.half[sl] = o.half_cycle_amb; s.cnr[sl] = o.cnr_dbhz
    end
    return collect(values(out))
end

"""Groups per-system MSM frames by timestamp and calls `onepoch(t, sats)` for a
timestamp once a newer one arrives."""
mutable struct Grouper
    week::Int
    groups::Dict{Int,Vector{GV.SatObs}}
    times::Dict{Int,GV.GpsTime}
    cur::Int
    onepoch::Function
    pref::Symbol
end
Grouper(week::Int, onepoch::Function; pref::Symbol = :civil) =
    Grouper(week, Dict{Int,Vector{GV.SatObs}}(), Dict{Int,GV.GpsTime}(), -1, onepoch, pref)

function push_frame!(g::Grouper, ep::GaiaSeis.MsmEpoch)
    t = epoch_time(ep.sys, ep.epoch_ms, g.week)
    isfinite(t.tow) || return
    key = round(Int, t.tow * 10)
    if key > g.cur
        g.cur >= 0 && haskey(g.groups, g.cur) && flush_group!(g, g.cur)
        g.cur = key
    elseif key < g.cur
        return                                     # late frame of an older epoch
    end
    append!(get!(g.groups, key, GV.SatObs[]), to_satobs(ep, g.pref))
    g.times[key] = t
end

function flush_group!(g::Grouper, key::Int)
    sats = pop!(g.groups, key); t = pop!(g.times, key)
    g.onepoch(t, sats)
end

# One CSV line per solution, same columns as compare_live.jl
function write_solution(io, t::GV.GpsTime, sol, ms)
    rej = join((Printf.@sprintf("%c%02d", r.sys, r.prn) for r in sol.sats if r.reason == :rejected), ";")
    slp = join((Printf.@sprintf("%c%02d", r.sys, r.prn) for r in sol.sats if r.reason == :slip), ";")
    clk(c) = get(sol.clk, c, NaN)
    used = [r.resid_m for r in sol.sats if r.used && isfinite(r.resid_m)]
    rms_u = isempty(used) ? NaN : sqrt(sum(abs2, used) / length(used))
    Printf.@printf(io, "%.1f,%.1f,%.6f,%.6f,%.6f,%d,%d,%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.2f,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f\n",
            mod(t.tow, 86400.0), t.tow, sol.dneu[1], sol.dneu[2], sol.dneu[3], sol.nsat_in, sol.nsat_used, sol.nsys,
            sol.pdop, sol.rms_m, sol.sig_neu[1], sol.sig_neu[2], sol.sig_neu[3], ms, sol.status, rej, slp,
            clk('G'), clk('R'), clk('E'), clk('C'), rms_u)
end
const CSV_HEADER = "sod,tow,N,E,U,nsat_in,nsat_used,nsys,pdop,rms,sigN,sigE,sigU,solve_ms,status,rejected,slips,clkG,clkR,clkE,clkC,rms_u"
