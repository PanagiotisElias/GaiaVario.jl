# Replay a raw relay capture (logs/rtcm_raw_<stamp>.bin from compare_live.jl)
# through GaiaVario with one or more option variants; one CSV per variant.
#
#   cd NtripClient.jl && julia --project=. ../GaiaVario.jl/tools/replay_gaiavario.jl logs/rtcm_raw_<stamp>.bin OUTPREFIX \
#        "base" "sin:weights=sin" "cnr:weights=cnr" "gal7:gal=7" "bds26:bds1=2,bds2=6" ...
#
# Variant syntax: name[:key=value,key=value]; keys: weights (sin2|sin|sin2floor|cnr|sin2cnr|uniform),
# floor (deg), gal (slot2 candidates for Galileo, e.g. 7 or 7;8;5), bds1, bds2 (idem for BeiDou),
# gps2 (GPS slot2 candidates), mode (dual|single), k (reject_k), mask (deg).
# Output: OUTPREFIX_<name>.csv in the same columns as the live tool.

using GaiaSeis, Dates, Printf
include(joinpath(@__DIR__, "..", "src", "GaiaVario.jl"))
using .GaiaVario
const GV = GaiaVario
include(joinpath(@__DIR__, "live_common.jl"))

raw_path = ARGS[1]; prefix = ARGS[2]
specs = length(ARGS) >= 3 ? ARGS[3:end] : ["base"]

function parse_variant(spec::String)
    i = findfirst(':', spec)
    name = i === nothing ? spec : spec[1:i-1]
    rest = i === nothing ? "" : spec[i+1:end]
    opt = Options(mode = :dual, elev_mask_deg = 10.0)
    pref = :civil
    isempty(rest) && return name, opt, pref
    for kv in split(rest, ',')
        k, v = split(kv, '=')
        slots(v) = [parse(Int, x) for x in split(v, ';')]
        if k == "weights"; opt.weight_model = Symbol(v)
        elseif k == "floor"; opt.weight_floor_deg = parse(Float64, v)
        elseif k == "gal"; opt.slot2_by_sys['E'] = slots(v)
        elseif k == "bds1"; opt.slot1_by_sys['C'] = slots(v)
        elseif k == "bds2"; opt.slot2_by_sys['C'] = slots(v)
        elseif k == "gps2"; opt.slot2_by_sys['G'] = slots(v)
        elseif k == "mode"; opt.mode = Symbol(v)
        elseif k == "k"; opt.reject_k = parse(Float64, v)
        elseif k == "mask"; opt.elev_mask_deg = parse(Float64, v)
        elseif k == "pref"; pref = Symbol(v)
        else; error("unknown key $k")
        end
    end
    return name, opt, pref
end

# read all records once
records = Tuple{Float64,UInt8,Vector{UInt8}}[]
open(raw_path) do io
    try
        while !eof(io)
            tw = read(io, Float64); st = read(io, UInt8); n = read(io, UInt32)
            push!(records, (tw, st, read(io, n)))
        end
    catch err
        err isa EOFError || rethrow()            # a capture still being written
    end
end
week = floor(Int, (records[1][1] + 18.0 - Dates.value(DateTime(1980, 1, 6) - DateTime(1970, 1, 1)) / 1000) / 604800)
println("records $(length(records)), GPS week $week, span $(round(records[end][1] - records[1][1]; digits = 1)) s")

for spec in specs
    name, opt, pref = parse_variant(spec)
    out = "$(prefix)_$(name).csv"
    dec = GaiaSeis.RtcmDecoder()
    nav = NavStore()
    rec = Ref{Union{Nothing,NTuple{3,Float64}}}(nothing)
    eng = Ref{Union{Nothing,Engine}}(nothing)
    io = open(out, "w")
    println(io, "# GaiaVario replay of $raw_path  variant=$spec  week=$week")
    println(io, CSV_HEADER)
    nsol = Ref(0); nep = Ref(0)
    function onepoch(t, sats)
        rec[] === nothing && return
        if eng[] === nothing
            eng[] = Engine(rec[]; nav = nav, opt = opt)
            println(io, "# station ECEF $(rec[][1]) $(rec[][2]) $(rec[][3])")
        end
        nep[] += 1
        t0 = time_ns(); sol = step!(eng[], Epoch(t, sats)); ms = (time_ns() - t0) / 1e6
        sol.status == :ok && (nsol[] += 1)
        write_solution(io, t, sol, ms)
    end
    g = Grouper(week, onepoch; pref = pref)
    for (tw, st, buf) in records
        stream = st == 0 ? :obs : :eph
        tnear = GpsTime(week, tw + 18.0 - Dates.value(DateTime(1980, 1, 6) - DateTime(1970, 1, 1)) / 1000 - week * 604800)
        for ev in GaiaSeis.feed_bytes!(dec, buf; stream = stream)
            if ev isa GaiaSeis.EvStation
                rec[] === nothing && (rec[] = (ev.coord.x, ev.coord.y, ev.coord.z))
            elseif ev isa GaiaSeis.EvEphemeris
                ev.record.valid && push!(nav, eph_from_row(ev.record.sys, ev.record.prn, ev.record.row, tnear))
            elseif ev isa GaiaSeis.EvEpoch
                push_frame!(g, ev.epoch)
            end
        end
    end
    close(io)
    @printf("%-28s epochs %5d solved %5d -> %s\n", spec, nep[], nsol[], out)
end
