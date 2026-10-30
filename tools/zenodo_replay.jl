# Offline replay of a recorded RTCM file (e.g. the Timisoara shaking-table
# data set, zenodo.org/records/18847292) through BOTH engines on identical
# epochs: the GaiaSeis bridge with PPPHVA (production code path) and GaiaVario.
# Broadcast ephemerides come from a daily RINEX navigation file when the
# stream carries none, served to both engines with the same schedule.
#
#   cd NtripClient.jl && julia --project=. ../GaiaVario.jl/tools/zenodo_replay.jl \
#       RAW.rtcm3 NAV.rnx OUTPREFIX 2025-10-21 X Y Z T0_UTC_SOD T1_UTC_SOD [variant ...]
#
# T0/T1: window in UTC seconds of day to solve (a lead of 60 s is added so both engines
# are warm). Variants (GaiaVario only) as in replay_gaiavario.jl, default "vario".
# Output: OUTPREFIX_va.txt (bridge log: SOD nsat dX dY dZ ...) and OUTPREFIX_<variant>.csv.

using GaiaSeis, Dates, Printf
include(joinpath(@__DIR__, "..", "src", "GaiaVario.jl"))
using .GaiaVario
const GV = GaiaVario
include(joinpath(@__DIR__, "live_common.jl"))
include(joinpath(@__DIR__, "rinex_nav_rows.jl"))

raw, navpath, prefix, datestr = ARGS[1:4]
X, Y, Z = parse.(Float64, ARGS[5:7])
t0u, t1u = parse.(Float64, ARGS[8:9])
specs = length(ARGS) >= 10 ? ARGS[10:end] : ["vario"]
const LEAP = 18.0
date = Date(datestr)
days = (date - Date(1980, 1, 6)).value
week = div(days, 7)
day0 = GV.GpsTime(week, mod(days, 7) * 86400.0)
tw0 = day0 + (t0u + LEAP - 60.0)          # first epoch handled (lead)
tw1 = day0 + (t1u + LEAP + 10.0)
tsolve0 = day0 + (t0u + LEAP)
println("week $week, window GPS tow $(tw0.tow) .. $(tw1.tow)")

const PREF = Ref(:civil)
function parse_variant(spec::String)
    i = findfirst(':', spec)
    name = i === nothing ? spec : spec[1:i-1]
    rest = i === nothing ? "" : spec[i+1:end]
    opt = Options(mode = :dual, elev_mask_deg = 10.0)
    for kv in (isempty(rest) ? String[] : split(rest, ','))
        k, v = split(kv, '=')
        slots(v) = [parse(Int, x) for x in split(v, ';')]
        if k == "weights"; opt.weight_model = Symbol(v)
        elseif k == "mode"; opt.mode = Symbol(v)
        elseif k == "gal"; opt.slot2_by_sys['E'] = slots(v)
        elseif k == "bds1"; opt.slot1_by_sys['C'] = slots(v)
        elseif k == "bds2"; opt.slot2_by_sys['C'] = slots(v)
        elseif k == "gps2"; opt.slot2_by_sys['G'] = slots(v)
        elseif k == "mask"; opt.elev_mask_deg = parse(Float64, v)
        elseif k == "k"; opt.reject_k = parse(Float64, v)
        elseif k == "pref"; PREF[] = Symbol(v)
        elseif k == "sys"; opt.systems = collect(Char, v)
        elseif k == "robust"; opt.robust = Symbol(v)
        elseif k == "clock"; opt.clock_model = Symbol(v)
        elseif k == "gftau"; opt.gf_tau_s = parse(Float64, v)
        elseif k == "gfdb"; opt.gf_deadband_mm_s = parse(Float64, v)
        elseif k == "huber"; opt.huber_c = parse(Float64, v)
        elseif k == "hyst"; opt.hyst_epochs = parse(Int, v)
        elseif k == "taper"; opt.taper_deg = parse(Float64, v)
        elseif k == "admit"; opt.admit_epochs = parse(Int, v)
        else; error("unknown key $k")
        end
    end
    return name, opt
end

# ------------------------------------------------------------ navigation: RINEX schedule or the stream itself
const NAV_FROM_STREAM = (navpath == "stream")
recs = NAV_FROM_STREAM ? Tuple{Char,Int,GV.GpsTime,Vector{Float64}}[] : rinex_nav_rows(navpath)
NAV_FROM_STREAM || println("RINEX navigation records: $(length(recs))  (", join(["$s:$(count(r -> r[1] == s, recs))" for s in "GREC"], " "), ")")
NAV_FROM_STREAM && println("ephemerides from the stream's own messages")
sched = NavSchedule(recs)
const T_NEAR = day0 + ((t0u + t1u) / 2 + LEAP)        # reference time that places sod/sow references
const GLO_K = Dict{Int,Int}(r[2] => Int(round(r[4][2])) for r in recs if r[1] == 'R')

# ------------------------------------------------------------ PPPHVA through the bridge
prof = GaiaSeis.NtripProfile(name = "ZENODO_va", caster = "replay", mountpoint = "replay",
                             process_mode = get(ENV, "VA_MODE", "DFV"), elv_cutoff_deg = 10.0,
                             antex_file = joinpath(pwd(), "igs14_2163.atx"),
                             enable_gps = true, enable_glo = true, enable_gal = true, enable_bds = true,
                             use_ssr = false, use_eph = false, solve_incomplete_epochs = false,
                             max_gap_epochs = 5, outlier_mad_threshold = 6.0, outlier_window = 20,
                             miniseed_mseed_dir = "", miniseed_csv_path = "", datalink_host = "")
sess = GaiaSeis.NtripSession(prof)
sess.stats = GaiaSeis.SessionStats(started_at = Dates.now())
sess.bridge_state = nothing
va_log = prefix * "_va.txt"
sess.log_file = open(va_log, "w")
println(sess.log_file, "# GaiaSeis replay  profile=ZENODO_va  raw=$raw")
println(sess.log_file, "#    SOD nsat            dX            dY            dZ           DTG           DTR           DTE")
dec = GaiaSeis.RtcmDecoder()
dec.station = GaiaSeis.StationCoord(1, X, Y, Z, NaN)
dec.descriptor = GaiaSeis.StationDescriptor(1, "JAVRINGANT_G5T  JAVC", "", "JAVAD_DELTA_3S", "", "")

# ------------------------------------------------------------ GaiaVario variants
struct Variant
    name::String; opt::Options; nav::GV.NavStore
    eng::Base.RefValue{Union{Nothing,GV.Engine}}; io::IOStream; grp::Grouper
    nsol::Base.RefValue{Int}; nep::Base.RefValue{Int}
end
variants = Variant[]
for spec in specs
    PREF[] = :civil
    name, opt = parse_variant(spec)
    pref = PREF[]
    io = open("$(prefix)_$(name).csv", "w")
    println(io, "# GaiaVario zenodo replay of $raw  variant=$spec  week=$week")
    println(io, CSV_HEADER); println(io, "# station ECEF $X $Y $Z")
    nav = GV.NavStore(); eng = Ref{Union{Nothing,GV.Engine}}(nothing)
    nsol = Ref(0); nep = Ref(0)
    function onepoch(t, sats)
        t < tsolve0 && eng[] === nothing && false
        if eng[] === nothing
            eng[] = GV.Engine((X, Y, Z); nav = nav, opt = opt)
        end
        nep[] += 1
        t0 = time_ns(); sol = GV.step!(eng[], GV.Epoch(t, sats)); ms = (time_ns() - t0) / 1e6
        t < tsolve0 && return                              # lead-in epochs only warm the engine
        sol.status == :ok && (nsol[] += 1)
        write_solution(io, t, sol, ms)
    end
    push!(variants, Variant(name, opt, nav, eng, io, Grouper(week, onepoch; pref = pref), nsol, nep))
end

# ------------------------------------------------------------ main loop
include(joinpath(@__DIR__, "obs_dump.jl"))
nva = 0; nstep = 0; last_t = GV.GpsTime(week, 0.0)
const FIXED_DESC = GaiaSeis.StationDescriptor(1, "JAVRINGANT_G5T  JAVC", "", "RECEIVER", "", "")

function process_epoch!(ep::GaiaSeis.MsmEpoch)
    global nva, nstep, last_t
    t = epoch_time(ep.sys, ep.epoch_ms, week)
    (isfinite(t.tow) && !(t < tw0) && !(tw1 < t)) || return
    if !NAV_FROM_STREAM && last_t < t
        for (sys, prn, toc, row) in advance!(sched, t)
            dec.eph[(sys, prn)] = GaiaSeis.EphemerisRecord(sys, prn, copy(row), true)
            for v in variants
                push!(v.nav, GV.eph_from_row(sys, prn, row, t))
            end
        end
        last_t = t
    end
    ep.sys == GaiaSeis._primary_sys(prof) && GaiaSeis._note_input_epoch!(sess.stats, ep.sod)
    sess.stats.last_msg_type = ep.msg_type
    try
        solved = GaiaSeis.bridge_step!(sess, dec, ep; stale = false)
        isnan(solved) || (nva += 1)
    catch err
        nstep += 1
        nstep <= 3 && @warn "bridge step failed" err
    end
    for v in variants
        push_frame!(v.grp, ep)
    end
end

if endswith(raw, ".dump")
    for ep in read_obs_dump(raw, week, GLO_K)
        process_epoch!(ep)
    end
else
    open(raw) do io
        while !eof(io)
            buf = read(io, 1 << 18)
            evs = GaiaSeis.feed_bytes!(dec, buf; stream = :obs)
            dec.station = GaiaSeis.StationCoord(1, X, Y, Z, NaN)       # identical position and antenna for every receiver
            dec.descriptor = FIXED_DESC
            for ev in evs
                if ev isa GaiaSeis.EvEphemeris && NAV_FROM_STREAM && ev.record.valid
                    for v in variants
                        try; push!(v.nav, GV.eph_from_row(ev.record.sys, ev.record.prn, ev.record.row, T_NEAR)); catch; end
                    end
                elseif ev isa GaiaSeis.EvEpoch
                    process_epoch!(ev.epoch)
                end
            end
        end
    end
end
close(sess.log_file)
for v in variants
    close(v.io)
    @printf("GaiaVario %-14s epochs %5d solved %5d -> %s_%s.csv\n", v.name, v.nep[], v.nsol[], prefix, v.name)
end
println("PPPHVA bridge: groups solved $nva, step errors $nstep -> $va_log")
