# Live dual-engine comparison on the same relay bytes: a GaiaSeis session
# (PPPHVA, production settings, all outputs blank) and GaiaVario fed from a
# second connection to the same local relays. Neither touches the receiver or
# the production outputs. The relay bytes are also recorded for offline replay.
#
#   cd NtripClient.jl && GV_MODE=single julia -t auto --project=. ../GaiaVario.jl/tools/compare_live.jl [minutes]
#   (GV_MODE: single | dual | mean | l1gf; default single)
#
# Writes logs/UPC10Z_cmp_<stamp>.txt (PPPHVA per-epoch log, by the session),
# logs/gaiavario_cmp_<stamp>.csv (GaiaVario per-epoch solutions) and
# logs/rtcm_raw_<stamp>.bin (records of [Float64 wall time][UInt8 stream 0=obs
# 1=eph][UInt32 length][bytes]; see replay_gaiavario.jl).

using GaiaSeis, Dates, Sockets, Printf
include(joinpath(@__DIR__, "..", "src", "GaiaVario.jl"))
using .GaiaVario
const GV = GaiaVario
include(joinpath(@__DIR__, "live_common.jl"))

mins = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 10.0
const OBS = ("127.0.0.1", 2199, "B")
const EPH = ("127.0.0.1", 2198, "BCEP00BKG0")
stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
mkpath("logs")
csv_path = joinpath("logs", "gaiavario_cmp_$(stamp).csv")
raw_path = joinpath("logs", "rtcm_raw_$(stamp).bin")

# ---------------- PPPHVA side: a GaiaSeis session with the production settings
p = GaiaSeis.NtripProfile(name = "UPC10Z_cmp", caster = OBS[1], port = OBS[2], mountpoint = OBS[3],
                          use_eph = true, eph_caster = EPH[1], eph_port = EPH[2], eph_mountpoint = EPH[3],
                          process_mode = "DFV", elv_cutoff_deg = 10.0,
                          antex_file = joinpath(pwd(), "STHCR3-G3.atx"),
                          enable_gps = true, enable_glo = true, enable_gal = true, enable_bds = true,
                          miniseed_mseed_dir = "", miniseed_csv_path = "", datalink_host = "",
                          output_filename = "", log_dir = "logs")
session = GaiaSeis.NtripSession(p)

# ---------------- GaiaVario side
function ntrip_open(host, port, mount)
    sock = connect(host, port)
    write(sock, "GET /$mount HTTP/1.1\r\nHost: $host\r\nNtrip-Version: Ntrip/2.0\r\n" *
                "User-Agent: NTRIP GaiaVario/0.1\r\nAuthorization: Basic eDp4\r\nConnection: close\r\n\r\n")
    hdr = UInt8[]
    while isopen(sock)
        push!(hdr, read(sock, UInt8))
        length(hdr) >= 4 && hdr[end-3:end] == UInt8['\r', '\n', '\r', '\n'] && break
    end
    return sock
end

# current GPS seconds from the wall clock (UTC + 18 s leap seconds)
unix_gps_s() = Dates.value(now(UTC) - DateTime(1980, 1, 6)) / 1000 + 18.0
const WEEK = floor(Int, unix_gps_s() / 604800)

dec = GaiaSeis.RtcmDecoder()
declock = ReentrantLock()
nav = NavStore()
rec = Ref{Union{Nothing,NTuple{3,Float64}}}(nothing)
eng = Ref{Union{Nothing,Engine}}(nothing)
csv = open(csv_path, "w")
println(csv, "# GaiaVario live comparison  obs=$(OBS[1]):$(OBS[2])/$(OBS[3])  eph=$(EPH[1]):$(EPH[2])/$(EPH[3])  started $(now(UTC)) UTC  week=$WEEK")
println(csv, CSV_HEADER)
raw = open(raw_path, "w")
rawlock = ReentrantLock()
function dump_raw(stream::Symbol, buf)
    lock(rawlock) do
        write(raw, time()); write(raw, UInt8(stream === :obs ? 0 : 1)); write(raw, UInt32(length(buf))); write(raw, buf)
    end
end
nsolved = Ref(0); nepochs = Ref(0)
running = Ref(true)

function onepoch(t, sats)
    rec[] === nothing && return
    if eng[] === nothing
        eng[] = Engine(rec[]; nav = nav, opt = Options(mode = Symbol(get(ENV, "GV_MODE", "single")), elev_mask_deg = 10.0))
        println(csv, "# GaiaVario mode ", eng[].opt.mode)
        println(csv, "# station ECEF $(rec[][1]) $(rec[][2]) $(rec[][3])")
    end
    nepochs[] += 1
    t0 = time_ns()
    sol = step!(eng[], Epoch(t, sats))
    ms = (time_ns() - t0) / 1e6
    sol.status == :ok && (nsolved[] += 1)
    write_solution(csv, t, sol, ms)
end
grouper = Grouper(WEEK, onepoch)

function handle_events(evs)
    for ev in evs
        if ev isa GaiaSeis.EvStation
            rec[] === nothing && (rec[] = (ev.coord.x, ev.coord.y, ev.coord.z))
        elseif ev isa GaiaSeis.EvEphemeris
            ev.record.valid || continue
            tnear = GpsTime(WEEK, unix_gps_s() - WEEK * 604800)
            try
                push!(nav, eph_from_row(ev.record.sys, ev.record.prn, ev.record.row, tnear))
            catch err
                @warn "eph adapter" sys = ev.record.sys prn = ev.record.prn err
            end
        elseif ev isa GaiaSeis.EvEpoch
            push_frame!(grouper, ev.epoch)
        end
    end
end

function reader(stream::Symbol, host, port, mount)
    while running[]
        try
            sock = ntrip_open(host, port, mount)
            while running[] && isopen(sock)
                buf = readavailable(sock)
                isempty(buf) && (sleep(0.02); continue)
                dump_raw(stream, buf)
                lock(declock) do
                    handle_events(GaiaSeis.feed_bytes!(dec, buf; stream = stream))
                end
            end
            close(sock)
        catch err
            running[] || break
            @warn "reader $stream reconnecting" err
            sleep(2)
        end
    end
end

println("comparison start $(now()) for $(mins) min; PPPHVA log by the session, GaiaVario -> $csv_path"); flush(stdout)
GaiaSeis.start!(session)
t_eph = @async reader(:eph, EPH...)
t_obs = @async reader(:obs, OBS...)
t_end = time() + mins * 60
while time() < t_end
    sleep(30)
    st = session.stats
    @printf("[%s] PPPHVA epochs=%d missing=%d nsat=%d | GaiaVario epochs=%d solved=%d nav=%d station=%s\n",
            Dates.format(now(), "HH:MM:SS"), st.input_epochs, st.input_missing, st.last_epoch_nsat,
            nepochs[], nsolved[], length(nav), rec[] === nothing ? "waiting" : "ok")
    flush(stdout); flush(csv)
end
running[] = false
try; GaiaSeis.stop!(session); catch; end
sleep(1)
close(csv)
close(raw)
println("comparison end $(now()); GaiaVario epochs=$(nepochs[]) solved=$(nsolved[])")
println("gaiavario_csv=", csv_path)
println("raw_rtcm=", raw_path)
