# Scan a raw RTCM 3 file with the GaiaSeis decoder: message counts, constellations, signal codes,
# ephemeris messages, station position and descriptor. Used to decide whether a recording needs a
# separate navigation file (the Javad stream of the example has MSM4 observations only).
#
#   cd NtripClient.jl && julia --project=. ../GaiaVario.jl/tools/scan_rtcm.jl FILE.rtcm3
#
# Data set of the GaiaVario examples (the data are not part of this repository, cite the record):
#   Kudlacik, I., Kaplon, J., Nastase, E. I., Tiganescu, A., Kaczmarek, A., Elias, P., Karakonstantis, A.,
#   Zuliani, D., Nistor, S. (2026). High-Rate GNSS and Multi-Sensor Benchmark Dataset from Earthquake
#   Simulation Experiments (v1). Zenodo. https://doi.org/10.5281/zenodo.18847292  (CC BY 4.0)
using GaiaSeis
path = ARGS[1]
dec = GaiaSeis.RtcmDecoder()
nep = Dict{Char,Int}(); neph = Dict{Char,Int}(); first_sod = Dict{Char,Float64}(); last_sod = Dict{Char,Float64}()
sta = nothing; desc = nothing; codes = Dict{Char,Set{String}}()
open(path) do io
    while !eof(io)
        buf = read(io, 1 << 20)
        for ev in GaiaSeis.feed_bytes!(dec, buf; stream = :obs)
            if ev isa GaiaSeis.EvEpoch
                s = ev.sys; nep[s] = get(nep, s, 0) + 1
                haskey(first_sod, s) || (first_sod[s] = ev.epoch.sod); last_sod[s] = ev.epoch.sod
                cs = get!(codes, s, Set{String}()); for o in ev.epoch.obs; push!(cs, o.rinex_code); end
            elseif ev isa GaiaSeis.EvEphemeris
                neph[ev.record.sys] = get(neph, ev.record.sys, 0) + 1
            elseif ev isa GaiaSeis.EvStation
                global sta = ev.coord
            elseif ev isa GaiaSeis.EvDescriptor
                global desc = ev.desc
            end
        end
    end
end
println("bytes seen ", dec.bytes_seen, "  crc failed ", dec.frames_crc_failed)
println("message counts: ", sort(collect(dec.msg_counts)))
println("MSM epochs per system: ", nep)
for s in keys(first_sod); println("  ", s, " sod ", round(first_sod[s]; digits=1), " .. ", round(last_sod[s]; digits=1), "  codes ", sort(collect(codes[s]))); end
println("ephemeris messages: ", neph, "  distinct sats in decoder: ", length(dec.eph))
println("station: ", sta)
println("descriptor: ", desc)
