# Record broadcast ephemerides from an NTRIP ephemeris stream into a TSV file
# that GaiaVario's tests and the synthetic-truth harness can replay offline.
#
#   julia --project=../NtripClient.jl tools/record_brdc.jl [seconds] [host] [port] [mount] [out.tsv]
#
# Defaults: 90 s from the local ephemeris relay 127.0.0.1:2198 mount BCEP00BKG0
# (the relay ignores credentials). The latest record per (system, PRN) is kept
# and written as: sys prn col1 ... col39 (the decoder's row layout).

using GaiaSeis
using Sockets
using Dates

secs = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 90.0
host = length(ARGS) >= 2 ? ARGS[2] : "127.0.0.1"
port = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2198
mount = length(ARGS) >= 4 ? ARGS[4] : "BCEP00BKG0"
out = length(ARGS) >= 5 ? ARGS[5] :
      joinpath(@__DIR__, "..", "test", "data", "brdc_" * Dates.format(now(UTC), "yyyymmdd_HHMM") * ".tsv")

sock = connect(host, port)
write(sock, "GET /$mount HTTP/1.1\r\nHost: $host\r\nNtrip-Version: Ntrip/2.0\r\n" *
            "User-Agent: NTRIP GaiaVario/0.1\r\nAuthorization: Basic eDp4\r\nConnection: close\r\n\r\n")
# skip the response header
hdr = UInt8[]
while isopen(sock)
    push!(hdr, read(sock, UInt8))
    length(hdr) >= 4 && hdr[end-3:end] == UInt8['\r', '\n', '\r', '\n'] && break
end
println(stderr, "header: ", strip(split(String(copy(hdr)), "\r\n")[1]))

dec = GaiaSeis.RtcmDecoder()
records = Dict{Tuple{Char,Int},Vector{Float64}}()
counts = Dict{Char,Int}()
t_end = time() + secs
while time() < t_end && isopen(sock)
    buf = readavailable(sock)
    isempty(buf) && (sleep(0.05); continue)
    for ev in GaiaSeis.feed_bytes!(dec, buf; stream = :eph)
        ev isa GaiaSeis.EvEphemeris || continue
        r = ev.record
        r.valid || continue
        records[(r.sys, r.prn)] = copy(r.row)
        counts[r.sys] = get(counts, r.sys, 0) + 1
    end
end
close(sock)

mkpath(dirname(out))
open(out, "w") do io
    println(io, "# broadcast ephemerides recorded ", Dates.format(now(UTC), "yyyy-mm-dd HH:MM:SS"), " UTC from ", host, ":", port, "/", mount)
    println(io, "# sys prn col1..col39 (GaiaSeis decoder row layout; see NtripClient.jl/src/RTCM/messages_ephemeris.jl)")
    for k in sort(collect(keys(records)))
        print(io, k[1], '\t', k[2])
        for v in records[k]
            print(io, '\t', isnan(v) ? "nan" : repr(v))
        end
        println(io)
    end
end
println(stderr, "messages per system: ", counts)
println(stderr, "records: ", length(records), " -> ", out)
