# RINEX 3 navigation file (mixed GNSS) -> records in the GaiaSeis decoder's
# 39-column layout (see NtripClient.jl/src/RTCM/messages_ephemeris.jl), so that
# both engines can be driven from a daily broadcast file when the RTCM stream
# carries no ephemeris messages. Written from the RINEX 3.05 specification.
# Expects `GV` (GaiaVario module) and Dates to be loaded by the includer.
#
#   recs = rinex_nav_rows(path)   # Vector of (sys, prn, toc::GpsTime, row::Vector{Float64}), sorted by toc

function _rn_num(l::AbstractString, k::Int)
    a = 4 + 19 * (k - 1) + 1; b = 4 + 19 * k
    a > lastindex(l) && return NaN
    s = strip(l[a:min(b, lastindex(l))])
    isempty(s) && return NaN
    return parse(Float64, replace(replace(s, 'D' => 'e'), 'd' => 'e'))
end

function _rn_gpstime(y, mo, d, h, mi, s)
    days = (Date(y, mo, d) - Date(1980, 1, 6)).value
    return GV.GpsTime(div(days, 7), mod(days, 7) * 86400.0 + h * 3600.0 + mi * 60.0 + s)
end

function rinex_nav_rows(path::AbstractString)
    lines = readlines(path)
    i = findfirst(l -> occursin("END OF HEADER", l), lines)
    i === nothing && error("no RINEX header end in $path")
    i += 1
    out = Tuple{Char,Int,GV.GpsTime,Vector{Float64}}[]
    while i <= length(lines)
        l = lines[i]
        (isempty(strip(l)) || !isletter(l[1])) && (i += 1; continue)
        sys = l[1]
        nl = sys in ('G', 'E', 'C', 'J', 'I') ? 8 : sys == 'R' ? 4 : sys == 'S' ? 4 : 8
        (i + nl - 1 > length(lines)) && break
        blk = lines[i:i+nl-1]; i += nl
        sys in ('G', 'E', 'C', 'R') || continue
        prn = parse(Int, l[2:3])
        y = parse(Int, l[5:8]); mo = parse(Int, l[10:11]); d = parse(Int, l[13:14])
        h = parse(Int, l[16:17]); mi = parse(Int, l[19:20]); s = parse(Int, l[22:23])
        # the clock fields start at column 24, i.e. the first value sits at offset 23 instead of 4
        f(k) = parse(Float64, replace(replace(strip(l[24+19*(k-1):min(23+19*k, lastindex(l))]), 'D' => 'e'), 'd' => 'e'))
        c1, c2, c3 = f(1), f(2), f(3)
        toc = _rn_gpstime(y, mo, d, h, mi, s)
        row = fill(NaN, 39); row[1] = prn
        if sys == 'R'
            health = _rn_num(blk[2], 4)
            (isnan(health) || health != 0) && continue
            toc = toc + 18.0                                   # UTC -> GPS
            row[2] = _rn_num(blk[3], 4)                        # FDMA channel
            row[8] = c1; row[9] = c2                           # -tau_n, +gamma_n
            row[11] = _rn_num(blk[2], 1); row[12] = _rn_num(blk[2], 2); row[13] = _rn_num(blk[2], 3)
            row[15] = _rn_num(blk[3], 1); row[16] = _rn_num(blk[3], 2); row[17] = _rn_num(blk[3], 3)
            row[19] = _rn_num(blk[4], 1); row[20] = _rn_num(blk[4], 2); row[21] = _rn_num(blk[4], 3)
            row[39] = mod(toc.tow, 86400.0)
        else
            health = _rn_num(blk[7], 2)
            (isnan(health) || health != 0) && continue
            row[8] = c1; row[9] = c2; row[10] = c3
            row[11] = _rn_num(blk[2], 1); row[12] = _rn_num(blk[2], 2); row[13] = _rn_num(blk[2], 3); row[14] = _rn_num(blk[2], 4)
            row[15] = _rn_num(blk[3], 1); row[16] = _rn_num(blk[3], 2); row[17] = _rn_num(blk[3], 3); row[18] = _rn_num(blk[3], 4)
            row[19] = _rn_num(blk[4], 1); row[20] = _rn_num(blk[4], 2); row[21] = _rn_num(blk[4], 3); row[22] = _rn_num(blk[4], 4)
            row[23] = _rn_num(blk[5], 1); row[24] = _rn_num(blk[5], 2); row[25] = _rn_num(blk[5], 3); row[26] = _rn_num(blk[5], 4)
            row[27] = _rn_num(blk[6], 1)
            if sys == 'G'
                row[2] = _rn_num(blk[6], 3)                    # GPS week
                row[39] = mod(toc.tow, 86400.0)
            elseif sys == 'E'
                row[39] = mod(toc.tow, 86400.0)                # row[2] stays NaN as in the decoder
            else                                               # BeiDou: toc is BDT
                row[2] = _rn_num(blk[6], 3) + 1356.0
                toc = toc + 14.0                               # BDT -> GPS
                row[39] = mod(toc.tow - 14.0, 86400.0) + 14.0  # decoder convention
            end
        end
        push!(out, (sys, prn, toc, row))
    end
    sort!(out; by = r -> (r[3].week, r[3].tow))
    return out
end

"""Per-satellite schedule: `advance!(sched, t)` returns the records whose toc has
been reached since the previous call (latest per satellite)."""
mutable struct NavSchedule
    recs::Vector{Tuple{Char,Int,GV.GpsTime,Vector{Float64}}}
    next::Int
end
NavSchedule(recs) = NavSchedule(recs, 1)

function advance!(s::NavSchedule, t::GV.GpsTime)
    fresh = Dict{Tuple{Char,Int},Tuple{Char,Int,GV.GpsTime,Vector{Float64}}}()
    while s.next <= length(s.recs) && !(t < s.recs[s.next][3])
        r = s.recs[s.next]; fresh[(r[1], r[2])] = r; s.next += 1
    end
    return collect(values(fresh))
end
