# Adapters from the GaiaSeis RTCM decoder's representations. The decoder stores
# each broadcast record in a 39-column row (see NtripClient.jl
# src/RTCM/messages_ephemeris.jl for the column assignment); these functions
# rebuild typed records with continuous GPS times.
#
#   Kepler rows: [1] prn, [2] week (GPS: absolute; BDS: GPS week; Galileo: NaN),
#   [8..10] af0 af1 af2, [11] iode, [12] crs, [13] dn, [14] M0, [15] cuc, [16] e,
#   [17] cus, [18] sqrtA, [19] toe (seconds of week; BDS in BDT), [20] cic,
#   [21] Omega0, [22] cis, [23] i0, [24] crc, [25] omega, [26] Omegadot,
#   [27] idot, [39] toc as GPS-scale seconds of day.
#   GLONASS rows: [2] channel k, [8] -tau_n, [9] gamma_n, [11,15,19] x y z km,
#   [12,16,20] velocities km/s, [13,17,21] accelerations km/s^2,
#   [39] tb as GPS-scale seconds of day.

"""Rebuild a continuous GPS time from a seconds-of-day reference and a GPS
time near it (the record is valid within half a day of its reference)."""
function _sod_to_gps(sod::Float64, t_near::GpsTime)
    cand = GpsTime(t_near.week, floor(t_near.tow / DAY_S) * DAY_S + sod)
    return cand + DAY_S * round((t_near - cand) / DAY_S)
end

"""Seconds-of-week reference to a continuous GPS time near `t_near`."""
function _sow_to_gps(sow::Float64, t_near::GpsTime)
    cand = GpsTime(t_near.week, sow)
    return cand + GPS_WEEK_S * round((t_near - cand) / GPS_WEEK_S)
end

"""Typed record from a decoder row. `t_near` is the current GPS time, used to
place seconds-of-day and seconds-of-week references on the continuous scale."""
function eph_from_row(sys::Char, prn::Int, row::AbstractVector{Float64}, t_near::GpsTime)
    if sys == 'R'
        tb = _sod_to_gps(row[39], t_near)
        return GlonassEph(prn, Int(round(row[2])), tb,
                          (row[11] * 1e3, row[15] * 1e3, row[19] * 1e3),
                          (row[12] * 1e3, row[16] * 1e3, row[20] * 1e3),
                          (row[13] * 1e3, row[17] * 1e3, row[21] * 1e3),
                          row[8], row[9])
    end
    toe_sow = row[19]
    # BeiDou: toe is BDT seconds of week; the GPS-scale instant is 14 s later.
    toe = sys == 'C' ? _sow_to_gps(toe_sow - BDT_MINUS_GPST_S, t_near) : _sow_to_gps(toe_sow, t_near)
    toc = _sod_to_gps(row[39], t_near)
    week = isnan(row[2]) ? gps_week_tow(toe)[1] : Int(round(row[2]))
    return KeplerEph(sys, prn, Int(round(row[11])), toe, toc,
                     row[18], row[16], row[23], row[21], row[25], row[14], row[13],
                     row[27], row[26], row[15], row[17], row[24], row[12], row[20], row[22],
                     row[8], row[9], row[10], toe_sow, week)
end

"""Load a TSV written by `tools/record_brdc.jl` into a NavStore. `t_near` is a
GPS time within half a day of the recording."""
function load_brdc_tsv(path::AbstractString, t_near::GpsTime)
    nav = NavStore()
    for line in eachline(path)
        (isempty(line) || startswith(line, '#')) && continue
        f = split(line, '\t')
        length(f) < 41 && continue
        sys = f[1][1]; prn = parse(Int, f[2])
        row = [v == "nan" ? NaN : parse(Float64, v) for v in f[3:41]]
        push!(nav, eph_from_row(sys, prn, row, t_near))
    end
    return nav
end

# In kepler_position the Earth-rotation term uses toe_sow of the constellation's
# own week count; tk is computed on the continuous scale. For BeiDou tk must be
# BDT-based too: toe stored is the GPS-scale instant, and t is GPS, so t - toe
# is the same in both scales.
