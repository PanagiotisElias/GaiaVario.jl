# Time. The engine runs on one continuous scale, GPS time, kept as a week
# number and seconds of week. Differences between two GpsTime values are exact
# Float64 seconds, and every model evaluates time only through such
# differences (t - toe, t - t_prev, transit time applied after the reference
# is removed). There is no seconds-of-day arithmetic, so a day boundary cannot
# wrap, and no absolute Float64 second count, whose 0.24 microsecond
# resolution at 1.45e9 s would quantise satellite positions at 0.1 mm.

struct GpsTime
    week::Int
    tow::Float64
    function GpsTime(week::Integer, tow::Real)
        w = Int(week); s = Float64(tow)
        if !(0.0 <= s < GPS_WEEK_S)
            dw = floor(Int, s / GPS_WEEK_S)
            w += dw; s -= dw * GPS_WEEK_S
        end
        new(w, s)
    end
end

Base.:-(a::GpsTime, b::GpsTime) = (a.week - b.week) * GPS_WEEK_S + (a.tow - b.tow)
Base.:+(a::GpsTime, s::Real) = GpsTime(a.week, a.tow + s)
Base.:-(a::GpsTime, s::Real) = GpsTime(a.week, a.tow - s)
Base.isless(a::GpsTime, b::GpsTime) = a.week < b.week || (a.week == b.week && a.tow < b.tow)
Base.:(==)(a::GpsTime, b::GpsTime) = a.week == b.week && a.tow == b.tow
Base.show(io::IO, t::GpsTime) = print(io, "GpsTime(", t.week, ", ", t.tow, ")")

"""GPS time from week number and seconds of week."""
gps_time(week::Integer, tow::Real) = GpsTime(week, tow)

"""(week, seconds of week)."""
gps_week_tow(t::GpsTime) = (t.week, t.tow)

"""Seconds since the GPS epoch as one Float64 (for low-precision uses only)."""
to_seconds(t::GpsTime) = t.week * GPS_WEEK_S + t.tow

"""dt wrapped into (-302400, 302400] s (IS-GPS-200 20.3.3.4.3)."""
@inline function wrap_week(dt::Real)
    dt > GPS_WEEK_S / 2 && return dt - GPS_WEEK_S
    dt <= -GPS_WEEK_S / 2 && return dt + GPS_WEEK_S
    return dt
end

"""dt wrapped into (-43200, 43200] s, for references given as seconds of day."""
@inline function wrap_day(dt::Real)
    dt > DAY_S / 2 && return dt - DAY_S
    dt <= -DAY_S / 2 && return dt + DAY_S
    return dt
end

"""Julian date (GPS time taken as TT: the 1-minute difference is irrelevant
for Sun/Moon positions and sidereal time at the centimetre level)."""
julian_date(t::GpsTime) = 2_444_244.5 + t.week * 7 + t.tow / DAY_S

"""Day of year (1-based, fractional), ignoring leap seconds."""
function day_of_year(t::GpsTime)
    jd = julian_date(t)
    # Fliegel–Van Flandern inverse
    j = floor(Int, jd + 0.5)
    f = jd + 0.5 - j
    l = j + 68569
    n = div(4 * l, 146097)
    l = l - div(146097 * n + 3, 4)
    i = div(4000 * (l + 1), 1461001)
    l = l - div(1461 * i, 4) + 31
    jm = div(80 * l, 2447)
    day = l - div(2447 * jm, 80)
    l = div(jm, 11)
    month = jm + 2 - 12 * l
    year = 100 * (n - 49) + i + l
    leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
    cum = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
    doy = cum[month] + day + (leap && month > 2 ? 1 : 0)
    return doy + f
end

"""Greenwich mean sidereal time [rad] (IAU 1982 expression)."""
function gmst(t::GpsTime)
    jd = julian_date(t)
    T = (jd - 2_451_545.0) / 36525.0
    s = 67310.54841 + (876600.0 * 3600 + 8640184.812866) * T + 0.093104 * T^2 - 6.2e-6 * T^3
    return mod(s / 240.0, 360.0) * π / 180
end
