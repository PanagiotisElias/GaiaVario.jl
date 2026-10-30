# Observation models other than the orbit. Each is written from the cited
# source; in a variometric pair only the change between two epochs matters,
# which relaxes the accuracy needed from every term here.
#   geodesy: standard closed forms
#   troposphere: Saastamoinen zenith delays with a standard atmosphere, Niell
#     mapping functions (Niell 1996, JGR 101, 3227). The pair difference of the
#     troposphere is < 1 mm per 0.1 s even at 10 deg elevation, so a 1 %
#     mapping-function error is far below the phase noise.
#   solid Earth tide: degree-2 term of IERS Conventions 2010 eq. 7.5 with
#     nominal h2, l2; Sun/Moon from low-precision series (Montenbruck & Gill,
#     Satellite Orbits, 3.3.2).
#   phase wind-up: Wu et al. 1993, Manuscripta Geodaetica 18, 91, with the
#     nominal yaw attitude (body z to Earth, y perpendicular to the Sun).

# ---- Geodesy -------------------------------------------------------------------

"""Geodetic latitude, longitude [rad] and height [m] of an ECEF point."""
function ecef2geodetic(p::NTuple{3,Float64})
    x, y, z = p
    lon = atan(y, x)
    r = hypot(x, y)
    lat = atan(z, r * (1 - WGS84_E2))
    h = 0.0
    for _ in 1:6
        N = WGS84_A / sqrt(1 - WGS84_E2 * sin(lat)^2)
        h = r / cos(lat) - N
        lat = atan(z, r * (1 - WGS84_E2 * N / (N + h)))
    end
    return lat, lon, h
end

"""Rows: unit vectors of north, east, up in ECEF."""
function neu_frame(lat::Float64, lon::Float64)
    sl, cl = sincos(lat); so, co = sincos(lon)
    n = (-sl * co, -sl * so, cl)
    e = (-so, co, 0.0)
    u = (cl * co, cl * so, sl)
    return n, e, u
end

@inline dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline norm3(a) = sqrt(dot3(a, a))
@inline unit3(a) = (n = norm3(a); (a[1] / n, a[2] / n, a[3] / n))
@inline sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline add3(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline scale3(a, s) = (a[1] * s, a[2] * s, a[3] * s)

"""Elevation and azimuth [rad] of a satellite seen from the receiver."""
function elev_azim(frame, rec::NTuple{3,Float64}, sat::NTuple{3,Float64})
    n, e, u = frame
    d = unit3(sub3(sat, rec))
    dn, de, du = dot3(d, n), dot3(d, e), dot3(d, u)
    return asin(clamp(du, -1.0, 1.0)), atan(de, dn)
end

# ---- Troposphere ---------------------------------------------------------------

"""Zenith hydrostatic and wet delay [m] from a standard atmosphere (Saastamoinen)."""
function zenith_delays(lat::Float64, h::Float64)
    hk = clamp(h, -500.0, 9000.0)
    P = 1013.25 * (1 - 2.2557e-5 * hk)^5.2568
    T = 15.0 - 0.0065 * hk + 273.15
    e = 0.5 * 6.108 * exp((17.15 * T - 4684.0) / (T - 38.45))   # 50 % humidity
    zhd = 0.0022768 * P / (1 - 0.00266 * cos(2lat) - 0.00028e-3 * hk)
    zwd = 0.002277 * (1255.0 / T + 0.05) * e
    return zhd, zwd
end

# Niell 1996 table coefficients (latitude 15, 30, 45, 60, 75 deg)
const _NMF_LAT = (15.0, 30.0, 45.0, 60.0, 75.0)
const _NMF_HA = (1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3)
const _NMF_HB = (2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3)
const _NMF_HC = (62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3)
const _NMF_HA_AMP = (0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5)
const _NMF_HB_AMP = (0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5)
const _NMF_HC_AMP = (0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5)
const _NMF_WA = (5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4)
const _NMF_WB = (1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3)
const _NMF_WC = (4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2)

@inline _nmf_frac(a, b, c, sinel) = (1 + a / (1 + b / (1 + c))) / (sinel + a / (sinel + b / (sinel + c)))

function _nmf_interp(tab, lat_deg)
    la = abs(lat_deg)
    la <= _NMF_LAT[1] && return tab[1]
    la >= _NMF_LAT[5] && return tab[5]
    i = findlast(<=(la), _NMF_LAT)
    f = (la - _NMF_LAT[i]) / 15.0
    return tab[i] + f * (tab[i + 1] - tab[i])
end

"""Niell hydrostatic and wet mapping functions at elevation `el` [rad]."""
function niell_mf(el::Float64, lat::Float64, h::Float64, doy::Float64)
    lat_deg = lat * 180 / π
    phase = cos(2π * (doy - 28.0 + (lat < 0 ? 182.625 : 0.0)) / 365.25)
    a = _nmf_interp(_NMF_HA, lat_deg) - _nmf_interp(_NMF_HA_AMP, lat_deg) * phase
    b = _nmf_interp(_NMF_HB, lat_deg) - _nmf_interp(_NMF_HB_AMP, lat_deg) * phase
    c = _nmf_interp(_NMF_HC, lat_deg) - _nmf_interp(_NMF_HC_AMP, lat_deg) * phase
    s = sin(el)
    mh = _nmf_frac(a, b, c, s)
    # height correction
    mh += (1 / s - _nmf_frac(2.53e-5, 5.49e-3, 1.14e-3, s)) * h / 1000
    mw = _nmf_frac(_nmf_interp(_NMF_WA, lat_deg), _nmf_interp(_NMF_WB, lat_deg),
                   _nmf_interp(_NMF_WC, lat_deg), s)
    return mh, mw
end

"""Slant tropospheric delay [m]."""
function tropo_delay(el::Float64, lat::Float64, h::Float64, doy::Float64, zhd::Float64, zwd::Float64)
    mh, mw = niell_mf(max(el, 0.5 * π / 180), lat, h, doy)
    return mh * zhd + mw * zwd
end

# ---- Sun and Moon ----------------------------------------------------------------

"""Sun position [m] in ECEF (mean equator and equinox, rotated by GMST)."""
function sun_ecef(t::GpsTime)
    T = (julian_date(t) - 2_451_545.0) / 36525.0
    M = deg2rad(mod(357.5256 + 35999.049 * T, 360.0))
    lam = deg2rad(mod(282.9400 + rad2deg(M) + (6892 / 3600) * sin(M) + (72 / 3600) * sin(2M), 360.0))
    r = (149.619 - 2.499 * cos(M) - 0.021 * cos(2M)) * 1e9
    eps = deg2rad(23.43929111)
    xi = r * cos(lam); yi = r * sin(lam) * cos(eps); zi = r * sin(lam) * sin(eps)
    return _eci2ecef((xi, yi, zi), t)
end

"""Moon position [m] in ECEF (Montenbruck & Gill low-precision series, ~1 arcmin)."""
function moon_ecef(t::GpsTime)
    T = (julian_date(t) - 2_451_545.0) / 36525.0
    L0 = mod(0.606433 + 1336.851344 * T, 1.0) * 2π
    l = mod(0.374897 + 1325.552410 * T, 1.0) * 2π
    lp = mod(0.993133 + 99.997361 * T, 1.0) * 2π
    D = mod(0.827361 + 1236.853086 * T, 1.0) * 2π
    F = mod(0.259086 + 1342.227825 * T, 1.0) * 2π
    dL = 22640 * sin(l) - 4586 * sin(l - 2D) + 2370 * sin(2D) + 769 * sin(2l) -
         668 * sin(lp) - 412 * sin(2F) - 212 * sin(2l - 2D) - 206 * sin(l + lp - 2D) +
         192 * sin(l + 2D) - 165 * sin(lp - 2D) - 125 * sin(D) - 110 * sin(l + lp) +
         148 * sin(l - lp) - 55 * sin(2F - 2D)
    L = L0 + deg2rad(dL / 3600)
    S = F + deg2rad((dL + 412 * sin(2F) + 541 * sin(lp)) / 3600)
    hh = F - 2D
    N = -526 * sin(hh) + 44 * sin(l + hh) - 31 * sin(-l + hh) - 23 * sin(lp + hh) +
        11 * sin(-lp + hh) - 25 * sin(-2l + F) + 21 * sin(-l + F)
    B = deg2rad((18520 * sin(S) + N) / 3600)
    r = (385000 - 20905 * cos(l) - 3699 * cos(2D - l) - 2956 * cos(2D) - 570 * cos(2l) +
         246 * cos(2l - 2D) - 205 * cos(lp - 2D) - 171 * cos(l + 2D) - 152 * cos(l + lp - 2D)) * 1e3
    eps = deg2rad(23.43929111)
    xe = r * cos(L) * cos(B); ye = r * sin(L) * cos(B); ze = r * sin(B)
    # ecliptic -> equatorial
    xi = xe; yi = cos(eps) * ye - sin(eps) * ze; zi = sin(eps) * ye + cos(eps) * ze
    return _eci2ecef((xi, yi, zi), t)
end

function _eci2ecef(p, t)
    th = gmst(t)
    s, c = sincos(th)
    return (c * p[1] + s * p[2], -s * p[1] + c * p[2], p[3])
end

# ---- Solid Earth tide ----------------------------------------------------------------

"""Degree-2 solid Earth tide displacement [m] in ECEF at receiver `rec`."""
function solid_tide(rec::NTuple{3,Float64}, t::GpsTime)
    rhat = unit3(rec)
    R = norm3(rec)
    d = (0.0, 0.0, 0.0)
    for (gm, body) in ((GM_MOON, moon_ecef(t)), (GM_SUN, sun_ecef(t)))
        Rb = norm3(body)
        bhat = unit3(body)
        cosz = dot3(bhat, rhat)
        f = gm / GM_EARTH * R^4 / Rb^3
        radial = scale3(rhat, LOVE_H2 * (1.5 * cosz^2 - 0.5))
        trans = scale3(sub3(bhat, scale3(rhat, cosz)), 3 * SHIDA_L2 * cosz)
        d = add3(d, scale3(add3(radial, trans), f))
    end
    return d
end

# ---- Phase wind-up ----------------------------------------------------------------

"""Phase wind-up [cycles] between satellite at `sat` and receiver at `rec`,
receiver dipole along local north/east. Only differences between epochs are
used by the engine, so no unwrapping state is needed here."""
function windup_cycles(sat::NTuple{3,Float64}, rec::NTuple{3,Float64}, sun::NTuple{3,Float64}, frame)
    n, e, u = frame
    k = unit3(sub3(rec, sat))                     # satellite -> receiver
    zs = unit3(scale3(sat, -1.0))                 # body z toward Earth centre
    es = unit3(sub3(sun, sat))
    ys = unit3(cross3(zs, es))
    xs = cross3(ys, zs)
    Ds = sub3(sub3(xs, scale3(k, dot3(k, xs))), cross3(k, ys))
    Dr = add3(sub3(n, scale3(k, dot3(k, n))), cross3(k, e))
    cosw = clamp(dot3(Ds, Dr) / (norm3(Ds) * norm3(Dr)), -1.0, 1.0)
    w = acos(cosw) / 2π
    return dot3(k, cross3(Ds, Dr)) < 0 ? -w : w
end

"""Difference of wind-up between two geometries, with the half-cycle wrap removed."""
@inline function windup_diff(w_now::Float64, w_prev::Float64)
    d = w_now - w_prev
    d > 0.5 && (d -= 1.0)
    d < -0.5 && (d += 1.0)
    return d
end
