# Satellite position and clock from broadcast records.
#   Kepler sets: IS-GPS-200 table 20-IV; Galileo OS SIS ICD 5.1.1; BeiDou SIS
#   ICD (B1I) 5.2.4.12 including the GEO rotation.
#   GLONASS: ICD 5.1 appendix 3, RK4 integration of the PZ-90 equations.
# Times are GPS seconds; the record times already are too (types.jl).

function _kepler_E(M::Float64, e::Float64)
    E = e < 0.8 ? M : π
    for _ in 1:30
        dE = (E - e * sin(E) - M) / (1 - e * cos(E))
        E -= dE
        abs(dE) < 1e-13 && break
    end
    return E
end

_sys_consts(sys::Char) = sys == 'G' ? (MU_GPS, OMEGA_E_GPS, F_REL_GPS) :
                        sys == 'E' ? (MU_GAL, OMEGA_E_GAL, F_REL_GAL) :
                                     (MU_BDS, OMEGA_E_BDS, F_REL_BDS)

is_bds_geo(prn::Int) = prn <= 5 || prn >= 59

"""ECEF position [m] and relativistic clock term [s] of a Kepler record at
GPS time t + dt (dt is applied after the reference time is removed, so a
transit time of 70 ms costs no precision)."""
function kepler_position(eph::KeplerEph, t::GpsTime, dt::Real = 0.0)
    mu, we, Frel = _sys_consts(eph.sys)
    A = eph.sqrtA^2
    n = sqrt(mu / A^3) + eph.dn
    tk = wrap_week(t - eph.toe) + dt
    M = eph.M0 + n * tk
    E = _kepler_E(M, eph.e)
    sE, cE = sincos(E)
    nu = atan(sqrt(1 - eph.e^2) * sE, cE - eph.e)
    phi = nu + eph.omega
    s2, c2 = sincos(2phi)
    u = phi + eph.cus * s2 + eph.cuc * c2
    r = A * (1 - eph.e * cE) + eph.crs * s2 + eph.crc * c2
    i = eph.i0 + eph.idot * tk + eph.cis * s2 + eph.cic * c2
    xo = r * cos(u); yo = r * sin(u)
    rel = Frel * eph.e * eph.sqrtA * sE
    if eph.sys == 'C' && is_bds_geo(eph.prn)
        # GEO: inertial-like frame with Omega excluding Earth rotation, then
        # rotate by Rz(we*tk) Rx(-5 deg) (BDS ICD 5.2.4.12)
        Om = eph.Omega0 + eph.Omegadot * tk - we * eph.toe_sow
        sO, cO = sincos(Om); si, ci = sincos(i)
        xg = xo * cO - yo * ci * sO
        yg = xo * sO + yo * ci * cO
        zg = yo * si
        a5 = -5.0 * π / 180
        s5, c5 = sincos(a5)
        # Rx(-5 deg)
        x1 = xg; y1 = c5 * yg + s5 * zg; z1 = -s5 * yg + c5 * zg
        # Rz(we*tk)
        sw, cw = sincos(we * tk)
        return (cw * x1 + sw * y1, -sw * x1 + cw * y1, z1), rel
    end
    Om = eph.Omega0 + (eph.Omegadot - we) * tk - we * eph.toe_sow
    sO, cO = sincos(Om); si, ci = sincos(i)
    return (xo * cO - yo * ci * sO, xo * sO + yo * ci * cO, yo * si), rel
end

"""Satellite clock offset [s] including the relativistic term (TGD omitted: it
cancels in the ionosphere-free pair and is constant in a single-frequency pair)."""
function kepler_clock(eph::KeplerEph, t::GpsTime, rel::Float64, dt::Real = 0.0)
    dtc = wrap_week(t - eph.toc) + dt
    return eph.af0 + eph.af1 * dtc + eph.af2 * dtc^2 + rel
end

# ---- GLONASS -----------------------------------------------------------------

@inline function _glo_accel(x, y, z, vx, vy, ax, ay, az)
    r2 = x * x + y * y + z * z
    r = sqrt(r2)
    mur3 = MU_GLO / (r2 * r)
    j = 1.5 * GLO_J2 * MU_GLO * GLO_AE^2 / (r2 * r2 * r)
    z2r2 = 5 * z * z / r2
    we = OMEGA_E_GLO
    fx = -mur3 * x - j * x * (1 - z2r2) + we^2 * x + 2we * vy + ax
    fy = -mur3 * y - j * y * (1 - z2r2) + we^2 * y - 2we * vx + ay
    fz = -mur3 * z - j * z * (3 - z2r2) + az
    return fx, fy, fz
end

"""PZ-90 position [m] at GPS time t by RK4 from tb (step <= 60 s, ICD 5.1 A3.1.2)."""
function glonass_position(eph::GlonassEph, t::GpsTime, dtau::Real = 0.0)
    x, y, z = eph.pos; vx, vy, vz = eph.vel; ax, ay, az = eph.acc
    dt = wrap_day(t - eph.tb) + dtau
    n = max(1, ceil(Int, abs(dt) / 60.0))
    h = dt / n
    for _ in 1:n
        k1 = _glo_accel(x, y, z, vx, vy, ax, ay, az)
        x2 = x + 0.5h * vx; y2 = y + 0.5h * vy; z2 = z + 0.5h * vz
        vx2 = vx + 0.5h * k1[1]; vy2 = vy + 0.5h * k1[2]; vz2 = vz + 0.5h * k1[3]
        k2 = _glo_accel(x2, y2, z2, vx2, vy2, ax, ay, az)
        x3 = x + 0.5h * vx2; y3 = y + 0.5h * vy2; z3 = z + 0.5h * vz2
        vx3 = vx + 0.5h * k2[1]; vy3 = vy + 0.5h * k2[2]; vz3 = vz + 0.5h * k2[3]
        k3 = _glo_accel(x3, y3, z3, vx3, vy3, ax, ay, az)
        x4 = x + h * vx3; y4 = y + h * vy3; z4 = z + h * vz3
        vx4 = vx + h * k3[1]; vy4 = vy + h * k3[2]; vz4 = vz + h * k3[3]
        k4 = _glo_accel(x4, y4, z4, vx4, vy4, ax, ay, az)
        x += h / 6 * (vx + 2vx2 + 2vx3 + vx4)
        y += h / 6 * (vy + 2vy2 + 2vy3 + vy4)
        z += h / 6 * (vz + 2vz2 + 2vz3 + vz4)
        vx += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1])
        vy += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2])
        vz += h / 6 * (k1[3] + 2k2[3] + 2k3[3] + k4[3])
    end
    return (x, y, z)
end

glonass_clock(eph::GlonassEph, t::GpsTime, dtau::Real = 0.0) = eph.tau_neg + eph.gamma * (wrap_day(t - eph.tb) + dtau)

# ---- Unified access ------------------------------------------------------------

"""Satellite ECEF position [m] and clock offset [s] at GPS time t + dt."""
function sat_state(eph::KeplerEph, t::GpsTime, dt::Real = 0.0)
    p, rel = kepler_position(eph, t, dt)
    return p, kepler_clock(eph, t, rel, dt)
end
sat_state(eph::GlonassEph, t::GpsTime, dt::Real = 0.0) = glonass_position(eph, t, dt), glonass_clock(eph, t, dt)

"""Satellite velocity [m/s] by central difference (used for nothing critical)."""
function sat_velocity(eph::AnyEph, t::GpsTime; h = 0.5)
    p1, _ = sat_state(eph, t, -h); p2, _ = sat_state(eph, t, h)
    return ((p2[1] - p1[1]) / 2h, (p2[2] - p1[2]) / 2h, (p2[3] - p1[3]) / 2h)
end

"""Geometric range with signal transit time and Earth rotation (Sagnac):
returns range [m], the rotated satellite position used, and the satellite
clock at transmit time [s]. `t_rx` is the receiver time of reception."""
function sagnac_range(eph::AnyEph, t_rx::GpsTime, rec::NTuple{3,Float64})
    tau = 0.075
    we = eph isa GlonassEph ? OMEGA_E_GLO : _sys_consts(eph.sys)[2]
    local ps, clk, rho
    for _ in 1:3
        ps, clk = sat_state(eph, t_rx, -tau)
        a = we * tau
        sa, ca = sincos(a)
        xr = ca * ps[1] + sa * ps[2]; yr = -sa * ps[1] + ca * ps[2]
        ps = (xr, yr, ps[3])
        rho = sqrt((ps[1] - rec[1])^2 + (ps[2] - rec[2])^2 + (ps[3] - rec[3])^2)
        tau = rho / C_LIGHT
    end
    return rho, ps, clk
end
# The satellite clock offset (up to 1 ms) is not applied to the transmit time:
# the resulting position error is below 1 mm and constant between the two
# epochs of a pair to far better than the phase noise.
