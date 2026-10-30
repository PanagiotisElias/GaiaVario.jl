# Physical and constellation constants. Sources: IS-GPS-200 (GPS), Galileo OS
# SIS ICD, BeiDou SIS ICD (B1I), GLONASS ICD 5.1 (PZ-90), IERS Conventions 2010.

const C_LIGHT = 299_792_458.0            # m/s

# Gravitational parameter and Earth rotation rate per constellation [m^3/s^2], [rad/s]
const MU_GPS = 3.986005e14
const MU_GAL = 3.986004418e14
const MU_BDS = 3.986004418e14
const MU_GLO = 3.986004418e14            # PZ-90.11 (ICD 5.1 uses 398600.4418 km^3/s^2)
const OMEGA_E_GPS = 7.2921151467e-5
const OMEGA_E_GAL = 7.2921151467e-5
const OMEGA_E_BDS = 7.292115e-5
const OMEGA_E_GLO = 7.292115e-5

# Relativistic clock constant F = -2 sqrt(mu)/c^2 [s/sqrt(m)]
const F_REL_GPS = -4.442807633e-10
const F_REL_GAL = -4.442807309e-10
const F_REL_BDS = -4.442807309e-10

# GLONASS PZ-90 second zonal harmonic and equatorial radius
const GLO_AE = 6_378_136.0
const GLO_J2 = 1.0826257e-3

# WGS-84 ellipsoid
const WGS84_A = 6_378_137.0
const WGS84_F = 1.0 / 298.257223563
const WGS84_E2 = WGS84_F * (2 - WGS84_F)

# Time scale offsets
const BDT_MINUS_GPST_S = -14.0           # BDT = GPST - 14 s
const GPS_WEEK_S = 604_800.0
const DAY_S = 86_400.0
const BDT_WEEK0_GPS_WEEK = 1356

# Carrier frequencies [Hz] by (system, frequency slot as used by the RTCM
# decoder: 1 = L1/E1/B1, 2 = L2/B1I(2I), 5 = L5/E5a/B2a, 6 = E6/B3I, 7 = E5b/B2b,
# 8 = E5ab). GLONASS slots take the FDMA channel number k.
const FREQ_GPS = Dict(1 => 1575.42e6, 2 => 1227.60e6, 5 => 1176.45e6)
const FREQ_GAL = Dict(1 => 1575.42e6, 5 => 1176.45e6, 6 => 1278.75e6, 7 => 1207.14e6, 8 => 1191.795e6)
const FREQ_BDS = Dict(1 => 1575.42e6, 2 => 1561.098e6, 5 => 1176.45e6, 6 => 1268.52e6, 7 => 1207.14e6)
glo_freq(slot::Int, k::Int) = slot == 1 ? 1602.0e6 + k * 562.5e3 :
                              slot == 2 ? 1246.0e6 + k * 437.5e3 : NaN

function carrier_freq(sys::Char, slot::Int, k::Int = 0)
    sys == 'G' && return get(FREQ_GPS, slot, NaN)
    sys == 'E' && return get(FREQ_GAL, slot, NaN)
    sys == 'C' && return get(FREQ_BDS, slot, NaN)
    sys == 'R' && return glo_freq(slot, k)
    return NaN
end

# Love and Shida numbers of degree 2 (IERS 2010, nominal values)
const LOVE_H2 = 0.6078
const SHIDA_L2 = 0.0847
const GM_SUN = 1.32712440018e20
const GM_MOON = 4.902800066e12
const GM_EARTH = 3.986004418e14
const AU_M = 149_597_870_700.0
