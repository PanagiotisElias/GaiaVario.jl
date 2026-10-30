# Synthetic-truth generator: carrier phases from broadcast records, a prescribed
# station motion, noise, multipath, cycle slips and ephemeris changeovers. The
# forward model is written independently of the engine's (same physics, own
# code path through the public functions) so a wrong sign in the engine shows
# up as a bias against the truth rather than cancelling.

"""Prescribed motions (NEU metres as a function of time since start)."""
motion_static(τ) = (0.0, 0.0, 0.0)
motion_step(amp = 0.05, t0 = 30.0) = τ -> τ >= t0 ? (amp, 0.0, 0.0) : (0.0, 0.0, 0.0)
motion_sine(amp = 0.02, period = 5.0) = τ -> (amp * sin(2π * τ / period), 0.0, 0.5amp * cos(2π * τ / period))
"""Ricker-like pulse: an earthquake-shaped displacement in N and U."""
motion_pulse(amp = 0.10, t0 = 60.0, width = 2.0) = τ -> begin
    a = (τ - t0) / width
    g = (1 - 2a^2) * exp(-a^2)
    (amp * g, -0.6amp * g, 0.4amp * g)
end

Base.@kwdef mutable struct SynthOptions
    rate_hz::Float64 = 10.0
    duration_s::Float64 = 120.0
    noise_m::Float64 = 0.002             # phase noise per frequency at zenith [m]
    multipath_m::Float64 = 0.0           # amplitude of slow sinusoidal multipath
    multipath_period_s::Float64 = 300.0
    elev_mask_deg::Float64 = 5.0
    receiver_clock::Bool = true          # random-walk receiver clock [m]
    slips::Vector{Tuple{Char,Int,Float64}} = []     # (sys, prn, time since start)
    changeovers::Vector{Tuple{Char,Int,Float64,Float64}} = []  # (sys, prn, time, radial offset m)
    slots::Dict{Char,Tuple{Int,Int}} = Dict('G' => (1, 2), 'R' => (1, 2), 'E' => (1, 5), 'C' => (1, 7))
    seed::Int = 1
end

"""Simple deterministic RNG (xorshift) so tests are reproducible without Random."""
mutable struct XorShift; s::UInt64; end
function randu(r::XorShift)
    x = r.s; x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17; r.s = x
    return (x >> 11) * (1.0 / 9007199254740992)
end
function randn1(r::XorShift)
    u1 = max(randu(r), 1e-300); u2 = randu(r)
    return sqrt(-2log(u1)) * cos(2π * u2)
end

"""Shift a Kepler record radially by `dr` metres (a changed sqrtA), with a new IOD."""
function perturb_eph(e::KeplerEph, dr::Float64)
    A = e.sqrtA^2 + dr
    KeplerEph(e.sys, e.prn, e.iod + 1, e.toe, e.toc, sqrt(A), e.e, e.i0, e.Omega0, e.omega, e.M0,
              e.dn, e.idot, e.Omegadot, e.cuc, e.cus, e.crc, e.crs, e.cic, e.cis,
              e.af0, e.af1, e.af2, e.toe_sow, e.week)
end

"""Generate epochs. Returns (epochs, truth) where truth[i] is the NEU
displacement of epoch i relative to the start, and a NavStore that receives
the changeover records at the right times through `nav_events`
(vector of (time, record)). `t0` is the GPS time of the first epoch."""
function synthesize(nav0::NavStore, rec::NTuple{3,Float64}, t0::GpsTime, motion;
                    so::SynthOptions = SynthOptions())
    rng = XorShift(UInt64(0x9E3779B97F4A7C15) ⊻ UInt64(so.seed))
    lat, lon, h = ecef2geodetic(rec)
    frame = neu_frame(lat, lon)
    nf, ef, uf = frame
    zhd, zwd = zenith_delays(lat, h)
    n = round(Int, so.duration_s * so.rate_hz)
    epochs = Vector{Epoch}(undef, n)
    truth = Vector{NTuple{3,Float64}}(undef, n)
    nav_events = Tuple{GpsTime,AnyEph}[]
    # records the generator uses: start with nav0, apply changeovers at their times
    current = Dict{Tuple{Char,Int},AnyEph}(k => v for (k, v) in nav0.eph)
    pending = sort(so.changeovers; by = x -> x[3])
    # epoch i has tau = (i - 1) / rate, so a slip at tau lands on index round(tau*rate) + 1
    slips = Set((s, p, round(Int, τ * so.rate_hz) + 1) for (s, p, τ) in so.slips)
    lockcount = Dict{Tuple{Char,Int},Int}()
    ambig = Dict{Tuple{Char,Int,Int},Float64}()   # integer ambiguity in metres per slot
    clk = 0.0
    mp_phase = Dict{Tuple{Char,Int,Int},Float64}()
    for i in 1:n
        τ = (i - 1) / so.rate_hz
        t = t0 + τ
        # A changeover is a NEW RECORD for the engine, not a new orbit: the
        # satellite keeps flying the truth orbit (`current` is not touched);
        # the engine receives a record that differs from the truth by `dr`.
        while !isempty(pending) && pending[1][3] <= τ
            sys, prn, _, dr = popfirst!(pending)
            e = current[(sys, prn)]
            e isa KeplerEph || continue
            push!(nav_events, (t, perturb_eph(e, dr)))
        end
        dneu = motion(τ)
        truth[i] = dneu
        disp = add3(add3(scale3(nf, dneu[1]), scale3(ef, dneu[2])), scale3(uf, dneu[3]))
        tide = solid_tide(rec, t)
        pos = add3(add3(rec, disp), tide)
        so.receiver_clock && (clk += 0.02 * randn1(rng))     # 2 cm random walk per epoch
        sun = sun_ecef(t)
        doy = day_of_year(t)
        sats = SatObs[]
        for ((sys, prn), e) in current
            rho, ps, sclk = sagnac_range(e, t, pos)
            el, az = elev_azim(frame, rec, ps)
            el < deg2rad(so.elev_mask_deg) && continue
            trop = tropo_delay(el, lat, h, doy, zhd, zwd)
            w = windup_cycles(ps, rec, sun, frame)
            k = e isa GlonassEph ? e.k : 0
            s = SatObs(sys, prn)
            key = (sys, prn)
            lc = get(lockcount, key, 0) + 1
            if (sys, prn, i) in slips
                lc = 0
                for sl in so.slots[sys]
                    delete!(ambig, (sys, prn, sl))
                end
            end
            lockcount[key] = lc
            for sl in so.slots[sys]
                f = carrier_freq(sys, sl, k)
                lam = C_LIGHT / f
                amb = get!(ambig, (sys, prn, sl)) do
                    lam * round(2e6 * (randu(rng) - 0.5))
                end
                mp = 0.0
                if so.multipath_m > 0
                    ph = get!(mp_phase, (sys, prn, sl), 2π * randu(rng))
                    mp = so.multipath_m / max(sin(el), 0.2) * sin(2π * τ / so.multipath_period_s + ph)
                end
                noise = so.noise_m / max(sin(el), 0.2) * randn1(rng)
                phase = rho + clk - C_LIGHT * sclk + trop + lam * w + amb + mp + noise
                code = rho + clk - C_LIGHT * sclk + trop + 0.3 * randn1(rng)
                s.phase[sl] = phase; s.code[sl] = code
                s.lock[sl] = min(lc, 1023); s.half[sl] = false; s.cnr[sl] = 45.0
            end
            push!(sats, s)
        end
        epochs[i] = Epoch(t, sats)
    end
    return epochs, truth, nav_events
end

"""GLONASS-like records: circular PZ-90 orbits (25 508 km, 64.8 deg), state
vectors at the quarter hour before `t0`, channels k = -7..6."""
function synthetic_glonass(t0::GpsTime; n::Int = 8)
    nav = NavStore()
    r = 25_508_000.0
    v = sqrt(MU_GLO / r)
    tb = GpsTime(t0.week, floor(t0.tow / 900) * 900)
    inc = deg2rad(64.8)
    for prn in 1:n
        Om = deg2rad(120.0 * (prn % 3)) + 0.7
        u = deg2rad(45.0 * prn)
        # position/velocity in the orbital plane, rotated to ECEF (inertial
        # orientation at tb; the RK4 propagation includes Earth rotation terms)
        sO, cO = sincos(Om); si, ci = sincos(inc); su, cu = sincos(u)
        p = (r * (cu * cO - su * ci * sO), r * (cu * sO + su * ci * cO), r * su * si)
        vi = (v * (-su * cO - cu * ci * sO), v * (-su * sO + cu * ci * cO), v * cu * si)
        # ECEF velocity = inertial velocity - omega x r
        ve = (vi[1] + OMEGA_E_GLO * p[2], vi[2] - OMEGA_E_GLO * p[1], vi[3])
        push!(nav, GlonassEph(prn, (prn % 14) - 7, tb, p, ve, (0.0, 0.0, 0.0), 1e-5 * prn, 1e-12))
    end
    return nav
end

"""A plausible GPS/Galileo-like constellation of `n` Kepler records around GPS
time `t0`, for tests that need no real broadcast data. Six planes, evenly
spaced satellites; orbital elements of the GPS design (26 560 km, 55 deg)."""
function synthetic_constellation(t0::GpsTime; n::Int = 24, sys::Char = 'G')
    nav = NavStore()
    w, tow = gps_week_tow(t0)
    toe_sow = floor(tow / 7200) * 7200
    toe = gps_time(w, toe_sow)
    per_plane = cld(n, 6)
    prn = 1
    for plane in 0:5, j in 0:(per_plane - 1)
        prn > n && break
        Om0 = deg2rad(60.0 * plane) + 0.3
        M0 = deg2rad(360.0 / per_plane * j + 30.0 * plane)
        push!(nav, KeplerEph(sys, prn, 100 + prn, toe, toe,
                             sqrt(26_560_000.0), 0.01 + 0.002 * (prn % 3), deg2rad(55.0),
                             Om0, deg2rad(40.0 + 10prn), M0, 4.5e-9, 1e-10, -8.0e-9,
                             1e-6, 7e-6, 200.0, 50.0, 1e-7, -1e-7,
                             1e-4 * (prn % 5), 1e-11, 0.0, toe_sow, w))
        prn += 1
    end
    return nav
end
