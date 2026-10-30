# The variometric engine.
#
# For every satellite seen at two consecutive epochs with continuous lock, the
# observation is the change of the (ionosphere-free or single-frequency)
# carrier phase minus the change of everything modelled: geometric range to the
# a-priori station, satellite clock, troposphere, solid tide, wind-up. What is
# left is the receiver displacement projected on the line of sight plus the
# receiver clock change of that constellation. Weighted least squares with
# iterative rejection gives the displacement over the pair.
#
# Pair rule: both epochs of a pair are modelled with the SAME broadcast record
# per satellite, the record that was current when the previous epoch was
# solved. A record that arrives between the two epochs is used from the next
# pair on. The same will hold for SSR issue numbers.

mutable struct SatState
    phase_prev::Float64          # combined phase at the previous epoch [m]
    lock_prev::Dict{Int,Int}
    half_prev::Dict{Int,Bool}
    gf_prev::Float64             # geometry-free phase at the previous epoch [m]
    pair::Tuple{Int,Int}         # frequency slots used (second is 0 in :single)
    eph_pinned::Union{AnyEph,Nothing}   # record current at the previous epoch
    t_prev::GpsTime
    valid::Int                   # consecutive epochs with a valid pair
    hold::Int                    # epochs left in the hysteresis hold after a rejection
end

# Geometry-free history of one satellite (mode :l1gf): it survives short gaps of the satellite,
# so a dropout does not restart the smoothing of the ionospheric rate.
mutable struct GfState
    gf_last::Float64
    t_last::GpsTime
    rate::Float64                # smoothed rate of the geometry-free phase [m/s]
    span::Float64                # duration of the continuous arc behind the rate [s]
end

mutable struct Engine
    opt::Options
    nav::NavStore
    rec::NTuple{3,Float64}       # a-priori station ECEF
    lat::Float64; lon::Float64; hgt::Float64
    frame::Tuple{NTuple{3,Float64},NTuple{3,Float64},NTuple{3,Float64}}
    zhd::Float64; zwd::Float64
    sats::Dict{Tuple{Char,Int},SatState}
    t_prev::Union{GpsTime,Nothing}
    rms_hist::Vector{Float64}
    epochs::Int
    solved::Int
    gfs::Dict{Tuple{Char,Int},GfState}
end

function Engine(rec::NTuple{3,Float64}; nav::NavStore = NavStore(), opt::Options = Options())
    lat, lon, h = ecef2geodetic(rec)
    zhd, zwd = zenith_delays(lat, h)
    Engine(opt, nav, rec, lat, lon, h, neu_frame(lat, lon), zhd, zwd,
           Dict{Tuple{Char,Int},SatState}(), nothing, Float64[], 0, 0, Dict{Tuple{Char,Int},GfState}())
end

"""Forget the previous epoch and every per-satellite history: the next epoch
starts a new pair (the platform calls this after a hole in the input)."""
function reset!(eng::Engine)
    empty!(eng.sats)
    isdefined(eng, :gfs) && empty!(eng.gfs)
    eng.t_prev = nothing
    return eng
end

"""Set a new a-priori position (the engine differences, so metres of error are harmless)."""
function set_position!(eng::Engine, rec::NTuple{3,Float64})
    eng.rec = rec
    eng.lat, eng.lon, eng.hgt = ecef2geodetic(rec)
    eng.frame = neu_frame(eng.lat, eng.lon)
    eng.zhd, eng.zwd = zenith_delays(eng.lat, eng.hgt)
    eng
end

# First slot of the candidate list that carries a finite phase, or 0.
function _pick_slot(cands::Vector{Int}, s::SatObs)
    for sl in cands
        isfinite(get(s.phase, sl, NaN)) && return sl
    end
    return 0
end

# Combined observable of a satellite: ionosphere-free pair or one frequency.
# Returns (phase [m], geometry-free [m] or NaN, wavelength factor for wind-up
# [m/cycle], slot pair, ok)
function _combine(opt::Options, s::SatObs, k::Int)
    s1 = _pick_slot(get(opt.slot1_by_sys, s.sys, [1]), s)
    s1 == 0 && return NaN, NaN, NaN, (0, 0), false
    f1 = carrier_freq(s.sys, s1, k)
    p1 = s.phase[s1]
    if opt.mode == :single
        isnan(f1) && return NaN, NaN, NaN, (s1, 0), false
        return p1, NaN, C_LIGHT / f1, (s1, 0), true
    end
    s2 = _pick_slot(get(opt.slot2_by_sys, s.sys, [2]), s)
    if s2 == 0
        opt.mode == :l1gf && return p1, NaN, C_LIGHT / f1, (s1, 0), true      # no second frequency: plain L1
        return NaN, NaN, NaN, (s1, 0), false
    end
    f2 = carrier_freq(s.sys, s2, k)
    p2 = s.phase[s2]
    (isnan(f1) || isnan(f2) || f1 == f2) && return NaN, NaN, NaN, (s1, s2), false
    if opt.mode == :l1gf
        # L1 phase; the second frequency only supplies the geometry-free phase for the ionospheric correction
        return p1, p1 - p2, C_LIGHT / f1, (s1, s2), true
    end
    if opt.mode == :mean
        # plain average of the two phases [m]: no noise amplification; the
        # ionospheric change between epochs is left in (it is far below the
        # noise at 10 Hz); wind-up enters with the mean wavelength
        return 0.5 * (p1 + p2), p1 - p2, 0.5 * C_LIGHT * (1 / f1 + 1 / f2), (s1, s2), true
    end
    a = f1^2 / (f1^2 - f2^2); b = -f2^2 / (f1^2 - f2^2)
    return a * p1 + b * p2, p1 - p2, C_LIGHT / (f1 + f2), (s1, s2), true
end

# Lock continuity between the previous and the current epoch on every slot used.
function _lock_ok(opt::Options, st::SatState, s::SatObs, pair::Tuple{Int,Int})
    pair == st.pair || return false
    for sl in pair
        sl == 0 && continue
        lp = get(st.lock_prev, sl, -1); lc = get(s.lock, sl, -1)
        (lp < 0 || lc < 0 || lc < lp) && return false
        get(st.half_prev, sl, false) != get(s.half, sl, false) && return false
    end
    return true
end

# Observation weight from elevation [rad] and the weakest C/N0 of the pair [dB-Hz].
function _weight(opt::Options, el::Float64, cnr::Float64)
    w = _weight0(opt, el, cnr)
    if opt.taper_deg > 0
        s = clamp((el - deg2rad(opt.elev_mask_deg)) / deg2rad(opt.taper_deg), 0.0, 1.0)
        w *= s * s * (3 - 2s)                      # smoothstep: a rising satellite fades in
    end
    return w
end
function _weight0(opt::Options, el::Float64, cnr::Float64)
    m = opt.weight_model
    m == :uniform && return 1.0
    m == :sin && return sin(el)
    m == :sin2floor && return sin(min(el, deg2rad(opt.weight_floor_deg)))^2
    c = isfinite(cnr) ? 10.0^((cnr - 45.0) / 10) : 1.0
    m == :cnr && return c
    m == :sin2cnr && return sin(el)^2 * c
    return sin(el)^2
end

# Weakest C/N0 over the slots of the pair (NaN when not reported).
function _pair_cnr(s::SatObs, pair::Tuple{Int,Int})
    c = Inf
    for sl in pair
        sl == 0 && continue
        c = min(c, get(s.cnr, sl, NaN))
    end
    return isfinite(c) ? c : NaN
end

# All model terms of one satellite at one epoch with one record [m]; also
# elevation/azimuth and the satellite position (for wind-up and weights).
function _model(eng::Engine, eph::AnyEph, t::GpsTime, sun, tide, doy)
    rho, ps, clk = sagnac_range(eph, t, eng.rec)
    el, az = elev_azim(eng.frame, eng.rec, ps)
    m = rho - C_LIGHT * clk
    if eng.opt.tropo
        m += tropo_delay(el, eng.lat, eng.hgt, doy, eng.zhd, eng.zwd)
    end
    if eng.opt.tide
        # tide moves the receiver toward the satellite by tide . e_los
        e = unit3(sub3(ps, eng.rec))
        m -= dot3(tide, e)
    end
    w = eng.opt.windup ? windup_cycles(ps, eng.rec, sun, eng.frame) : 0.0
    return m, w, el, az, ps
end

"""Advance the engine by one epoch and return the solution of the pair
(previous epoch -> this epoch). The first epoch returns status :first."""
function step!(eng::Engine, ep::Epoch)
    opt = eng.opt
    eng.epochs += 1
    t = ep.t
    dt = eng.t_prev === nothing ? NaN : t - eng.t_prev
    doy = day_of_year(t)
    sun_now = sun_ecef(t)
    sun_prev = isnan(dt) ? sun_now : sun_ecef(eng.t_prev)
    tide_now = opt.tide ? solid_tide(eng.rec, t) : (0.0, 0.0, 0.0)
    tide_prev = (opt.tide && !isnan(dt)) ? solid_tide(eng.rec, eng.t_prev) : (0.0, 0.0, 0.0)

    rows = Vector{NTuple{3,Float64}}()      # -e_los
    sysv = Char[]
    obs = Float64[]
    wts = Float64[]
    results = SatResult[]
    idx = Int[]                               # index into results for rows
    seen = Set{Tuple{Char,Int}}()

    for s in ep.sats
        key = (s.sys, s.prn)
        s.sys in opt.systems || continue
        push!(seen, key)
        eph_now = get(eng.nav, s.sys, s.prn)
        k = eph_now isa GlonassEph ? eph_now.k : 0
        phase, gf, lam, pair, ok = _combine(opt, s, k)
        st = get(eng.sats, key, nothing)
        reason = :ok
        el_deg = NaN; az_deg = NaN
        if eph_now === nothing
            reason = :noeph
        elseif !ok
            reason = :nofreq
        elseif st === nothing || st.eph_pinned === nothing || isnan(dt)
            reason = :new
        elseif !_lock_ok(opt, st, s, pair)
            reason = :slip
        elseif !isnan(gf) && !isnan(st.gf_prev) && abs(gf - st.gf_prev) > opt.slip_gf_m
            reason = :slip
        elseif st.valid < opt.admit_epochs
            reason = :wait                          # admission control: not yet stable
        elseif st.hold > 0
            reason = :held                          # hysteresis hold after a rejection
        else
            eph = st.eph_pinned
            m_now, w_now, el, az, ps = _model(eng, eph, t, sun_now, tide_now, doy)
            m_prev, w_prev, _, _, _ = _model(eng, eph, st.t_prev, sun_prev, tide_prev, doy)
            el_deg = rad2deg(el); az_deg = rad2deg(az)
            if el < deg2rad(opt.elev_mask_deg)
                reason = :elev
            else
                dmodel = (m_now - m_prev) + lam * windup_diff(w_now, w_prev)
                l = (phase - st.phase_prev) - dmodel
                gs = get(eng.gfs, key, nothing)
                if opt.mode == :l1gf && pair[2] != 0 && !isnan(gf) && gs !== nothing && gs.span >= 0.5 * opt.gf_tau_s
                    # remove the ionospheric change of the L1 phase: p1 = rho - I1 and gf = (g - 1) I1, g = (f1/f2)^2
                    f1 = carrier_freq(s.sys, pair[1], k); f2 = carrier_freq(s.sys, pair[2], k)
                    r = gs.rate; db = 1e-3 * opt.gf_deadband_mm_s
                    r = sign(r) * max(abs(r) - db, 0.0)
                    l += r * dt / ((f1 / f2)^2 - 1)
                end
                e = unit3(sub3(ps, eng.rec))
                push!(rows, (-e[1], -e[2], -e[3]))
                push!(sysv, s.sys)
                push!(obs, l)
                push!(wts, _weight(opt, el, _pair_cnr(s, pair)))
                push!(idx, length(results) + 1)
            end
        end
        push!(results, SatResult(s.sys, s.prn, el_deg, az_deg, NaN, NaN, false, reason))
        # Update the per-satellite state for the next pair: pin the record that
        # is current NOW (it will be used for both epochs of the next pair).
        if opt.mode == :l1gf && !isnan(gf) && eph_now !== nothing
            _gf_update!(eng, key, gf, t, opt)
        end
        if ok && eph_now !== nothing
            if st === nothing
                eng.sats[key] = SatState(phase, copy(s.lock), copy(s.half), gf, pair, eph_now, t, 0, 0)
            else
                if reason in (:ok, :elev, :wait, :held)
                    st.valid += 1
                else
                    st.valid = 0
                end
                st.hold > 0 && (st.hold -= 1)
                st.phase_prev = phase; st.gf_prev = gf; st.t_prev = t; st.pair = pair
                st.lock_prev = copy(s.lock); st.half_prev = copy(s.half)
                st.eph_pinned = eph_now
            end
        elseif st !== nothing
            delete!(eng.sats, key)
        end
    end
    # Satellites that vanished lose their state (a reappearance starts a new pair).
    for key in collect(keys(eng.sats))
        key in seen || delete!(eng.sats, key)
    end
    eng.t_prev = t
    if opt.mode == :l1gf
        for key in collect(keys(eng.gfs))
            (t - eng.gfs[key].t_last) > opt.gf_gap_s && delete!(eng.gfs, key)
        end
    end
    isnan(dt) && return nofit(t, dt, :first, results)

    sol = _solve(eng, t, dt, rows, sysv, obs, wts, results, idx)
    sol.status == :ok && (eng.solved += 1)
    return sol
end

# Update of the geometry-free history of one satellite with this epoch's value.
function _gf_update!(eng::Engine, key, gf::Float64, t::GpsTime, opt::Options)
    gs = get(eng.gfs, key, nothing)
    if gs !== nothing
        dts = t - gs.t_last
        if 0 < dts <= opt.gf_gap_s && abs(gf - gs.gf_last - gs.rate * dts) < opt.slip_gf_m
            α = min(1.0, max(dts / opt.gf_tau_s, 1 / (1 + gs.span / max(dts, 1e-3))))
            gs.rate += α * ((gf - gs.gf_last) / dts - gs.rate)
            gs.span += dts; gs.gf_last = gf; gs.t_last = t
            return
        end
    end
    eng.gfs[key] = GfState(gf, t, 0.0, 0.0)          # new arc (first sight, a slip, or a gap that was too long)
end

# Weighted least squares with iterative rejection on standardised residuals.
function _solve(eng::Engine, t, dt, rows, sysv, obs, wts, results, idx)
    opt = eng.opt
    n = length(obs)
    systems = unique(sysv)
    nsys = length(systems)
    common = opt.clock_model == :common
    npar = 3 + (common ? 1 : nsys)
    n < max(opt.min_sats, npar + 1) && return nofit(t, dt, :few_sats, results)
    syscol = Dict(s => (common ? 4 : 3 + i) for (i, s) in enumerate(systems))
    A = zeros(n, npar)
    for i in 1:n
        A[i, 1], A[i, 2], A[i, 3] = rows[i]
        A[i, syscol[sysv[i]]] = 1.0
    end
    active = trues(n)
    v = zeros(n)
    x = zeros(npar)
    rms = NaN
    max_rej = floor(Int, opt.max_reject_frac * n)
    nrej = 0
    if opt.robust == :none
        max_rej = 0
    elseif opt.robust == :huber
        max_rej = 0
    end
    hw = ones(n)                               # Huber factors (1 = untouched)
    for it in 1:(opt.robust == :huber ? 4 : max_rej + 1)
        ia = findall(active)
        length(ia) < max(opt.min_sats, npar + 1) && return nofit(t, dt, :few_sats, results)
        Aa = A[ia, :]; W = wts[ia] .* hw[ia]; la = obs[ia]
        N = Aa' * (W .* Aa)
        rank(N) < npar && return nofit(t, dt, :singular, results)
        x = N \ (Aa' * (W .* la))
        v .= obs .- A * x
        rms = sqrt(sum(W .* v[ia] .^ 2) / max(length(ia) - npar, 1))
        # robust scale of the SIGNED weighted residuals (the MAD of absolute
        # values would be ~0.6 sigma and reject inliers), floored by sigma0_m
        z = v[ia] .* sqrt.(W)
        med = _median(z)
        mad = max(1.4826 * _median(abs.(z .- med)), opt.sigma0_m)
        dev = abs.(z .- med)
        if opt.robust == :huber
            # Huber weights from the unweighted-by-hw standardised residuals: nothing is removed
            zz = v .* sqrt.(wts)
            for i in 1:n
                a = abs(zz[i] - med) / mad
                hw[i] = a <= opt.huber_c ? 1.0 : opt.huber_c / a
            end
            continue
        end
        worst = argmax(dev)
        if dev[worst] > opt.reject_k * mad && nrej < max_rej
            active[ia[worst]] = false
            nrej += 1
            if opt.robust == :hyst
                r = results[idx[ia[worst]]]
                st = get(eng.sats, (r.sys, r.prn), nothing)
                st === nothing || (st.hold = opt.hyst_epochs)
            end
        else
            break
        end
    end
    ia = findall(active)
    # PDOP from the unweighted geometry of the used satellites
    G = A[ia, :]
    Q = inv(G' * G)
    pdop = sqrt(Q[1, 1] + Q[2, 2] + Q[3, 3])
    # formal covariance of the position from the weighted fit, in NEU
    Qw = inv(G' * (wts[ia] .* G))
    nf, ef, uf = eng.frame
    R = [nf[1] nf[2] nf[3]; ef[1] ef[2] ef[3]; uf[1] uf[2] uf[3]]
    Qneu = R * Qw[1:3, 1:3] * R'
    sig_neu = (rms * sqrt(max(Qneu[1, 1], 0.0)), rms * sqrt(max(Qneu[2, 2], 0.0)), rms * sqrt(max(Qneu[3, 3], 0.0)))
    for i in 1:n
        r = results[idx[i]]
        results[idx[i]] = SatResult(r.sys, r.prn, r.elev_deg, r.azim_deg, v[i], wts[i],
                                    active[i], active[i] ? :ok : :rejected)
    end
    dxyz = (x[1], x[2], x[3])
    dneu = (dot3(dxyz, nf), dot3(dxyz, ef), dot3(dxyz, uf))
    clk = Dict(s => x[syscol[s]] for s in systems)
    push!(eng.rms_hist, rms)
    length(eng.rms_hist) > opt.sigma_window && popfirst!(eng.rms_hist)
    sigma = _median(eng.rms_hist) * 1.0
    return Solution(t, dt, dxyz, dneu, clk, n, length(ia), nsys, pdop, rms, sig_neu, sigma, results, :ok)
end

function _median(v::AbstractVector{<:Real})
    isempty(v) && return NaN
    s = sort(v)
    m = length(s)
    return isodd(m) ? Float64(s[(m + 1) ÷ 2]) : 0.5 * (s[m ÷ 2] + s[m ÷ 2 + 1])
end
