using Test
using GaiaVario
const GV = GaiaVario

# Patras-like station (ECEF, metres)
const REC = (4.65456e6, 1.86062e6, 3.93057e6)   # lat 38.29 N, lon 21.79 E, h ~100 m
# An instant mid-week: GPS week 2400, Wednesday 10:00
const T0 = gps_time(2400, 3 * 86400 + 36000.0)

# Feed a synthetic run through an engine; returns per-epoch NEU error (sol - truth increment)
function run_engine(nav, epochs, truth, nav_events; opt = Options())
    eng = Engine(REC; nav = nav, opt = opt)
    ev = copy(nav_events)
    sols = Solution[]
    errs = NTuple{3,Float64}[]
    for (i, ep) in enumerate(epochs)
        while !isempty(ev) && ev[1][1] <= ep.t
            push!(eng.nav, popfirst!(ev)[2])
        end
        sol = step!(eng, ep)
        push!(sols, sol)
        if sol.status == :ok
            d = (truth[i][1] - truth[i - 1][1], truth[i][2] - truth[i - 1][2], truth[i][3] - truth[i - 1][3])
            push!(errs, (sol.dneu[1] - d[1], sol.dneu[2] - d[2], sol.dneu[3] - d[3]))
        end
    end
    return eng, sols, errs
end
maxabs(errs) = maximum(e -> max(abs(e[1]), abs(e[2]), abs(e[3])), errs)
rms(v) = sqrt(sum(abs2, v) / length(v))

@testset "time" begin
    w, tow = gps_week_tow(T0)
    @test w == 2400 && tow ≈ 3 * 86400 + 36000.0
    @test GV.wrap_week(400_000.0) ≈ 400_000.0 - 604_800
    @test GV.wrap_day(-50_000.0) ≈ 36_400.0
    @test GV.julian_date(GpsTime(0, 0.0)) == 2_444_244.5
    d0 = GV.day_of_year(T0)
    @test GV.day_of_year(T0 + 86400) ≈ d0 + 1 atol = 1e-9
    # 1980-01-06 is day 6 of 1980
    @test floor(GV.day_of_year(GpsTime(0, 0.0))) == 6
end

@testset "orbits" begin
    nav = synthetic_constellation(T0; n = 24)
    for prn in 1:24
        e = get(nav, 'G', prn)
        p, clk = sat_state(e, T0 + 1234.5)
        r = sqrt(sum(abs2, p))
        @test 26_100e3 < r < 27_000e3        # A(1 ± e), e up to 0.014
        # ECEF velocity: inertial 3.87 km/s minus the Earth-rotation share
        v = GV.sat_velocity(e, T0 + 1234.5)
        @test 2.4e3 < sqrt(sum(abs2, v)) < 4.0e3
        @test abs(clk) < 1e-3
    end
    rho, ps, clk = sagnac_range(get(nav, 'G', 1), T0, REC)
    @test 19_000e3 < rho < 27_000e3
    # GLONASS: radius conserved by the RK4 integration over 30 min
    gnav = GV.synthetic_glonass(T0)
    g = get(gnav, 'R', 3)
    r0 = sqrt(sum(abs2, GV.glonass_position(g, g.tb)))
    r1 = sqrt(sum(abs2, GV.glonass_position(g, g.tb + 1800)))
    @test abs(r1 - r0) < 5e3
    @test abs(r0 - 25_508e3) < 1.0
    # BeiDou GEO stays near its longitude over an hour. The broadcast i0 of a
    # GEO carries the 5 deg tilt that the Rx(-5 deg) rotation removes; the
    # tilt is about the x axis of the frame fixed at toe; with the ICD's
    # R_X(-5 deg) an equatorial orbit appears in that frame with i = 5 deg and
    # its ascending node at Omega_k = pi (checked by hand: the point of
    # maximum latitude maps to z = 0 only for that node).
    w, tow = gps_week_tow(T0)
    om0 = mod(GV.OMEGA_E_BDS * tow + π, 2π)
    geo = KeplerEph('C', 1, 1, T0, T0, sqrt(42_164_000.0), 0.0005, deg2rad(5.0), om0, 0.5, 2.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, tow, w)
    p1, _ = GV.kepler_position(geo, T0)
    p2, _ = GV.kepler_position(geo, T0 + 3600)
    @test sqrt(sum(abs2, p1 .- p2)) < 150e3
    @test abs(sqrt(sum(abs2, p1)) - 42_164e3) < 50e3
end

@testset "models" begin
    lat, lon, h = ecef2geodetic(REC)
    @test 0.5 < lat < 0.8 && 0.3 < lon < 0.5 && -50 < h < 500
    zhd, zwd = GV.zenith_delays(lat, h)
    @test 1.8 < zhd < 2.4 && 0.0 < zwd < 0.4
    mh, mw = GV.niell_mf(deg2rad(10.0), lat, h, 120.0)
    @test 5.4 < mh < 5.9 && 5.4 < mw < 6.0
    mh90, _ = GV.niell_mf(deg2rad(90.0), lat, h, 120.0)
    @test mh90 ≈ 1.0 atol = 0.01
    s = GV.sun_ecef(T0); m = GV.moon_ecef(T0)
    @test 1.45e11 < sqrt(sum(abs2, s)) < 1.53e11
    @test 3.5e8 < sqrt(sum(abs2, m)) < 4.1e8
    tide = GV.solid_tide(REC, T0)
    @test sqrt(sum(abs2, tide)) < 0.5
    # the tide changes by micrometres per 0.1 s
    t2 = GV.solid_tide(REC, T0 + 0.1)
    @test sqrt(sum(abs2, tide .- t2)) < 2e-5
end

@testset "engine: static, no noise" begin
    nav = synthetic_constellation(T0)
    so = SynthOptions(duration_s = 30.0, noise_m = 0.0, receiver_clock = true)
    epochs, truth, ev = synthesize(nav, REC, T0, motion_static; so = so)
    eng, sols, errs = run_engine(nav, epochs, truth, ev)
    @test sols[1].status == :first
    @test all(s -> s.status == :ok, sols[2:end])
    @test maxabs(errs) < 2e-5
    @test all(s -> s.nsat_used == s.nsat_in, sols[2:end])
end

@testset "engine: step and pulse, no noise" begin
    nav = synthetic_constellation(T0)
    so = SynthOptions(duration_s = 60.0, noise_m = 0.0)
    epochs, truth, ev = synthesize(nav, REC, T0, motion_step(0.05, 30.0); so = so)
    _, sols, errs = run_engine(nav, epochs, truth, ev)
    @test maxabs(errs) < 5e-5
    i = findfirst(s -> abs(s.dneu[1]) > 0.04, sols)
    @test i !== nothing
    @test sols[i].dneu[1] ≈ 0.05 atol = 1e-4
    epochs, truth, ev = synthesize(nav, REC, T0, motion_pulse(0.10, 30.0, 2.0); so = so)
    _, sols, errs = run_engine(nav, epochs, truth, ev)
    @test maxabs(errs) < 5e-5
end

# Real broadcast records recorded by tools/record_brdc.jl (newest file), or nothing
using Dates
function newest_brdc()
    dir = joinpath(@__DIR__, "data")
    isdir(dir) || return nothing
    files = sort(filter(f -> startswith(f, "brdc_") && endswith(f, ".tsv"), readdir(dir)))
    isempty(files) && return nothing
    path = joinpath(dir, last(files))
    stamp = match(r"brdc_(\d{8})_(\d{4})", basename(path))
    dtm = DateTime(stamp[1] * stamp[2], dateformat"yyyymmddHHMM")
    t_rec = GpsTime(0, Dates.value(dtm - DateTime(1980, 1, 6)) / 1000 + 18.0)   # UTC -> GPS
    return GaiaVario.load_brdc_tsv(path, t_rec), t_rec, basename(path)
end

@testset "engine: noise 2 mm" begin
    # With the synthetic constellation only ~6 satellites are above 10 deg at
    # Patras, so the test of the noise propagation uses the real constellation
    # when a recording exists: ~30 satellites of four systems, PDOP < 2.
    rb = newest_brdc()
    nav, t0 = rb === nothing ? (synthetic_constellation(T0), T0) : (rb[1], rb[2])
    so = SynthOptions(duration_s = 60.0, noise_m = 0.002, seed = 7)
    epochs, truth, ev = synthesize(nav, REC, t0, motion_sine(0.02, 5.0); so = so)
    _, sols, errs = run_engine(nav, epochs, truth, ev)
    pd = maximum(s -> s.pdop, sols[2:end])
    # the formal sigma of the fit must predict the empirical error per component
    ok = filter(s -> s.status == :ok, sols)
    for c in 1:3
        emp = rms([e[c] for e in errs])
        formal = sqrt(sum(s -> s.sig_neu[c]^2, ok) / length(ok))
        @test 0.7 < emp / formal < 1.4
    end
    if rb !== nothing
        @test pd < 2.5
        # 2 mm per frequency -> ~12 mm per ionosphere-free pair observation
        @test rms([e[1] for e in errs]) < 0.012
        @test rms([e[3] for e in errs]) < 0.030
    end
end

@testset "engine: cycle slip is isolated" begin
    nav = synthetic_constellation(T0)
    # pick a satellite the engine actually uses at Patras (the synthetic
    # constellation has only ~6 above the mask)
    so0 = SynthOptions(duration_s = 1.0, noise_m = 0.0)
    ep0, tr0, ev0 = synthesize(nav, REC, T0, motion_static; so = so0)
    _, sols0, _ = run_engine(nav, ep0, tr0, ev0)
    prn = first(r.prn for r in sols0[end].sats if r.used)
    so = SynthOptions(duration_s = 20.0, noise_m = 0.0, slips = [('G', prn, 10.0)])
    epochs, truth, ev = synthesize(nav, REC, T0, motion_static; so = so)
    _, sols, errs = run_engine(nav, epochs, truth, ev)
    @test maxabs(errs) < 2e-5
    i = round(Int, 10.0 * so.rate_hz) + 1
    r = only(filter(x -> x.sys == 'G' && x.prn == prn, sols[i].sats))
    @test r.reason == :slip && !r.used
    r2 = only(filter(x -> x.sys == 'G' && x.prn == prn, sols[i + 1].sats))
    @test r2.used
end

@testset "engine: ephemeris changeover produces no step" begin
    nav = synthetic_constellation(T0)
    # a 2 m radial change of the record for G07 at 10 s, and another for G12 at 12 s
    so = SynthOptions(duration_s = 20.0, noise_m = 0.0,
                      changeovers = [('G', 7, 10.0, 2.0), ('G', 12, 12.0, -3.0)])
    epochs, truth, ev = synthesize(nav, REC, T0, motion_static; so = so)
    @test length(ev) == 2
    eng, sols, errs = run_engine(nav, epochs, truth, ev)
    @test eng.nav.updates == length(nav) + 2      # initial records + the two changeovers
    @test maxabs(errs) < 1e-4
    @test all(s -> s.nsat_used == s.nsat_in, sols[2:end])
end

@testset "engine: day and week boundaries" begin
    for t0 in (gps_time(2400, 3 * 86400 - 5.0), gps_time(2400, 604_800 - 5.0))
        nav = synthetic_constellation(t0)
        so = SynthOptions(duration_s = 10.0, noise_m = 0.0)
        epochs, truth, ev = synthesize(nav, REC, t0, motion_static; so = so)
        _, sols, errs = run_engine(nav, epochs, truth, ev)
        @test all(s -> s.status == :ok, sols[2:end])
        @test maxabs(errs) < 2e-5
    end
end

@testset "engine: four systems, single frequency" begin
    nav = synthetic_constellation(T0; n = 10, sys = 'G')
    for e in values(synthetic_constellation(T0 + 1.0; n = 8, sys = 'E').eph)
        push!(nav, e)
    end
    for e in values(synthetic_constellation(T0 + 2.0; n = 8, sys = 'C').eph)
        # keep away from the GEO PRNs
        push!(nav, KeplerEph('C', e.prn + 10, e.iod, e.toe, e.toc, e.sqrtA, e.e, e.i0, e.Omega0 + 0.4,
                             e.omega, e.M0 + 0.3, e.dn, e.idot, e.Omegadot, e.cuc, e.cus, e.crc, e.crs,
                             e.cic, e.cis, e.af0, e.af1, e.af2, e.toe_sow, e.week))
    end
    for e in values(GV.synthetic_glonass(T0).eph)
        push!(nav, e)
    end
    so = SynthOptions(duration_s = 20.0, noise_m = 0.0)
    epochs, truth, ev = synthesize(nav, REC, T0, motion_sine(0.01, 4.0); so = so)
    _, sols, errs = run_engine(nav, epochs, truth, ev)
    @test all(s -> s.status == :ok, sols[2:end])
    @test maximum(s -> s.nsys, sols[2:end]) == 4
    @test maxabs(errs) < 5e-5
    _, sols1, errs1 = run_engine(nav, epochs, truth, ev; opt = Options(mode = :single))
    @test all(s -> s.status == :ok, sols1[2:end])
    @test maxabs(errs1) < 5e-5
end

@testset "adapter: decoder rows round-trip" begin
    nav = synthetic_constellation(T0)
    e = get(nav, 'G', 3)
    row = fill(NaN, 39)
    row[1] = 3; row[2] = e.week; row[8] = e.af0; row[9] = e.af1; row[10] = e.af2; row[11] = e.iod
    row[12] = e.crs; row[13] = e.dn; row[14] = e.M0; row[15] = e.cuc; row[16] = e.e; row[17] = e.cus
    row[18] = e.sqrtA; row[19] = e.toe_sow; row[20] = e.cic; row[21] = e.Omega0; row[22] = e.cis
    row[23] = e.i0; row[24] = e.crc; row[25] = e.omega; row[26] = e.Omegadot; row[27] = e.idot
    row[39] = mod(e.toe_sow, 86400.0)
    e2 = eph_from_row('G', 3, row, T0)
    @test e2.toe == e.toe && e2.toc == e.toc
    p1, _ = sat_state(e, T0 + 100); p2, _ = sat_state(e2, T0 + 100)
    @test p1 == p2
    # Galileo row without a week column
    row[2] = NaN
    e3 = eph_from_row('E', 3, row, T0 + 86400 * 2)
    @test e3.toe == e.toe
end

# ---- Real broadcast records (recorded by tools/record_brdc.jl) -------------------
using Dates
const BRDC_FILES = filter(f -> startswith(f, "brdc_") && endswith(f, ".tsv"),
                          isdir(joinpath(@__DIR__, "data")) ? readdir(joinpath(@__DIR__, "data")) : String[])

if !isempty(BRDC_FILES)
    @testset "real broadcast records: $(last(sort(BRDC_FILES)))" begin
        path = joinpath(@__DIR__, "data", last(sort(BRDC_FILES)))
        # recording time from the file name, UTC -> GPS (+18 s leap seconds since 2017)
        stamp = match(r"brdc_(\d{8})_(\d{4})", basename(path))
        dtm = DateTime(stamp[1] * stamp[2], dateformat"yyyymmddHHMM")
        t_rec = GpsTime(0, Dates.value(dtm - DateTime(1980, 1, 6)) / 1000 + 18.0)
        nav = GaiaVario.load_brdc_tsv(path, t_rec)
        @test length(nav) > 80
        ranges = Dict('G' => (26.0e6, 27.2e6), 'E' => (29.0e6, 30.2e6), 'R' => (25.0e6, 26.0e6))
        n_by = Dict{Char,Int}()
        for ((sys, prn), e) in nav.eph
            n_by[sys] = get(n_by, sys, 0) + 1
            p, clk = sat_state(e, t_rec)
            r = sqrt(sum(abs2, p))
            @test abs(clk) < 1e-2            # Galileo IOV clocks reach several ms
            if sys == 'C'
                if GV.is_bds_geo(prn)
                    @test 41.5e6 < r < 42.8e6
                    lon = rad2deg(atan(p[2], p[1]))
                    @test 55.0 < lon < 165.0           # BDS GEO slots 58.75E .. 160E
                    p2, _ = sat_state(e, t_rec + 3600)
                    @test sqrt(sum(abs2, p .- p2)) < 400e3
                else
                    @test (27.4e6 < r < 28.5e6) || (41.5e6 < r < 42.8e6)   # MEO or IGSO
                end
            elseif sys == 'E' && e.e > 0.05
                @test 23.0e6 < r < 33.0e6          # E14/E18 on eccentric orbits (e = 0.17)
            else
                lo, hi = ranges[sys]
                @test lo < r < hi
            end
        end
        @test all(haskey(n_by, s) for s in "GREC")
        # the real constellation seen from Patras drives a four-system synthetic run
        so = SynthOptions(duration_s = 15.0, noise_m = 0.0)
        epochs, truth, ev = synthesize(nav, REC, t_rec, motion_sine(0.01, 3.0); so = so)
        eng = Engine(REC; nav = nav)
        errs = NTuple{3,Float64}[]
        nsys = 0; nused = 0
        for (i, ep) in enumerate(epochs)
            sol = step!(eng, ep)
            sol.status == :ok || continue
            nsys = max(nsys, sol.nsys); nused = max(nused, sol.nsat_used)
            d = truth[i] .- truth[i - 1]
            push!(errs, sol.dneu .- d)
        end
        @test nsys == 4
        @test nused >= 25
        @test maxabs(errs) < 5e-5
    end
end

# ---- The observables respond to the ionosphere as the combinations predict ---------------------
if !isempty(BRDC_FILES)
    @testset "observables: ionosphere ramp (IF cancels it, L1 sees it, the mean sees 1.32 x L1)" begin
        path = joinpath(@__DIR__, "data", last(sort(BRDC_FILES)))
        stamp = match(r"brdc_(\d{8})_(\d{4})", basename(path))
        dtm = DateTime(stamp[1] * stamp[2], dateformat"yyyymmddHHMM")
        t_rec = GpsTime(0, Dates.value(dtm - DateTime(1980, 1, 6)) / 1000 + 18.0)
        nav = GaiaVario.load_brdc_tsv(path, t_rec)
        rec = (4157049.9942, 1614214.9074, 4544925.5357)
        lat, lon, _ = ecef2geodetic(rec)
        epochs, _, _ = synthesize(nav, rec, t_rec, motion_static; so = SynthOptions(duration_s = 30.0, noise_m = 0.0))
        function with_iono(which)
            map(epochs) do ep
                sats = map(ep.sats) do s
                    ns = SatObs(s.sys, s.prn, copy(s.phase), copy(s.code), copy(s.lock), copy(s.half), copy(s.cnr))
                    e = get(nav, s.sys, s.prn)
                    ps, _ = sat_state(e, ep.t)
                    el, _ = GaiaVario.elev_azim(GaiaVario.neu_frame(lat, lon), rec, ps)
                    F = 1 / sqrt(1 - (6371e3 * cos(el) / (6371e3 + 350e3))^2)
                    tec = 0.05 * (ep.t - epochs[1].t) * F * 1e16
                    for (sl, ph) in s.phase
                        f = GaiaVario.carrier_freq(s.sys, sl, e isa GlonassEph ? e.k : 0)
                        (which == :both || (which == :L2 && sl != 1)) && (ns.phase[sl] = ph - 40.3 * tec / f^2)
                    end
                    ns
                end
                Epoch(ep.t, sats)
            end
        end
        function err_u(eps, opt)
            eng = Engine(rec; nav = NavStore(copy(nav.eph), 0), opt = opt)
            u = 0.0
            for ep in eps
                s = step!(eng, ep)
                s.status == :ok && (u += s.dneu[3])
            end
            return 1000u
        end

        @testset "l1gf: L1 corrected with the smoothed geometry-free rate" begin
            epochs2, _, _ = synthesize(nav, rec, t_rec, motion_static; so = SynthOptions(duration_s = 90.0, noise_m = 0.0))
            function iono2(rate)
                map(epochs2) do ep
                    sats = map(ep.sats) do s
                        ns = SatObs(s.sys, s.prn, copy(s.phase), copy(s.code), copy(s.lock), copy(s.half), copy(s.cnr))
                        e = get(nav, s.sys, s.prn); ps, _ = sat_state(e, ep.t)
                        el, _ = GaiaVario.elev_azim(GaiaVario.neu_frame(lat, lon), rec, ps)
                        F = 1 / sqrt(1 - (6371e3 * cos(el) / (6371e3 + 350e3))^2)
                        tec = rate * (ep.t - epochs2[1].t) * F * 1e16
                        for (sl, ph) in s.phase
                            f = GaiaVario.carrier_freq(s.sys, sl, e isa GlonassEph ? e.k : 0)
                            ns.phase[sl] = ph - 40.3 * tec / f^2
                        end
                        ns
                    end
                    Epoch(ep.t, sats)
                end
            end
            function u_end(eps, opt)
                eng = Engine(rec; nav = NavStore(copy(nav.eph), 0), opt = opt)
                u = 0.0
                for ep in eps
                    s = step!(eng, ep); s.status == :ok && (u += s.dneu[3])
                end
                1000u
            end
            ramp = iono2(0.05)
            @test abs(u_end(ramp, Options(mode = :l1gf, gf_tau_s = 30.0))) < 0.35 * abs(u_end(ramp, Options(mode = :single)))
            # no ionosphere: no cost beyond the wind-up difference between the two wavelengths (0.6 mm in 90 s here)
            @test abs(u_end(epochs2, Options(mode = :l1gf)) - u_end(epochs2, Options(mode = :single))) < 2.0
        end
        both = with_iono(:both); l2 = with_iono(:L2)
        e_if = err_u(both, Options(mode = :dual)); e_l1 = err_u(both, Options(mode = :single)); e_mn = err_u(both, Options(mode = :mean))
        @test abs(e_if) < 0.5                                   # ionosphere-free: nothing
        @test abs(e_l1) > 100                                   # L1: a large vertical drift
        @test 1.2 < e_mn / e_l1 < 1.5                           # mean: (1 + 1.647)/2 = 1.32 times L1
        @test abs(err_u(l2, Options(mode = :single))) < 0.5     # L2-only ionosphere: L1 does not react
        @test 0.7 < err_u(l2, Options(mode = :mean)) / e_l1 < 1.0   # but the mean does, with weight 1/2 of 1.647
    end
end

@testset "default observable is l1gf; every observable recovers the motion" begin
    @test Options().mode == :l1gf
    nav = synthetic_constellation(T0)
    so = SynthOptions(duration_s = 40.0, noise_m = 0.0)
    epochs, truth, ev = synthesize(nav, REC, T0, motion_pulse(0.05, 20.0, 2.0); so = so)
    for opt in (Options(), Options(mode = :dual), Options(mode = :single), Options(mode = :mean))
        _, sols, errs = run_engine(nav, epochs, truth, ev; opt = opt)
        @test maxabs(errs) < 1e-4                  # every observable recovers the prescribed motion
    end
end
