# Data types shared by the engine, the adapters and the synthetic generator.

"""One satellite's observations at one epoch. Phases and codes are in metres,
indexed by frequency slot (1, 2, 5, 6, 7, 8; GLONASS 1, 2). `lock` is the RTCM
lock-time indicator per slot (a decrease means a loss of lock)."""
struct SatObs
    sys::Char
    prn::Int
    phase::Dict{Int,Float64}
    code::Dict{Int,Float64}
    lock::Dict{Int,Int}
    half::Dict{Int,Bool}      # half-cycle ambiguity flag per slot
    cnr::Dict{Int,Float64}
end
SatObs(sys, prn) = SatObs(sys, prn, Dict{Int,Float64}(), Dict{Int,Float64}(),
                          Dict{Int,Int}(), Dict{Int,Bool}(), Dict{Int,Float64}())

"""One receiver epoch: time in GPS seconds (all constellations already on the
GPS scale) and the satellites observed."""
struct Epoch
    t::GpsTime
    sats::Vector{SatObs}
end

# ---- Broadcast ephemerides -------------------------------------------------

"""Keplerian broadcast record (GPS LNAV, Galileo F/NAV-I/NAV, BeiDou D1/D2).
`toe` and `toc` are GPS seconds (continuous); for BeiDou they are already
shifted from BDT to GPST. `iod` is IODE / IODnav / AODE."""
struct KeplerEph
    sys::Char
    prn::Int
    iod::Int
    toe::GpsTime
    toc::GpsTime
    sqrtA::Float64
    e::Float64
    i0::Float64
    Omega0::Float64
    omega::Float64
    M0::Float64
    dn::Float64
    idot::Float64
    Omegadot::Float64
    cuc::Float64; cus::Float64
    crc::Float64; crs::Float64
    cic::Float64; cis::Float64
    af0::Float64; af1::Float64; af2::Float64
    toe_sow::Float64           # seconds of week as broadcast (for the -omega_e*toe term)
    week::Int                  # week of toe in the constellation's own count
end

"""GLONASS record: PZ-90 state vector at `tb` (GPS seconds), km units
converted to metres; clock as -tau_n (bias) and gamma_n (relative drift)."""
struct GlonassEph
    prn::Int
    k::Int                     # FDMA channel number
    tb::GpsTime
    pos::NTuple{3,Float64}
    vel::NTuple{3,Float64}
    acc::NTuple{3,Float64}
    tau_neg::Float64           # -tau_n [s]
    gamma::Float64             # gamma_n
end

const AnyEph = Union{KeplerEph,GlonassEph}

"""Latest broadcast record per (sys, prn). Records are kept by reference, so a
pointer taken at one epoch keeps pointing at that record after a newer one
arrives (the engine relies on this for the variometric pair)."""
mutable struct NavStore
    eph::Dict{Tuple{Char,Int},AnyEph}
    updates::Int
end
NavStore() = NavStore(Dict{Tuple{Char,Int},AnyEph}(), 0)

function Base.push!(nav::NavStore, e::KeplerEph)
    nav.eph[(e.sys, e.prn)] = e; nav.updates += 1; nav
end
function Base.push!(nav::NavStore, e::GlonassEph)
    nav.eph[('R', e.prn)] = e; nav.updates += 1; nav
end
Base.get(nav::NavStore, sys::Char, prn::Int) = get(nav.eph, (sys, prn), nothing)
Base.length(nav::NavStore) = length(nav.eph)

# ---- SSR (phase 2; the types exist so the pair bookkeeping covers them) ----

struct SsrCorr
    sys::Char
    prn::Int
    t::GpsTime                 # correction epoch
    iod::Int
    orb::NTuple{3,Float64}     # radial, along, cross [m]
    orb_rate::NTuple{3,Float64}
    clk::NTuple{3,Float64}     # c0 [m], c1 [m/s], c2 [m/s^2]
end

# ---- Engine configuration and output -----------------------------------------

Base.@kwdef mutable struct Options
    mode::Symbol = :l1gf              # DEFAULT :l1gf = L1 corrected with a smoothed geometry-free (ionospheric) rate per
                                      #         satellite (equal to :single in a quiet ionosphere, protects against ramps);
                                      #         user decision 2026-10-30. :single = L1 only, :dual = ionosphere-free pair,
                                      #         :mean = average of two phases.
    # frequency slots tried in order per system (first available is used;
    # the pair is kept per satellite and a change counts as a slip)
    # BeiDou: B1I + B3I first (best agreement with PPPHVA on live data, 2026-10-29)
    slot1_by_sys::Dict{Char,Vector{Int}} = Dict('G' => [1], 'R' => [1], 'E' => [1], 'C' => [2, 1])
    slot2_by_sys::Dict{Char,Vector{Int}} = Dict('G' => [2, 5], 'R' => [2], 'E' => [5, 7, 8, 6], 'C' => [6, 7, 5, 2])
    elev_mask_deg::Float64 = 10.0
    min_sats::Int = 5
    sigma0_m::Float64 = 0.002         # floor of the robust residual scale [m] (zenith-equivalent)
    # :sin2      w = sin^2(el)                     (noise ~ 1/sin el)
    # :sin       w = sin(el)
    # :sin2floor w = sin^2(min(el, weight_floor_deg))  (no extra weight above the floor)
    # :cnr       w = 10^((cnr - 45)/10), cnr = weakest signal of the pair
    # :sin2cnr   product of :sin2 and :cnr
    # :uniform   w = 1
    weight_model::Symbol = :sin2
    weight_floor_deg::Float64 = 30.0
    reject_k::Float64 = 4.0           # rejection threshold in robust sigmas (MAD-based scale)
    # Options that keep the satellite set and the weights from changing between
    # consecutive epochs (the increments of the pair solutions then telescope when summed):
    #   robust     :reject  hard iterative rejection at reject_k (default)
    #              :none    no rejection
    #              :huber   IRLS with Huber weights (continuous down-weighting, nothing removed)
    #              :hyst    hard rejection, and a rejected satellite stays out for hyst_epochs epochs
    #   taper_deg  weights fade in over this many degrees above the elevation mask (0 = off)
    #   admit_epochs  a satellite joins the fit only after this many consecutive valid epochs
    # receiver clock: one parameter per constellation (default) or one common to all
    # constellations (their epoch-to-epoch differences are about 1 mm, far below the clock noise)
    clock_model::Symbol = :per_system
    gf_tau_s::Float64 = 30.0          # time constant of the geometry-free rate smoother of mode :l1gf [s]
    gf_gap_s::Float64 = 10.0          # the geometry-free history of a satellite survives gaps up to this long [s]
    gf_deadband_mm_s::Float64 = 0.0   # :l1gf corrects only the part of the smoothed rate above this (soft threshold) [mm/s]
    robust::Symbol = :reject
    huber_c::Float64 = 1.5
    hyst_epochs::Int = 5
    taper_deg::Float64 = 0.0
    admit_epochs::Int = 0
    max_reject_frac::Float64 = 0.3    # never reject more than this fraction of satellites
    slip_gf_m::Float64 = 0.05         # geometry-free jump that marks a slip (dual)
    sigma_window::Int = 300           # epochs for the running robust sigma
    tide::Bool = true
    windup::Bool = true
    tropo::Bool = true
    systems::Vector{Char} = ['G', 'R', 'E', 'C']
end

"""Per-satellite line of the solution."""
struct SatResult
    sys::Char
    prn::Int
    elev_deg::Float64
    azim_deg::Float64
    resid_m::Float64           # post-fit residual of the pair observation
    weight::Float64
    used::Bool
    reason::Symbol             # :ok :new :slip :elev :noeph :nofreq :rejected
end

"""Solution of one variometric pair (previous epoch -> this epoch)."""
struct Solution
    t::GpsTime
    dt::Float64                # pair interval [s]
    dxyz::NTuple{3,Float64}    # ECEF displacement over the pair [m] (NaN when no fit)
    dneu::NTuple{3,Float64}
    clk::Dict{Char,Float64}    # receiver clock change per system [m]
    nsat_in::Int
    nsat_used::Int
    nsys::Int
    pdop::Float64
    rms_m::Float64             # weighted post-fit RMS (unit-weight sigma)
    sig_neu::NTuple{3,Float64} # formal 1-sigma of dneu from the fit [m]
    sigma_m::Float64           # windowed robust sigma of the fit RMS history
    sats::Vector{SatResult}
    status::Symbol             # :ok :first :few_sats :no_position :singular
end

nofit(t, dt, status, sats = SatResult[]) =
    Solution(t, dt, (NaN, NaN, NaN), (NaN, NaN, NaN), Dict{Char,Float64}(),
             length(sats), 0, 0, NaN, NaN, (NaN, NaN, NaN), NaN, sats, status)
