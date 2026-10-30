# Reader for the observation dump written by crx_to_dump.py: RINEX observations
# (e.g. a receiver's own RINEX log) turned into GaiaSeis MsmEpoch objects, so
# both engines get the observations through exactly the adapters used for RTCM.
#
#   E <gps sow>
#   S G07 C1C=21234567.123 L1C=111111111.123:0:7 S1C=45.250 ...
#
# Expects GaiaSeis, GV (GaiaVario) and live_common.jl's helpers to be loaded.

const _OBS_SYS_LABEL = Dict('G' => 1077, 'R' => 1087, 'E' => 1097, 'C' => 1127)

function read_obs_dump(path::AbstractString, week::Int, glo_k::Dict{Int,Int})
    eps = GaiaSeis.MsmEpoch[]
    lock = Dict{Tuple{Char,Int,String},Int}()
    cur_sow = NaN
    seen = Set{Tuple{Char,Int,String}}()
    bysys = Dict{Char,Vector{GaiaSeis.MsmObservation}}()
    flush_epoch() = begin
        isnan(cur_sow) && return
        for k in collect(keys(lock)); k in seen || delete!(lock, k); end     # signals that vanished lose their lock count
        empty!(seen)
        for sys in ('G', 'R', 'E', 'C')
            obs = get(bysys, sys, nothing)
            (obs === nothing || isempty(obs)) && continue
            if sys == 'G' || sys == 'E'
                ems = cur_sow * 1000.0; sod = mod(cur_sow, 86400.0)
            elseif sys == 'C'
                b = cur_sow - 14.0; b < 0 && (b += 604800.0)
                ems = b * 1000.0; sod = mod(b, 86400.0)
            else                                  # GLONASS: day of week and time of day in UTC+3h
                u = cur_sow - 18.0 + 3 * 3600.0
                u < 0 && (u += 604800.0)
                dow = floor(u / 86400.0); tod = u - dow * 86400.0
                ems = dow * 2.0^27 + tod * 1000.0
                sod = mod(tod - 3 * 3600.0, 86400.0)
            end
            push!(eps, GaiaSeis.MsmEpoch(sys, 1, ems, sod, false, obs, UInt64(0), UInt32(0), _OBS_SYS_LABEL[sys]))
        end
        empty!(bysys)
    end
    for line in eachline(path)
        isempty(line) && continue
        if line[1] == 'E'
            flush_epoch()
            cur_sow = parse(Float64, split(line)[2])
        elseif line[1] == 'S'
            f = split(line)
            sat = f[2]; sys = sat[1]; prn = parse(Int, sat[2:3])
            sys in ('G', 'R', 'E', 'C') || continue
            vals = Dict{String,Tuple{Float64,Int,Float64}}()      # code -> (value, lli, snr)
            for tok in f[3:end]
                code, rest = split(tok, '='; limit = 2)
                parts = split(rest, ':')
                vals[code] = (parse(Float64, parts[1]), length(parts) >= 2 ? parse(Int, parts[2]) : 0, NaN)
            end
            for code in keys(vals)
                (code[1] == 'L') || continue
                suf = code[2:3]
                slot = parse(Int, suf[1])
                k = sys == 'R' ? get(glo_k, prn, 0) : 0
                fq = GV.carrier_freq(sys, slot, k)
                isfinite(fq) || continue
                ph, lli, _ = vals[code]
                key = (sys, prn, suf); push!(seen, key)
                lock[key] = (lli & 1 != 0 || !haskey(lock, key)) ? 0 : min(lock[key] + 1, 1023)
                pr = haskey(vals, "C" * suf) ? vals["C" * suf][1] : NaN
                snr = haskey(vals, "S" * suf) ? vals["S" * suf][1] : 45.0
                push!(get!(bysys, sys, GaiaSeis.MsmObservation[]),
                      GaiaSeis.MsmObservation(sys, prn, suf, slot, pr, ph * GV.C_LIGHT / fq, lock[key], false, snr, NaN))
            end
            # a satellite-code that vanished from this epoch loses its lock counter
        end
    end
    flush_epoch()
    return eps
end
