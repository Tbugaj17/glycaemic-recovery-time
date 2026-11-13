using JLD2
using DataFrames
using Statistics
using Random
using HypothesisTests      # if missing:  ] add HypothesisTests

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------
const DIR      = "<path to data>"
const WIN      = 28        # days in development window
const PER_DAY  = 288       # 5-min slots per day
const MAX_LAG  = 288       # autocorrelation up to 24 h
const STEP_MIN = 5
const MAX_VALID_MGDL = 401.5   # above Dexcom's 401 "HIGH" code = artifact

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
# Mean across days at each time-of-day slot, ignoring NaN
function profile_mean(M::AbstractMatrix)
    mu = fill(NaN, size(M, 1))
    for t in 1:size(M, 1)
        v = filter(isfinite, @view M[t, :])
        isempty(v) || (mu[t] = mean(v))
    end
    return mu
end

# Autocorrelation using only pairs where both values are observed.
# out[k+1] is the autocorrelation at lag k.
function acf_nan(r::AbstractVector, maxlag::Int)
    ok = isfinite.(r)
    m = mean(r[ok])
    v = mean((r[ok] .- m) .^ 2)
    out = fill(NaN, maxlag + 1)
    out[1] = 1.0
    n = length(r)
    for k in 1:maxlag
        s = 0.0; c = 0
        @inbounds for i in 1:(n - k)
            a = r[i]; b = r[i + k]
            if isfinite(a) && isfinite(b)
                s += (a - m) * (b - m); c += 1
            end
        end
        c > 30 && (out[k + 1] = s / c / v)
    end
    return out
end

# Empirical half-life: first lag where ACF drops below 0.5 (interpolated), in minutes
function halflife_minutes(ac::AbstractVector; level = 0.5)
    for j in 2:length(ac)
        if isfinite(ac[j]) && ac[j] < level
            a0, a1 = ac[j - 1], ac[j]
            return STEP_MIN * ((j - 2) + (a0 - level) / (a0 - a1))
        end
    end
    return NaN   # never crossed within MAX_LAG
end

# ------------------------------------------------------------------
# Per-participant persistence metrics
# ------------------------------------------------------------------
eligible = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]

res = DataFrame(Subject = Int[], Regimen = String[],
                HalfLife = Float64[], HalfLifeAR1 = Float64[], Rho1 = Float64[],
                MinACF6h = Float64[], SigmaR = Float64[], RoutineStrength = Float64[],
                MeanGlucose = Float64[], CV = Float64[],
                TIR = Float64[], TAR = Float64[], TBR = Float64[], PctObserved = Float64[])
acfs = Dict{Int,Vector{Float64}}()

for r in eachrow(eligible)
    g = subjects[r.Subject]
    rows = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
    @assert count(rows) == WIN * PER_DAY
    @assert g.Measurement[findfirst(rows)] == 1       # window starts at 00:00

    G = copy(g.Glucose[rows])                         # mg/dL
    G[G .> MAX_VALID_MGDL] .= NaN

    M  = reshape(G, PER_DAY, WIN)                     # time-of-day × day
    mu = profile_mean(M)                              # daily profile
    R  = vec(M .- mu)                                 # residuals, NaN kept

    ac = acf_nan(R, MAX_LAG)
    acfs[r.Subject] = ac
    rho1 = ac[2]

    obs  = filter(isfinite, G)
    robs = filter(isfinite, R)

    push!(res, (r.Subject,
                r.BasalOnly ? "basal" : "basal-bolus",
                halflife_minutes(ac),
                0 < rho1 < 1 ? STEP_MIN * log(0.5) / log(rho1) : NaN,
                rho1,
                minimum(filter(isfinite, ac[2:73])),            # lags 5 min – 6 h
                std(robs),
                var(mu) / var(obs),
                mean(obs),
                100 * std(obs) / mean(obs),
                100 * count(x -> 70 <= x <= 180, obs) / length(obs),
                100 * count(>(180), obs) / length(obs),
                100 * count(<(70), obs) / length(obs),
                100 * length(obs) / length(G)))
end

@save joinpath(DIR, "persistence_v1.jld2") res acfs
println("Participants analysed: ", nrow(res),
        "   (half-life not reached within 24 h: ", count(isnan, res.HalfLife), ")")

# ------------------------------------------------------------------
# Group summaries
# ------------------------------------------------------------------
iqr_str(x; digits = 1) = begin
    v = filter(isfinite, x)
    q = quantile(v, [0.25, 0.5, 0.75])
    string(round(q[2]; digits = digits), " [", round(q[1]; digits = digits), "–",
           round(q[3]; digits = digits), "]")
end

# Columns on a 0–1 scale need more decimals to be readable
const FINE_COLS = (:Rho1, :MinACF6h, :RoutineStrength)

println("\n================ MEDIAN [IQR] BY REGIMEN ================")
for col in [:HalfLife, :HalfLifeAR1, :Rho1, :MinACF6h, :SigmaR, :RoutineStrength,
            :MeanGlucose, :CV, :TIR, :TAR, :TBR]
    b  = res[res.Regimen .== "basal", col]
    bb = res[res.Regimen .== "basal-bolus", col]
    dg = col in FINE_COLS ? 3 : 1
    println(rpad(string(col), 16), rpad("basal (n=$(length(b))): " * iqr_str(b; digits = dg), 40),
            "basal-bolus (n=$(length(bb))): ", iqr_str(bb; digits = dg))
end

# ------------------------------------------------------------------
# Test 1: half-life, basal-only vs basal-bolus
# ------------------------------------------------------------------
hl_b  = filter(isfinite, res.HalfLife[res.Regimen .== "basal"])
hl_bb = filter(isfinite, res.HalfLife[res.Regimen .== "basal-bolus"])

mw = MannWhitneyUTest(hl_bb, hl_b)

rng = MersenneTwister(2026)
boot = [median(rand(rng, hl_bb, length(hl_bb))) - median(rand(rng, hl_b, length(hl_b)))
        for _ in 1:10_000]
ci = quantile(boot, [0.025, 0.975])

println("\n================ TEST 1: HALF-LIFE BY REGIMEN ================")
println("Median half-life, basal-only:  ", round(median(hl_b);  digits = 1), " min (n=", length(hl_b), ")")
println("Median half-life, basal-bolus: ", round(median(hl_bb); digits = 1), " min (n=", length(hl_bb), ")")
println("Difference (bolus − basal):    ", round(median(hl_bb) - median(hl_b); digits = 1),
        " min, 95% bootstrap CI [", round(ci[1]; digits = 1), ", ", round(ci[2]; digits = 1), "]")
println("Mann–Whitney p-value:          ", round(pvalue(mw); digits = 4))
