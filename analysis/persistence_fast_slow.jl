using JLD2
using DataFrames
using Statistics
using StatsBase

# ------------------------------------------------------------------
# Fast vs slow persistence
#   The residual ACF decays quickly over ~4 h, then stays at ~0.15
#   out to 24 h. This script separates the two components:
#     (1) filtering: remove a 24-h centred moving average (slow part),
#         then compute the half-life of the remainder (fast part);
#     (2) fitting:   ACF(h) ≈ w·exp(-h/τ) + c on lags 5 min – 24 h;
#     (3) drift:     trend in daily mean residuals across the window.
# ------------------------------------------------------------------
const DIR       = "<path to data>"
const WIN       = 28
const PER_DAY   = 288
const MAX_LAG   = 288
const STEP_MIN  = 5
const MA_HALF   = 144          # ±12 h → 24-h centred moving average
const MAX_VALID = 401.5

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
@load joinpath(DIR, "cohort_v2.jld2") cohort
@load joinpath(DIR, "persistence_v1.jld2") res acfs

# ---------------- Helpers (same as persistence_test1.jl) ----------------
function profile_mean(M::AbstractMatrix)
    mu = fill(NaN, size(M, 1))
    for t in 1:size(M, 1)
        v = filter(isfinite, @view M[t, :])
        isempty(v) || (mu[t] = mean(v))
    end
    return mu
end

function acf_nan(r::AbstractVector, maxlag::Int)
    ok = isfinite.(r)
    m = mean(r[ok]); v = mean((r[ok] .- m) .^ 2)
    out = fill(NaN, maxlag + 1); out[1] = 1.0
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

function halflife_minutes(ac::AbstractVector; level = 0.5)
    for j in 2:length(ac)
        if isfinite(ac[j]) && ac[j] < level
            a0, a1 = ac[j - 1], ac[j]
            return STEP_MIN * ((j - 2) + (a0 - level) / (a0 - a1))
        end
    end
    return NaN
end

# ---------------- New helpers ----------------
# Centred moving average ignoring NaN; requires >= 50% of the window observed
function centred_mean(r::AbstractVector, half::Int)
    n = length(r)
    ok = isfinite.(r)
    cs = vcat(0.0, cumsum(ifelse.(ok, r, 0.0)))
    cc = vcat(0, cumsum(ok))
    out = fill(NaN, n)
    for i in 1:n
        a = max(1, i - half); b = min(n, i + half)
        c = cc[b + 1] - cc[a]
        c >= 0.5 * (b - a + 1) && (out[i] = (cs[b + 1] - cs[a]) / c)
    end
    return out
end

# Fit ACF(h) ≈ w·exp(-h/τ) + c by grid search over τ, least squares for (w, c)
function fit_exp_plateau(ac::AbstractVector)
    h = STEP_MIN .* collect(1:MAX_LAG)
    y = ac[2:end]
    ok = isfinite.(y); h = h[ok]; y = y[ok]
    best_sse = Inf; best = (NaN, NaN, NaN)
    for τ in 10.0:5.0:720.0
        X = hcat(exp.(-h ./ τ), ones(length(h)))
        β = X \ y
        sse = sum(abs2, y .- X * β)
        if sse < best_sse
            best_sse = sse; best = (τ, β[1], β[2])
        end
    end
    τ, w, c = best
    return τ * log(2), w, c        # fast half-life (min), fast weight, plateau
end

# ---------------- Per participant ----------------
eligible = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]

out = DataFrame(Subject = Int[], Regimen = String[],
                HalfLife = Float64[], HalfLifeFast = Float64[], SlowShare = Float64[],
                ExpHalfLife = Float64[], ExpWeight = Float64[], Plateau = Float64[],
                ACF24h = Float64[], DayMeanShare = Float64[],
                TrendSlope = Float64[], TrendRho = Float64[])

for r in eachrow(eligible)
    g = subjects[r.Subject]
    rows = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
    G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
    M = reshape(G, PER_DAY, WIN)
    R = vec(M .- profile_mean(M))

    # (1) Filtering
    slow = centred_mean(R, MA_HALF)
    fast = R .- slow
    ac_fast = acf_nan(fast, MAX_LAG)
    okR = isfinite.(R) .& isfinite.(slow)
    slow_share = var(slow[okR]) / var(R[okR])

    # (2) Curve fit on the original residual ACF
    ac = acfs[r.Subject]
    hl_exp, w, c = fit_exp_plateau(ac)

    # (3) Drift: daily mean residuals over the 28 days (days with no data skipped)
    Rm = reshape(R, PER_DAY, WIN)
    daymean = [any(isfinite, Rm[:, d]) ? mean(filter(isfinite, Rm[:, d])) : NaN for d in 1:WIN]
    okday = isfinite.(daymean)
    days_i = collect(1:WIN)[okday]; dm = daymean[okday]
    slope = cov(days_i, dm) / var(days_i)                  # mg/dL per day
    rho = corspearman(days_i, dm)
    okd = isfinite.(R) .& repeat(okday', PER_DAY)[:]
    daymean_share = var(repeat(daymean', PER_DAY)[:][okd]) / var(R[okd])

    push!(out, (r.Subject, r.BasalOnly ? "basal" : "basal-bolus",
                res.HalfLife[res.Subject .== r.Subject][1],
                halflife_minutes(ac_fast), slow_share,
                hl_exp, w, c, ac[end], daymean_share, slope, rho))
end

@save joinpath(DIR, "persistence_fast_slow.jld2") out

# ---------------- Summary, per cohort (never pooled) ----------------
med(x) = begin
    v = filter(isfinite, x); q = quantile(v, [0.25, 0.5, 0.75])
    string(round(q[2]; digits = 3), " [", round(q[1]; digits = 3), "–", round(q[3]; digits = 3), "]")
end

for grp in ["basal", "basal-bolus"]
    d = out[out.Regimen .== grp, :]
    println("\n================ FAST vs SLOW: ", grp, " (n=", nrow(d), ") ================")
    for col in [:HalfLife, :HalfLifeFast, :ExpHalfLife, :ExpWeight, :Plateau,
                :SlowShare, :DayMeanShare, :ACF24h, :TrendSlope, :TrendRho]
        println(rpad(string(col), 14), med(d[!, col]))
    end

    sp(x, y) = begin
        ok = isfinite.(d[!, x]) .& isfinite.(d[!, y])
        round(corspearman(d[ok, x], d[ok, y]); digits = 2)
    end
    println("\nSpearman correlations:")
    println("  HalfLife vs HalfLifeFast:    ", sp(:HalfLife, :HalfLifeFast))
    println("  HalfLife vs ExpHalfLife:     ", sp(:HalfLife, :ExpHalfLife))
    println("  HalfLifeFast vs ExpHalfLife: ", sp(:HalfLifeFast, :ExpHalfLife))
    println("  HalfLife vs SlowShare:       ", sp(:HalfLife, :SlowShare))
    println("  HalfLifeFast vs SlowShare:   ", sp(:HalfLifeFast, :SlowShare))

    println("\nDrift across the 28-day window:")
    println("  Negative trend (rho < -0.3): ", count(<(-0.3), d.TrendRho), " participants")
    println("  Positive trend (rho >  0.3): ", count(>(0.3),  d.TrendRho), " participants")
    println("  Median slope: ", round(median(filter(isfinite, d.TrendSlope)); digits = 2), " mg/dL per day")
end
