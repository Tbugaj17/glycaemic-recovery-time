using JLD2, DataFrames, Statistics, StatsBase

# ------------------------------------------------------------------
# Reliability of alternative persistence estimators
#   Same two-window design as halflife_reliability.jl:
#     window A = run days 1 … W,  window B = run days 29 … 28+W
#   Estimators (all from the same residuals):
#     HL          crossing of ACF = 0.5 (primary definition)
#     HL_fast     crossing after removing 24-h centred moving average
#     HL_exp      half-life from fit ACF(h) ≈ w·exp(-h/τ) + c
#     AUC4h       area under ACF, 5 min – 4 h (minutes; "integrated time")
#     ACF2h       ACF at lag 120 min
#     Plateau     c from the fit (slow component)
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const PER_DAY      = 288
const MIN_COMPLETE = 70.0
const MAX_LAG      = 288
const STEP_MIN     = 5
const MAX_VALID    = 401.5
const WINDOWS      = [14, 21, 28]
const B_START      = 29
const NEED_DAYS    = 56
const MA_HALF      = 144

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ---------------- Helpers ----------------
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
            x1 = r[i]; x2 = r[i + k]
            if isfinite(x1) && isfinite(x2)
                s += (x1 - m) * (x2 - m); c += 1
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

function centred_mean(r::AbstractVector, half::Int)
    n = length(r); ok = isfinite.(r)
    cs = vcat(0.0, cumsum(ifelse.(ok, r, 0.0))); cc = vcat(0, cumsum(ok))
    out = fill(NaN, n)
    for i in 1:n
        lo = max(1, i - half); hi = min(n, i + half)
        c = cc[hi + 1] - cc[lo]
        c >= 0.5 * (hi - lo + 1) && (out[i] = (cs[hi + 1] - cs[lo]) / c)
    end
    return out
end

function fit_exp_plateau(ac::AbstractVector)
    h = STEP_MIN .* collect(1:MAX_LAG); y = ac[2:end]
    ok = isfinite.(y); h = h[ok]; y = y[ok]
    best_sse = Inf; best = (NaN, NaN, NaN)
    for τ in 10.0:5.0:720.0
        X = hcat(exp.(-h ./ τ), ones(length(h)))
        β = X \ y
        sse = sum(abs2, y .- X * β)
        sse < best_sse && (best_sse = sse; best = (τ, β[1], β[2]))
    end
    τ, w, c = best
    return τ * log(2), c
end

function longest_run(flags)
    best = 0; best_start = 0; cur = 0; cur_start = 0
    for (i, f) in enumerate(flags)
        if f
            cur == 0 && (cur_start = i)
            cur += 1
            cur > best && (best = cur; best_start = cur_start)
        else
            cur = 0
        end
    end
    return best_start, best
end

function estimators(M::AbstractMatrix)
    R = vec(M .- profile_mean(M))
    ac = acf_nan(R, MAX_LAG)
    fast = R .- centred_mean(R, MA_HALF)
    hl_exp, plateau = fit_exp_plateau(ac)
    return (HL      = halflife_minutes(ac),
            HL_fast = halflife_minutes(acf_nan(fast, MAX_LAG)),
            HL_exp  = hl_exp,
            AUC4h   = STEP_MIN * sum(ac[2:49]),        # lags 5–240 min
            ACF2h   = ac[25],                          # lag 120 min
            Plateau = plateau)
end

function icc21(x::AbstractVector, y::AbstractVector)
    Y = hcat(x, y); n, k = size(Y); gm = mean(Y)
    ssr = k * sum(abs2, mean(Y; dims = 2) .- gm)
    ssc = n * sum(abs2, mean(Y; dims = 1) .- gm)
    sse = sum(abs2, Y .- gm) - ssr - ssc
    msr = ssr / (n - 1); msc = ssc / (k - 1); mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
end

const EST = [:HL, :HL_fast, :HL_exp, :AUC4h, :ACF2h, :Plateau]
const LOG_SCALE = Set([:HL, :HL_fast, :HL_exp, :AUC4h])      # ICC on log for time-type measures

# ---------------- Compute in both windows ----------------
base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
rows_out = DataFrame(Subject = Int[], Regimen = String[], W = Int[], Estimator = Symbol[],
                     A = Float64[], B = Float64[])

for r in eachrow(base)
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
    L < NEED_DAYS && continue
    run_days = d.Day[s:s+NEED_DAYS-1]
    g = subjects[r.Subject]
    sel = (g.Day .>= run_days[1]) .& (g.Day .<= run_days[end])
    G = copy(g.Glucose[sel]); G[G .> MAX_VALID] .= NaN
    M = reshape(G, PER_DAY, NEED_DAYS)
    reg = r.BasalOnly ? "basal" : "basal-bolus"

    for W in WINDOWS
        eA = estimators(M[:, 1:W])
        eB = estimators(M[:, B_START:B_START+W-1])
        for e in EST
            push!(rows_out, (r.Subject, reg, W, e, getfield(eA, e), getfield(eB, e)))
        end
    end
end

@save joinpath(DIR, "estimator_reliability.jld2") rows_out

# ---------------- Summary ----------------
for grp in ["basal", "basal-bolus"]
    dg = rows_out[rows_out.Regimen .== grp, :]
    println("\n================ COHORT: ", grp, " (n=", length(unique(dg.Subject)), ") ================")
    println(rpad("Estimator", 10), rpad("W", 5), rpad("ICC(2,1)", 10), rpad("Spearman", 10), "Median A / B")
    for e in EST, W in WINDOWS
        x = dg[(dg.Estimator .== e) .& (dg.W .== W), :]
        ok = isfinite.(x.A) .& isfinite.(x.B) .& (e in LOG_SCALE ? (x.A .> 0) .& (x.B .> 0) : trues(nrow(x)))
        va, vb = x.A[ok], x.B[ok]
        ia, ib = e in LOG_SCALE ? (log.(va), log.(vb)) : (va, vb)
        println(rpad(string(e), 10), rpad(W, 5),
                rpad(round(icc21(ia, ib); digits = 2), 10),
                rpad(round(corspearman(va, vb); digits = 2), 10),
                round(median(va); digits = 2), " / ", round(median(vb); digits = 2),
                count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
    end
end
