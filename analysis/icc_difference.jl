using JLD2, DataFrames, Statistics, Random

# ------------------------------------------------------------------
# Is the fast half-life more reliable than the autocorrelation indices?
#   Within-participant comparison of ICC(2,1) at 28 days, same participants
#   and windows (window A: days 1-28, window B: days 29-56 of the longest
#   run of usable days). Paired bootstrap of the ICC difference.
#   Computed directly from the CGM data (no saved results needed).
#   Set DATASET and run; repeat for the other dataset.
# ------------------------------------------------------------------
const DATASET = "diamont"          # "diamont" or "replacebg"
const DIR = DATASET == "diamont" ? "<path to data>" :
    "<path to data>"
const OUT = "<path to data>"   # results saved locally

const PER_DAY      = 288
const STEP_MIN     = 5
const MAX_LAG      = 288
const MA_HALF      = 144
const MAX_VALID    = 401.5
const MIN_COMPLETE = 70.0
const W            = 28
const NEED_DAYS    = 56
const N_BOOT       = 2000
const AC_LAGS      = 3:3:30          # 15-150 min

# ======================= Helpers (identical definitions) =======================
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

function longest_run(flags)
    best = 0; bs = 0; cur = 0; cs = 0
    for (i, f) in enumerate(flags)
        if f
            cur == 0 && (cs = i)
            cur += 1
            cur > best && (best = cur; bs = cs)
        else
            cur = 0
        end
    end
    return bs, best
end

function icc21(x::AbstractVector, y::AbstractVector)
    Y = hcat(x, y); n, k = size(Y); gm = mean(Y)
    ssr = k * sum(abs2, mean(Y; dims = 2) .- gm)
    ssc = n * sum(abs2, mean(Y; dims = 1) .- gm)
    sse = sum(abs2, Y .- gm) - ssr - ssc
    msr = ssr / (n - 1); msc = ssc / (k - 1); mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
end

function measures(M::AbstractMatrix)
    acG = acf_nan(vec(M), maximum(AC_LAGS))
    acv = [acG[k + 1] for k in AC_LAGS]
    R = vec(M .- profile_mean(M))
    F = R .- centred_mean(R, MA_HALF)
    return (HL_fast = halflife_minutes(acf_nan(F, MAX_LAG)),
            HL      = halflife_minutes(acf_nan(R, MAX_LAG)),
            AC_Mean = mean(acv), AC_Var = var(acv))
end

r2(x) = round(x; digits = 2)
label(g) = DATASET == "replacebg" ? "T1D (pump)" : g

# ======================= Main =======================
function main()
    @load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    elig = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]

    # ---------------- Compute measures in windows A and B ----------------
    rows = DataFrame(Subject = Int[], Regimen = String[],
                     fA = Float64[], fB = Float64[], uA = Float64[], uB = Float64[],
                     mA = Float64[], mB = Float64[], vA = Float64[], vB = Float64[])
    for r in eachrow(elig)
        d = sort(days[days.Subject .== r.Subject, :], :Day)
        s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
        L < NEED_DAYS && continue
        g = subjects[r.Subject]
        sel = (g.Day .>= d.Day[s]) .& (g.Day .<= d.Day[s] + NEED_DAYS - 1)
        G = copy(g.Glucose[sel]); G[G .> MAX_VALID] .= NaN
        M = reshape(G, PER_DAY, NEED_DAYS)
        a = measures(M[:, 1:W]); b = measures(M[:, W+1:2W])
        push!(rows, (r.Subject, r.BasalOnly ? "basal" : "basal-bolus",
                     a.HL_fast, b.HL_fast, a.HL, b.HL, a.AC_Mean, b.AC_Mean, a.AC_Var, b.AC_Var))
    end
    @save joinpath(OUT, "icc_difference_" * DATASET * ".jld2") rows

    # ---------------- ICCs and paired bootstrap differences ----------------
    for grp in unique(rows.Regimen)
        d = rows[rows.Regimen .== grp, :]
        d = d[isfinite.(d.fA) .& isfinite.(d.fB) .& (d.fA .> 0) .& (d.fB .> 0) .&
              isfinite.(d.mA) .& isfinite.(d.mB) .& isfinite.(d.vA) .& isfinite.(d.vB), :]
        n = nrow(d)
        stats(i) = (fast = icc21(log.(d.fA[i]), log.(d.fB[i])),
                    unf  = all(isfinite.(d.uA[i]) .& isfinite.(d.uB[i]) .& (d.uA[i] .> 0) .& (d.uB[i] .> 0)) ?
                           icc21(log.(d.uA[i]), log.(d.uB[i])) : NaN,
                    acm  = icc21(d.mA[i], d.mB[i]),
                    acv  = icc21(d.vA[i], d.vB[i]))
        est = stats(1:n)
        rng = MersenneTwister(2026)
        boot = [stats(rand(rng, 1:n, n)) for _ in 1:N_BOOT]
        ci(f) = quantile(filter(isfinite, [f(b) for b in boot]), [0.025, 0.975])

        println("\n================ ", uppercase(DATASET), " – ", label(grp), " (n = ", n, ", 28-day windows) ================")
        for (lab, key) in [("Fast half-life", :fast), ("Unfiltered half-life", :unf),
                           ("AC_Mean", :acm), ("AC_Var", :acv)]
            c = ci(b -> getfield(b, key))
            println("  ICC ", rpad(lab, 22), r2(getfield(est, key)), " [", r2(c[1]), ", ", r2(c[2]), "]")
        end
        for (lab, key) in [("AC_Mean", :acm), ("AC_Var", :acv), ("unfiltered HL", :unf)]
            diffs = [b.fast - getfield(b, key) for b in boot]
            ok = filter(isfinite, diffs)
            c = quantile(ok, [0.025, 0.975])
            println("  Difference, fast HL − ", rpad(lab, 14), r2(est.fast - getfield(est, key)),
                    " [", r2(c[1]), ", ", r2(c[2]), "]   share of bootstrap differences > 0: ",
                    round(mean(ok .> 0); digits = 3))
        end
    end
    println("\nSaved to: ", joinpath(OUT, "icc_difference_" * DATASET * ".jld2"))
end

main()
