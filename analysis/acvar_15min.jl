using JLD2, DataFrames, Statistics, StatsBase, Random

# ------------------------------------------------------------------
# Fast half-life vs AC_Mean and AC_Var (Sugimoto et al.) at 5- and 15-min sampling (DiaMonT)
#   The 5-min series is thinned to every third reading (phase 1).
#   AC_Mean, AC_Var: mean and variance of the ACF of glucose itself
#   (no profile removal, no filtering) at lags 15, 30, …, 150 min
#   (= readings 3, 6, …, 30 at 5 min and 1, 2, …, 10 at 15 min).
#   Fast half-life: as in the paper (28-day profile, 24-h moving average).
#
#   A. Analysis windows (n = 50 and 80): agreement between 5 and 15 min
#      for each measure (median, Spearman, ICC).
#   B. Test–retest reliability (>= 56 consecutive usable days, W = 28):
#      ICC(2,1) at 5 and 15 min for each measure (log for the half-life,
#      original scale for AC_Mean and AC_Var, as in the paper), and the
#      within-participant ICC differences, fast half-life − AC_Mean and
#      fast half-life − AC_Var, with paired bootstrap 95% CIs
#      (2,000 resamples of participants).
#   Only aggregated results are printed. Run in a fresh Julia session.
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const WIN          = 28
const MIN_COMPLETE = 70.0
const MAX_VALID    = 401.5
const B_START      = 29
const NEED_DAYS    = 56
const PER_DAY_5    = 288
const N_BOOT       = 2000
const RATES        = [("5 min", 1), ("15 min", 3)]      # (label, thinning factor); phase 1

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ---------------- Helpers (same definitions as the other scripts) ----------------
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

function halflife_minutes(ac::AbstractVector, step_min::Real; level = 0.5)
    for j in 2:length(ac)
        if isfinite(ac[j]) && ac[j] < level
            a0, a1 = ac[j - 1], ac[j]
            return step_min * ((j - 2) + (a0 - level) / (a0 - a1))
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

function icc21(x::AbstractVector, y::AbstractVector)
    Y = hcat(x, y); n, k = size(Y); gm = mean(Y)
    ssr = k * sum(abs2, mean(Y; dims = 2) .- gm)
    ssc = n * sum(abs2, mean(Y; dims = 1) .- gm)
    sse = sum(abs2, Y .- gm) - ssr - ssc
    msr = ssr / (n - 1); msc = ssc / (k - 1); mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
end

# Fast half-life, AC_Mean and AC_Var of a 5-min glucose vector covering whole days,
# after thinning by factor f (phase 1)
function measures(G5::AbstractVector, ndays::Int, f::Int)
    per_day = PER_DAY_5 ÷ f                             # 288 or 96
    G = G5[1:f:end]
    M = reshape(G, per_day, ndays)
    R = vec(M .- profile_mean(M))
    F = R .- centred_mean(R, per_day ÷ 2)               # ±12 h
    lags = (3 ÷ f):(3 ÷ f):(30 ÷ f)                     # 15–150 min
    acG = acf_nan(G, last(lags))
    acv = [acG[k + 1] for k in lags]
    return (HL_fast = halflife_minutes(acf_nan(F, per_day), 5 * f),
            AC_Mean = mean(acv), AC_Var = var(acv))
end

r2(x) = round(x; digits = 2); r1(x) = round(x; digits = 1)
regimen(r) = r.BasalOnly ? "basal" : "basal-bolus"
fmt(x) = abs(x) < 0.1 ? string(round(x; digits = 3)) : string(r1(x))
base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]

# ---------------- A. Analysis windows ----------------
prim = DataFrame(Subject = Int[], Regimen = String[], Rate = String[],
                 HL_fast = Float64[], AC_Mean = Float64[], AC_Var = Float64[])
for r in eachrow(base)
    g = subjects[r.Subject]
    rows = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
    G5 = copy(g.Glucose[rows]); G5[G5 .> MAX_VALID] .= NaN
    for (lab, f) in RATES
        m = measures(G5, WIN, f)
        push!(prim, (r.Subject, regimen(r), lab, m.HL_fast, m.AC_Mean, m.AC_Var))
    end
end

# ---------------- B. Test–retest reliability, W = 28 ----------------
rel = DataFrame(Subject = Int[], Regimen = String[], Rate = String[],
                HL_A = Float64[], HL_B = Float64[], ACM_A = Float64[], ACM_B = Float64[],
                ACV_A = Float64[], ACV_B = Float64[])
for r in eachrow(base)
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
    L < NEED_DAYS && continue
    run_days = d.Day[s:s+NEED_DAYS-1]
    g = subjects[r.Subject]
    sel = (g.Day .>= run_days[1]) .& (g.Day .<= run_days[end])
    G5 = copy(g.Glucose[sel]); G5[G5 .> MAX_VALID] .= NaN
    M5 = reshape(G5, PER_DAY_5, NEED_DAYS)
    for (lab, f) in RATES
        a = measures(vec(M5[:, 1:WIN]), WIN, f)
        b = measures(vec(M5[:, B_START:B_START+WIN-1]), WIN, f)
        push!(rel, (r.Subject, regimen(r), lab, a.HL_fast, b.HL_fast, a.AC_Mean, b.AC_Mean, a.AC_Var, b.AC_Var))
    end
end

@save joinpath(DIR, "acvar_15min.jld2") prim rel

# ---------------- Summary per cohort ----------------
for grp in ["basal", "basal-bolus"]
    p = prim[prim.Regimen .== grp, :]
    println("\n================ COHORT: ", grp, " (n = ", length(unique(p.Subject)), ") ================")

    println("A. Analysis windows: 15 min vs 5 min")
    p5  = sort(p[p.Rate .== "5 min", :], :Subject)
    p15 = sort(p[p.Rate .== "15 min", :], :Subject)
    println(rpad("Measure", 10), rpad("Median 5 min", 14), rpad("Median 15 min", 15), rpad("Spearman", 10), "ICC(2,1)")
    for (col, logs) in ((:HL_fast, true), (:AC_Mean, false), (:AC_Var, false))
        a, b = p5[!, col], p15[!, col]
        ok = isfinite.(a) .& isfinite.(b) .& (logs ? (a .> 0) .& (b .> 0) : trues(length(a)))
        ia, ib = logs ? (log.(a[ok]), log.(b[ok])) : (a[ok], b[ok])
        println(rpad(string(col), 10), rpad(fmt(median(a[ok])), 14), rpad(fmt(median(b[ok])), 15),
                rpad(r2(corspearman(a[ok], b[ok])), 10), r2(icc21(ia, ib)),
                count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
    end

    q = rel[rel.Regimen .== grp, :]
    println("\nB. Test–retest reliability, W = 28 (n = ", length(unique(q.Subject)),
            " with >= ", NEED_DAYS, " consecutive usable days)")
    println(rpad("Sampling", 10), rpad("ICC fast HL", 13), rpad("ICC AC_Mean", 13), rpad("ICC AC_Var", 12),
            rpad("Fast HL − AC_Mean [95% CI]", 30), "Fast HL − AC_Var [95% CI]")
    for (lab, _) in RATES
        x = q[q.Rate .== lab, :]
        ok = isfinite.(x.HL_A) .& isfinite.(x.HL_B) .& (x.HL_A .> 0) .& (x.HL_B .> 0) .&
             isfinite.(x.ACM_A) .& isfinite.(x.ACM_B) .& isfinite.(x.ACV_A) .& isfinite.(x.ACV_B)
        la, lb = log.(x.HL_A[ok]), log.(x.HL_B[ok])
        ma, mb = x.ACM_A[ok], x.ACM_B[ok]
        va, vb = x.ACV_A[ok], x.ACV_B[ok]
        n = length(la)
        iccs(i) = (icc21(la[i], lb[i]), icc21(ma[i], mb[i]), icc21(va[i], vb[i]))
        rng = MersenneTwister(2026)
        boot = [iccs(rand(rng, 1:n, n)) for _ in 1:N_BOOT]
        ci(f) = (d = filter(isfinite, [f(b) for b in boot]); quantile(d, [0.025, 0.975]))
        cm = ci(b -> b[1] - b[2]); cv = ci(b -> b[1] - b[3])
        h, m, v = iccs(1:n)
        println(rpad(lab, 10), rpad(r2(h), 13), rpad(r2(m), 13), rpad(r2(v), 12),
                rpad(string(r2(h - m), " [", r2(cm[1]), ", ", r2(cm[2]), "]"), 30),
                r2(h - v), " [", r2(cv[1]), ", ", r2(cv[2]), "]",
                count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
    end
end
