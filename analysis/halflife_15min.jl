using JLD2, DataFrames, Statistics, StatsBase

# ------------------------------------------------------------------
# Fast half-life at 15-min sampling compared with 5-min sampling (DiaMonT)
#   The 5-min series is thinned to every third reading. Starting point
#   (phase) 1, 2 or 3 is used to check that the choice does not matter.
#   Everything else is unchanged: 28-day daily profile, 24-h centred
#   moving average (±12 h), ACF to 24 h, linear interpolation at 0.5.
#
#   A. Analysis windows (basal-only n = 50, basal-bolus n = 80):
#      5-min vs 15-min fast half-life (median, median difference,
#      Spearman, ICC on log), and agreement between phases.
#   B. Test–retest reliability (>= 56 consecutive usable days):
#      window A = run days 1…W, window B = run days 29…28+W, at 5 and 15 min.
#   Participants and days are selected with the 5-min usable-day rule,
#   so the results are directly comparable with the paper.
#   Only aggregated results are printed. Run in a fresh Julia session.
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const WIN          = 28
const MIN_COMPLETE = 70.0
const MAX_VALID    = 401.5
const WINDOWS      = [14, 21, 28]
const B_START      = 29
const NEED_DAYS    = 56
const PER_DAY_5    = 288

# Sampling settings: (label, thinning factor, phase)
const RATES = [("5 min", 1, 1), ("15 min, phase 1", 3, 1), ("15 min, phase 2", 3, 2), ("15 min, phase 3", 3, 3)]

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

# Step size in minutes is an argument here (5 or 15)
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

# Fast half-life of a 5-min glucose vector covering whole days,
# after thinning by `f` from starting point `phase`
function fast_halflife(G5::AbstractVector, ndays::Int, f::Int, phase::Int)
    per_day = PER_DAY_5 ÷ f                       # 288 or 96
    G = G5[phase:f:end]
    M = reshape(G, per_day, ndays)
    R = vec(M .- profile_mean(M))
    fast = R .- centred_mean(R, per_day ÷ 2)      # ±12 h
    return halflife_minutes(acf_nan(fast, per_day), 5 * f)   # lags to 24 h
end

ok2(a, b) = isfinite.(a) .& isfinite.(b) .& (a .> 0) .& (b .> 0)
r2(x) = round(x; digits = 2); r1(x) = round(x; digits = 1)
regimen(r) = r.BasalOnly ? "basal" : "basal-bolus"

base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]

# ---------------- A. Analysis windows ----------------
prim = DataFrame(Subject = Int[], Regimen = String[], Rate = String[], HL = Float64[])
for r in eachrow(base)
    g = subjects[r.Subject]
    rows = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
    G5 = copy(g.Glucose[rows]); G5[G5 .> MAX_VALID] .= NaN
    for (lab, f, ph) in RATES
        push!(prim, (r.Subject, regimen(r), lab, fast_halflife(G5, WIN, f, ph)))
    end
end

# ---------------- B. Test–retest reliability ----------------
rel = DataFrame(Subject = Int[], Regimen = String[], Rate = String[], W = Int[], A = Float64[], B = Float64[])
for r in eachrow(base)
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
    L < NEED_DAYS && continue
    run_days = d.Day[s:s+NEED_DAYS-1]
    g = subjects[r.Subject]
    sel = (g.Day .>= run_days[1]) .& (g.Day .<= run_days[end])
    G5 = copy(g.Glucose[sel]); G5[G5 .> MAX_VALID] .= NaN
    M5 = reshape(G5, PER_DAY_5, NEED_DAYS)
    for (lab, f, ph) in RATES, W in WINDOWS
        push!(rel, (r.Subject, regimen(r), lab, W,
                    fast_halflife(vec(M5[:, 1:W]), W, f, ph),
                    fast_halflife(vec(M5[:, B_START:B_START+W-1]), W, f, ph)))
    end
end

@save joinpath(DIR, "halflife_15min.jld2") prim rel

# ---------------- Summary per cohort ----------------
hl_of(df, lab) = sort(df[df.Rate .== lab, :], :Subject).HL

for grp in ["basal", "basal-bolus"]
    p = prim[prim.Regimen .== grp, :]
    println("\n================ COHORT: ", grp, " (n = ", length(unique(p.Subject)), ") ================")

    println("A. Analysis windows: 15 min vs 5 min")
    h5 = hl_of(p, "5 min")
    println(rpad("Sampling", 18), rpad("Median HL", 11), rpad("Median diff", 13),
            rpad("Spearman", 10), "ICC(2,1) log")
    for (lab, _, _) in RATES
        h = hl_of(p, lab); ok = ok2(h5, h)
        println(rpad(lab, 18), rpad(r1(median(h[ok])), 11),
                rpad(r1(median(h[ok] .- h5[ok])), 13),
                rpad(r2(corspearman(h5[ok], h[ok])), 10),
                r2(icc21(log.(h5[ok]), log.(h[ok]))),
                count(.!ok) > 0 ? "   ($(count(.!ok)) without half-life)" : "")
    end
    println("Agreement between 15-min phases (Spearman):")
    for (a, b) in [(2, 3), (2, 4), (3, 4)]
        ha, hb = hl_of(p, RATES[a][1]), hl_of(p, RATES[b][1]); ok = ok2(ha, hb)
        println("  ", RATES[a][1], " vs ", RATES[b][1], ": ", r2(corspearman(ha[ok], hb[ok])))
    end

    q = rel[rel.Regimen .== grp, :]
    println("\nB. Test–retest reliability (n = ", length(unique(q.Subject)), " with >= ", NEED_DAYS,
            " consecutive usable days)")
    println(rpad("Sampling", 18), rpad("W", 5), rpad("ICC(2,1) log", 14), rpad("Spearman", 10), "Median A / B (min)")
    for (lab, _, _) in RATES, W in WINDOWS
        x = q[(q.Rate .== lab) .& (q.W .== W), :]
        ok = ok2(x.A, x.B); va, vb = x.A[ok], x.B[ok]
        println(rpad(lab, 18), rpad(W, 5), rpad(r2(icc21(log.(va), log.(vb))), 14),
                rpad(r2(corspearman(va, vb)), 10), r1(median(va)), " / ", r1(median(vb)),
                count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
    end
end
