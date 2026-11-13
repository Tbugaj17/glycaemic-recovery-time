using JLD2, DataFrames, Statistics, StatsBase

# ------------------------------------------------------------------
# Filter-width sensitivity of the fast half-life
#   The slow component is removed with a centred moving average.
#   This script repeats the reliability test for three widths:
#     12 h (MA_HALF = 72), 24 h (144, primary), 48 h (288)
#   and checks whether participants keep their ranking across widths.
#   Two-window design as before: A = run days 1…W, B = run days 29…28+W
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
const MA_WIDTHS    = [72, 144, 288]             # half-widths in 5-min steps: ±6 h, ±12 h, ±24 h
# Label shows total window and half-width, e.g. "24h (±12h)"
width_label(h) = string(2 * h * STEP_MIN ÷ 60, "h (±", h * STEP_MIN ÷ 60, "h)")

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

function fast_halflife(M::AbstractMatrix, half::Int)
    R = vec(M .- profile_mean(M))
    return halflife_minutes(acf_nan(R .- centred_mean(R, half), MAX_LAG))
end

function icc21(x::AbstractVector, y::AbstractVector)
    Y = hcat(x, y); n, k = size(Y); gm = mean(Y)
    ssr = k * sum(abs2, mean(Y; dims = 2) .- gm)
    ssc = n * sum(abs2, mean(Y; dims = 1) .- gm)
    sse = sum(abs2, Y .- gm) - ssr - ssc
    msr = ssr / (n - 1); msc = ssc / (k - 1); mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
end

# ---------------- Compute ----------------
base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
sens = DataFrame(Subject = Int[], Regimen = String[], MA = Int[], W = Int[],
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

    for h in MA_WIDTHS, W in WINDOWS
        push!(sens, (r.Subject, reg, h, W,
                     fast_halflife(M[:, 1:W], h),
                     fast_halflife(M[:, B_START:B_START+W-1], h)))
    end
end

@save joinpath(DIR, "filter_sensitivity.jld2") sens

# ---------------- Summary per cohort ----------------
for grp in ["basal", "basal-bolus"]
    dg = sens[sens.Regimen .== grp, :]
    println("\n================ COHORT: ", grp, " (n=", length(unique(dg.Subject)), ") ================")

    println("Reliability of the fast half-life by filter width")
    println(rpad("Filter", 14), rpad("W", 5), rpad("ICC(2,1) log", 14), rpad("Spearman", 10), "Median A / B (min)")
    for h in MA_WIDTHS, W in WINDOWS
        x = dg[(dg.MA .== h) .& (dg.W .== W), :]
        ok = isfinite.(x.A) .& isfinite.(x.B) .& (x.A .> 0) .& (x.B .> 0)
        va, vb = x.A[ok], x.B[ok]
        println(rpad(width_label(h), 14), rpad(W, 5),
                rpad(round(icc21(log.(va), log.(vb)); digits = 2), 14),
                rpad(round(corspearman(va, vb); digits = 2), 10),
                round(median(va); digits = 1), " / ", round(median(vb); digits = 1),
                count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
    end

    # Do participants keep their ranking across filter widths? (window A, W = 28)
    println("\nRank agreement across filter widths (window A, W = 28, Spearman):")
    hl = Dict(h => sort(dg[(dg.MA .== h) .& (dg.W .== 28), :], :Subject).A for h in MA_WIDTHS)
    for (h1, h2) in [(72, 144), (144, 288), (72, 288)]
        ok = isfinite.(hl[h1]) .& isfinite.(hl[h2])
        println("  ", width_label(h1), " vs ", width_label(h2), ": ",
                round(corspearman(hl[h1][ok], hl[h2][ok]); digits = 2))
    end
end
