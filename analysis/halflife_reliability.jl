using JLD2, DataFrames, Statistics, StatsBase

# ------------------------------------------------------------------
# Reliability of the half-life (test–retest), per cohort
#   For participants with >= 56 consecutive usable days:
#     window A = run days 1 … W
#     window B = run days 29 … 28+W
#   Each window gets its own daily profile and its own half-life.
#   W = 7, 14, 21, 28 shows how many days the half-life needs.
#   Windows never overlap; for W = 28 they are back to back.
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const PER_DAY      = 288
const MIN_COMPLETE = 70.0
const MAX_LAG      = 288
const STEP_MIN     = 5
const MAX_VALID    = 401.5
const WINDOWS      = [7, 14, 21, 28]
const B_START      = 29
const NEED_DAYS    = 56

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ---------------- Helpers (same definitions as persistence_test1.jl) ----------------
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

function longest_run(flags)
    best = 0; best_start = 0; cur = 0; cur_start = 0
    for (i, f) in enumerate(flags)
        if f
            cur == 0 && (cur_start = i)
            cur += 1
            if cur > best
                best = cur; best_start = cur_start
            end
        else
            cur = 0
        end
    end
    return best_start, best
end

halflife_of(M) = halflife_minutes(acf_nan(vec(M .- profile_mean(M)), MAX_LAG))

# ICC(2,1): two-way random effects, absolute agreement, single measurement
function icc21(x::AbstractVector, y::AbstractVector)
    Y = hcat(x, y); n, k = size(Y)
    gm = mean(Y)
    ssr = k * sum(abs2, mean(Y; dims = 2) .- gm)
    ssc = n * sum(abs2, mean(Y; dims = 1) .- gm)
    sse = sum(abs2, Y .- gm) - ssr - ssc
    msr = ssr / (n - 1); msc = ssc / (k - 1); mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
end

# ---------------- Half-lives in both windows ----------------
base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
rel = DataFrame(Subject = Int[], Regimen = String[], W = Int[], HL_A = Float64[], HL_B = Float64[])

for r in eachrow(base)
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
    L < NEED_DAYS && continue

    run_days = d.Day[s:s+NEED_DAYS-1]
    g = subjects[r.Subject]
    rows = (g.Day .>= run_days[1]) .& (g.Day .<= run_days[end])
    G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
    M = reshape(G, PER_DAY, NEED_DAYS)

    for W in WINDOWS
        hlA = halflife_of(M[:, 1:W])
        hlB = halflife_of(M[:, B_START:B_START+W-1])
        push!(rel, (r.Subject, r.BasalOnly ? "basal" : "basal-bolus", W, hlA, hlB))
    end
end

@save joinpath(DIR, "halflife_reliability.jld2") rel

# ---------------- Summary per cohort ----------------
for grp in ["basal", "basal-bolus"]
    dg = rel[rel.Regimen .== grp, :]
    println("\n================ COHORT: ", grp, " (n=", length(unique(dg.Subject)),
            " with >= ", NEED_DAYS, " consecutive usable days) ================")
    println(rpad("W (days)", 10), rpad("ICC(2,1) log", 14), rpad("Spearman", 10),
            rpad("Median HL A", 13), rpad("Median HL B", 13), "Median |A−B| (min)")
    for W in WINDOWS
        x = dg[dg.W .== W, :]
        ok = isfinite.(x.HL_A) .& isfinite.(x.HL_B)
        a, b = x.HL_A[ok], x.HL_B[ok]
        println(rpad(W, 10),
                rpad(round(icc21(log.(a), log.(b)); digits = 2), 14),
                rpad(round(corspearman(a, b); digits = 2), 10),
                rpad(round(median(a); digits = 1), 13),
                rpad(round(median(b); digits = 1), 13),
                round(median(abs.(a .- b)); digits = 1),
                count(.!ok) > 0 ? "   ($(count(.!ok)) without half-life)" : "")
    end
end
