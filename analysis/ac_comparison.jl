using JLD2, DataFrames, Statistics, StatsBase, Random, GLM

# ------------------------------------------------------------------
# Comparison with autocorrelation indices of glucose (Sugimoto et al.)
#   AC_Mean, AC_Var: mean and variance of the ACF of glucose ITSELF
#   (no profile removal, no filtering) at lags 15, 30, ..., 150 min
#   (= 3, 6, ..., 30 readings in the 5-min data).
#
#   1. Primary window: Spearman correlations with the fast and unfiltered
#      half-lives and with routine strength
#   2. Reliability: two non-overlapping windows (W = 14, 21, 28 days),
#      ICC(2,1) with bootstrap 95% CI and Spearman, for AC_Mean, AC_Var,
#      fast and unfiltered half-life (this also gives the CIs for the
#      main reliability table)
#   3. Non-redundancy: does AC_Var / AC_Mean explain the fast half-life
#      beyond mean, CV and TIR? And how much of AC_Var do they explain?
#   4. DiaMonT only: Spearman with baseline HbA1c (descriptive)
#
#   Set DATASET below and run. Cohorts are analysed separately.
# ------------------------------------------------------------------
const DATASET = "diamont"          # "diamont" or "replacebg"
const DIR = DATASET == "diamont" ? "<path to data>" :
    "<path to data>"
const DO_HBA1C = DATASET == "diamont"

const WIN          = 28
const PER_DAY      = 288
const STEP_MIN     = 5
const MAX_LAG      = 288
const MA_HALF      = 144
const MAX_VALID    = 401.5
const MIN_COMPLETE = 70.0
const REL_WINDOWS  = [14, 21, 28]
const B_START      = 29
const NEED_DAYS    = 56
const N_BOOT       = 2000
const AC_LAGS      = 3:3:30          # 15–150 min in 5-min steps

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

function icc21_ci(x, y; rng = MersenneTwister(2026))
    n = length(x)
    b = [icc21(x[i], y[i]) for i in (rand(rng, 1:n, n) for _ in 1:N_BOOT)]
    return quantile(filter(isfinite, b), [0.025, 0.975])
end

# All measures for one window matrix (time-of-day × days)
function measures(M::AbstractMatrix)
    G  = vec(M)
    acG = acf_nan(G, maximum(AC_LAGS))                  # glucose itself
    acv = [acG[k + 1] for k in AC_LAGS]
    mu = profile_mean(M)
    R  = vec(M .- mu)
    F  = R .- centred_mean(R, MA_HALF)
    obs = filter(isfinite, G)
    return (AC_Mean = mean(acv), AC_Var = var(acv),
            HL_fast = halflife_minutes(acf_nan(F, MAX_LAG)),
            HL      = halflife_minutes(acf_nan(R, MAX_LAG)),
            RoutineStrength = var(filter(isfinite, mu)) / var(obs),
            Mean = mean(obs), CV = 100 * std(obs) / mean(obs),
            TIR  = 100 * count(x -> 70 <= x <= 180, obs) / length(obs))
end

r2s(x) = round(x; digits = 2)
grp_label(g) = DATASET == "replacebg" ? "T1D (pump)" : g

# ======================= Main =======================
function main()
    @load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    elig = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
    regimen(r) = r.BasalOnly ? "basal" : "basal-bolus"

    function window(sid, startday, ndays)
        g = subjects[sid]
        rows = (g.Day .>= startday) .& (g.Day .<= startday + ndays - 1)
        @assert count(rows) == ndays * PER_DAY
        G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
        return reshape(G, PER_DAY, ndays)
    end

    # ---------------- Baseline HbA1c (DiaMonT only) ----------------
    base_hba = Dict{Int,Float64}()
    if DO_HBA1C
        @eval using ReadStatTables
        lab = DataFrame(Base.invokelatest(readstat, joinpath(DIR, "adapt_ano_labka_a.sas7bdat")))
        lab = lab[lab.Analysekode .== "HbA1ci", :]
        lab.OriginalID = Int.(lab.id_Patient)
        lab.StudyDay = fld.(Float64.(lab.MaalingTidFraInkl), 86400) .+ 1
        for sub in groupby(lab[lab.StudyDay .<= 14, :], :OriginalID)
            base_hba[sub.OriginalID[1]] = sub.InterntResultatNumerisk[argmin(sub.StudyDay)]
        end
    end

    # ---------------- 1. Primary window ----------------
    prim = DataFrame(Subject = Int[], Regimen = String[], AC_Mean = Float64[], AC_Var = Float64[],
                     HL_fast = Float64[], HL = Float64[], RoutineStrength = Float64[],
                     Mean = Float64[], CV = Float64[], TIR = Float64[], HbA1c = Float64[])
    for r in eachrow(elig)
        m = measures(window(r.Subject, r.WindowStartDay, WIN))
        push!(prim, (r.Subject, regimen(r), m.AC_Mean, m.AC_Var, m.HL_fast, m.HL,
                     m.RoutineStrength, m.Mean, m.CV, m.TIR, get(base_hba, r.OriginalID, NaN)))
    end

    # ---------------- 2. Reliability windows ----------------
    rel = DataFrame(Subject = Int[], Regimen = String[], W = Int[], Measure = Symbol[],
                    A = Float64[], B = Float64[])
    for r in eachrow(elig)
        d = sort(days[days.Subject .== r.Subject, :], :Day)
        s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
        L < NEED_DAYS && continue
        M56 = window(r.Subject, d.Day[s], NEED_DAYS)
        for W in REL_WINDOWS
            mA = measures(M56[:, 1:W]); mB = measures(M56[:, B_START:B_START+W-1])
            for e in (:HL_fast, :HL, :AC_Mean, :AC_Var)
                push!(rel, (r.Subject, regimen(r), W, e, getfield(mA, e), getfield(mB, e)))
            end
        end
    end

    @save joinpath(DIR, "ac_comparison.jld2") prim rel

    # ---------------- Output per cohort ----------------
    for grp in unique(prim.Regimen)
        p = prim[prim.Regimen .== grp, :]
        println("\n================ ", uppercase(DATASET), " – ", grp_label(grp), " (n = ", nrow(p), ") ================")

        println("\n---- 1. Primary window: medians [IQR] ----")
        for col in (:AC_Mean, :AC_Var, :HL_fast, :HL)
            q = quantile(filter(isfinite, p[!, col]), [0.25, 0.5, 0.75])
            println("  ", rpad(string(col), 10), round(q[2]; sigdigits = 3), " [",
                    round(q[1]; sigdigits = 3), "–", round(q[3]; sigdigits = 3), "]")
        end
        sp(a, b) = begin
            ok = isfinite.(p[!, a]) .& isfinite.(p[!, b])
            r2s(corspearman(p[ok, a], p[ok, b]))
        end
        println("\n  Spearman correlations")
        println("  ", rpad("", 18), rpad("AC_Mean", 10), rpad("AC_Var", 10), "HL_fast")
        for (lab, col) in [("Fast half-life", :HL_fast), ("Unfiltered HL", :HL),
                           ("Routine strength", :RoutineStrength), ("Mean glucose", :Mean),
                           ("CV", :CV), ("TIR", :TIR)]
            println("  ", rpad(lab, 18), rpad(sp(:AC_Mean, col), 10), rpad(sp(:AC_Var, col), 10),
                    col == :HL_fast ? "–" : sp(:HL_fast, col))
        end
        println("  AC_Mean vs AC_Var: ", sp(:AC_Mean, :AC_Var))

        println("\n---- 2. Reliability (n = ", length(unique(rel.Subject[rel.Regimen .== grp])),
                " with ≥ ", NEED_DAYS, " consecutive usable days) ----")
        println("  ", rpad("Measure", 10), rpad("W", 5), rpad("ICC(2,1) [95% CI]", 24), "Spearman")
        for e in (:HL_fast, :HL, :AC_Mean, :AC_Var), W in REL_WINDOWS
            x = rel[(rel.Regimen .== grp) .& (rel.Measure .== e) .& (rel.W .== W), :]
            logs = e in (:HL_fast, :HL)
            ok = isfinite.(x.A) .& isfinite.(x.B) .& (logs ? (x.A .> 0) .& (x.B .> 0) : trues(nrow(x)))
            a, b = x.A[ok], x.B[ok]
            ia, ib = logs ? (log.(a), log.(b)) : (a, b)
            ci = icc21_ci(ia, ib)
            println("  ", rpad(string(e), 10), rpad(W, 5),
                    rpad(string(r2s(icc21(ia, ib)), " [", r2s(ci[1]), "–", r2s(ci[2]), "]"), 24),
                    r2s(corspearman(a, b)), count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
        end

        println("\n---- 3. Non-redundancy (adjusted rank R²) ----")
        ok = isfinite.(p.HL_fast) .& isfinite.(p.AC_Var) .& isfinite.(p.AC_Mean)
        rk = DataFrame(HL = tiedrank(p.HL_fast[ok]), Mean = tiedrank(p.Mean[ok]),
                       CV = tiedrank(p.CV[ok]), TIR = tiedrank(p.TIR[ok]),
                       ACV = tiedrank(p.AC_Var[ok]), ACM = tiedrank(p.AC_Mean[ok]))
        println("  Fast HL ~ mean + CV + TIR:                    ", r2s(adjr2(lm(@formula(HL ~ Mean + CV + TIR), rk))))
        println("  Fast HL ~ mean + CV + TIR + AC_Var:           ", r2s(adjr2(lm(@formula(HL ~ Mean + CV + TIR + ACV), rk))))
        println("  Fast HL ~ mean + CV + TIR + AC_Var + AC_Mean: ", r2s(adjr2(lm(@formula(HL ~ Mean + CV + TIR + ACV + ACM), rk))))
        println("  AC_Var  ~ mean + CV + TIR:                    ", r2s(adjr2(lm(@formula(ACV ~ Mean + CV + TIR), rk))))
        println("  AC_Mean ~ mean + CV + TIR:                    ", r2s(adjr2(lm(@formula(ACM ~ Mean + CV + TIR), rk))))

        if DO_HBA1C
            println("\n---- 4. Spearman with baseline HbA1c (descriptive) ----")
            for col in (:HL_fast, :AC_Mean, :AC_Var)
                ok = isfinite.(p[!, col]) .& isfinite.(p.HbA1c)
                println("  ", rpad(string(col), 10), r2s(corspearman(p[ok, col], p.HbA1c[ok])), "  (n = ", count(ok), ")")
            end
        end
    end
    println("\nSaved to: ", joinpath(DIR, "ac_comparison.jld2"))
end

main()
