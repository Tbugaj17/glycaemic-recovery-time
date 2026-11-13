using JLD2, DataFrames, Statistics, StatsBase, Random, HypothesisTests, GLM

# ------------------------------------------------------------------
# T1D external validation (REPLACE-BG), one cohort of pump users
#   Same definitions as the DiaMonT scripts:
#     persistence_test1.jl, persistence_fast_slow.jl, estimator_reliability.jl,
#     filter_sensitivity.jl, fast_nonredundancy.jl, persistence_profile.jl
#   Input : cgm_subjects_v2.jld2 and cohort_v2.jld2 from build_replacebg.jl
#   Output: t1d_results.jld2 in the same folder
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const WIN          = 28
const PER_DAY      = 288
const STEP_MIN     = 5
const MAX_LAG      = 288          # 24 h
const MA_HALF      = 144          # 24-h centred moving average (primary)
const MA_WIDTHS    = [72, 144, 288]   # 12 h, 24 h, 48 h
const MAX_VALID    = 401.5
const MIN_COMPLETE = 70.0
const REL_WINDOWS  = [14, 21, 28]
const B_START      = 29
const NEED_DAYS    = 56
const BIN_SLOTS    = 24           # 2-hour blocks
const NBINS        = PER_DAY ÷ BIN_SLOTS
const MAX_LAG_BLK  = 96           # 8 h
const MIN_PAIRS    = 50
const N_BOOT       = 2000
const ACF_LAGS     = [(6, "30 min"), (12, "1 h"), (24, "2 h"), (48, "4 h"),
                      (72, "6 h"), (144, "12 h"), (288, "24 h")]

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
    return τ * log(2), w, c
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

block_of(meas) = (meas - 1) ÷ BIN_SLOTS + 1
block_label(b) = string(lpad(2 * (b - 1), 2, '0'), "–", lpad(2 * b, 2, '0'))

function block_acf(F::AbstractVector, meas::AbstractVector, b::Int)
    n = length(F)
    origins = [i for i in 1:n if isfinite(F[i]) && block_of(meas[i]) == b]
    ac = fill(NaN, MAX_LAG_BLK + 1); ac[1] = 1.0
    xs = Vector{Float64}(undef, length(origins)); ys = similar(xs)
    for k in 1:MAX_LAG_BLK
        m = 0
        for i in origins
            j = i + k
            if j <= n && isfinite(F[j])
                m += 1; xs[m] = F[i]; ys[m] = F[j]
            end
        end
        m >= MIN_PAIRS && (ac[k + 1] = cor(view(xs, 1:m), view(ys, 1:m)))
    end
    return ac
end

# All persistence summaries for one window (used for reliability and sensitivity)
function estimators(M::AbstractMatrix)
    R  = vec(M .- profile_mean(M))
    ac = acf_nan(R, MAX_LAG)
    hl_exp, _, plateau = fit_exp_plateau(ac)
    fast(h) = halflife_minutes(acf_nan(R .- centred_mean(R, h), MAX_LAG))
    return (HL = halflife_minutes(ac), HL_fast = fast(MA_HALF), HL_exp = hl_exp,
            AUC4h = STEP_MIN * sum(ac[2:49]), ACF2h = ac[25], Plateau = plateau,
            HL_fast12 = fast(72), HL_fast48 = fast(288))
end

med_iqr(x; d = 1) = begin
    v = filter(isfinite, x); q = quantile(v, [0.25, 0.5, 0.75])
    string(round(q[2]; digits = d), " [", round(q[1]; digits = d), "–", round(q[3]; digits = d), "]")
end
r2s(x) = round(x; digits = 2)

# ======================= Main =======================
function main()
    @load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    elig = cohort[cohort.HasWindow, :]

    function window(sid, startday, ndays)
        g = subjects[sid]
        rows = (g.Day .>= startday) .& (g.Day .<= startday + ndays - 1)
        @assert count(rows) == ndays * PER_DAY
        G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
        return reshape(G, PER_DAY, ndays), g.Measurement[rows]
    end

    # ---------------- 1. Primary window: persistence and standard metrics ----------------
    prim = DataFrame(Subject = Int[], HL = Float64[], HL_fast = Float64[], HL_AR1 = Float64[],
                     Rho1 = Float64[], HL_exp = Float64[], ExpW = Float64[], Plateau = Float64[],
                     SlowShare = Float64[], DayMeanShare = Float64[], TrendSlope = Float64[],
                     TrendRho = Float64[], SigmaR = Float64[], RoutineStrength = Float64[],
                     MeanGlucose = Float64[], SD = Float64[], CV = Float64[], TIR = Float64[],
                     TAR = Float64[], TBR = Float64[], Active = Float64[])
    acf_lag = Dict{Int,Vector{Float64}}()
    Fs = Dict{Int,Tuple{Vector{Float64},Vector{Int}}}()

    for r in eachrow(elig)
        M, meas = window(r.Subject, r.WindowStartDay, WIN)
        mu = profile_mean(M)
        R  = vec(M .- mu)
        ac = acf_nan(R, MAX_LAG)
        slow = centred_mean(R, MA_HALF); F = R .- slow
        Fs[r.Subject] = (F, collect(meas))
        okR = isfinite.(R) .& isfinite.(slow)
        hl_exp, w, c = fit_exp_plateau(ac)
        Rm = reshape(R, PER_DAY, WIN)
        dm = [mean(filter(isfinite, Rm[:, d])) for d in 1:WIN]
        di = collect(1:WIN)
        okd = isfinite.(R)
        obs = filter(isfinite, M); rho1 = ac[2]
        acf_lag[r.Subject] = [ac[k + 1] for (k, _) in ACF_LAGS]
        push!(prim, (r.Subject, halflife_minutes(ac), halflife_minutes(acf_nan(F, MAX_LAG)),
                     0 < rho1 < 1 ? STEP_MIN * log(0.5) / log(rho1) : NaN, rho1,
                     hl_exp, w, c,
                     var(slow[okR]) / var(R[okR]),
                     var(repeat(dm', PER_DAY)[okd]) / var(R[okd]),
                     cov(di, dm) / var(di), corspearman(di, dm),
                     std(filter(isfinite, R)), var(filter(isfinite, mu)) / var(obs),
                     mean(obs), std(obs), 100 * std(obs) / mean(obs),
                     100 * count(x -> 70 <= x <= 180, obs) / length(obs),
                     100 * count(>(180), obs) / length(obs),
                     100 * count(<(70), obs) / length(obs),
                     100 * length(obs) / length(M)))
    end

    println("================ T1D VALIDATION: REPLACE-BG (n = ", nrow(prim), ") ================")
    println("\n---- Cohort characteristics (primary 28-day window), median [IQR] ----")
    for (lab, col, d) in [("Wear time, %", :Active, 1), ("Mean glucose, mg/dL", :MeanGlucose, 0),
                          ("SD, mg/dL", :SD, 0), ("CV, %", :CV, 1), ("TIR, %", :TIR, 1),
                          ("TAR (>180), %", :TAR, 1), ("TBR (<70), %", :TBR, 1)]
        println("  ", rpad(lab, 22), med_iqr(prim[!, col]; d = d))
    end

    println("\n---- Residual autocorrelation (mean across participants) ----")
    L = reduce(hcat, values(acf_lag))
    for (j, (_, lab)) in enumerate(ACF_LAGS)
        println("  ", rpad(lab, 8), round(mean(filter(isfinite, L[j, :])); digits = 2))
    end

    println("\n---- Persistence and slow component, median [IQR] ----")
    for (lab, col, d) in [("Unfiltered half-life, min", :HL, 1), ("Fast half-life, min", :HL_fast, 1),
                          ("AR(1)-implied half-life, min", :HL_AR1, 1), ("Lag-1 autocorrelation", :Rho1, 3),
                          ("Exp-plateau half-life, min", :HL_exp, 1), ("Plateau c", :Plateau, 3),
                          ("Slow share", :SlowShare, 3), ("Day-mean share", :DayMeanShare, 3),
                          ("Trend, mg/dL per day", :TrendSlope, 2)]
        println("  ", rpad(lab, 30), med_iqr(prim[!, col]; d = d))
    end
    println("  Half-life not reached within 24 h: unfiltered ", count(isnan, prim.HL),
            ", fast ", count(isnan, prim.HL_fast))
    println("  Trends: downward (rho < -0.3) ", count(<(-0.3), prim.TrendRho),
            ", upward (rho > 0.3) ", count(>(0.3), prim.TrendRho))
    sp(a, b) = r2s(corspearman(prim[!, a], prim[!, b]))
    println("  Spearman: unfiltered vs slow share ", sp(:HL, :SlowShare),
            "; fast vs slow share ", sp(:HL_fast, :SlowShare),
            "; unfiltered vs fast ", sp(:HL, :HL_fast))

    # ---------------- 2. Relation to standard metrics (non-redundancy) ----------------
    println("\n---- Relation to standard metrics (fast half-life) ----")
    for (lab, col) in [("Mean glucose", :MeanGlucose), ("CV", :CV), ("TIR", :TIR), ("TAR", :TAR),
                       ("Residual SD", :SigmaR), ("Routine strength", :RoutineStrength),
                       ("Slow share", :SlowShare)]
        println("  Spearman with ", rpad(lab, 18), sp(:HL_fast, col))
    end
    rk = DataFrame(HL = tiedrank(prim.HL_fast), Mean = tiedrank(prim.MeanGlucose),
                   CV = tiedrank(prim.CV), TIR = tiedrank(prim.TIR),
                   SDR = tiedrank(prim.SigmaR), RS = tiedrank(prim.RoutineStrength))
    m1 = lm(@formula(HL ~ Mean + CV + TIR), rk)
    m2 = lm(@formula(HL ~ Mean + CV + TIR + SDR + RS), rk)
    println("  Adjusted rank R², standard metrics:                    ", r2s(adjr2(m1)))
    println("  Adjusted rank R², + residual SD and routine strength:  ", r2s(adjr2(m2)))

    # ---------------- 3. Reliability (two non-overlapping windows) ----------------
    rel = DataFrame(Subject = Int[], W = Int[], Estimator = Symbol[], A = Float64[], B = Float64[])
    for r in eachrow(elig)
        d = sort(days[days.Subject .== r.Subject, :], :Day)
        s, Lrun = longest_run(d.Completeness .>= MIN_COMPLETE)
        Lrun < NEED_DAYS && continue
        M56, _ = window(r.Subject, d.Day[s], NEED_DAYS)
        for W in REL_WINDOWS
            eA = estimators(M56[:, 1:W]); eB = estimators(M56[:, B_START:B_START+W-1])
            for e in keys(eA)
                push!(rel, (r.Subject, W, e, getfield(eA, e), getfield(eB, e)))
            end
        end
    end
    nrel = length(unique(rel.Subject))
    logscale = Set([:HL, :HL_fast, :HL_exp, :AUC4h, :HL_fast12, :HL_fast48])

    println("\n---- Reliability (n = ", nrel, " with ≥ ", NEED_DAYS, " consecutive usable days) ----")
    println("  ", rpad("Summary", 12), rpad("W", 5), rpad("ICC(2,1) [95% CI]", 24), rpad("Spearman", 10), "Median A / B")
    for e in [:HL_fast, :HL, :HL_exp, :AUC4h, :ACF2h, :Plateau], W in REL_WINDOWS
        x = rel[(rel.Estimator .== e) .& (rel.W .== W), :]
        ok = isfinite.(x.A) .& isfinite.(x.B) .& (e in logscale ? (x.A .> 0) .& (x.B .> 0) : trues(nrow(x)))
        a, b = x.A[ok], x.B[ok]
        ia, ib = e in logscale ? (log.(a), log.(b)) : (a, b)
        ci = icc21_ci(ia, ib)
        println("  ", rpad(string(e), 12), rpad(W, 5),
                rpad(string(r2s(icc21(ia, ib)), " [", r2s(ci[1]), "–", r2s(ci[2]), "]"), 24),
                rpad(r2s(corspearman(a, b)), 10),
                round(median(a); digits = 2), " / ", round(median(b); digits = 2),
                count(.!ok) > 0 ? "   ($(count(.!ok)) excluded)" : "")
    end

    println("\n---- Filter-width sensitivity of the fast half-life (W = 28) ----")
    for (e, lab) in [(:HL_fast12, "12 h"), (:HL_fast, "24 h"), (:HL_fast48, "48 h")]
        x = rel[(rel.Estimator .== e) .& (rel.W .== 28), :]
        ok = isfinite.(x.A) .& isfinite.(x.B) .& (x.A .> 0) .& (x.B .> 0)
        println("  ", rpad(lab, 6), "ICC ", r2s(icc21(log.(x.A[ok]), log.(x.B[ok]))),
                "   Spearman ", r2s(corspearman(x.A[ok], x.B[ok])),
                "   median A ", round(median(x.A[ok]); digits = 1), " min")
    end
    hlA(e) = sort(rel[(rel.Estimator .== e) .& (rel.W .== 28), :], :Subject).A
    for (e1, e2, lab) in [(:HL_fast12, :HL_fast, "12 h vs 24 h"), (:HL_fast, :HL_fast48, "24 h vs 48 h"),
                          (:HL_fast12, :HL_fast48, "12 h vs 48 h")]
        a, b = hlA(e1), hlA(e2); ok = isfinite.(a) .& isfinite.(b)
        println("  Rank agreement ", lab, ": ", r2s(corspearman(a[ok], b[ok])))
    end

    # ---------------- 4. 24-hour persistence profile ----------------
    prof = DataFrame(Subject = Int[], Block = Int[], HalfLife = Float64[])
    for (sid, (F, meas)) in Fs, b in 1:NBINS
        push!(prof, (sid, b, halflife_minutes(block_acf(F, meas, b))))
    end
    rng = MersenneTwister(2026)
    println("\n---- 24-hour persistence profile (fast half-life by block of onset) ----")
    println("  ", rpad("Start", 8), rpad("n", 5), rpad("Median", 9), rpad("IQR", 15), "95% CI of median")
    meds = Float64[]
    for b in 1:NBINS
        v = filter(isfinite, prof.HalfLife[prof.Block .== b])
        q = quantile(v, [0.25, 0.5, 0.75]); push!(meds, q[2])
        ci = quantile([median(rand(rng, v, length(v))) for _ in 1:N_BOOT], [0.025, 0.975])
        println("  ", rpad(block_label(b), 8), rpad(length(v), 5), rpad(round(q[2]; digits = 1), 9),
                rpad(string(round(q[1]; digits = 1), "–", round(q[3]; digits = 1)), 15),
                round(ci[1]; digits = 1), "–", round(ci[2]; digits = 1))
    end
    println("  Slowest block: ", block_label(argmax(meds)), "  fastest block: ", block_label(argmin(meds)),
            "  ratio ", r2s(maximum(meds) / minimum(meds)))
    night = Float64[]; day = Float64[]
    for sub in groupby(prof, :Subject)
        nv = filter(isfinite, sub.HalfLife[sub.Block .<= 3]); dv = filter(isfinite, sub.HalfLife[sub.Block .> 3])
        if !isempty(nv) && !isempty(dv)
            push!(night, mean(nv)); push!(day, mean(dv))
        end
    end
    println("  Night (00–06) vs day (06–24), within participant: median ",
            round(median(night); digits = 1), " vs ", round(median(day); digits = 1),
            " min, Wilcoxon p = ", round(pvalue(SignedRankTest(night, day)); sigdigits = 3),
            " (n = ", length(night), ")")

    @save joinpath(DIR, "t1d_results.jld2") prim rel prof
    println("\nSaved to: ", joinpath(DIR, "t1d_results.jld2"))
end

main()
