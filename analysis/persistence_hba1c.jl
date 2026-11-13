using JLD2, DataFrames, Statistics, StatsBase, GLM, ReadStatTables

# ------------------------------------------------------------------
# Persistence and HbA1c, per cohort (never pooled)
#   First window : the primary 28-day window (from cohort_v2.jld2)
#   Last window  : latest 28-day window ending on or before the day of
#                  the end-of-trial HbA1c, same eligibility rule
#                  (>= 70% overall, >= 26 of 28 days individually >= 70%),
#                  not overlapping the first window.
#   Analyses:
#     A. Baseline:  first-window persistence vs baseline HbA1c
#     B. End:       last-window persistence  vs end-of-trial HbA1c
#     C. Change:    Δ persistence (last − first) vs Δ HbA1c
#   Persistence measures: fast half-life (all day) and overnight
#   fast half-life (deviations starting 00–06).
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const WIN          = 28
const PER_DAY      = 288
const STEP_MIN     = 5
const MA_HALF      = 144
const MAX_LAG      = 288
const MAX_LAG_NT   = 96
const MAX_VALID    = 401.5
const MIN_COMPLETE = 70.0
const MIN_USABLE   = 26
const NIGHT_SLOTS  = 72            # 00:00–06:00

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ---------------- Helpers (same definitions as earlier scripts) ----------------
function profile_mean(M::AbstractMatrix)
    mu = fill(NaN, size(M, 1))
    for t in 1:size(M, 1)
        v = filter(isfinite, @view M[t, :])
        isempty(v) || (mu[t] = mean(v))
    end
    return mu
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

# Correlation between F at selected origins and F k steps later
function origin_acf(F::AbstractVector, origins::Vector{Int}, maxlag::Int)
    n = length(F)
    ac = fill(NaN, maxlag + 1); ac[1] = 1.0
    xs = Vector{Float64}(undef, length(origins)); ys = similar(xs)
    for k in 1:maxlag
        m = 0
        for i in origins
            j = i + k
            if j <= n && isfinite(F[j])
                m += 1; xs[m] = F[i]; ys[m] = F[j]
            end
        end
        m >= 50 && (ac[k + 1] = cor(view(xs, 1:m), view(ys, 1:m)))
    end
    return ac
end

function window_ok(d::DataFrame, s::Int)
    r = s:s+WIN-1
    c = 100 * sum(d.N_observed[r]) / (WIN * PER_DAY)
    return c >= MIN_COMPLETE && count(d.Completeness[r] .>= MIN_COMPLETE) >= MIN_USABLE
end

# Persistence and standard metrics for one 28-day window starting at grid day `startday`
function window_metrics(g::DataFrame, startday::Int)
    rows = (g.Day .>= startday) .& (g.Day .<= startday + WIN - 1)
    G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
    meas = g.Measurement[rows]
    M = reshape(G, PER_DAY, WIN)
    R = vec(M .- profile_mean(M))
    F = R .- centred_mean(R, MA_HALF)
    night = [i for i in eachindex(F) if isfinite(F[i]) && meas[i] <= NIGHT_SLOTS]
    obs = filter(isfinite, G)
    return (HLfast  = halflife_minutes(acf_nan(F, MAX_LAG)),
            HLnight = halflife_minutes(origin_acf(F, night, MAX_LAG_NT)),
            Mean    = mean(obs),
            CV      = 100 * std(obs) / mean(obs),
            TIR     = 100 * count(x -> 70 <= x <= 180, obs) / length(obs))
end

# ---------------- HbA1c ----------------
lab = DataFrame(readstat(joinpath(DIR, "adapt_ano_labka_a.sas7bdat")))
lab = lab[lab.Analysekode .== "HbA1ci", :]
lab.OriginalID = Int.(lab.id_Patient)
lab.StudyDay = fld.(Float64.(lab.MaalingTidFraInkl), 86400) .+ 1

hba = DataFrame(OriginalID = Int[], HbA1c_base = Float64[], HbA1c_end = Float64[], EndDay = Int[])
for sub in groupby(lab, :OriginalID)
    base = sub[sub.StudyDay .<= 14, :]
    fin  = sub[sub.StudyDay .>= 60, :]
    (nrow(base) == 0 || nrow(fin) == 0) && continue
    b = base[argmin(base.StudyDay), :]; e = fin[argmax(fin.StudyDay), :]
    push!(hba, (sub.OriginalID[1], b.InterntResultatNumerisk, e.InterntResultatNumerisk, e.StudyDay))
end

# ---------------- Per participant: first and last window ----------------
eligible = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
res = DataFrame(Subject = Int[], Regimen = String[], HbA1c_base = Float64[], HbA1c_end = Float64[],
                HLfast_1 = Float64[], HLnight_1 = Float64[], Mean_1 = Float64[], CV_1 = Float64[], TIR_1 = Float64[],
                HLfast_2 = Float64[], HLnight_2 = Float64[], Mean_2 = Float64[], CV_2 = Float64[], TIR_2 = Float64[],
                LastStartStudyDay = Int[])

for r in eachrow(eligible)
    h = hba[hba.OriginalID .== r.OriginalID, :]
    nrow(h) == 0 && continue
    g = subjects[r.Subject]
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    m1 = window_metrics(g, r.WindowStartDay)

    # Last window: latest start whose window ends by the end-HbA1c day and after the first window
    first_end = r.WindowStartDay + WIN - 1
    last_s = 0
    for s in (nrow(d) - WIN + 1):-1:1
        d.Day[s] <= first_end && break
        d.StudyDay[s + WIN - 1] <= h.EndDay[1] || continue
        if window_ok(d, s)
            last_s = s; break
        end
    end
    m2 = last_s == 0 ? (HLfast = NaN, HLnight = NaN, Mean = NaN, CV = NaN, TIR = NaN) :
                       window_metrics(g, d.Day[last_s])

    push!(res, (r.Subject, r.BasalOnly ? "basal" : "basal-bolus", h.HbA1c_base[1], h.HbA1c_end[1],
                m1.HLfast, m1.HLnight, m1.Mean, m1.CV, m1.TIR,
                m2.HLfast, m2.HLnight, m2.Mean, m2.CV, m2.TIR,
                last_s == 0 ? -1 : d.StudyDay[last_s]))
end

@save joinpath(DIR, "persistence_hba1c.jld2") res hba

# ---------------- Analysis helpers ----------------
function sp(x, y)
    ok = isfinite.(x) .& isfinite.(y)
    return count(ok) < 10 ? (NaN, count(ok)) : (corspearman(x[ok], y[ok]), count(ok))
end

# Does a persistence measure add to standard metrics in explaining HbA1c? (rank-based)
function incremental(y, hl, mean_g, cv, tir)
    ok = isfinite.(y) .& isfinite.(hl) .& isfinite.(mean_g) .& isfinite.(cv) .& isfinite.(tir)
    rk = DataFrame(Y = tiedrank(y[ok]), HL = tiedrank(hl[ok]), Mean = tiedrank(mean_g[ok]),
                   CV = tiedrank(cv[ok]), TIR = tiedrank(tir[ok]))
    m0 = lm(@formula(Y ~ Mean + CV + TIR), rk)
    m1 = lm(@formula(Y ~ Mean + CV + TIR + HL), rk)
    ct = coeftable(m1)
    return r2(m0), r2(m1), ct.cols[ct.pvalcol][end], count(ok)
end

fmt(x) = isnan(x) ? "  –  " : string(round(x; digits = 2))

# ---------------- Results per cohort ----------------
for grp in ["basal", "basal-bolus"]
    d = res[res.Regimen .== grp, :]
    println("\n================ COHORT: ", grp, " ================")
    println("With baseline and end HbA1c: ", nrow(d),
            "   with eligible last window: ", count(d.LastStartStudyDay .> 0))
    println("HbA1c (mmol/mol): baseline median ", median(d.HbA1c_base),
            ", end median ", median(d.HbA1c_end),
            ", median change ", median(d.HbA1c_end .- d.HbA1c_base))
    ls = d.LastStartStudyDay[d.LastStartStudyDay .> 0]
    isempty(ls) || println("Last window start (study day): median ", median(ls), ", range ", extrema(ls))

    println("\nSpearman correlations with HbA1c (n in brackets):")
    println(rpad("", 26), rpad("A. baseline", 16), rpad("B. end", 16), "C. change vs ΔHbA1c")
    for (lab_, x1, x2) in [("Fast half-life", d.HLfast_1, d.HLfast_2),
                           ("Overnight half-life", d.HLnight_1, d.HLnight_2),
                           ("Mean glucose (ref.)", d.Mean_1, d.Mean_2),
                           ("TIR (ref.)", d.TIR_1, d.TIR_2)]
        a = sp(x1, d.HbA1c_base); b = sp(x2, d.HbA1c_end)
        c = sp(x2 .- x1, d.HbA1c_end .- d.HbA1c_base)
        println(rpad(lab_, 26), rpad(string(fmt(a[1]), " (", a[2], ")"), 16),
                rpad(string(fmt(b[1]), " (", b[2], ")"), 16), fmt(c[1]), " (", c[2], ")")
    end

    println("\nIncremental value beyond standard metrics (rank-based R², HbA1c ~ mean + CV + TIR [+ persistence]):")
    for (lab_, y, hl, mg, cv, tir) in [
            ("A. Baseline, fast",      d.HbA1c_base, d.HLfast_1,  d.Mean_1, d.CV_1, d.TIR_1),
            ("A. Baseline, overnight", d.HbA1c_base, d.HLnight_1, d.Mean_1, d.CV_1, d.TIR_1),
            ("B. End, fast",           d.HbA1c_end,  d.HLfast_2,  d.Mean_2, d.CV_2, d.TIR_2),
            ("B. End, overnight",      d.HbA1c_end,  d.HLnight_2, d.Mean_2, d.CV_2, d.TIR_2)]
        r0, r1, p, n = incremental(y, hl, mg, cv, tir)
        println("  ", rpad(lab_, 24), "R² ", fmt(r0), " → ", fmt(r1),
                "   p(persistence) = ", round(p; digits = 3), "   n = ", n)
    end
end
