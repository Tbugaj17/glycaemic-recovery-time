using JLD2, DataFrames, Statistics, StatsBase, Random, HypothesisTests

# ------------------------------------------------------------------
# 24-hour persistence profile
#   For deviations STARTING in each 2-hour block of the day, how long
#   until they have halved?  Uses the same fast residuals as the
#   primary measure (residual minus 24-h centred moving average).
#
#   For block b and lag k:
#     ρ_b(k) = correlation between F(s) and F(s+k),
#              over all origins s whose clock time falls in block b
#   Half-life of block b = first lag where ρ_b(k) < 0.5
#
#   Each cohort is summarised separately (never pooled).
# ------------------------------------------------------------------
const DIR       = "<path to data>"
const WIN       = 28
const PER_DAY   = 288
const STEP_MIN  = 5
const MA_HALF   = 144          # 24-h centred moving average (±12 h)
const MAX_VALID = 401.5
const BIN_SLOTS = 24           # 2-hour blocks
const NBINS     = PER_DAY ÷ BIN_SLOTS
const MAX_LAG   = 96           # follow deviations up to 8 h ahead
const MIN_PAIRS = 50
const N_BOOT    = 2000

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
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

function halflife_minutes(ac::AbstractVector; level = 0.5)
    for j in 2:length(ac)
        if isfinite(ac[j]) && ac[j] < level
            a0, a1 = ac[j - 1], ac[j]
            return STEP_MIN * ((j - 2) + (a0 - level) / (a0 - a1))
        end
    end
    return NaN
end

block_of(meas) = (meas - 1) ÷ BIN_SLOTS + 1
block_label(b) = string(lpad(2 * (b - 1), 2, '0'), "–", lpad(2 * b, 2, '0'))

# Correlation between F at origins in block b and F k steps later
function block_acf(F::AbstractVector, meas::AbstractVector, b::Int)
    n = length(F)
    origins = [i for i in 1:n if isfinite(F[i]) && block_of(meas[i]) == b]
    ac = fill(NaN, MAX_LAG + 1); ac[1] = 1.0
    xs = Vector{Float64}(undef, length(origins)); ys = similar(xs)
    for k in 1:MAX_LAG
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

# ---------------- Per participant ----------------
eligible = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
prof = DataFrame(Subject = Int[], Regimen = String[], Block = Int[], HalfLife = Float64[])

for r in eachrow(eligible)
    g = subjects[r.Subject]
    rows = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
    G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
    meas = g.Measurement[rows]
    M = reshape(G, PER_DAY, WIN)
    R = vec(M .- profile_mean(M))
    F = R .- centred_mean(R, MA_HALF)                  # fast residuals
    reg = r.BasalOnly ? "basal" : "basal-bolus"
    for b in 1:NBINS
        push!(prof, (r.Subject, reg, b, halflife_minutes(block_acf(F, meas, b))))
    end
end

# ---------------- Cohort summaries ----------------
rng = MersenneTwister(2026)
summary = DataFrame(Regimen = String[], Block = Int[], Label = String[], N = Int[],
                    Median = Float64[], Q1 = Float64[], Q3 = Float64[],
                    CI_lo = Float64[], CI_hi = Float64[])

for grp in ["basal", "basal-bolus"], b in 1:NBINS
    v = filter(isfinite, prof.HalfLife[(prof.Regimen .== grp) .& (prof.Block .== b)])
    q = quantile(v, [0.25, 0.5, 0.75])
    boot = [median(rand(rng, v, length(v))) for _ in 1:N_BOOT]
    ci = quantile(boot, [0.025, 0.975])
    push!(summary, (grp, b, block_label(b), length(v), q[2], q[1], q[3], ci[1], ci[2]))
end

@save joinpath(DIR, "persistence_profile.jld2") prof summary

for grp in ["basal", "basal-bolus"]
    s = summary[summary.Regimen .== grp, :]
    println("\n================ 24-h PERSISTENCE PROFILE: ", grp, " ================")
    println(rpad("Start", 8), rpad("n", 5), rpad("Median HL (min)", 17), rpad("IQR", 16), "95% CI of median")
    for x in eachrow(s)
        println(rpad(x.Label, 8), rpad(x.N, 5), rpad(round(x.Median; digits = 1), 17),
                rpad(string(round(x.Q1; digits = 1), "–", round(x.Q3; digits = 1)), 16),
                round(x.CI_lo; digits = 1), "–", round(x.CI_hi; digits = 1))
    end
    println("Ratio of slowest to fastest block (medians): ",
            round(maximum(s.Median) / minimum(s.Median); digits = 2))

    # Within-participant contrast: night (00–06) vs day (06–24)
    d = prof[prof.Regimen .== grp, :]
    night = Float64[]; day = Float64[]
    for sub in groupby(d, :Subject)
        n_v = filter(isfinite, sub.HalfLife[sub.Block .<= 3])
        d_v = filter(isfinite, sub.HalfLife[sub.Block .> 3])
        if !isempty(n_v) && !isempty(d_v)
            push!(night, mean(n_v)); push!(day, mean(d_v))
        end
    end
    t = SignedRankTest(night, day)
    println("Night (00–06) vs day (06–24), within participant: median ",
            round(median(night); digits = 1), " vs ", round(median(day); digits = 1),
            " min, Wilcoxon p = ", round(pvalue(t); digits = 4), " (n=", length(night), ")")
end
