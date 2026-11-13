using JLD2, DataFrames, Statistics, StatsBase, GLM   # ] add StatsBase GLM  if missing

# ------------------------------------------------------------------
# Non-redundancy of the FAST half-life (primary measure), per cohort
#   Requires:
#     persistence_v1.jld2         (from persistence_test1.jl)
#     persistence_fast_slow.jld2  (from persistence_fast_slow.jl)
#     filter_sensitivity.jld2     (from filter_sensitivity.jl, for reliability)
#   Cohorts are never pooled.
# ------------------------------------------------------------------
const DIR = "<path to data>"

@load joinpath(DIR, "persistence_v1.jld2") res
@load joinpath(DIR, "persistence_fast_slow.jld2") out
@load joinpath(DIR, "filter_sensitivity.jld2") sens

df = innerjoin(res, out[:, [:Subject, :HalfLifeFast, :SlowShare]], on = :Subject)

const METRICS = [:CV, :MeanGlucose, :TIR, :TAR, :SigmaR, :RoutineStrength, :SlowShare]

for grp in ["basal", "basal-bolus"]
    d = df[df.Regimen .== grp, :]
    n = nrow(d)
    adj(m, p) = 1 - (1 - r2(m)) * (n - 1) / (n - p - 1)

    println("\n================ COHORT: ", grp, " (n=", n, ") ================")
    q = quantile(d.HalfLifeFast, [0.25, 0.5, 0.75])
    println("Fast half-life: median ", round(q[2]; digits = 1), " min [",
            round(q[1]; digits = 1), "–", round(q[3]; digits = 1), "]")

    # 1. Spearman correlations
    println("\nSpearman correlation with fast half-life:")
    for v in METRICS
        println("  ", rpad(string(v), 16), round(corspearman(d.HalfLifeFast, d[!, v]); digits = 2))
    end

    # 2. Variance explained (rank-based, comparable to Spearman)
    rk = DataFrame(HL = tiedrank(d.HalfLifeFast))
    for v in [:CV, :MeanGlucose, :TIR, :SigmaR, :RoutineStrength]
        rk[!, v] = tiedrank(d[!, v])
    end
    m_std = lm(@formula(HL ~ CV + MeanGlucose + TIR), rk)
    m_all = lm(@formula(HL ~ CV + MeanGlucose + TIR + SigmaR + RoutineStrength), rk)

    println("\nVariance of fast half-life explained (rank-based):")
    println("  Standard metrics (CV, mean, TIR):  R² = ", round(r2(m_std); digits = 2),
            "   adjusted = ", round(adj(m_std, 3); digits = 2))
    println("  All five metrics:                  R² = ", round(r2(m_all); digits = 2),
            "   adjusted = ", round(adj(m_all, 5); digits = 2))

    # 3. Rough "new, reproducible information": rank reliability minus rank R²
    s = sens[(sens.Regimen .== grp) .& (sens.MA .== 144) .& (sens.W .== 28), :]
    ok = isfinite.(s.A) .& isfinite.(s.B)
    rel = corspearman(s.A[ok], s.B[ok])
    println("\nRank reliability (Spearman, 28-day windows, n=", count(ok), "): ", round(rel; digits = 2))
    println("Approximate new reproducible share (reliability − adjusted R² standard): ",
            round(rel - adj(m_std, 3); digits = 2))
end
