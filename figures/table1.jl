using JLD2, DataFrames, Statistics, StatsBase, HypothesisTests, ReadStatTables, Printf

# ------------------------------------------------------------------
# Table 1: cohort characteristics following the international
# consensus on CGM metrics (Battelino et al., Diabetes Care 2019)
#   Computed in each participant's primary 28-day window.
#   Cohorts are described separately; the p-value column is a
#   descriptive comparison (Mann–Whitney), not a pooled analysis.
# ------------------------------------------------------------------
const DIR       = "<path to data>"
const WIN       = 28
const PER_DAY   = 288
const MAX_VALID = 401.5

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ---------------- Consensus CGM metrics for one window ----------------
function cgm_metrics(G::AbstractVector)
    obs = filter(isfinite, G)
    n = length(obs); m = mean(obs); sd = std(obs)
    pct(f) = 100 * count(f, obs) / n
    return (Active = 100 * n / length(G),
            Mean   = m,
            SD     = sd,
            GMI    = 3.31 + 0.02392 * m,                 # GMI (%) from mean glucose in mg/dL
            CV     = 100 * sd / m,
            TAR2   = pct(>(250)),                         # level 2 hyperglycaemia
            TAR1   = pct(x -> 180 < x <= 250),            # level 1 hyperglycaemia
            TIR    = pct(x -> 70 <= x <= 180),
            TBR1   = pct(x -> 54 <= x < 70),              # level 1 hypoglycaemia
            TBR2   = pct(<(54)))                          # level 2 hypoglycaemia
end

# ---------------- Baseline HbA1c ----------------
lab = DataFrame(readstat(joinpath(DIR, "adapt_ano_labka_a.sas7bdat")))
lab = lab[lab.Analysekode .== "HbA1ci", :]
lab.OriginalID = Int.(lab.id_Patient)
lab.StudyDay = fld.(Float64.(lab.MaalingTidFraInkl), 86400) .+ 1
base_hba = Dict{Int,Float64}()
for sub in groupby(lab[lab.StudyDay .<= 14, :], :OriginalID)
    base_hba[sub.OriginalID[1]] = sub.InterntResultatNumerisk[argmin(sub.StudyDay)]
end

# ---------------- Per participant ----------------
eligible = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
tab = DataFrame()
for r in eachrow(eligible)
    g = subjects[r.Subject]
    rows = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
    G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
    m = cgm_metrics(G)
    push!(tab, merge((Subject = r.Subject,
                      Regimen = r.BasalOnly ? "basal" : "basal-bolus",
                      HbA1c = get(base_hba, r.OriginalID, NaN)), m); cols = :union)
end

@save joinpath(DIR, "table1.jld2") tab

# ---------------- Formatting helpers ----------------
function med_iqr(x; d = 1)
    v = filter(isfinite, x)
    q = quantile(v, [0.25, 0.5, 0.75])
    fmt(z) = d == 0 ? string(round(Int, z)) : string(round(z; digits = d))
    return string(fmt(q[2]), " [", fmt(q[1]), "–", fmt(q[3]), "]")
end
n_pct(b) = string(count(b), " (", round(Int, 100 * count(b) / length(b)), ")")
pfmt(p) = p < 0.001 ? "<0.001" : string(round(p; digits = 3))

B  = tab[tab.Regimen .== "basal", :]
BB = tab[tab.Regimen .== "basal-bolus", :]

rows_cont = [
    ("CGM wear time, % of possible readings", :Active, 1),
    ("Mean glucose, mg/dL",                   :Mean,   0),
    ("Glucose SD, mg/dL",                     :SD,     0),
    ("GMI, %",                                :GMI,    1),
    ("Coefficient of variation, %",           :CV,     1),
    ("TAR level 2 (>250 mg/dL), %",           :TAR2,   1),
    ("TAR level 1 (181–250 mg/dL), %",        :TAR1,   1),
    ("TIR (70–180 mg/dL), %",                 :TIR,    1),
    ("TBR level 1 (54–69 mg/dL), %",          :TBR1,   1),
    ("TBR level 2 (<54 mg/dL), %",            :TBR2,   1),
    ("HbA1c at baseline, mmol/mol",           :HbA1c,  0)]

rows_target = [
    ("TIR >70%",     d -> d.TIR .> 70),
    ("TBR <4%",      d -> (d.TBR1 .+ d.TBR2) .< 4),
    ("TBR level 2 <1%", d -> d.TBR2 .< 1),
    ("TAR <25%",     d -> (d.TAR1 .+ d.TAR2) .< 25),
    ("TAR level 2 <5%", d -> d.TAR2 .< 5),
    ("CV ≤36%",      d -> d.CV .<= 36)]

# ---------------- Print a Quarto/Markdown table ----------------
println("\nTable: Cohort characteristics in the primary 28-day window. Values are median [interquartile range] or n (%). ",
        "CGM metrics follow the international consensus on time in range. {#tbl-cohort}\n")
println("| Characteristic | Basal-only (n = ", nrow(B), ") | Basal-bolus (n = ", nrow(BB), ") | p |")
println("|---|---|---|---|")
for (lab_, col, d) in rows_cont
    x = filter(isfinite, B[!, col]); y = filter(isfinite, BB[!, col])
    p = pvalue(MannWhitneyUTest(x, y))
    println("| ", lab_, " | ", med_iqr(B[!, col]; d = d), " | ", med_iqr(BB[!, col]; d = d), " | ", pfmt(p), " |")
end
println("| **Consensus targets met** | | | |")
for (lab_, f) in rows_target
    println("| ", lab_, " | ", n_pct(f(B)), " | ", n_pct(f(BB)), " | |")
end
println("\nMissing baseline HbA1c: basal ", count(isnan, B.HbA1c), ", basal-bolus ", count(isnan, BB.HbA1c))
