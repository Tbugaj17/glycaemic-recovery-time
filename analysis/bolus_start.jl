using ReadStatTables, JLD2, DataFrames, Statistics, HypothesisTests

# Did starting bolus insulin change recovery time?  DiaMonT, exploratory.
#   Starters: no short-acting doses in the first window, short-acting doses on
#             at least 14 days of the late window.
#   Controls: basal-only (no short-acting doses in either window) and
#             basal-bolus throughout (short-acting doses in both windows).
#   First window: the analysis window of the paper. Late window: the latest
#   eligible 28-day window (>= 70% readings, >= 26 days >= 70%) after it.
const DIR = "<path to data>"

fin(x) = filter(isfinite, x)
prof(M) = [let v = fin(M[t, :]); isempty(v) ? NaN : mean(v) end for t in 1:288]
function slow(r)
    n = length(r); ok = isfinite.(r)
    cs = [0.0; cumsum(ifelse.(ok, r, 0.0))]; cc = [0; cumsum(ok)]
    [begin lo = max(1, i - 144); hi = min(n, i + 144); c = cc[hi+1] - cc[lo]
         c >= (hi - lo + 1) / 2 ? (cs[hi+1] - cs[lo]) / c : NaN end for i in 1:n]
end
function halflife(r)
    x = fin(r); m = mean(x); v = mean((x .- m) .^ 2); prev = 1.0
    for k in 1:288
        p = [(r[i] - m) * (r[i+k] - m) for i in 1:length(r)-k if isfinite(r[i]) && isfinite(r[i+k])]
        ac = mean(p) / v
        ac < 0.5 && return 5 * ((k - 1) + (prev - 0.5) / (prev - ac))
        prev = ac
    end
    NaN
end
function window(g, start)                           # glucose of a 28-day window, or nothing
    G = Float64.(g.Glucose[(g.Day .>= start) .& (g.Day .< start + 28)])
    length(G) == 8064 || return nothing
    G[G .> 401.5] .= NaN
    G
end
fasthl(G) = (M = reshape(G, 288, 28); R = vec(M .- prof(M)); halflife(R .- slow(R)))

subjects = load(joinpath(DIR, "cgm_subjects_v2.jld2"), "subjects")
cohort = load(joinpath(DIR, "cohort_v2.jld2"), "cohort")
ins = DataFrame(readstat(joinpath(DIR, "adapt_ano_insulin.sas7bdat")))
println("loaded"); flush(stdout)

t = Float64.(coalesce.(ins.MaalingTidFraInkl, NaN))
ins.Day = floor.(Int, (maximum(fin(t)) > 1000 ? t ./ 86400 : t)) .+ 1
short = ins[ins.InsulinType .== "short", [:id_Patient, :Day]]
shortdays(id, s) = length(unique(short.Day[(short.id_Patient .== id) .& (short.Day .>= s) .& (short.Day .< s + 28)]))

c = cohort[cohort.HasWindow .& (cohort.N_long .> 0), :]
res = DataFrame(Group = String[], HL1 = Float64[], HL2 = Float64[], Mean1 = Float64[], Mean2 = Float64[])
for r in eachrow(c)
    g = subjects[r.Subject]
    days = sort(unique(g.Day))
    comp = Dict(d => mean(isfinite.(g.Glucose[g.Day .== d]) .& (g.Glucose[g.Day .== d] .<= 401.5)) for d in days)
    ok(s) = (v = [get(comp, d, 0.0) for d in s:s+27]; mean(v) >= 0.7 && sum(v .>= 0.7) >= 26)
    late = [s for s in (r.WindowStartDay + 28):(maximum(days) - 27) if ok(s)]
    isempty(late) && continue
    s2 = maximum(late)
    n1 = shortdays(r.OriginalID, r.WindowStartDay); n2 = shortdays(r.OriginalID, s2)
    grp = n1 == 0 && n2 >= 14 ? "Started bolus" :
          n1 == 0 && n2 == 0  ? "Basal-only throughout" :
          n1 > 0  && n2 > 0   ? "Basal-bolus throughout" : "Other"
    G1 = window(g, r.WindowStartDay); G2 = window(g, s2)
    (G1 === nothing || G2 === nothing) && continue
    push!(res, (grp, fasthl(G1), fasthl(G2), mean(fin(G1)), mean(fin(G2))))
end

md(x) = round(median(fin(x)); digits = 1)
for grp in ("Started bolus", "Basal-only throughout", "Basal-bolus throughout", "Other")
    d = res[res.Group .== grp, :]
    nrow(d) == 0 && continue
    ch = d.HL2 .- d.HL1; okc = isfinite.(ch)
    p = sum(okc) >= 5 ? round(pvalue(ApproximateSignedRankTest(d.HL2[okc], d.HL1[okc])); digits = 3) : NaN
    println("\n$grp (n = $(nrow(d)))")
    println("  Fast half-life, min: first ", md(d.HL1), "  late ", md(d.HL2), "  median change ", md(ch), "  p = ", p)
    println("  Mean glucose, mg/dL: first ", md(d.Mean1), "  late ", md(d.Mean2))
end

st = res[res.Group .== "Started bolus", :]; ch_st = fin(st.HL2 .- st.HL1)
for ref in ("Basal-only throughout", "Basal-bolus throughout")
    rf = res[res.Group .== ref, :]; ch_rf = fin(rf.HL2 .- rf.HL1)
    (length(ch_st) >= 5 && length(ch_rf) >= 5) || continue
    println("Change in starters vs ", ref, ": p = ",
            round(pvalue(ApproximateMannWhitneyUTest(ch_st, ch_rf)); digits = 3))
end
