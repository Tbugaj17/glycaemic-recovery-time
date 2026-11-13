using ReadStatTables, JLD2, DataFrames, Statistics, StatsBase, HypothesisTests

# ------------------------------------------------------------------
# Representativeness of the DiaMonT analysis cohorts (Supplement).
#   A. Included vs excluded, within each regimen, among participants with
#      CGM and insulin-pen records: sex, age >= 60, BMI, diabetes duration,
#      baseline HbA1c.
#      Basal-only:  included (eligible window) vs no eligible window.
#      Basal-bolus: included vs no eligible window vs no short-acting
#                   insulin in the window (the 17).
#   B. Reliability subset vs the rest of each analysis cohort: fast
#      half-life, mean glucose, CV, TIR in the primary window.
# Groups are aggregated (all n >= 5), so nothing individual is printed.
# Run in a fresh Julia session.
# ------------------------------------------------------------------
const DIR = "<path to data>"

cohort = load(joinpath(DIR, "cohort_v2.jld2"), "cohort")
c = cohort[cohort.N_long .> 0, :]                                   # CGM and pen records
c.Group = [r.BasalOnly ? (r.HasWindow ? "Basal-only, included" : "Basal-only, no eligible window") :
           r.HasWindow ? "Basal-bolus, included" :
           r.NoBolusInWindow ? "Basal-bolus, no short-acting in window" : "Basal-bolus, no eligible window"
           for r in eachrow(c)]

# ---- characteristics at inclusion ----
st = DataFrame(readstat(joinpath(DIR, "adapt_ano_redcap_stamdata.sas7bdat")))
st = st[:, [:id_Patient, :Koen, :AlderInklusionInterval]]
st.Women = string.(st.Koen) .== "Kvinde"
st.Age60 = [something(tryparse(Int, first(split(string(a), "-"))), -1) >= 60 for a in st.AlderInklusionInterval]

q = DataFrame(readstat(joinpath(DIR, "adapt_ano_redcap_spg_ink.sas7bdat")))
answer(question, name) = (b = q[(q.Spoergsmaal .== question) .& .!ismissing.(q.SvarKode), [:id_Patient, :SvarKode]];
                          combine(groupby(b, :id_Patient), :SvarKode => (x -> Float64(first(x))) => name))

lab = DataFrame(readstat(joinpath(DIR, "adapt_ano_labka_a.sas7bdat")))
lab = lab[lab.Analysekode .== "HbA1ci", :]
lab.StudyDay = fld.(Float64.(lab.MaalingTidFraInkl), 86400) .+ 1
lab = sort(lab[lab.StudyDay .<= 14, :], :StudyDay)
hba = combine(groupby(lab, :id_Patient), :InterntResultatNumerisk => (x -> Float64(first(x))) => :HbA1c)

d = leftjoin(c[:, [:OriginalID, :Subject, :Group]], st[:, [:id_Patient, :Women, :Age60]], on = :OriginalID => :id_Patient)
d = leftjoin(d, answer("BMI", :BMI), on = :OriginalID => :id_Patient)
d = leftjoin(d, answer("Antal år med diabetes", :Duration), on = :OriginalID => :id_Patient)
d = leftjoin(d, hba, on = :OriginalID => :id_Patient)

fin(x) = Float64.(filter(v -> !ismissing(v) && isfinite(v), x))
mq(x) = (x = fin(x); isempty(x) ? "–" : string(round(median(x); digits = 1), " [", round(quantile(x, 0.25); digits = 1),
                                              "–", round(quantile(x, 0.75); digits = 1), "]"))
pc(x) = (x = collect(skipmissing(x)); string(sum(x), " (", round(Int, 100 * mean(x)), ")"))
pmw(a, b) = (a = fin(a); b = fin(b); (length(a) < 5 || length(b) < 5) ? "–" : string(round(pvalue(MannWhitneyUTest(a, b)); digits = 3)))
pfx(a, b) = (a = collect(skipmissing(a)); b = collect(skipmissing(b));
             string(round(pvalue(FisherExactTest(sum(a), sum(b), length(a) - sum(a), length(b) - sum(b))); digits = 3)))

println("\n================ A. Included vs excluded (CGM and insulin-pen records) ================")
for (ref, others) in (("Basal-only, included", ["Basal-only, no eligible window"]),
                      ("Basal-bolus, included", ["Basal-bolus, no eligible window", "Basal-bolus, no short-acting in window"]))
    a = d[d.Group .== ref, :]
    for o in others
        b = d[d.Group .== o, :]
        println("\n", ref, " (n = ", nrow(a), ")  vs  ", o, " (n = ", nrow(b), ")")
        println("  Women, n (%):            ", pc(a.Women), "  vs  ", pc(b.Women), "   p = ", pfx(a.Women, b.Women))
        println("  Age >= 60, n (%):        ", pc(a.Age60), "  vs  ", pc(b.Age60), "   p = ", pfx(a.Age60, b.Age60))
        for (col, lab_) in ((:BMI, "BMI, kg/m2"), (:Duration, "Diabetes duration, y"), (:HbA1c, "Baseline HbA1c, mmol/mol"))
            println("  ", rpad(lab_ * ":", 25), mq(a[!, col]), "  vs  ", mq(b[!, col]), "   p = ", pmw(a[!, col], b[!, col]),
                    "   (available ", length(fin(a[!, col])), " / ", length(fin(b[!, col])), ")")
        end
    end
end

println("\n================ B. Reliability subset vs rest of the analysis cohort ================")
prim = load(joinpath(DIR, "ac_comparison.jld2"), "prim")
reliable = load(joinpath(DIR, "reliability_windows.jld2"), "reliable")
for grp in ("basal", "basal-bolus")
    p = prim[prim.Regimen .== grp, :]
    inr = in.(p.Subject, Ref(reliable))
    a, b = p[inr, :], p[.!inr, :]
    println("\n", grp, ": reliability subset (n = ", nrow(a), ")  vs  rest (n = ", nrow(b), ")")
    for (col, lab_) in ((:HL_fast, "Fast half-life, min"), (:Mean, "Mean glucose, mg/dL"), (:CV, "CV, %"), (:TIR, "TIR, %"))
        println("  ", rpad(lab_ * ":", 22), mq(a[!, col]), "  vs  ", mq(b[!, col]), "   p = ", pmw(a[!, col], b[!, col]))
    end
end
