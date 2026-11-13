using ReadStatTables, JLD2, DataFrames, Statistics

# BMI and diabetes duration at inclusion for the two DiaMonT cohorts, median [IQR].
const DIR = "<path to data>"

q = DataFrame(readstat(joinpath(DIR, "adapt_ano_redcap_spg_ink.sas7bdat")))
cohort = load(joinpath(DIR, "cohort_v2.jld2"), "cohort")
c = cohort[cohort.HasWindow .& (cohort.N_long .> 0), [:OriginalID, :BasalOnly]]

function variable(question, name)
    b = q[(q.Spoergsmaal .== question) .& .!ismissing.(q.SvarKode), [:id_Patient, :SvarKode]]
    b = combine(groupby(b, :id_Patient), :SvarKode => (x -> Float64(first(x))) => name)
    innerjoin(c, b, on = :OriginalID => :id_Patient)
end

for (question, name, unit) in (("BMI", :BMI, "kg/m²"), ("Antal år med diabetes", :Duration, "years"))
    d = variable(question, name)
    for (grp_name, grp) in (("Basal-only", true), ("Basal-bolus", false))
        x = d[d.BasalOnly .== grp, name]
        println(grp_name, ": ", name, " ", round(median(x); digits = 1), " [", round(quantile(x, 0.25); digits = 1),
                "–", round(quantile(x, 0.75); digits = 1), "] ", unit, ", available for ", length(x),
                " of ", sum(c.BasalOnly .== grp))
    end
end
