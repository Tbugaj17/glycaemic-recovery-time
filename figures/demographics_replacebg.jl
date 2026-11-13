using CSV, JLD2, DataFrames, Statistics, StatsBase

# Sex, age groups, BMI and diabetes duration for the REPLACE-BG cohort (n = 216).
const BASE = "<path to data>"
rd(f) = CSV.read(joinpath(BASE, "Data Tables", f), DataFrame; delim = '|')

scr = rd("HScreening.txt")[:, [:PtID, :Gender, :DiagAge, :Weight, :Height]]
ros = rd("HPtRoster.txt")[:, [:PtID, :AgeAsOfEnrollDt]]
cohort = load(joinpath(BASE, "derived", "cohort_v2.jld2"), "cohort")
ids = DataFrame(PtID = cohort.OriginalID[cohort.HasWindow])

d = innerjoin(innerjoin(ids, scr, on = :PtID), ros, on = :PtID)
num(x) = something(tryparse(Float64, string(x)), NaN)
age = num.(d.AgeAsOfEnrollDt); w = num.(d.Weight); h = num.(d.Height); dx = num.(d.DiagAge)
bmi = w ./ (h ./ 100) .^ 2
dur = age .- dx

qs(x) = (x = filter(isfinite, x); string(round(median(x); digits = 1), " [", round(quantile(x, 0.25); digits = 1),
                                       "–", round(quantile(x, 0.75); digits = 1), "] (n = ", length(x), ")"))
pct(k) = string(k, " (", round(Int, 100k / nrow(d)), ")")

println("REPLACE-BG cohort: n = ", nrow(d))
println("Sex: ", join(["$k " * pct(v) for (k, v) in sort(collect(countmap(string.(d.Gender))))], ",  "))
println("Age < 50: ", pct(sum(age .< 50)), ",  50-59: ", pct(sum(50 .<= age .< 60)),
        ",  60-69: ", pct(sum(60 .<= age .< 70)), ",  >= 70: ", pct(sum(age .>= 70)))
println("Weight: ", qs(w), ";  height: ", qs(h), "   (check units: kg and cm expected)")
println("BMI, kg/m2: ", qs(bmi))
println("Diabetes duration, years: ", qs(dur))
