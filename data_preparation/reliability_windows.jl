using JLD2, DataFrames, ReadStatTables

# ------------------------------------------------------------------
# Which DiaMonT participants qualify for the reliability analysis?
#   Same run selection as the reliability scripts: the longest run of
#   usable days (>= 70% complete), at least 56 days; half A = run days
#   1-28, half B = run days 29-56.
#   Basal-only (no short-acting records in the trial): qualify as before.
#   Basal-bolus: additionally need short-acting doses in BOTH halves.
# Saves the qualifying Subject IDs to reliability_windows.jld2, which the
# five reliability scripts load.
# ------------------------------------------------------------------
const DIR = "<path to data>"
const MIN_COMPLETE, NEED_DAYS = 70.0, 56

@load joinpath(DIR, "cgm_subjects_v2.jld2") days
@load joinpath(DIR, "cohort_v2.jld2") cohort

ins = DataFrame(readstat(joinpath(DIR, "adapt_ano_insulin.sas7bdat")))
t = Float64.(coalesce.(ins.MaalingTidFraInkl, NaN))
ins.StudyDay = floor.(Int, (maximum(filter(isfinite, t)) > 1000 ? t ./ 86400 : t)) .+ 1
sh = ins[ins.InsulinType .== "short", :]
short_days = Dict(Int(k.id_Patient) => g.StudyDay for (k, g) in pairs(groupby(sh, :id_Patient)))
short_in(oid, a, b) = count(d -> a <= d <= b, get(short_days, oid, Int[]))

function longest_run(flags)                      # identical to the reliability scripts
    best = 0; best_start = 0; cur = 0; cur_start = 0
    for (i, f) in enumerate(flags)
        if f
            cur == 0 && (cur_start = i)
            cur += 1
            cur > best && (best = cur; best_start = cur_start)
        else
            cur = 0
        end
    end
    best_start, best
end

base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
out = DataFrame(Subject = Int[], BasalOnly = Bool[], Qualifies = Bool[], ShortA = Int[], ShortB = Int[])
for r in eachrow(base)
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
    L < NEED_DAYS && continue
    sd = d.StudyDay[s:s+NEED_DAYS-1]
    nA = short_in(r.OriginalID, sd[1], sd[28]); nB = short_in(r.OriginalID, sd[29], sd[56])
    ok = r.BasalOnly || (nA > 0 && nB > 0)
    push!(out, (r.Subject, r.BasalOnly, ok, nA, nB))
end

for (name, grp) in (("Basal-only", true), ("Basal-bolus", false))
    x = out[out.BasalOnly .== grp, :]
    println(name, ": ", nrow(x), " with >= 56 consecutive usable days, ",
            count(x.Qualifies), " qualify for the reliability analysis")
end
reliable = Set(out.Subject[out.Qualifies])
@save joinpath(DIR, "reliability_windows.jld2") reliable
println("Saved to: ", joinpath(DIR, "reliability_windows.jld2"))
