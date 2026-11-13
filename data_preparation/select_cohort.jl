using JLD2
using DataFrames
using Statistics
using ReadStatTables

# ------------------------------------------------------------------
# Pre-specified cohort rules (state these in the Methods)
# ------------------------------------------------------------------
const DIR = "<path to data>"
const WIN = 28      # development window length (days)
const MIN_COMPLETE = 70.0    # % observed CGM within the window (consensus threshold)
const MIN_USABLE_DAYS = 26   # of 28 days, each >= 70% complete (max 2 bad days);
# rules out the control arm's 20-day blinded periods
const PER_DAY = 288

# Cohort definition (revised):
#   basal-only:  no short-acting insulin records at any time in the trial
#   basal-bolus: short-acting doses recorded WITHIN the analysis window
#   participants with short-acting records only outside the window are
#   excluded from both cohorts (HasWindow = false, NoBolusInWindow = true)

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects id_mapping report days

# ---------------- Back up the previous cohort file (once) ----------------
old = joinpath(DIR, "cohort_v2.jld2")
bak = joinpath(DIR, "cohort_v2_before_window_bolus_rule.jld2")
isfile(old) && !isfile(bak) && cp(old, bak)

# ---------------- Insulin regimen per participant ----------------
ins = DataFrame(readstat(joinpath(DIR, "adapt_ano_insulin.sas7bdat")))
ins.id_Patient = Int.(ins.id_Patient)
ins_summary = combine(groupby(ins, :id_Patient),
    :InsulinType => (x -> count(==("long"), x)) => :N_long,
    :InsulinType => (x -> count(==("short"), x)) => :N_short)
ins_counts = Dict(r.id_Patient => (r.N_long, r.N_short) for r in eachrow(ins_summary))

# Study day of each short-acting dose (time from inclusion; seconds or days)
t = Float64.(coalesce.(ins.MaalingTidFraInkl, NaN))
ins.StudyDay = floor.(Int, (maximum(filter(isfinite, t)) > 1000 ? t ./ 86400 : t)) .+ 1
short_days = Dict(k.id_Patient => g.StudyDay
                  for (k, g) in pairs(groupby(ins[ins.InsulinType .== "short", :], :id_Patient)))
short_in(oid, s) = count(d -> s <= d <= s + WIN - 1, get(short_days, oid, Int[]))

# ---------------- Helpers ----------------
function first_window(d::DataFrame)
    obs = d.N_observed
    usable = d.Completeness .>= MIN_COMPLETE
    n = length(obs)
    n < WIN && return nothing
    for s in 1:(n-WIN+1)
        r = s:(s+WIN-1)
        c = 100 * sum(@view obs[r]) / (WIN * PER_DAY)
        if c >= MIN_COMPLETE && count(@view usable[r]) >= MIN_USABLE_DAYS
            return (s, c)
        end
    end
    return nothing
end

function longest_run(flags)
    best = 0
    cur = 0
    for f in flags
        cur = f ? cur + 1 : 0
        best = max(best, cur)
    end
    return best
end

# ---------------- Per-participant eligibility ----------------
cohort = DataFrame(Subject=Int[], OriginalID=Int[], N_long=Int[], N_short=Int[],
    N_short_window=Int[], NoBolusInWindow=Bool[],
    BasalOnly=Bool[], HasWindow=Bool[],
    WindowStartDay=Int[], WindowStartStudyDay=Int[],
    WindowCompleteness=Float64[], UsableDaysTotal=Int[],
    LongestUsableRun=Int[], PctCeilingWindow=Float64[])

for sid in sort(collect(keys(subjects)))
    oid = id_mapping.OriginalID[id_mapping.Subject .== sid][1]
    nl, ns = get(ins_counts, oid, (0, 0))

    d = sort(days[days.Subject .== sid, :], :Day)
    usable = d.Completeness .>= MIN_COMPLETE
    w = first_window(d)

    if w === nothing
        push!(cohort, (sid, oid, nl, ns, 0, false, ns == 0, false, -1, -1, NaN,
            count(usable), longest_run(usable), NaN))
    else
        s, c = w
        g = subjects[sid]
        in_win = (g.Day .>= d.Day[s]) .& (g.Day .<= d.Day[s] + WIN - 1) .& isfinite.(g.Glucose)
        pct_ceiling = 100 * count(g.Censored[in_win] .== 1) / max(count(in_win), 1)
        nsw = short_in(oid, d.StudyDay[s])
        nobolus = ns > 0 && nsw == 0                 # bolus during the trial, but not in the window
        push!(cohort, (sid, oid, nl, ns, nsw, nobolus, ns == 0, !nobolus, d.Day[s], d.StudyDay[s], c,
            count(usable), longest_run(usable), pct_ceiling))
    end
end

# ---------------- Flow diagram counts ----------------
has_long = cohort.N_long .> 0
basal = has_long .& cohort.BasalOnly
primary = basal .& cohort.HasWindow
bolus_w = has_long .& .!cohort.BasalOnly .& cohort.HasWindow
excluded_nobolus = has_long .& cohort.NoBolusInWindow

println("================ COHORT FLOW ================")
println("Participants with CGM data:                 ", nrow(cohort))
println("  with long-acting insulin records:         ", count(has_long))
println("  basal-only (no short-acting records):     ", count(basal))
println("  + first eligible 28-day window:           ", count(primary), "   <-- BASAL-ONLY COHORT")
println("    (>=70% overall and >=", MIN_USABLE_DAYS, " of 28 days individually >=70%)")
println()
println("Basal-bolus, short-acting doses in window:  ", count(bolus_w), "   <-- BASAL-BOLUS COHORT")
println("Excluded: short-acting only outside window: ", count(excluded_nobolus))
println("No insulin records at all:                  ", count(.!has_long))
println("Check: windows where Day != StudyDay:       ",
        count(cohort.HasWindow .& (cohort.WindowStartDay .!= cohort.WindowStartStudyDay)))

p = cohort[primary, :]
println("\n---------- Basal-only cohort ----------")
println("Window start (study day): median ", median(p.WindowStartStudyDay),
    ", range ", extrema(p.WindowStartStudyDay))
b = cohort[bolus_w, :]
println("---------- Basal-bolus cohort ----------")
println("Window start (study day): median ", median(b.WindowStartStudyDay),
    ", range ", extrema(b.WindowStartStudyDay))
println("Participants with >2% readings at ceiling: basal-only ", count(>(2), p.PctCeilingWindow),
        ", basal-bolus ", count(>(2), b.PctCeilingWindow))

@save joinpath(DIR, "cohort_v2.jld2") cohort
println("\nSaved to: ", joinpath(DIR, "cohort_v2.jld2"), "  (backup: ", bak, ")")
