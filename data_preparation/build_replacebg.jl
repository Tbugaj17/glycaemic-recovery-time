using CSV, DataFrames, Dates, JLD2, Statistics

# ------------------------------------------------------------------
# REPLACE-BG: build the CGM dataset and select the cohort
#   Same rules as the DiaMonT pipeline (build_dataset_v2.jl + select_cohort.jl):
#     - absolute 5-min grid, round(t / 300), clock days 00:00-23:55
#     - duplicates < 2.5 min: identical kept once, conflicting set missing
#     - collisions moved to next free slot if error <= 5 min, else dropped
#     - LOW (39 mg/dL) and HIGH (400/401 mg/dL) kept and flagged as censored
#     - first 28-day window with >= 70% overall and >= 26 of 28 days >= 70%
#   REPLACE-BG specifics:
#     - display clock (DeviceDtTmDaysFromEnroll + DeviceTm)
#     - study period = randomization day to last recorded visit
#       (removes the blinded run-in and wrongly dated readings)
#     - glucose already in mg/dL (no conversion)
#     - windows overlapping a pump interruption (injection period) are skipped
#   Output is written to a SEPARATE folder, with the same file names and
#   structure as DiaMonT, so the analysis scripts only need DIR changed.
# ------------------------------------------------------------------
const DATA_DIR = "<path to data>"
const OUT_DIR  = "<path to data>"
mkpath(OUT_DIR)

const SLOT         = 300       # seconds per CGM slot
const PER_DAY      = 288
const LOW_CODE     = 39.5      # Dexcom LOW exported as 39 mg/dL
const HIGH_LIMIT   = 400.0     # Dexcom HIGH exported as 400 or 401 mg/dL
const MAX_VALID    = 401.5     # anything above is an artefact
const WIN          = 28
const MIN_COMPLETE = 70.0
const MIN_USABLE   = 26
const NEED_REL     = 56

# ---------------- Grid snapping (identical to build_dataset_v2.jl) ----------------
function snap_to_grid(t::AbstractVector{Float64}, g::AbstractVector{Float64})
    slots = Int[]; vals = Float64[]; times = Float64[]
    dup_same = 0; dup_conflict = 0; shifted = 0; dropped = 0
    for (ti, gi) in zip(t, g)
        s = round(Int, ti / SLOT)
        if !isempty(slots) && s <= slots[end]
            if ti - times[end] < SLOT / 2
                if isapprox(gi, vals[end]; atol = 1e-6)
                    dup_same += 1
                else
                    vals[end] = NaN
                    dup_conflict += 1
                end
                continue
            elseif abs(ti / SLOT - (slots[end] + 1)) <= 1.0
                s = slots[end] + 1
                shifted += 1
            else
                dropped += 1
                continue
            end
        end
        push!(slots, s); push!(vals, gi); push!(times, ti)
    end
    max_err_min = maximum(abs.(times ./ SLOT .- slots)) * SLOT / 60
    return slots, vals, (dup_same = dup_same, dup_conflict = dup_conflict,
                         shifted = shifted, dropped = dropped, max_err_min = max_err_min)
end

# ---------------- Study period per participant ----------------
vis = CSV.read(joinpath(DATA_DIR, "HVisitInfo.txt"), DataFrame; delim = '|')
vis = dropmissing(vis, :VisitDtDaysFromEnroll)
rand_day = Dict(r.PtID => r.VisitDtDaysFromEnroll for r in eachrow(vis) if r.Visit == "Randomization")
last_day = Dict{Int,Int}()
for r in eachrow(vis)
    last_day[r.PtID] = max(get(last_day, r.PtID, typemin(Int)), r.VisitDtDaysFromEnroll)
end

# ---------------- Pump interruptions (injection periods after enrolment) ----------------
ins = CSV.read(joinpath(DATA_DIR, "HInsulin.txt"), DataFrame; delim = '|')
interrupts = Dict{Int,Vector{Tuple{Int,Int}}}()
inj_unknown_start = Int[]; inj_open_end = Int[]
for r in eachrow(ins)
    coalesce(r.InsRoute, "") == "Injection" || continue
    a = r.InsTypeStartDtDaysFromEnroll
    b = r.InsTypeStopDtDaysFromEnroll
    if ismissing(a)
        push!(inj_unknown_start, r.PtID); continue
    end
    if ismissing(b)
        push!(inj_open_end, r.PtID); b = a          # stop unknown: mark the start day only
    end
    push!(get!(interrupts, r.PtID, Tuple{Int,Int}[]), (Int(a), Int(b)))
end

# ---------------- Load CGM (display clock, study period only) ----------------
df = CSV.read(joinpath(DATA_DIR, "HDeviceCGM.txt"), DataFrame; delim = '|',
              select = [:PtID, :DeviceDtTmDaysFromEnroll, :DeviceTm, :RecordType, :GlucoseValue])
df = df[df.RecordType .== "CGM", :]
keep = [haskey(rand_day, p) && haskey(last_day, p) && rand_day[p] <= d <= last_day[p]
        for (p, d) in zip(df.PtID, df.DeviceDtTmDaysFromEnroll)]
df = df[keep, :]
df.t = Float64.(df.DeviceDtTmDaysFromEnroll .* 86400 .+ hour.(df.DeviceTm) .* 3600 .+
                minute.(df.DeviceTm) .* 60 .+ second.(df.DeviceTm))
df.g = Float64.(df.GlucoseValue)
select!(df, [:PtID, :t, :g])
sort!(df, [:PtID, :t])

original_ids = sort(unique(df.PtID))
id_mapping = DataFrame(Subject = 1:length(original_ids), OriginalID = original_ids)
id_map = Dict(o => s for (s, o) in enumerate(original_ids))

# ---------------- Build grid per participant ----------------
subjects = Dict{Int,DataFrame}()
report = DataFrame(Subject = Int[], OriginalID = Int[], N_readings = Int[],
                   Dup_same = Int[], Dup_conflict = Int[], Shifted = Int[], Dropped = Int[],
                   Max_align_err_min = Float64[], Censored_low = Int[], Censored_high = Int[],
                   Above_valid = Int[])
days = DataFrame(Subject = Int[], Day = Int[], StudyDay = Int[],
                 N_observed = Int[], Completeness = Float64[])

for sub in groupby(df, :PtID)
    oid = sub.PtID[1]; sid = id_map[oid]
    slots, vals, st = snap_to_grid(sub.t, sub.g)

    first_slot = fld(slots[1], PER_DAY) * PER_DAY
    last_slot  = (fld(slots[end], PER_DAY) + 1) * PER_DAY - 1
    n = last_slot - first_slot + 1
    abs_slot = first_slot:last_slot

    glu = fill(NaN, n)
    glu[slots .- first_slot .+ 1] = vals
    above = count(x -> isfinite(x) && x > MAX_VALID, glu)
    glu[glu .> MAX_VALID] .= NaN

    censored = zeros(Int8, n)                       # -1 LOW, +1 HIGH
    censored[glu .<= LOW_CODE]   .= -1
    censored[glu .>= HIGH_LIMIT] .=  1

    full = DataFrame(
        Subject     = fill(sid, n),
        OriginalID  = fill(oid, n),
        GridIndex   = 1:n,
        Day         = fld.(0:n-1, PER_DAY) .+ 1,     # day within participant grid
        StudyDay    = fld.(abs_slot, PER_DAY) .+ 1,  # day 1 = enrolment day (day 0 from enrolment)
        Measurement = mod.(abs_slot, PER_DAY) .+ 1,  # 1 = 00:00 clock time
        Glucose     = glu,                           # mg/dL, NaN = missing
        Censored    = censored)
    full.Time = Time(0) .+ Minute.(5 .* (full.Measurement .- 1))
    full.InsulinLong = zeros(n); full.InsulinShort = zeros(n)      # kept for compatibility
    full.InsulinUntyped = zeros(n); full.Insulin = zeros(n)
    subjects[sid] = full

    for d in groupby(full, :Day)
        nobs = count(isfinite, d.Glucose)
        push!(days, (sid, d.Day[1], d.StudyDay[1], nobs, 100 * nobs / PER_DAY))
    end
    push!(report, (sid, oid, nrow(sub), st.dup_same, st.dup_conflict, st.shifted, st.dropped,
                   st.max_err_min, count(==(-1), censored), count(==(1), censored), above))
end
sort!(report, :Subject)

@save joinpath(OUT_DIR, "cgm_subjects_v2.jld2") subjects id_mapping report days

println("================ BUILD ================")
println("Participants with CGM in the study period: ", length(subjects))
println("Readings: ", sum(report.N_readings))
println("Duplicates kept once (identical):       ", sum(report.Dup_same))
println("Duplicates set missing (conflicting):   ", sum(report.Dup_conflict))
println("Readings shifted one slot (jitter):     ", sum(report.Shifted),
        "  (", round(100 * sum(report.Shifted) / sum(report.N_readings); digits = 2), "%)")
println("Readings dropped (collision > 1 slot):  ", sum(report.Dropped))
println("Max alignment error (min):              ", maximum(report.Max_align_err_min))
println("LOW-coded readings:  ", sum(report.Censored_low))
println("HIGH-coded readings: ", sum(report.Censored_high))
println("Readings above ", MAX_VALID, " mg/dL set missing: ", sum(report.Above_valid))

# ---------------- Cohort selection (same rule as select_cohort.jl) ----------------
function window_ok(d::DataFrame, s::Int)
    r = s:s+WIN-1
    c = 100 * sum(d.N_observed[r]) / (WIN * PER_DAY)
    return c >= MIN_COMPLETE && count(d.Completeness[r] .>= MIN_COMPLETE) >= MIN_USABLE, c
end

# Window covers enrolment days (StudyDay - 1) ... (StudyDay - 1 + 27)
function overlaps_interruption(oid::Int, start_study_day::Int)
    haskey(interrupts, oid) || return false
    lo = start_study_day - 1; hi = lo + WIN - 1
    return any(a <= hi && b >= lo for (a, b) in interrupts[oid])
end

function longest_run(flags)
    best = 0; cur = 0
    for f in flags
        cur = f ? cur + 1 : 0
        best = max(best, cur)
    end
    return best
end

cohort = DataFrame(Subject = Int[], OriginalID = Int[], N_long = Int[], N_short = Int[],
                   BasalOnly = Bool[], HasWindow = Bool[],
                   WindowStartDay = Int[], WindowStartStudyDay = Int[],
                   WindowCompleteness = Float64[], UsableDaysTotal = Int[],
                   LongestUsableRun = Int[], PctCeilingWindow = Float64[],
                   RandDay = Int[], StartDaysAfterRand = Int[], SkippedForInterruption = Int[])

for sid in sort(collect(keys(subjects)))
    oid = id_mapping.OriginalID[sid]
    d = sort(days[days.Subject .== sid, :], :Day)
    usable = d.Completeness .>= MIN_COMPLETE
    found = 0; comp = NaN; skipped = 0
    for s in 1:(nrow(d) - WIN + 1)
        ok, c = window_ok(d, s)
        ok || continue
        if overlaps_interruption(oid, d.StudyDay[s])
            skipped += 1; continue
        end
        found = s; comp = c; break
    end
    # All pump users: coded as basal-bolus so the DiaMonT scripts run unchanged
    if found == 0
        push!(cohort, (sid, oid, 1, 1, false, false, -1, -1, NaN,
                       count(usable), longest_run(usable), NaN, rand_day[oid], -1, skipped))
    else
        g = subjects[sid]
        in_win = (g.Day .>= d.Day[found]) .& (g.Day .<= d.Day[found] + WIN - 1) .& isfinite.(g.Glucose)
        pct_ceiling = 100 * count(g.Censored[in_win] .== 1) / max(count(in_win), 1)
        push!(cohort, (sid, oid, 1, 1, false, true, d.Day[found], d.StudyDay[found], comp,
                       count(usable), longest_run(usable), pct_ceiling,
                       rand_day[oid], (d.StudyDay[found] - 1) - rand_day[oid], skipped))
    end
end

@save joinpath(OUT_DIR, "cohort_v2.jld2") cohort

w = cohort[cohort.HasWindow, :]
println("\n================ COHORT ================")
println("Randomized participants with CGM after randomization: ", nrow(cohort))
println("  with an eligible 28-day window:                     ", nrow(w), "   <-- VALIDATION COHORT")
println("  without an eligible window:                         ", count(.!cohort.HasWindow))
println("  with >= ", NEED_REL, " consecutive usable days (reliability):  ",
        count(>=(NEED_REL), w.LongestUsableRun))
println("Window start, days after randomization: median ", median(w.StartDaysAfterRand),
        ", range ", extrema(w.StartDaysAfterRand))
println("Window completeness (%): median ", round(median(w.WindowCompleteness); digits = 1),
        ", range ", round.(extrema(w.WindowCompleteness); digits = 1))
println("Participants with >2% of window readings at the upper limit: ", count(>(2), w.PctCeilingWindow))
println("Windows skipped because of a pump interruption: ", sum(cohort.SkippedForInterruption),
        " (in ", count(>(0), cohort.SkippedForInterruption), " participants)")
isempty(inj_unknown_start) || println("Injection records with unknown start (not excluded), PtID: ", unique(inj_unknown_start))
isempty(inj_open_end)      || println("Injection records with unknown stop (start day only excluded), PtID: ", unique(inj_open_end))
println("\nSaved to: ", OUT_DIR)
