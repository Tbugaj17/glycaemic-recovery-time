using DataFrames
using Dates
using ReadStatTables
using JLD2
using Statistics

# ------------------------------------------------------------------
# Constants (state these in the Methods)
# ------------------------------------------------------------------
const MMOL_TO_MGDL    = 18.016   # single conversion point for the whole pipeline
const SLOT            = 300      # seconds per CGM slot
const PER_DAY         = 288      # slots per day
const LOW_CODE_MMOL   = 2.17     # Dexcom "LOW" is exported as 39 mg/dL = 2.165 mmol/L
const HIGH_LIMIT_MMOL = 22.2     # Dexcom upper limit, 400 mg/dL

# ------------------------------------------------------------------
# Snap readings to an ABSOLUTE 5-min grid (no cumulative drift).
# MaalingTidFraInkl is anchored at midnight of the inclusion day
# (confirmed by insulin timing), so slot mod 288 = real clock slot.
# ------------------------------------------------------------------
function snap_to_grid(t::AbstractVector{Float64}, g::AbstractVector{Float64})
    slots = Int[]; vals = Float64[]; times = Float64[]
    dup_same = 0; dup_conflict = 0; shifted = 0; dropped = 0

    for (ti, gi) in zip(t, g)
        s = round(Int, ti / SLOT)
        if !isempty(slots) && s <= slots[end]
            if ti - times[end] < SLOT / 2
                # True duplicate (< 2.5 min apart)
                if isapprox(gi, vals[end]; atol = 1e-6)
                    dup_same += 1                     # identical: keep one copy
                else
                    vals[end] = NaN                   # conflicting: invalid (as in paper)
                    dup_conflict += 1
                end
                continue
            elseif abs(ti / SLOT - (slots[end] + 1)) <= 1.0
                s = slots[end] + 1                    # timing jitter: next free slot
                shifted += 1
            else
                dropped += 1                          # would exceed 1-slot error
                continue
            end
        end
        push!(slots, s); push!(vals, gi); push!(times, ti)
    end

    # Alignment check: assigned slot vs true time, in minutes
    max_err_min = maximum(abs.(times ./ SLOT .- slots)) * SLOT / 60

    return slots, vals,
           (dup_same = dup_same, dup_conflict = dup_conflict,
            shifted = shifted, dropped = dropped, max_err_min = max_err_min)
end

function build_dataset(;
        path_cgm = "<path to data>",
        path_ins = "<path to data>",
        out      = "<path to data>")

    # ---------------- Load CGM ----------------
    df = DataFrame(readstat(path_cgm))
    df.id_Patient = Int.(df.id_Patient)
    df.t = Float64.(df.MaalingTidFraInkl)
    df.g = Float64.(df.GlukoseVaerdi)            # mmol/L
    sort!(df, [:id_Patient, :t])

    original_ids = sort(unique(df.id_Patient))
    @assert length(original_ids) == 300
    id_mapping = DataFrame(Subject = 1:length(original_ids), OriginalID = original_ids)
    id_map = Dict(o => s for (s, o) in enumerate(original_ids))

    # ---------------- Load insulin ----------------
    ins = DataFrame(readstat(path_ins))
    ins.id_Patient = Int.(ins.id_Patient)
    ins.t = Float64.(ins.MaalingTidFraInkl)
    ins_by_id = Dict(k.id_Patient => DataFrame(v)
                     for (k, v) in pairs(groupby(ins, :id_Patient)))

    subjects = Dict{Int,DataFrame}()
    report = DataFrame(Subject = Int[], OriginalID = Int[], N_readings = Int[],
                       Dup_same = Int[], Dup_conflict = Int[], Shifted = Int[],
                       Dropped = Int[], Max_align_err_min = Float64[],
                       Censored_low = Int[], Above_limit = Int[], Max_mmol = Float64[],
                       Insulin_outside_grid = Int[], Insulin_untyped = Int[])
    days = DataFrame(Subject = Int[], Day = Int[], StudyDay = Int[],
                     N_observed = Int[], Completeness = Float64[])

    for sub in groupby(df, :id_Patient)
        oid = sub.id_Patient[1]
        sid = id_map[oid]

        slots, vals, st = snap_to_grid(sub.t, sub.g)

        # Grid from midnight of first day to 23:55 of last day
        first_slot = fld(slots[1], PER_DAY) * PER_DAY
        last_slot  = (fld(slots[end], PER_DAY) + 1) * PER_DAY - 1
        n = last_slot - first_slot + 1
        abs_slot = first_slot:last_slot

        mmol = fill(NaN, n)
        mmol[slots .- first_slot .+ 1] = vals

        censored = zeros(Int8, n)                 # -1 LOW, +1 at/above upper limit
        censored[mmol .<= LOW_CODE_MMOL]   .= -1
        censored[mmol .>= HIGH_LIMIT_MMOL] .=  1

        full = DataFrame(
            Subject     = fill(sid, n),
            OriginalID  = fill(oid, n),
            GridIndex   = 1:n,
            Day         = fld.(0:n-1, PER_DAY) .+ 1,         # day within participant grid
            StudyDay    = fld.(abs_slot, PER_DAY) .+ 1,      # day 1 = inclusion day
            Measurement = mod.(abs_slot, PER_DAY) .+ 1,      # 1 = 00:00 clock time
            Glucose     = mmol .* MMOL_TO_MGDL,              # mg/dL, NaN = missing
            Censored    = censored,
        )
        full.Time = Time(0) .+ Minute.(5 .* (full.Measurement .- 1))

        # ---------------- Insulin on the same absolute grid ----------------
        full.InsulinLong    = zeros(n)
        full.InsulinShort   = zeros(n)
        full.InsulinUntyped = zeros(n)
        outside = 0; untyped = 0
        if haskey(ins_by_id, oid)
            for row in eachrow(ins_by_id[oid])
                k = round(Int, row.t / SLOT) - first_slot + 1
                if !(1 <= k <= n)
                    outside += 1
                    continue
                end
                if row.InsulinType == "long"
                    full.InsulinLong[k] += row.InsulinVaerdi
                elseif row.InsulinType == "short"
                    full.InsulinShort[k] += row.InsulinVaerdi
                else
                    full.InsulinUntyped[k] += row.InsulinVaerdi
                    untyped += 1
                end
            end
        end
        full.Insulin = full.InsulinLong .+ full.InsulinShort .+ full.InsulinUntyped

        subjects[sid] = full

        # ---------------- Day-level completeness (clock days) ----------------
        for d in groupby(full, :Day)
            nobs = count(isfinite, d.Glucose)
            push!(days, (sid, d.Day[1], d.StudyDay[1], nobs, 100 * nobs / PER_DAY))
        end

        push!(report, (sid, oid, nrow(sub), st.dup_same, st.dup_conflict, st.shifted,
                       st.dropped, st.max_err_min,
                       count(==(-1), censored), count(==(1), censored),
                       maximum(filter(isfinite, mmol)), outside, untyped))
    end

    sort!(report, :Subject)

    # ---------------- Build summary ----------------
    println("Subjects: ", length(subjects))
    println("Duplicates kept once (identical): ", sum(report.Dup_same))
    println("Duplicates set missing (conflicting): ", sum(report.Dup_conflict))
    println("Readings shifted one slot (jitter): ", sum(report.Shifted))
    println("Readings dropped (collision > 1 slot): ", sum(report.Dropped))
    println("Max alignment error (min): ", maximum(report.Max_align_err_min))
    println("LOW-coded readings: ", sum(report.Censored_low))
    println("Readings >= 22.2 mmol/L: ", sum(report.Above_limit),
            " in ", count(>(0), report.Above_limit), " participants")
    println("Insulin records outside CGM grid: ", sum(report.Insulin_outside_grid))
    println("Insulin records without type: ", sum(report.Insulin_untyped))

    @save out subjects id_mapping report days
    println("Saved to: ", out)

    return subjects, id_mapping, report, days
end

subjects, id_mapping, report, days = build_dataset()
