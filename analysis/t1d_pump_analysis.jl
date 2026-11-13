using CSV, DataFrames, Dates, JLD2, Statistics, StatsBase

# ------------------------------------------------------------------
# REPLACE-BG pump data vs persistence (in each participant's 28-day window)
#   From HDeviceBasal + HDeviceBolus:
#     - total daily insulin (basal + bolus), basal fraction, per kg if weight found
#     - dawn ratio: programmed basal rate 04-08 / 00-03
#     - night interventions (00-06): correction boluses and temp basal/suspend
#       starts per night
#   Related (Spearman) to the fast half-life, the overnight half-life
#   (mean of blocks 00-06), the night/day ratio, and time below range.
#   Inputs: Data Tables (pump, screening) and derived (cohort, t1d_results).
# ------------------------------------------------------------------
const DATA_DIR = "<path to data>"
const DIR      = "<path to data>"
const WIN          = 28
const PER_DAY      = 288
const MIN_COVERAGE = 0.70        # minimum share of the window with basal records

tosec(x::Time) = hour(x) * 3600 + minute(x) * 60 + second(x)
tosec(::Missing) = missing
tofloat(x::Missing) = NaN
tofloat(x::Number) = Float64(x)
tofloat(x) = something(tryparse(Float64, strip(string(x))), NaN)

function main()
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    t1 = load(joinpath(DIR, "t1d_results.jld2"))
    prim, prof = t1["prim"], t1["prof"]
    w = cohort[cohort.HasWindow, :]

    hl_night = Dict(k.Subject => mean(filter(isfinite, g.HalfLife))
                    for (k, g) in pairs(groupby(prof[prof.Block .<= 3, :], :Subject)))
    hl_day   = Dict(k.Subject => mean(filter(isfinite, g.HalfLife))
                    for (k, g) in pairs(groupby(prof[prof.Block .> 3, :], :Subject)))

    # ---------------- Pump files (exact duplicates from repeated uploads removed) ----------------
    basal = CSV.read(joinpath(DATA_DIR, "HDeviceBasal.txt"), DataFrame; delim = '|')
    bolus = CSV.read(joinpath(DATA_DIR, "HDeviceBolus.txt"), DataFrame; delim = '|')
    nb0, nbo0 = nrow(basal), nrow(bolus)
    unique!(basal, Not([:RecID, :ParentHDeviceUploadsID]))
    unique!(bolus, Not([:RecID, :ParentHDeviceUploadsID]))
    println("Basal rows: ", nb0, " -> ", nrow(basal), " after removing duplicates")
    println("Bolus rows: ", nbo0, " -> ", nrow(bolus), " after removing duplicates")
    println("Basal types: ", countmap(basal.BasalType))
    println("Bolus types: ", countmap(bolus.BolusType))
    basal_by = Dict(k.PtID => g for (k, g) in pairs(groupby(basal, :PtID)))
    bolus_by = Dict(k.PtID => g for (k, g) in pairs(groupby(bolus, :PtID)))

    # ---------------- Weight (if available) ----------------
    weight = Dict{Int,Float64}()
    scr = CSV.read(joinpath(DATA_DIR, "HScreening.txt"), DataFrame; delim = '|')
    wcols = filter(c -> occursin("weight", lowercase(c)) || occursin("wt", lowercase(c)), names(scr))
    println("\nHScreening columns that may hold weight: ", wcols)
    if !isempty(wcols)
        wc = wcols[1]
        v = tofloat.(scr[!, wc])
        lbs = median(filter(isfinite, v)) > 130          # US data: likely pounds
        for (p, x) in zip(scr.PtID, v)
            isfinite(x) && (weight[p] = lbs ? x * 0.45359 : x)
        end
        println("Using column '", wc, "'", lbs ? " (interpreted as pounds, converted to kg)" : " (kg)",
                "; median ", round(median(collect(values(weight))); digits = 1), " kg, n = ", length(weight))
    else
        println("No weight column found: insulin per kg will be missing. Columns: ", names(scr))
    end

    # ---------------- Per participant, in the analysis window ----------------
    res = DataFrame(Subject = Int[], PtID = Int[], Coverage = Float64[], TDD = Float64[],
                    BasalPerDay = Float64[], BolusPerDay = Float64[], BasalFraction = Float64[],
                    TDDperKg = Float64[], DawnRatio = Float64[],
                    NightBolusPerNight = Float64[], NightTempPerNight = Float64[])
    n = WIN * PER_DAY
    for r in eachrow(w)
        oid = r.OriginalID
        d0 = r.WindowStartStudyDay - 1                  # first window day, days from enrolment
        d1 = d0 + WIN - 1
        haskey(basal_by, oid) || continue
        b = basal_by[oid]
        b = b[(b.DeviceDtTmDaysFromEnroll .>= d0 - 1) .& (b.DeviceDtTmDaysFromEnroll .<= d1), :]
        b = dropmissing(b, [:DeviceTm, :Duration])
        b.start = (b.DeviceDtTmDaysFromEnroll .- d0) .* 86400 .+ tosec.(b.DeviceTm)
        sort!(b, :start)

        actual = fill(NaN, n); programmed = fill(NaN, n)
        for s in eachrow(b)
            dur = tofloat(s.Duration) / 1000
            isfinite(dur) && dur > 0 || continue
            act  = s.BasalType == "suspend" ? 0.0 : tofloat(s.Rate)
            prog = s.BasalType == "scheduled" ? tofloat(s.Rate) : tofloat(s.SuprRate)
            i0 = max(1, floor(Int, s.start / 300) + 1)
            i1 = min(n, floor(Int, (s.start + dur - 1) / 300) + 1)
            i1 < i0 && continue
            isfinite(act)  && (actual[i0:i1] .= act)       # later segments overwrite earlier ones
            isfinite(prog) && (programmed[i0:i1] .= prog)
        end
        cov = count(isfinite, actual) / n
        basal_day = cov > 0 ? mean(filter(isfinite, actual)) * 24 : NaN

        bo = haskey(bolus_by, oid) ? bolus_by[oid] : bolus[1:0, :]
        bo = dropmissing(bo[(bo.DeviceDtTmDaysFromEnroll .>= d0) .& (bo.DeviceDtTmDaysFromEnroll .<= d1), :], :DeviceTm)
        units = tofloat.(bo.Normal) .+ map(x -> isfinite(x) ? x : 0.0, tofloat.(bo.Extended))
        units = map(x -> isfinite(x) ? x : 0.0, units)
        covered_days = cov * WIN
        bolus_day = covered_days > 0 ? sum(units) / covered_days : NaN

        P = reshape(programmed, PER_DAY, WIN)
        early = filter(isfinite, vec(P[1:36, :]))         # 00:00-02:59
        dawn  = filter(isfinite, vec(P[49:96, :]))        # 04:00-07:59
        dawn_ratio = (isempty(early) || isempty(dawn)) ? NaN : median(dawn) / median(early)

        night_bol = count((tosec.(bo.DeviceTm) .< 6 * 3600) .& (units .> 0))
        bw = b[(b.DeviceDtTmDaysFromEnroll .>= d0) .& in.(b.BasalType, Ref(["temp", "suspend"])), :]
        night_tmp = count(tosec.(bw.DeviceTm) .< 6 * 3600)

        ok = cov >= MIN_COVERAGE
        tdd = ok ? basal_day + bolus_day : NaN
        push!(res, (r.Subject, oid, cov, tdd, ok ? basal_day : NaN, ok ? bolus_day : NaN,
                    ok ? basal_day / tdd : NaN,
                    ok && haskey(weight, oid) ? tdd / weight[oid] : NaN,
                    ok ? dawn_ratio : NaN,
                    ok ? night_bol / covered_days : NaN, ok ? night_tmp / covered_days : NaN))
    end

    # ---------------- Join with persistence ----------------
    res.HL_fast   = [get(Dict(zip(prim.Subject, prim.HL_fast)), s, NaN) for s in res.Subject]
    res.TBR       = [get(Dict(zip(prim.Subject, prim.TBR)), s, NaN) for s in res.Subject]
    res.HL_night  = [get(hl_night, s, NaN) for s in res.Subject]
    res.NightDayRatio = [get(hl_night, s, NaN) / get(hl_day, s, NaN) for s in res.Subject]

    println("\n================ PUMP DATA IN THE ANALYSIS WINDOWS ================")
    println("Participants with pump records: ", nrow(res), ";  with basal coverage ≥ ",
            round(Int, 100 * MIN_COVERAGE), "%: ", count(>=(MIN_COVERAGE), res.Coverage))
    mi(x) = (v = filter(isfinite, x); isempty(v) ? "–" :
             (q = quantile(v, [0.25, 0.5, 0.75]);
              string(round(q[2]; sigdigits = 3), " [", round(q[1]; sigdigits = 3), "–", round(q[3]; sigdigits = 3), "]  n=", length(v))))
    for (lab, col) in [("Total daily insulin, U", :TDD), ("Basal, U/day", :BasalPerDay),
                       ("Bolus, U/day", :BolusPerDay), ("Basal fraction", :BasalFraction),
                       ("Total daily insulin, U/kg", :TDDperKg), ("Dawn ratio (04-08 / 00-03)", :DawnRatio),
                       ("Night boluses per night", :NightBolusPerNight),
                       ("Night temp/suspend starts per night", :NightTempPerNight)]
        println("  ", rpad(lab, 38), mi(res[!, col]))
    end

    println("\n---- Spearman correlations ----")
    targets = [(:HL_fast, "Fast HL"), (:HL_night, "Night HL"), (:NightDayRatio, "Night/day"), (:TBR, "TBR")]
    print("  ", rpad("", 38)); foreach(t -> print(rpad(t[2], 12)), targets); println()
    for (lab, col) in [("Total daily insulin, U", :TDD), ("Total daily insulin, U/kg", :TDDperKg),
                       ("Basal fraction", :BasalFraction), ("Dawn ratio", :DawnRatio),
                       ("Night boluses per night", :NightBolusPerNight),
                       ("Night temp/suspend per night", :NightTempPerNight), ("TBR", :TBR)]
        print("  ", rpad(lab, 38))
        for (t, _) in targets
            ok = isfinite.(res[!, col]) .& isfinite.(res[!, t])
            print(rpad(count(ok) > 10 && col != t ? string(round(corspearman(res[ok, col], res[ok, t]); digits = 2)) : "–", 12))
        end
        println()
    end

    @save joinpath(DIR, "pump_metrics.jld2") res
    println("\nSaved to: ", joinpath(DIR, "pump_metrics.jld2"))
end

main()
