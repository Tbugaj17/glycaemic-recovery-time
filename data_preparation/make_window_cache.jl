using JLD2, DataFrames

# ------------------------------------------------------------------
# Run ONCE per dataset. Loads the full CGM data (the slow part) and saves
# only what the window analyses need: each participant's 28-day window as
# a 288 x 28 matrix (clock time x day), with readings above the reporting
# limit set to NaN, plus the cohort label.
#
# The cache is ~150-220 participants x 8064 values, a few MB, and loads in
# well under a second. Saved LOCALLY (not on the Z: drive).
# ------------------------------------------------------------------
const DATASET = "diamont"          # "diamont" or "replacebg"
const DIR = DATASET == "diamont" ? "<path to data>" :
    "<path to data>"
const OUT = joinpath("<path to data>",
                     "window_cache_$(DATASET).jld2")
const WIN, PER_DAY, MAX_VALID = 28, 288, 401.5

function main()
    t0 = time()
    @load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    println("Full data loaded in ", round(time() - t0; digits = 1), " s"); flush(stdout)

    c = cohort[cohort.HasWindow, :]
    if DATASET == "diamont"
        c = c[c.N_long .> 0, :]
        groups = ifelse.(c.BasalOnly, "T2D, basal-only", "T2D, basal-bolus")
    else
        groups = fill("T1D, insulin pump", nrow(c))
    end

    windows = Dict{Int,Matrix{Float64}}()
    group = Dict{Int,String}()
    for (r, grp) in zip(eachrow(c), groups)
        g = subjects[r.Subject]
        sel = (g.Day .>= r.WindowStartDay) .& (g.Day .<= r.WindowStartDay + WIN - 1)
        G = Float64.(g.Glucose[sel])
        length(G) == WIN * PER_DAY || continue
        G[G .> MAX_VALID] .= NaN
        windows[r.Subject] = reshape(G, PER_DAY, WIN)
        group[r.Subject] = grp
    end
    subjects = nothing; GC.gc()                     # free the memory of the full data

    jldsave(OUT; windows, group)
    println("Saved ", length(windows), " windows to ", OUT, " in ", round(time() - t0; digits = 1), " s total")
end

main()
