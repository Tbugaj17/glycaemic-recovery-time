using JLD2
using DataFrames
using Statistics
using HypothesisTests      # if missing:  ] add HypothesisTests

# ------------------------------------------------------------------
# Profile window experiment
#   Question: how many training days give the best out-of-sample
#   daily profile, and does a weekly (day-of-week) profile do better?
#
#   Design: for each participant with a long run of usable days,
#   fixed 7-day test blocks are placed after the longest training
#   window, so every window length and profile type is evaluated
#   on exactly the same test data (paired comparison).
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const PER_DAY      = 288
const MIN_COMPLETE = 70.0
const TEST_DAYS    = 7
const WINDOWS      = [7, 14, 21, 28, 35, 42, 49, 56]   # multiples of 7 keep weekdays balanced
const LAMBDAS      = [0.0, 0.25, 0.5, 0.75, 1.0]       # 0 = daily, 1 = weekly, between = shrunken weekly
const REF_W        = 28
const MAX_VALID    = 401.5

@load joinpath(DIR, "cgm_subjects_v2.jld2") subjects days
@load joinpath(DIR, "cohort_v2.jld2") cohort

# ---------------- Helpers ----------------
function profile_mean(M::AbstractMatrix)
    mu = fill(NaN, size(M, 1))
    for t in 1:size(M, 1)
        v = filter(isfinite, @view M[t, :])
        isempty(v) || (mu[t] = mean(v))
    end
    return mu
end

function longest_run(flags)
    best = 0; best_start = 0; cur = 0; cur_start = 0
    for (i, f) in enumerate(flags)
        if f
            cur == 0 && (cur_start = i)
            cur += 1
            if cur > best
                best = cur; best_start = cur_start
            end
        else
            cur = 0
        end
    end
    return best_start, best
end

label(λ) = λ == 0.0 ? "daily" : λ == 1.0 ? "weekly" : "shrunk_$(λ)"

# ---------------- Participants ----------------
# Telemonitoring participants (eligible window), both regimens
base = cohort[(cohort.N_long .> 0) .& cohort.HasWindow, :]
need = maximum(WINDOWS) + TEST_DAYS

results = DataFrame(Subject = Int[], Regimen = String[], Origin = Int[],
                    W = Int[], Profile = String[], SSE = Float64[], N = Int[])

n_used = 0
for r in eachrow(base)
    d = sort(days[days.Subject .== r.Subject, :], :Day)
    s, L = longest_run(d.Completeness .>= MIN_COMPLETE)
    L < need && continue
    global n_used += 1

    run_days = d.Day[s:s+L-1]
    g = subjects[r.Subject]
    rows = (g.Day .>= run_days[1]) .& (g.Day .<= run_days[end])
    G = copy(g.Glucose[rows])
    G[G .> MAX_VALID] .= NaN
    M = reshape(G, PER_DAY, L)                 # time-of-day × day
    dow = mod.(run_days, 7)                    # day-of-cycle (weekday unknown, cycle is enough)

    regimen = r.BasalOnly ? "basal" : "basal-bolus"

    # Test blocks: non-overlapping 7-day blocks after the longest window
    for o in (maximum(WINDOWS) + 1):TEST_DAYS:(L - TEST_DAYS + 1)
        test = o:(o + TEST_DAYS - 1)

        for W in WINDOWS
            train = (o - W):(o - 1)
            daily = profile_mean(M[:, train])
            flat  = mean(filter(isfinite, M[:, train]))
            daily[isnan.(daily)] .= flat

            # Weekly profile: one curve per day-of-cycle
            weekly = Dict{Int,Vector{Float64}}()
            for k in 0:6
                cols = collect(train)[dow[train] .== k]
                wk = profile_mean(M[:, cols])
                wk[isnan.(wk)] .= daily[isnan.(wk)]
                weekly[k] = wk
            end

            # Flat 24-h mean (preprint's reference model)
            sse = 0.0; n = 0
            for j in test, t in 1:PER_DAY
                y = M[t, j]
                isfinite(y) || continue
                sse += (y - flat)^2; n += 1
            end
            push!(results, (r.Subject, regimen, o, W, "flat", sse, n))

            # Daily, shrunken weekly, weekly
            for λ in LAMBDAS
                sse = 0.0; n = 0
                for j in test
                    f = daily .+ λ .* (weekly[dow[j]] .- daily)
                    for t in 1:PER_DAY
                        y = M[t, j]
                        isfinite(y) || continue
                        sse += (y - f[t])^2; n += 1
                    end
                end
                push!(results, (r.Subject, regimen, o, W, label(λ), sse, n))
            end
        end
    end
end

println("Participants with >= ", need, " consecutive usable days: ", n_used,
        "   test blocks: ", length(unique(results[:, [:Subject, :Origin]]).Subject))

# ---------------- Per-participant MSE ----------------
pp = combine(groupby(results, [:Subject, :Regimen, :W, :Profile]),
             [:SSE, :N] => ((s, n) -> sum(s) / sum(n)) => :MSE)

ref = pp[(pp.W .== REF_W) .& (pp.Profile .== "daily"), [:Subject, :MSE]]
rename!(ref, :MSE => :MSE_ref)
pp = innerjoin(pp, ref, on = :Subject)
pp.Ratio = pp.MSE ./ pp.MSE_ref

# ---------------- Summary table ----------------
rows_out = DataFrame(Profile = String[], W = Int[], MeanMSE = Float64[],
                     MedianRMSE = Float64[], MedianRatio = Float64[], p_vs_daily28 = Float64[])
refv = sort(ref, :Subject)
for gdf in groupby(pp, [:Profile, :W])
    x = sort(gdf, :Subject)
    p = (x.Profile[1] == "daily" && x.W[1] == REF_W) ? NaN :
        pvalue(SignedRankTest(x.MSE, refv.MSE_ref))
    push!(rows_out, (x.Profile[1], x.W[1], mean(x.MSE), median(sqrt.(x.MSE)),
                     median(x.Ratio), p))
end
order = Dict("flat" => 0, "daily" => 1, "shrunk_0.25" => 2, "shrunk_0.5" => 3,
             "shrunk_0.75" => 4, "weekly" => 5)
rows_out.ord = [get(order, p, 9) for p in rows_out.Profile]
sort!(rows_out, [:ord, :W])
select!(rows_out, Not(:ord))

println("\nOut-of-sample error by profile and training window")
println("MSE in mg²/dL², RMSE in mg/dL; Ratio = MSE / MSE(daily, 28 days), per participant")
show(rows_out; allrows = true, eltypes = false)
println()

# ---------------- Best window per participant (daily profile) ----------------
dly = pp[pp.Profile .== "daily", :]
best = combine(groupby(dly, :Subject), [:W, :MSE] => ((w, m) -> w[argmin(m)]) => :BestW)
println("\nBest training window (daily profile), number of participants:")
for W in WINDOWS
    println("  ", lpad(W, 2), " days: ", count(==(W), best.BestW))
end

@save joinpath(DIR, "profile_window_experiment.jld2") results pp rows_out best
println("\nSaved to: ", joinpath(DIR, "profile_window_experiment.jld2"))
