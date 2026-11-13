using JLD2, DataFrames, Statistics, StatsBase, LinearAlgebra

# ------------------------------------------------------------------
# Fast half-life vs established variability metrics (MAGE, CONGA1, MAG),
# in the same 28-day windows as the paper. Run once per dataset.
#
#   MAGE:   per day, turning points of the observed readings; excursions
#           smaller than that day's SD are merged iteratively; MAGE is the
#           mean amplitude of the remaining excursions (both directions),
#           averaged over days
#   CONGA1: SD of differences between readings 1 hour (12 readings) apart
#   MAG:    mean absolute change between consecutive readings, per hour
#
# Output per cohort: Spearman correlations with the fast half-life, and
# adjusted rank R^2 of the fast half-life on standard metrics (mean glucose,
# CV, time in range) without and with each variability metric.
# ------------------------------------------------------------------
const DATASET = "replacebg"          # "diamont" or "replacebg"
const DIR = DATASET == "diamont" ? "<path to data>" :
    "<path to data>"

fin(x) = filter(isfinite, x)
prof(M) = [let v = fin(M[t, :]); isempty(v) ? NaN : mean(v) end for t in 1:288]

function slow(r)
    n = length(r); ok = isfinite.(r)
    cs = [0.0; cumsum(ifelse.(ok, r, 0.0))]; cc = [0; cumsum(ok)]
    [begin lo = max(1, i - 144); hi = min(n, i + 144); c = cc[hi+1] - cc[lo]
         c >= (hi - lo + 1) / 2 ? (cs[hi+1] - cs[lo]) / c : NaN end for i in 1:n]
end

function halflife(r)
    x = fin(r); m = mean(x); v = mean((x .- m) .^ 2); prev = 1.0
    for k in 1:288
        p = [(r[i] - m) * (r[i+k] - m) for i in 1:length(r)-k if isfinite(r[i]) && isfinite(r[i+k])]
        ac = mean(p) / v
        ac < 0.5 && return 5 * ((k - 1) + (prev - 0.5) / (prev - ac))
        prev = ac
    end
    NaN
end

function mage_day(x)
    x = fin(x); length(x) < 24 && return NaN
    sd = std(x)
    tp = [x[1]]                                          # turning points
    for i in 2:length(x)-1
        (x[i] - x[i-1]) * (x[i+1] - x[i]) < 0 && push!(tp, x[i])
    end
    push!(tp, x[end])
    while length(tp) > 2                                 # merge excursions smaller than 1 SD
        amp = abs.(diff(tp)); j = argmin(amp)
        amp[j] >= sd && break
        if j == 1;                deleteat!(tp, 2)
        elseif j == length(tp)-1; deleteat!(tp, j)
        else                      deleteat!(tp, j:j+1)
        end
    end
    amp = abs.(diff(tp)); amp = amp[amp .>= sd]
    isempty(amp) ? NaN : mean(amp)
end
mage(M) = mean(fin([mage_day(M[:, d]) for d in 1:size(M, 2)]))

function conga(G, k)
    d = [G[i+k] - G[i] for i in 1:length(G)-k if isfinite(G[i]) && isfinite(G[i+k])]
    std(d)
end
mag(G) = 12 * mean(fin(abs.(diff(G))))

"Adjusted R^2 of rank regression of y on the columns of X."
function adjr2(y, X)
    ok = isfinite.(y) .& vec(all(isfinite.(X); dims = 2))
    yr = tiedrank(y[ok]); Xr = hcat(ones(sum(ok)), [tiedrank(X[ok, j]) for j in 1:size(X, 2)]...)
    b = Xr \ yr; res = yr .- Xr * b
    n, p = size(Xr)
    1 - (sum(res .^ 2) / (n - p)) / var(yr)
end

function main()
    @load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    println("loaded"); flush(stdout)
    c = cohort[cohort.HasWindow, :]
    if DATASET == "diamont"
        c = c[c.N_long .> 0, :]
        grp = ifelse.(c.BasalOnly, "T2D, basal-only", "T2D, basal-bolus")
    else
        grp = fill("T1D, insulin pump", nrow(c))
    end

    res = DataFrame(Group = String[], HL = Float64[], Mean = Float64[], CV = Float64[], TIR = Float64[],
                    MAGE = Float64[], CONGA1 = Float64[], MAG = Float64[])
    for (r, gname) in zip(eachrow(c), grp)
        g = subjects[r.Subject]
        G = Float64.(g.Glucose[(g.Day .>= r.WindowStartDay) .& (g.Day .< r.WindowStartDay + 28)])
        length(G) == 8064 || continue
        G[G .> 401.5] .= NaN
        M = reshape(G, 288, 28)
        R = vec(M .- prof(M)); F = R .- slow(R)
        x = fin(G)
        push!(res, (gname, halflife(F), mean(x), 100 * std(x) / mean(x), 100 * mean(70 .<= x .<= 180),
                    mage(M), conga(G, 12), mag(G)))
    end

    for gname in unique(res.Group)
        d = res[res.Group .== gname, :]
        y = d.HL; std3 = hcat(d.Mean, d.CV, d.TIR)
        sp(v) = (ok = isfinite.(y) .& isfinite.(v); round(corspearman(y[ok], v[ok]); digits = 2))
        println("\n==== ", gname, " (n = ", nrow(d), ") ====")
        qs(v) = (v = fin(v); string(round(median(v); digits = 1), " [", round(quantile(v, 0.25); digits = 1), "–", round(quantile(v, 0.75); digits = 1), "]"))
        println("  MAGE ", qs(d.MAGE), " mg/dL;  CONGA1 ", qs(d.CONGA1), " mg/dL;  MAG ", qs(d.MAG), " mg/dL/h")
        println("  Spearman with fast half-life:  MAGE ", sp(d.MAGE), ",  CONGA1 ", sp(d.CONGA1), ",  MAG ", sp(d.MAG))
        println("  Adjusted rank R^2 of the fast half-life:")
        println("    standard metrics:                 ", round(adjr2(y, std3); digits = 2))
        println("    + MAGE:                           ", round(adjr2(y, hcat(std3, d.MAGE)); digits = 2))
        println("    + CONGA1:                         ", round(adjr2(y, hcat(std3, d.CONGA1)); digits = 2))
        println("    + MAG:                            ", round(adjr2(y, hcat(std3, d.MAG)); digits = 2))
        println("    + MAGE, CONGA1, and MAG:          ", round(adjr2(y, hcat(std3, d.MAGE, d.CONGA1, d.MAG)); digits = 2))
    end
end

main()
