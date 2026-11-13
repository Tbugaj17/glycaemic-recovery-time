using JLD2, DataFrames, Statistics, Plots

# ------------------------------------------------------------------
# Figure 3: test-retest reliability of the fast and unfiltered half-life
#   Each point is one participant: half-life in the first 28 days (A)
#   against the second 28 days (B), log scale, with the identity line.
#   Top row: fast half-life; bottom row: unfiltered half-life.
#   Columns: T2D basal-only, T2D basal-bolus, T1D.
#   Uses the results saved by icc_difference.jl (no data access needed).
# ------------------------------------------------------------------
const OUT = "<path to data>"

dia = load(joinpath(OUT, "icc_difference_diamont.jld2"), "rows")
rbg = load(joinpath(OUT, "icc_difference_replacebg.jld2"), "rows")
rbg.Regimen .= "T1D"

function icc21(x, y)
    Y = hcat(x, y); n, k = size(Y); gm = mean(Y)
    ssr = k * sum(abs2, mean(Y; dims = 2) .- gm); ssc = n * sum(abs2, mean(Y; dims = 1) .- gm)
    sse = sum(abs2, Y .- gm) - ssr - ssc
    msr = ssr / (n - 1); msc = ssc / (k - 1); mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
end

cohorts = [("T2D, basal-only", dia[dia.Regimen .== "basal", :]),
           ("T2D, basal-bolus", dia[dia.Regimen .== "basal-bolus", :]),
           ("T1D, insulin pump", rbg)]
lims = (25, 400)
ticks = ([30, 60, 120, 240], ["30", "60", "120", "240"])

panels = Plots.Plot[]
for (row, (a, b, lab)) in enumerate([(:fA, :fB, "Fast half-life"), (:uA, :uB, "Unfiltered half-life")])
    for (title, d) in cohorts
        ok = isfinite.(d[!, a]) .& isfinite.(d[!, b]) .& (d[!, a] .> 0) .& (d[!, b] .> 0)
        x, y = d[ok, a], d[ok, b]
        icc = icc21(log.(x), log.(y))
        p = scatter(x, y; xscale = :log10, yscale = :log10, xlims = lims, ylims = lims,
                    xticks = ticks, yticks = ticks, markersize = 3, markerstrokewidth = 0,
                    color = row == 1 ? :steelblue : :gray45, alpha = 0.7, legend = false,
                    aspect_ratio = :equal, framestyle = :box,
                    title = row == 1 ? title : "", titlefontsize = 9,
                    xlabel = row == 2 ? "Days 1–28 (min)" : "",
                    ylabel = title == "T2D, basal-only" ? "$lab\nDays 29–56 (min)" : "")
        plot!(p, [lims[1], lims[2]], [lims[1], lims[2]]; color = :black, ls = :dash, lw = 0.8)
        annotate!(p, 32, 300, text("ICC = $(round(icc; digits = 2))\nn = $(length(x))", 8, :left))
        push!(panels, p)
    end
end

fig = plot(panels...; layout = (2, 3), size = (900, 620), dpi = 300,
           left_margin = 6Plots.mm, bottom_margin = 4Plots.mm, guidefontsize = 8, tickfontsize = 7)
savefig(fig, joinpath(OUT, "figure_reliability.png"))
println("Saved: ", joinpath(OUT, "figure_reliability.png"))
