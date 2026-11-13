using Plots

# ------------------------------------------------------------------
# Figure 2: 24-hour persistence profile, all three cohorts
#   Values: median fast half-life (min) by 2-hour block of deviation onset,
#   with 95% bootstrap CI of the median.
#   DiaMonT values from Supplementary Table S6 (persistence_profile.jl);
#   T1D values from t1d_validation.jl. No data files needed.
# ------------------------------------------------------------------
const OUT_FILE = "persistence_profile.png"     # saved next to this script

x = 2 .* (1:12) .- 1                             # block midpoints (hours)

cohorts = [
    ("T2D, basal-only (n = 50)", :steelblue,
     [103.0, 121.2, 93.1, 60.2, 49.1, 58.5, 56.6, 64.7, 68.6, 56.2, 73.9, 84.4],
     [90.4,  94.0,  76.2, 48.1, 44.1, 50.7, 50.9, 56.1, 57.4, 53.0, 68.5, 68.6],
     [117.7, 130.9, 108.1, 76.7, 54.2, 67.0, 66.2, 73.9, 74.2, 66.5, 85.2, 97.7]),
    ("T2D, basal-bolus (n = 80)", :darkorange,
     [122.0, 138.0, 120.9, 71.5, 58.0, 71.3, 67.2, 68.6, 71.0, 73.6, 86.7, 96.3],
     [111.3, 123.3, 107.4, 64.6, 56.4, 64.8, 60.6, 63.1, 64.0, 62.0, 75.7, 84.6],
     [135.7, 153.8, 134.2, 82.7, 62.1, 78.0, 73.8, 74.1, 79.5, 78.0, 96.5, 109.4]),
    ("T1D, insulin pump (n = 216)", :seagreen,
     [129.5, 154.6, 116.7, 78.9, 75.9, 77.4, 73.6, 81.2, 78.7, 78.5, 84.8, 95.7],
     [121.9, 144.2, 110.8, 74.5, 72.2, 74.0, 70.8, 78.2, 74.6, 73.8, 81.1, 91.2],
     [135.9, 163.8, 124.6, 82.3, 79.0, 80.5, 76.2, 83.4, 82.9, 81.7, 88.6, 101.6]),
]

p = plot(xlabel = "Time of deviation onset (h)", ylabel = "Fast half-life (min)",
         xticks = 0:2:24, xlims = (0, 24), ylims = (0, 180), legend = :topright,
         size = (900, 520), dpi = 300, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
         guidefontsize = 10, legendfontsize = 8)
vspan!(p, [0, 6]; color = :gray, alpha = 0.10, label = "Night (00–06)")
for (lab, col, med, lo, hi) in cohorts
    plot!(p, x, med; ribbon = (med .- lo, hi .- med), fillalpha = 0.18,
          lw = 2, color = col, marker = :circle, markersize = 3, label = lab)
end
savefig(p, OUT_FILE)
println("Saved: ", abspath(OUT_FILE))
