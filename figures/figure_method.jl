using JLD2, DataFrames, Statistics, Plots

# ------------------------------------------------------------------
# Figure 1: how glycaemic recovery time is calculated, on one participant
#   Uses REPLACE-BG (public, de-identified data), so no individual
#   DiaMonT data are shown.
#   Participant chosen automatically: fast half-life close to the cohort
#   median AND slow share in the upper quarter of the cohort, so that the
#   day-level shifts are visible; the most complete window among those.
#   a. Three days of glucose with the average daily profile
#   b. Fourteen days of deviations with the slow component
#   c. ACF of the deviations (unfiltered) and of the fast component
#   Panels carry only a lowercase letter (no titles), and all text is
#   sans serif, 8-12 pt, as required by Communications Medicine.
#   Set PICK_SUBJECT to override.
# ------------------------------------------------------------------
const DIR          = "<path to data>"
const OUT_FILE     = "figure_method.png"
const PICK_SUBJECT = nothing
const WIN          = 28
const PER_DAY      = 288
const STEP_MIN     = 5
const MAX_LAG      = 288
const MA_HALF      = 144
const MAX_VALID    = 401.5
const SHOW_DAYS_A  = 3
const SHOW_DAYS_B  = 14
const SHOW_LAG_H   = 24

# ---------------- Helpers (identical definitions) ----------------
function profile_mean(M::AbstractMatrix)
    mu = fill(NaN, size(M, 1))
    for t in 1:size(M, 1)
        v = filter(isfinite, @view M[t, :])
        isempty(v) || (mu[t] = mean(v))
    end
    return mu
end

function acf_nan(r::AbstractVector, maxlag::Int)
    ok = isfinite.(r)
    m = mean(r[ok]); v = mean((r[ok] .- m) .^ 2)
    out = fill(NaN, maxlag + 1); out[1] = 1.0
    n = length(r)
    for k in 1:maxlag
        s = 0.0; c = 0
        @inbounds for i in 1:(n - k)
            x1 = r[i]; x2 = r[i + k]
            if isfinite(x1) && isfinite(x2)
                s += (x1 - m) * (x2 - m); c += 1
            end
        end
        c > 30 && (out[k + 1] = s / c / v)
    end
    return out
end

function halflife_minutes(ac::AbstractVector; level = 0.5)
    for j in 2:length(ac)
        if isfinite(ac[j]) && ac[j] < level
            a0, a1 = ac[j - 1], ac[j]
            return STEP_MIN * ((j - 2) + (a0 - level) / (a0 - a1))
        end
    end
    return NaN
end

function centred_mean(r::AbstractVector, half::Int)
    n = length(r); ok = isfinite.(r)
    cs = vcat(0.0, cumsum(ifelse.(ok, r, 0.0))); cc = vcat(0, cumsum(ok))
    out = fill(NaN, n)
    for i in 1:n
        lo = max(1, i - half); hi = min(n, i + half)
        c = cc[hi + 1] - cc[lo]
        c >= 0.5 * (hi - lo + 1) && (out[i] = (cs[hi + 1] - cs[lo]) / c)
    end
    return out
end

function window_glucose(subjects, sid, startday)
    g = subjects[sid]
    rows = (g.Day .>= startday) .& (g.Day .<= startday + WIN - 1)
    G = copy(g.Glucose[rows]); G[G .> MAX_VALID] .= NaN
    return G
end

function main()
    @load joinpath(DIR, "cgm_subjects_v2.jld2") subjects
    @load joinpath(DIR, "cohort_v2.jld2") cohort
    prim = load(joinpath(DIR, "t1d_results.jld2"), "prim")

    # ---------------- Choose a representative participant with visible day-level shifts ----------------
    d = innerjoin(select(prim, :Subject, :HL_fast, :SlowShare),
                  select(cohort[cohort.HasWindow, :], :Subject, :WindowStartDay, :WindowCompleteness),
                  on = :Subject)
    sid = if PICK_SUBJECT === nothing
        mh = median(d.HL_fast); q75 = quantile(d.SlowShare, 0.75)
        near = d[(abs.(d.HL_fast .- mh) .<= 10) .& (d.SlowShare .>= q75), :]
        isempty(near) && (near = d[d.SlowShare .>= q75, :])
        near[argmax(near.WindowCompleteness), :Subject]
    else
        PICK_SUBJECT
    end
    r = d[d.Subject .== sid, :][1, :]
    println("Participant: Subject ", sid, "  fast half-life ", round(r.HL_fast; digits = 1),
            " min (cohort median ", round(median(d.HL_fast); digits = 1), "), slow share ",
            round(r.SlowShare; digits = 2), " (cohort 75th percentile ", round(quantile(d.SlowShare, 0.75); digits = 2), ")")

    # ---------------- Compute everything for the window ----------------
    G  = window_glucose(subjects, sid, r.WindowStartDay)
    M  = reshape(G, PER_DAY, WIN)
    mu = profile_mean(M)
    R  = vec(M .- mu)
    S  = centred_mean(R, MA_HALF)
    F  = R .- S
    acR = acf_nan(R, MAX_LAG); acF = acf_nan(F, MAX_LAG)
    hlR = halflife_minutes(acR); hlF = halflife_minutes(acF)
    println("Unfiltered half-life ", round(hlR; digits = 1), " min, fast half-life ", round(hlF; digits = 1), " min")

    # ---------------- Choose the most complete stretch to show ----------------
    dayfrac = [count(isfinite, M[:, k]) / PER_DAY for k in 1:WIN]
    best_start(len) = argmax([mean(dayfrac[k:k+len-1]) for k in 1:WIN-len+1])
    dA = best_start(SHOW_DAYS_A)
    dB = best_start(SHOW_DAYS_B)
    println("Panel A shows days ", dA, "–", dA + SHOW_DAYS_A - 1,
            ", panel B days ", dB, "–", dB + SHOW_DAYS_B - 1)

    # ---------------- Panel A: glucose and average daily profile ----------------
    iA = ((dA - 1) * PER_DAY + 1):((dA - 1 + SHOW_DAYS_A) * PER_DAY)
    tA = (0:length(iA)-1) ./ PER_DAY
    gA = G[iA]
    ytop = ceil(maximum(filter(isfinite, vcat(gA, mu))) / 50) * 50 + 30
    pA = plot(tA, gA; color = :gray40, lw = 1.2, label = "CGM glucose",
              ylabel = "Glucose (mg/dL)", xlabel = "Day", xlims = (0, SHOW_DAYS_A),
              ylims = (40, ytop), legend = :outertopright,
              title = "a", titleloc = :left, titlefontsize = 12)
    plot!(pA, tA, repeat(mu, SHOW_DAYS_A); color = :steelblue, lw = 2.2, label = "Average daily profile")
    hline!(pA, [70, 180]; color = :black, ls = :dot, lw = 0.8, label = "")

    # ---------------- Panel B: deviations and slow component ----------------
    iB = ((dB - 1) * PER_DAY + 1):((dB - 1 + SHOW_DAYS_B) * PER_DAY)
    tB = (0:length(iB)-1) ./ PER_DAY
    lim = ceil(maximum(abs.(filter(isfinite, R[iB]))) / 25) * 25 + 10
    pB = plot(tB, R[iB]; color = :gray60, lw = 0.8, label = "Deviation from profile",
              ylabel = "Deviation (mg/dL)", xlabel = "Day", xlims = (0, SHOW_DAYS_B),
              ylims = (-lim, lim), legend = :outertopright,
              title = "b", titleloc = :left, titlefontsize = 12)
    plot!(pB, tB, S[iB]; color = :firebrick, lw = 2.2, label = "Slow component\n(24-h moving average)")
    hline!(pB, [0]; color = :black, lw = 0.8, label = "")

    # ---------------- Panel C: autocorrelation and half-lives ----------------
    nC = SHOW_LAG_H * 12 + 1
    lag_h = (0:nC-1) .* STEP_MIN ./ 60
    lowest = minimum(filter(isfinite, vcat(acR[1:nC], acF[1:nC])))
    ymin   = min(-0.1, floor(lowest * 10) / 10 - 0.05)
    pC = plot(lag_h, acR[1:nC]; color = :gray40, lw = 2, label = "Deviations (unfiltered)",
              ylabel = "Autocorrelation", xlabel = "Lag (hours)", xlims = (0, SHOW_LAG_H),
              ylims = (ymin, 1.02), legend = :outertopright, xticks = 0:2:SHOW_LAG_H,
              title = "c", titleloc = :left, titlefontsize = 12)
    plot!(pC, lag_h, acF[1:nC]; color = :steelblue, lw = 2.2, label = "Fast component")
    hline!(pC, [0.5]; color = :black, ls = :dash, lw = 0.8, label = "")
    hline!(pC, [0]; color = :black, lw = 0.5, label = "")
    vline!(pC, [hlF / 60]; color = :steelblue, ls = :dot, lw = 1.2,
           label = "Fast half-life\n= $(round(Int, hlF)) min")
    vline!(pC, [hlR / 60]; color = :gray40, ls = :dot, lw = 1.2,
           label = "Unfiltered half-life\n= $(round(Int, hlR)) min")

    fig = plot(pA, pB, pC; layout = (3, 1), size = (1050, 1000), dpi = 300,
               left_margin = 5Plots.mm, bottom_margin = 3Plots.mm, right_margin = 3Plots.mm,
               guidefontsize = 10, tickfontsize = 9, legendfontsize = 9, fontfamily = "Helvetica")
    savefig(fig, OUT_FILE)
    println("Saved: ", OUT_FILE)
end

main()