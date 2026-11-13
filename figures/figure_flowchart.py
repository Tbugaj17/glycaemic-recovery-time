import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# Figure 1: participant flow, DiaMonT (T2D) and REPLACE-BG (T1D).
WHITE, GREY, COHORT, SUBSET = "white", "#d9d9d9", "#a9cbea", "#e3eef9"
FS = 13

fig, ax = plt.subplots(figsize=(20, 17))
ax.set_xlim(0, 20); ax.set_ylim(0.6, 17.2); ax.axis("off")

def box(x, y, w, h, text, color=WHITE, bold=False):
    # Border style as a second cue besides colour (colour-blind readers):
    # exclusions dashed, analysis cohorts thick, all other boxes thin and solid.
    ls, lw = "-", 1.6
    if color == GREY:
        ls = (0, (4, 2.5))
    elif color == COHORT:
        lw = 3.2
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0,rounding_size=0.12",
                                fc=color, ec="black", lw=lw, ls=ls))
    ax.text(x, y, text, ha="center", va="center", fontsize=FS + (1 if bold else 0),
            fontweight="bold" if bold else "normal", linespacing=1.3)

def arrow(x, y0, y1):
    ax.annotate("", xy=(x, y1), xytext=(x, y0),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="black", mutation_scale=16))

def line(x0, y0, x1, y1):
    ax.plot([x0, x1], [y0, y1], color="black", lw=1.6)

H = 1.0           # height of main boxes
y_rand, y_cgm, y_pen, y_reg = 15.4, 13.3, 11.2, 9.4
y_elig, y_short, y_coh, y_rel = 7.9, 6.4, 4.9, 3.2

# ---------------- DiaMonT ----------------
ax.text(4.9, 16.6, "DiaMonT (T2D)", ha="center", fontsize=18, fontweight="bold")
xm, xb, xbb = 4.9, 2.3, 7.45
box(xm, y_rand, 4.6, H, "Randomised\nn = 331")
box(9.7, 14.35, 3.9, 0.85, "No CGM data\nn = 31", GREY)
arrow(xm, y_rand - H/2, y_cgm + H/2); line(xm, 14.35, 9.7 - 1.95, 14.35)
box(xm, y_cgm, 4.6, H, "CGM data\nn = 300")
box(9.7, 12.25, 3.9, 0.85, "No insulin-pen records\nn = 8", GREY)
arrow(xm, y_cgm - H/2, y_pen + H/2); line(xm, 12.25, 9.7 - 1.95, 12.25)
box(xm, y_pen, 5.3, H, "CGM and insulin-pen records\nn = 292")
y_split = y_pen - H/2 - 0.4
line(xm, y_pen - H/2, xm, y_split); line(xb, y_split, xbb, y_split)
arrow(xb, y_split, y_reg + H/2); arrow(xbb, y_split, y_reg + H/2)
box(xb, y_reg, 3.2, H, "Basal-only\nn = 137")
box(xbb, y_reg, 3.2, H, "Basal-bolus\nn = 155")

# basal-only branch
arrow(xb, y_reg - H/2, y_coh + H/2)
box(4.6, y_elig, 2.3, 1.2, "No eligible\nwindow$^{\\mathrm{a}}$\nn = 87", GREY); line(xb, y_elig, 4.6 - 1.15, y_elig)
box(xb, y_coh, 3.2, H, "Analysis cohort\nn = 50", COHORT, bold=True)
arrow(xb, y_coh - H/2, y_rel + 0.58)
box(xb, y_rel, 3.4, 1.15, "Reliability subset$^{\\mathrm{b}}$\nn = 35", SUBSET)

# basal-bolus branch
arrow(xbb, y_reg - H/2, y_coh + H/2)
box(9.75, y_elig, 2.3, 1.2, "No eligible\nwindow$^{\\mathrm{a}}$\nn = 58", GREY); line(xbb, y_elig, 9.75 - 1.15, y_elig)
box(9.75, y_short, 2.3, 1.2, "No short-acting\ninsulin in\nwindow$^{\\mathrm{d}}$, n = 17", GREY); line(xbb, y_short, 9.75 - 1.15, y_short)
box(xbb, y_coh, 3.2, H, "Analysis cohort\nn = 80", COHORT, bold=True)
arrow(xbb, y_coh - H/2, y_rel + 0.58)
box(xbb, y_rel, 3.4, 1.15, "Reliability subset$^{\\mathrm{b}}$\nn = 56", SUBSET)

# divider
line(11.85, 16.9, 11.85, 2.2); ax.lines[-1].set_linestyle("--"); ax.lines[-1].set_color("grey"); ax.lines[-1].set_lw(1)

# ---------------- REPLACE-BG ----------------
xr = 15.35
ax.text(xr, 16.6, "REPLACE-BG (T1D)", ha="center", fontsize=18, fontweight="bold")
box(xr, y_rand, 4.6, H, "Randomised\nn = 226")
arrow(xr, y_rand - H/2, y_coh + H/2)
box(18.0, y_pen - 1.1, 2.7, 1.2, "No eligible\nwindow$^{\\mathrm{a}}$\nn = 10", GREY); line(xr, y_pen - 1.1, 18.0 - 1.35, y_pen - 1.1)
box(xr, y_coh, 3.9, H, "Analysis cohort\nn = 216", COHORT, bold=True)
y_s2 = y_coh - H/2 - 0.4
line(xr, y_coh - H/2, xr, y_s2); line(13.8, y_s2, 16.9, y_s2)
arrow(13.8, y_s2, y_rel + 0.58); arrow(16.9, y_s2, y_rel + 0.58)
box(13.8, y_rel, 2.9, 1.15, "Reliability\nsubset$^{\\mathrm{b}}$\nn = 131", SUBSET)
box(16.9, y_rel, 2.9, 1.15, "Insulin-pump\nanalyses$^{\\mathrm{c}}$\nn = 197", SUBSET)

# ---------------- footnotes ----------------
foot = ("$^{\\mathrm{a}}$ No period of 28 consecutive days with ≥70% of CGM readings observed and ≥26 of the 28 days individually ≥70%; in DiaMonT, this included all\n"
        "   standard-of-care participants (blinded CGM during the first and last 20 days only); in REPLACE-BG, windows were required to lie between\n"
        "   randomisation and the last visit, without interruption of pump therapy.   $^{\\mathrm{b}}$ At least 56 consecutive usable days.   "
        "$^{\\mathrm{c}}$ Basal records covering ≥70% of the analysis window.\n"
        "$^{\\mathrm{d}}$ Short-acting insulin recorded during the study but not during the 28-day analysis window.")
ax.text(0.2, 1.9, foot, ha="left", va="top", fontsize=11.5, color="#333333", linespacing=1.5)

fig.savefig("figure_flowchart.png", dpi=150, bbox_inches="tight")
print("Saved figure_flowchart.png")
