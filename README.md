# glycaemic-recovery-time

Julia code for glycaemic recovery time: the half-life of deviations from a person's average daily CGM profile.

## The measure

Glycaemic recovery time describes how quickly glucose returns to a person's usual pattern after departing from it. It is calculated from continuous glucose monitoring (CGM) data in three steps:

1. **Average daily profile:** the mean glucose at each clock time (5-minute resolution) across a 28-day window.
2. **Deviations:** glucose minus the profile. A centred 24-hour moving average of the deviations (the slow component, mainly whole-day shifts) is removed, leaving the fast component.
3. **Half-life:** the lag at which the autocorrelation of the fast component first falls below 0.5, obtained by linear interpolation. This is the time after which a deviation is expected to have halved.

The repository also contains the analyses of test–retest reliability, the relation to standard CGM metrics and autocorrelation indices, the 24-hour recovery profile, HbA1c, insulin delivery, and sensitivity analyses.

## Data

No data are included.

- **DiaMonT** (type 2 diabetes): the data are confidential and cannot be shared. The DiaMonT scripts are provided for transparency.
- **REPLACE-BG** (type 1 diabetes): the data are publicly available from the Jaeb Center for Health Research (https://public.jaeb.org/dataset/546). Download the dataset and point the scripts to its `Data Tables` folder.

## Requirements

- Julia 1.12.1 with the packages CSV, DataFrames, Distributions, GLM, HypothesisTests, JLD2, Plots, ReadStatTables, and StatsBase
- Python 3 with matplotlib (participant flowchart only)

## Structure

- `data_preparation/`: reads the raw CGM data, places readings on a 5-minute grid, and selects cohorts, analysis windows and reliability subsets
- `analysis/`: recovery time and all analyses
- `figures/`: tables and figures

## Usage

1. Set the paths at the top of each script. They are marked as placeholders, for example `"<path to REPLACE-BG derived data>"`.
2. Run each script in a fresh Julia session.
3. Run the scripts in this order:

**Data preparation**
| Script | Purpose |
|---|---|
| `build_replacebg.jl` | REPLACE-BG: CGM on a 5-minute grid, cohort and analysis windows |
| `build_dataset_v2.jl` | DiaMonT: CGM and insulin records on a 5-minute grid |
| `select_cohort.jl` | DiaMonT: cohorts and analysis windows |
| `reliability_windows.jl` | Reliability subsets (at least 56 consecutive usable days) |
| `make_window_cache.jl` | Cache of analysis windows used by later scripts |

**Analyses**
| Script | Purpose |
|---|---|
| `t1d_validation.jl` | REPLACE-BG: recovery time, reliability, relation to CGM metrics, 24-hour profile |
| `t1d_pump_analysis.jl` | REPLACE-BG: insulin delivery and recovery time |
| `profile_window_experiment.jl` | Choice of training length and type of the daily profile |
| `persistence_test1.jl` | Unfiltered half-life and autocorrelation of the deviations |
| `persistence_fast_slow.jl` | Fast and slow components |
| `halflife_reliability.jl` | Test–retest reliability of the half-life |
| `estimator_reliability.jl` | Reliability of alternative summaries of the deviations |
| `filter_sensitivity.jl` | Sensitivity to the width of the moving average |
| `icc_difference.jl` | ICCs with bootstrap confidence intervals and within-participant differences |
| `fast_nonredundancy.jl` | Relation to standard CGM metrics |
| `variability_metrics.jl` | Relation to MAGE, CONGA1 and MAG |
| `ac_comparison.jl` | Comparison with the autocorrelation indices AC_Mean and AC_Var |
| `persistence_profile.jl` | 24-hour recovery profile |
| `persistence_hba1c.jl` | Relation to HbA1c |
| `selection_comparison.jl` | Analysed versus excluded participants |
| `bolus_start.jl`, `bolus_start_adjusted.jl` | Change in recovery time after starting short-acting insulin |
| `halflife_15min.jl`, `acvar_15min.jl`, `sampling_15min_t1d.jl` | Sensitivity to 15-minute sampling |

**Tables and figures**
| Script | Purpose |
|---|---|
| `table1.jl`, `demographics_replacebg.jl`, `bmi_diamont.jl` | Cohort characteristics |
| `figure_method.jl` | Calculation of recovery time in one participant |
| `figure_flowchart.py` | Participant flowchart |
| `figure_reliability.jl` | Test–retest reliability |
| `figure_profile.jl` | 24-hour recovery profile |

Scripts print only aggregated results. Results are saved as `.jld2` files in the data or results folder set at the top of each script.

## Use and citation

The code may be used for research purposes. If you use the code or the method, please cite:

[Citation will be added when the paper is published]

For questions, contact Tanja Kortsen Bugajski, Department of Mathematical Sciences, Aalborg University.