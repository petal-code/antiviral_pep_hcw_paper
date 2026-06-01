# =============================================================================
# run_single_simulation.R
#
# A minimal, heavily-commented walk-through of the full fiber simulation
# workflow: turn a scenario into model parameters, run ONE stochastic outbreak
# with branching_process_main(), then summarise the result as (1) a table of
# key metrics and (2) an epidemic-incidence plot.
#
# It deliberately defines NO new functions -- it only wires together the
# project's existing helpers (functions/setup_model_parameters.R) and the fiber
# model with plain base R, so it doubles as a readable example of how the pieces
# fit together.
#
# Run from the repo root with:  Rscript inst/run_single_simulation.R
# =============================================================================


# ---------------------------------------------------------------------------
# 0. Setup: locate the project, load fiber, source the parameter helpers
# ---------------------------------------------------------------------------
# here::here() finds the repo root by locating obv_hcw_paper.Rproj, so every
# path below resolves the same regardless of the directory this is launched from.
library(fiber)                                               # the model engine
source(here::here("functions", "setup_model_parameters.R"))  # parameter builders

SCENARIO_CSV <- here::here("data-processed", "final_four_scenario_values.csv")
OUTPUT_DIR   <- here::here("outputs", "single_simulation_demo")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# The two knobs a reader is most likely to want to change:
SCENARIO_ID <- "Worst_WestAfrica"   # which scenario's column-set to simulate
SEED        <- 1                    # fix the RNG so the run is reproducible

# Fail fast (with a clear message) if the installed fiber predates the NPI
# parameter interface this workflow targets.
check_model_function_version()


# ---------------------------------------------------------------------------
# 1. Read the scenario matrix
# ---------------------------------------------------------------------------
# The CSV holds one row per time-step per scenario: time-varying probabilities,
# delay factors, and the ETU / IPC coverage curves. read_scenario_matrix()
# loads it and checks the required columns are present.
scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)


# ---------------------------------------------------------------------------
# 2. Build the model parameters for this scenario
# ---------------------------------------------------------------------------
# make_model_parameters() ties three steps together:
#   * make_base_args()          -- scalar parameters + the natural-history delay
#                                  samplers, from DEFAULT_SCALAR_INPUTS with any
#                                  `overrides` applied on top;
#   * build_time_varying_args() -- turns the scenario rows into time-varying
#                                  curves (prob_hospitalised, prop_etu,
#                                  ppe_coverage_hcw, ...);
#   * returns them combined in `$args`, ready to hand straight to the model.
#
# `overrides` only touches scalars that exist in DEFAULT_SCALAR_INPUTS. Here we
# shrink the run for a fast, self-contained demo: more seeding cases (so the
# outbreak reliably takes off) and a smaller final-size cap. Drop `overrides`
# entirely to use the literature-informed defaults.
mp <- make_model_parameters(
  scenario_id     = SCENARIO_ID,
  scenario_matrix = scenario_matrix,
  overrides       = list(
    seeding_cases    = 20,      # initial community infections that seed the outbreak
    check_final_size = 10000    # stop once this many cases exist (demo cap for speed)
  )
)

# `mp$args` is the complete argument list for the model. The RNG seed is a
# direct argument to branching_process_main() (not a DEFAULT_SCALAR_INPUTS
# scalar), so we attach it to the args list here.
args <- mp$args
args$seed <- SEED


# ---------------------------------------------------------------------------
# 3. Run a single stochastic outbreak
# ---------------------------------------------------------------------------
# do.call() unpacks the named `args` list into branching_process_main(), which
# iteratively generates community, hospital and funeral offspring until the
# outbreak dies out or reaches check_final_size.
out <- do.call(branching_process_main, args)

# `out$tdf` is the transmission tree: one row per individual. Rows that were
# pre-allocated but never used have a missing infection time, so we keep only
# the realised cases.
tree  <- out$tdf
cases <- tree[!is.na(tree$time_infection_absolute), ]


# ---------------------------------------------------------------------------
# 4. Key summary metrics (computed inline from the transmission tree)
# ---------------------------------------------------------------------------
# Useful columns: `outcome` is TRUE for a death (FALSE = recovery); `class` is
# "HCW" or "genPop"; the *_absolute columns are calendar-day timings measured
# from the start of the outbreak (day 0).
died <- !is.na(cases$outcome) & cases$outcome
hcw  <- cases$class == "HCW"

# Weekly bins span from day 0 to just past the last event. We size the grid from
# BOTH infection and outcome times so deaths (which happen after infection)
# always fall inside the range when we reuse these breaks for the plot below.
bin_width   <- 7                                                  # days per bin (weekly)
last_day    <- max(c(cases$time_infection_absolute,
                     cases$time_outcome_absolute), na.rm = TRUE)
week_breaks <- seq(0, ceiling(last_day / bin_width) * bin_width, by = bin_width)

# hist(..., plot = FALSE) is a convenient way to bin a vector into fixed windows.
weekly_infections <- hist(cases$time_infection_absolute,
                          breaks = week_breaks, plot = FALSE)$counts

summary_table <- data.frame(
  metric = c("Scenario", "Total infections", "Total deaths",
             "Case fatality ratio", "HCW infections", "HCW deaths",
             "Outbreak duration (days)", "Generations", "Peak weekly infections"),
  value  = c(SCENARIO_ID,
             nrow(cases),
             sum(died),
             round(sum(died) / nrow(cases), 3),
             sum(hcw),
             sum(died & hcw),
             round(max(cases$time_outcome_absolute, na.rm = TRUE), 1),
             max(cases$generation),
             max(weekly_infections)),
  stringsAsFactors = FALSE
)

print(summary_table, row.names = FALSE)
write.csv(summary_table,
          file.path(OUTPUT_DIR, "single_simulation_summary.csv"),
          row.names = FALSE)


# ---------------------------------------------------------------------------
# 5. Epidemic-incidence plot (new infections per week, deaths overlaid)
# ---------------------------------------------------------------------------
# Bin deaths onto the SAME weekly grid, then draw with base graphics. The bar
# midpoints are the left edge of each week plus half a bin width.
weekly_deaths <- hist(cases$time_outcome_absolute[died],
                      breaks = week_breaks, plot = FALSE)$counts
week_mid <- week_breaks[-length(week_breaks)] + bin_width / 2

png(file.path(OUTPUT_DIR, "single_simulation_incidence.png"),
    width = 1800, height = 1100, res = 200)

plot(week_mid, weekly_infections, type = "h", lwd = 4, lend = 1,
     col = "steelblue",
     xlab = "Time since outbreak start (days)", ylab = "New cases per week",
     main = sprintf("Single simulated outbreak: %s", SCENARIO_ID))
lines(week_mid, weekly_deaths, type = "o", pch = 16, cex = 0.6, col = "firebrick")
legend("topright", bty = "n",
       legend = c("New infections", "Deaths"),
       col = c("steelblue", "firebrick"), lwd = c(4, 1), pch = c(NA, 16))

invisible(dev.off())

message("Wrote summary table and incidence plot to: ", OUTPUT_DIR)
