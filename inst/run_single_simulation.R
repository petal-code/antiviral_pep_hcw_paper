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
library(ggplot2)
source(here::here("functions", "setup_model_parameters.R"))  # parameter builders
source(here::here("functions", "calculate_model_approx_r0.R"))  # R0 <-> offspring-means solver

SCENARIO_CSV <- here::here("data-processed", "final_four_scenario_values.csv")

# The two knobs a reader is most likely to want to change:
SCENARIO_ID <- "Worst_WestAfrica"   # which scenario's column-set to simulate
SEED        <- 23                    # fix the RNG so the run is reproducible

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
  overrides       = list(       # You can use overrides to change any of the model parameters from the defaults
    seeding_cases    = 20,      # initial community infections that seed the outbreak
    check_final_size = 20000,   # stop the model once this many cases exist (demo cap for speed)
    mn_offspring_genPop = 1.2
  )
)

# --- Derive the offspring means from a target baseline R0 + funeral share ------
# Rather than hand-pick mn_offspring_genPop, we specify the epidemiology we want
# at t = 0 and let the solver back out the means. solve_offspring_means_for_R0()
# resolves prop_etu(0) and the efficacy scalars to get the direct (D) and funeral
# (F) R0 multipliers at t = 0, then inverts:
#     mn_offspring_genPop  = (1 - funeral_share) * R0 / D
#     mn_offspring_funeral =      funeral_share  * R0 / F
# R0 here is the BASELINE (t = 0) value; the West Africa scenario's control
# measures (rising ETU coverage, hospitalisation, IPC) pull Rt below 1 over time,
# so the outbreak peaks then declines instead of running to the cap.
TARGET_R0    <- 1.3    # baseline reproduction number (fitted WA was ~1.45 -> near the cap; 1.3 leaves headroom)
FUNERAL_FRAC <- 0.25   # share of t = 0 transmission via unsafe funerals

sol <- solve_offspring_means_for_R0(
  R0   = TARGET_R0,
  args = mp$args,                                       # carries prop_etu(0), efficacies, death/funeral probs
  proportion_transmission_from_funerals = FUNERAL_FRAC,
  seed = SEED                                           # reproduces the Monte-Carlo step inside the solver
)

# Hand the solved means to the model (mn_offspring_hcw stays at its default —
# it is not part of the genPop-dominant single-type R0 inversion).
# `mp$args` is the complete argument list for the model. The RNG seed is a
# direct argument to branching_process_main() (not a DEFAULT_SCALAR_INPUTS
# scalar), so we attach it to the args list here.
args <- mp$args
args$mn_offspring_genPop  <- sol$mn_offspring_genPop_required
args$mn_offspring_funeral <- sol$mn_offspring_funeral_required
args$seed <- SEED

## Plotting the time-varying parameters
# The time varying probabilities/hospitalisation delay factor - plotted below
days <- 0:max(mp$scenario_matrix$relative_day)   # daily grid spanning the scenario
curves <- c("prob_hospitalised_genPop", "hospitalisation_delay_factor",
            "prop_etu", "ppe_coverage_hcw",
            "p_unsafe_funeral_comm_genPop", "p_unsafe_funeral_hosp_genPop")

# resolve_time_varying() (from fiber) evaluates a scalar-or-function input at the
# given times; here it returns each curve's values on `days`. Stack to long form.
tv_values <- sapply(curves, function(nm) fiber:::resolve_time_varying(mp$tv_args[[nm]], days, nm))
tv_long <- data.frame(
  day   = rep(days, times = length(curves)),
  input = rep(curves, each  = length(days)),
  value = as.vector(tv_values)
)

ggplot(tv_long, aes(day, value)) +
  geom_line(colour = "steelblue", linewidth = 0.8) +
  facet_wrap(~ input, scales = "free_y") +   # free y: each input on its own scale
  labs(x = "Time since outbreak start (days)", y = "Input value",
       title = sprintf("Time-varying model inputs: %s", SCENARIO_ID)) +
  theme_bw()

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
# note that if the cap is reached, some of the infections in tdf won't have had offspring generated
# important to bear this in mind


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


# ---------------------------------------------------------------------------
# 5. Epidemic-incidence plot (new infections per week, deaths overlaid)
# ---------------------------------------------------------------------------
# Bin deaths onto the SAME weekly grid, then draw with base graphics. The bar
# midpoints are the left edge of each week plus half a bin width.
weekly_deaths <- hist(cases$time_outcome_absolute[died],
                      breaks = week_breaks, plot = FALSE)$counts
week_mid <- week_breaks[-length(week_breaks)] + bin_width / 2

plot(week_mid, weekly_infections, type = "h", lwd = 4, lend = 1,
     col = "steelblue",
     xlab = "Time since outbreak start (days)", ylab = "New cases per week",
     main = sprintf("Single simulated outbreak: %s", SCENARIO_ID))
lines(week_mid, weekly_deaths, type = "o", pch = 16, cex = 0.6, col = "firebrick")
legend("topright", bty = "n",
       legend = c("New infections", "Deaths"),
       col = c("steelblue", "firebrick"), lwd = c(4, 1), pch = c(NA, 16))

# ---------------------------------------------------------------------------
# 6. Repeat the SAME run with different seeds to see stochastic variability
# ---------------------------------------------------------------------------
# Identical parameters (`args`) -- only the RNG seed changes between runs. This
# shows how much two outbreaks can differ from pure chance alone. We collect
# each run's weekly incidence into one long data frame tagged by iteration, then
# facet so all 12 realisations can be compared side by side.
N_ITER <- 12

# Run all 12 first, keeping each run's realised infection-day vector, so we can
# then bin every realisation onto ONE shared weekly grid (comparable x-axes).
# (The anonymous function(i) is just an inline loop body, not a new helper.)
inf_times <- lapply(seq_len(N_ITER), function(i) {
  a <- args; a$seed <- i                       # same args, different seed each time
  t <- do.call(branching_process_main, a)$tdf$time_infection_absolute
  t[!is.na(t)]                                  # keep realised infections only
})

# Common weekly grid spanning the longest of the 12 epidemics (reuses bin_width
# = 7 days from step 4 so these curves match the single-run incidence plot).
breaks12 <- seq(0, ceiling(max(unlist(inf_times)) / bin_width) * bin_width, by = bin_width)
mids12   <- breaks12[-length(breaks12)] + bin_width / 2

# Bin each run onto that grid and stack into long form, tagged by iteration.
epi_long <- do.call(rbind, lapply(seq_len(N_ITER), function(i) {
  data.frame(iteration = i,
             day       = mids12,
             incidence = hist(inf_times[[i]], breaks = breaks12, plot = FALSE)$counts)
}))

# One panel per realisation, each coloured by its iteration number.
ggplot(epi_long, aes(day, incidence, colour = factor(iteration))) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ iteration) +
  labs(x = "Time since outbreak start (days)", y = "New cases per week",
       title = sprintf("12 stochastic realisations, identical parameters: %s", SCENARIO_ID)) +
  theme_bw() +
  theme(legend.position = "none")   # the facet strip already labels each iteration
