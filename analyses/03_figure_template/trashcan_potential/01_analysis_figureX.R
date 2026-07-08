# Computation step. Reads inputs, runs analysis, saves an intermediate .rds for
# the plot step to read. All paths are resolved from the repo root with
# here::here(), so the script behaves identically regardless of the working
# directory it is launched from.

ANALYSIS_DIR     <- here::here("analyses", "03_figure_template")
INTERMEDIATE_DIR <- file.path(ANALYSIS_DIR, "_intermediate")  # gitignored: heavy + regenerable
dir.create(INTERMEDIATE_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("functions", "<shared-helpers>.R"))            # if needed
source(file.path(ANALYSIS_DIR, "helper_functions_figureX.R"))

# ... analysis ...

saveRDS(results, file.path(INTERMEDIATE_DIR, "fig_X_results.rds"))
