# Plot step. Reads the intermediate saved by 01_analysis_figureX.R and renders
# the figure. Paths are resolved with here::here() so the working directory the
# script is launched from does not matter. The intermediate is read from the
# gitignored _intermediate/ folder; the finished figure is written to the
# per-analysis outputs/ folder.

ANALYSIS_DIR     <- here::here("analyses", "03_figure_template")
INTERMEDIATE_DIR <- file.path(ANALYSIS_DIR, "_intermediate")
OUTPUT_DIR       <- here::here("outputs", "03_figure_template")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(ANALYSIS_DIR, "helper_functions_figureX.R"))

results <- readRDS(file.path(INTERMEDIATE_DIR, "fig_X_results.rds"))

# ... plotting ...

ggsave(file.path(OUTPUT_DIR, "fig_X.png"), plot = p, width = 7, height = 5)
