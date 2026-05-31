# Computation step. Reads inputs, runs analysis, saves intermediate as .rds.

source(here::here("functions", "<shared-helpers>.R"))  # if needed
source("helper_functions_figureX.R")

# ... analysis ...

saveRDS(results, "fig_X_results.rds")
