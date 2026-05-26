# Computation step. Reads inputs, runs analysis, saves intermediate as .rds.

source(here::here("R", "<top-level-helpers>.R"))  # if needed
source("helper_functions.R")

# ... analysis ...

saveRDS(results, "fig_X_results.rds")
