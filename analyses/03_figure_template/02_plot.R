# Plot step. Reads intermediate from 01_compute.R, renders the figure.

source("helper_functions.R")

results <- readRDS("fig_X_results.rds")

# ... plotting ...

ggsave(here::here("outputs", "fig_X.png"), plot = p, width = 7, height = 5)
