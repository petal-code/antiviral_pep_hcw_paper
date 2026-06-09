# ============================================================
# 13_print_final_14_ipcQscaled_WAplusplusFixed_plots.R
#
# Plot-only script for the final 14 scenarios after the WA conflict++
# direction fix.
#
# Reads:
#   final_original_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv
#   final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv
#
# Does not rebuild, refit, or overwrite matrices.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

base_dir <- getwd()
if (!dir.exists(file.path(base_dir, "v6_run_outputs", "final_730day_outputs"))) {
  base_dir <- "C:/Users/jnstapley/Documents/Efficacy_curves/filovirus_14scenario_v8"
}

out_dir <- file.path(base_dir, "v6_run_outputs", "final_730day_outputs")

orig_matrix_path <- file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv")
rev_matrix_path  <- file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv")
audit_path       <- file.path(out_dir, "scenario_manifest_output_resolution_audit.csv")

if (!file.exists(orig_matrix_path)) stop("Cannot find fixed original matrix: ", orig_matrix_path)
if (!file.exists(rev_matrix_path))  stop("Cannot find fixed revised matrix: ", rev_matrix_path)
if (!file.exists(audit_path))       stop("Cannot find audit file: ", audit_path)

plot_parameters <- c(
  "prob_hosp",
  "delay_hosp",
  "prob_unsafe_funeral_comm",
  "prob_unsafe_funeral_hosp",
  "prop_etu",
  "ipc_helper"
)

param_alias <- c(
  "p_hosp" = "prob_hosp",
  "prob_hosp" = "prob_hosp",
  "prob_hospitalised_genPop; prob_hospitalised_hcw" = "prob_hosp",
  "delay_hosp" = "delay_hosp",
  "hospitalisation_delay_factor" = "delay_hosp",
  "p_unsafe_funeral_comm" = "prob_unsafe_funeral_comm",
  "prob_unsafe_funeral_comm" = "prob_unsafe_funeral_comm",
  "p_unsafe_funeral_comm_genPop; p_unsafe_funeral_comm_hcw" = "prob_unsafe_funeral_comm",
  "p_unsafe_funeral_hosp" = "prob_unsafe_funeral_hosp",
  "prob_unsafe_funeral_hosp" = "prob_unsafe_funeral_hosp",
  "p_unsafe_funeral_hosp_genPop; p_unsafe_funeral_hcw" = "prob_unsafe_funeral_hosp",
  "p_unsafe_funeral_hosp_genPop; p_unsafe_funeral_hosp_hcw" = "prob_unsafe_funeral_hosp",
  "p_ETU" = "prop_etu",
  "prop_etu" = "prop_etu",
  "latent_IPC" = "ipc_helper",
  "ipc_helper" = "ipc_helper",
  "ipc" = "ipc_helper",
  "ipc_index" = "ipc_helper",
  "ipc_ppe" = "ipc_helper",
  "ipc_helper; ppe_efficacy_hcw" = "ipc_helper"
)

map_parameter <- function(x) {
  x <- as.character(x)
  out <- unname(param_alias[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

panel_lookup <- c(
  "prob_hosp" = "Probability of hospitalisation",
  "delay_hosp" = "Delay to hospitalisation",
  "prob_unsafe_funeral_comm" = "Unsafe funeral after community death",
  "prob_unsafe_funeral_hosp" = "Unsafe funeral after hospital death",
  "prop_etu" = "Proportion in ETU / ETC",
  "ipc_helper" = "Latent IPC / PPE index"
)

panel_levels <- unname(panel_lookup[plot_parameters])

as_panel <- function(parameter) {
  factor(unname(panel_lookup[parameter]), levels = panel_levels)
}

clean_path <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return(NA_character_)
  x <- as.character(x)[1]
  if (!nzchar(x) || x == "NA") return(NA_character_)
  x
}

safe_read_csv <- function(path) {
  path <- clean_path(path)
  if (is.na(path) || !file.exists(path)) return(tibble())
  suppressMessages(readr::read_csv(path, show_col_types = FALSE))
}

orig <- readr::read_csv(orig_matrix_path, show_col_types = FALSE)
rev  <- readr::read_csv(rev_matrix_path, show_col_types = FALSE)
audit <- readr::read_csv(audit_path, show_col_types = FALSE)
all_final <- bind_rows(orig, rev)

# ---------------- helpers ----------------
get_audit_row <- function(method, scen_key) {
  audit %>% filter(methodology == method, scenario_key == scen_key) %>% slice(1)
}

resolve_anchor_path <- function(method, scen_key) {
  row <- get_audit_row(method, scen_key)
  p <- if (nrow(row) > 0 && "anchor_path" %in% names(row)) clean_path(row$anchor_path[[1]]) else NA_character_
  if (!is.na(p) && file.exists(p)) return(p)

  if (scen_key %in% c("worst_west_africa_conflict", "worst_west_africa_conflict_plusplus")) {
    base_row <- get_audit_row(method, "worst_west_africa")
    p2 <- if (nrow(base_row) > 0 && "anchor_path" %in% names(base_row)) clean_path(base_row$anchor_path[[1]]) else NA_character_
    if (!is.na(p2) && file.exists(p2)) return(p2)
  }

  NA_character_
}

empty_anchor_tbl <- function() {
  tibble(anchor_id = character(), parameter = character(), relative_day = numeric(),
         value_used = numeric(), value_low = numeric(), value_high = numeric(),
         panel_title = factor(levels = panel_levels))
}

read_anchors <- function(path) {
  a <- safe_read_csv(path)
  if (nrow(a) == 0) return(empty_anchor_tbl())

  if (!"anchor_id" %in% names(a)) {
    if ("Parameter" %in% names(a)) a$anchor_id <- as.character(a$Parameter)
    else a$anchor_id <- paste0("anchor_", seq_len(nrow(a)))
  }

  if (!"parameter" %in% names(a)) {
    if ("raw_parameter" %in% names(a)) a$parameter <- a$raw_parameter
    else return(empty_anchor_tbl())
  }

  if (!"relative_day" %in% names(a)) return(empty_anchor_tbl())

  if (!"value_used" %in% names(a)) {
    if ("y_obs" %in% names(a)) a$value_used <- a$y_obs
    else if ("mean" %in% names(a)) a$value_used <- a$mean
    else if ("value" %in% names(a)) a$value_used <- a$value
    else return(empty_anchor_tbl())
  }

  if (!"value_low" %in% names(a)) a$value_low <- a$value_used
  if (!"value_high" %in% names(a)) a$value_high <- a$value_used

  a %>%
    transmute(
      anchor_id = as.character(anchor_id),
      parameter = map_parameter(parameter),
      relative_day = as.numeric(relative_day),
      value_used = as.numeric(value_used),
      value_low = as.numeric(value_low),
      value_high = as.numeric(value_high)
    ) %>%
    filter(parameter %in% plot_parameters, is.finite(relative_day), is.finite(value_used)) %>%
    mutate(
      value_low = if_else(is.finite(value_low), value_low, value_used),
      value_high = if_else(is.finite(value_high), value_high, value_used),
      panel_title = as_panel(parameter)
    )
}

empty_sdb_tbl <- function() {
  tibble(relative_day = numeric(), value = numeric(), n_eligible_sum = numeric(),
         is_drc_plusplus_collapse_point = logical(),
         panel_title = factor(levels = panel_levels))
}

read_sdb <- function(path) {
  s <- safe_read_csv(path)
  if (nrow(s) == 0 || !"relative_day" %in% names(s)) return(empty_sdb_tbl())

  val <- NULL
  if ("p_unsafe_funeral_comm_proxy" %in% names(s)) val <- s$p_unsafe_funeral_comm_proxy
  else if ("success_avg_smoothed" %in% names(s)) val <- 1 - s$success_avg_smoothed
  else if ("prop_success_avg" %in% names(s)) val <- 1 - s$prop_success_avg
  else if ("prop_success_clean" %in% names(s)) val <- 1 - s$prop_success_clean
  else if ("value" %in% names(s)) val <- s$value
  else return(empty_sdb_tbl())

  n_eligible <- if ("n_eligible_sum" %in% names(s)) as.numeric(s$n_eligible_sum)
  else if ("n_eligible" %in% names(s)) as.numeric(s$n_eligible)
  else rep(NA_real_, nrow(s))

  collapse_flag <- if ("is_drc_plusplus_collapse_point" %in% names(s)) as.logical(s$is_drc_plusplus_collapse_point)
  else rep(FALSE, nrow(s))

  tibble(relative_day = as.numeric(s$relative_day), value = as.numeric(val),
         n_eligible_sum = n_eligible,
         is_drc_plusplus_collapse_point = collapse_flag,
         panel_title = factor("Unsafe funeral after community death", levels = panel_levels)) %>%
    filter(is.finite(relative_day), is.finite(value))
}

read_curve_ribbon <- function(path, final_long_for_scenario) {
  x <- safe_read_csv(path)

  fallback <- final_long_for_scenario %>%
    transmute(parameter, panel_title, relative_day, mean, q5 = mean, q95 = mean)

  if (nrow(x) == 0 || !"relative_day" %in% names(x)) return(fallback)

  out <- tibble()

  if ("parameter" %in% names(x)) {
    mean_col <- case_when(
      "mean" %in% names(x) ~ "mean",
      "median" %in% names(x) ~ "median",
      "value" %in% names(x) ~ "value",
      TRUE ~ NA_character_
    )

    if (!is.na(mean_col)) {
      q5_col <- case_when("q5" %in% names(x) ~ "q5", "lower" %in% names(x) ~ "lower", TRUE ~ NA_character_)
      q95_col <- case_when("q95" %in% names(x) ~ "q95", "upper" %in% names(x) ~ "upper", TRUE ~ NA_character_)

      out <- x %>%
        transmute(
          parameter = map_parameter(parameter),
          relative_day = as.numeric(relative_day),
          mean = as.numeric(.data[[mean_col]]),
          q5 = if (!is.na(q5_col)) as.numeric(.data[[q5_col]]) else as.numeric(.data[[mean_col]]),
          q95 = if (!is.na(q95_col)) as.numeric(.data[[q95_col]]) else as.numeric(.data[[mean_col]])
        )
    }
  } else {
    cols <- intersect(c(names(param_alias), plot_parameters), names(x))
    cols <- unique(cols)
    if (length(cols) > 0) {
      out <- x %>%
        select(relative_day, all_of(cols)) %>%
        pivot_longer(cols = all_of(cols), names_to = "parameter", values_to = "mean") %>%
        transmute(parameter = map_parameter(parameter),
                  relative_day = as.numeric(relative_day),
                  mean = as.numeric(mean),
                  q5 = mean,
                  q95 = mean)
    }
  }

  if (nrow(out) == 0) return(fallback)

  out <- out %>%
    filter(parameter %in% plot_parameters, is.finite(relative_day),
           is.finite(mean), is.finite(q5), is.finite(q95)) %>%
    mutate(q_low = pmin(q5, q95), q_high = pmax(q5, q95)) %>%
    group_by(parameter, relative_day) %>%
    summarise(mean = median(mean, na.rm = TRUE),
              q5 = median(q_low, na.rm = TRUE),
              q95 = median(q_high, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(parameter, relative_day) %>%
    mutate(panel_title = as_panel(parameter))

  if (nrow(out) == 0) return(fallback)

  end_ext <- out %>%
    group_by(parameter) %>%
    arrange(relative_day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    filter(relative_day < 730) %>%
    mutate(relative_day = 730)

  bind_rows(out, end_ext) %>% arrange(parameter, relative_day)
}

find_native_shared_points <- function(audit_row) {
  candidates <- character()
  for (p in c("matrix_path", "curve_path", "sdb_path", "anchor_path", "q_path")) {
    if (p %in% names(audit_row)) {
      val <- clean_path(audit_row[[p]][[1]])
      if (!is.na(val)) {
        d <- dirname(val)
        candidates <- c(candidates, list.files(d, pattern = "shared_conflict_Q_points\\.csv$", full.names = TRUE))
      }
    }
  }
  candidates <- c(
    candidates,
    list.files(file.path(base_dir, "v6_run_outputs"), pattern = "shared_conflict_Q_points\\.csv$",
               full.names = TRUE, recursive = TRUE)
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0) return(NA_character_)
  scen <- audit_row$scenario_key[[1]]
  if (grepl("plusplus", scen)) {
    pp <- candidates[grepl("plusplus|directCollapse|PlusPlus", basename(candidates), ignore.case = TRUE)]
    if (length(pp) > 0) return(pp[1])
  }
  cc <- candidates[grepl("conflict|warsame", basename(candidates), ignore.case = TRUE)]
  if (length(cc) > 0) return(cc[1])
  candidates[1]
}

plot_one <- function(method, scen_key) {
  scen_df <- all_final %>% filter(methodology == method, scenario_key == scen_key)
  scen_audit <- get_audit_row(method, scen_key)

  df_long <- scen_df %>%
    select(methodology, scenario_key, scenario, relative_day, all_of(plot_parameters)) %>%
    pivot_longer(cols = all_of(plot_parameters), names_to = "parameter", values_to = "mean") %>%
    mutate(parameter = as.character(parameter),
           relative_day = as.numeric(relative_day),
           mean = as.numeric(mean),
           panel_title = as_panel(parameter)) %>%
    arrange(parameter, relative_day)

  ribbon <- read_curve_ribbon(scen_audit$curve_path[[1]], df_long)

  if (method == "revised" && scen_key %in% c("middle_drc_conflict", "middle_drc_conflict_plusplus")) {
    ribbon <- ribbon %>%
      filter(parameter != "ipc_helper") %>%
      bind_rows(
        df_long %>%
          filter(parameter == "ipc_helper") %>%
          transmute(parameter, panel_title, relative_day, mean, q5 = mean, q95 = mean)
      )
  }

  anchor_path <- resolve_anchor_path(method, scen_key)
  anchors <- read_anchors(anchor_path) %>% filter(relative_day <= 730)

  if (method == "original" && grepl("middle_drc_conflict", scen_key)) {
    native_path <- find_native_shared_points(scen_audit)
    sdb <- read_sdb(native_path)
    if (nrow(sdb) == 0) sdb <- read_sdb(scen_audit$sdb_path[[1]])
  } else {
    sdb <- read_sdb(scen_audit$sdb_path[[1]])
  }

  sdb_regular <- sdb %>% filter(!is_drc_plusplus_collapse_point)
  sdb_collapse <- sdb %>% filter(is_drc_plusplus_collapse_point)

  # IPC 0-1 support only.
  ipc_axis_support <- tibble(
    panel_title = factor("Latent IPC / PPE index", levels = panel_levels),
    relative_day = 0,
    y = c(0, 1)
  )

  scenario_label <- if ("scenario" %in% names(scen_df)) unique(scen_df$scenario)[1] else scen_key
  title <- paste0(stringr::str_to_title(method), " methodology: ", scenario_label)

  p <- ggplot() +
    geom_blank(data = ipc_axis_support, aes(x = relative_day, y = y), inherit.aes = FALSE) +
    geom_ribbon(data = ribbon, aes(x = relative_day, ymin = q5, ymax = q95, group = parameter),
                alpha = 0.20, fill = "lightblue", inherit.aes = FALSE) +
    geom_line(data = df_long, aes(x = relative_day, y = mean, group = parameter),
              linewidth = 0.9, colour = "black", inherit.aes = FALSE) +
    geom_point(data = sdb_regular, aes(x = relative_day, y = value, size = n_eligible_sum),
               inherit.aes = FALSE, colour = "grey45", alpha = 0.75) +
    geom_point(data = sdb_collapse, aes(x = relative_day, y = value),
               inherit.aes = FALSE, shape = 21, colour = "red3", fill = NA, stroke = 1.0, size = 3.2) +
    geom_errorbar(data = anchors, aes(x = relative_day, ymin = value_low, ymax = value_high),
                  inherit.aes = FALSE, width = 0, colour = "orange3", linewidth = 0.55) +
    geom_point(data = anchors, aes(x = relative_day, y = value_used),
               inherit.aes = FALSE, colour = "orange3", size = 2.1) +
    geom_text(data = anchors, aes(x = relative_day, y = value_used, label = anchor_id),
              inherit.aes = FALSE, colour = "orange4", size = 2.4,
              hjust = -0.05, vjust = -0.6, check_overlap = TRUE) +
    scale_size_continuous(range = c(1.2, 4.0), guide = "none") +
    facet_wrap(~panel_title, scales = "free_y", ncol = 2) +
    labs(title = title,
         subtitle = "730-day final matrix plot; black line from WA++-fixed final matrix; native q5–q95 ribbons and orange anchors where available",
         x = "Relative day",
         y = "Parameter value") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  attr(p, "n_anchor_points") <- nrow(anchors)
  attr(p, "n_ribbon_rows") <- nrow(ribbon)
  attr(p, "n_ribbon_nonzero_width") <- sum(abs(ribbon$q95 - ribbon$q5) > 1e-10, na.rm = TRUE)

  p
}

# Safety check.
row_check <- all_final %>%
  group_by(methodology, scenario_key) %>%
  summarise(n_rows = n(), n_days = n_distinct(relative_day),
            min_day = min(relative_day), max_day = max(relative_day),
            duplicate_days = n_rows - n_days,
            .groups = "drop")

bad <- row_check %>% filter(n_rows != 731 | n_days != 731 | min_day != 0 | max_day != 730 | duplicate_days != 0)
if (nrow(bad) > 0) {
  print(bad)
  stop("Final matrices failed row-count check.")
}

scenario_order <- all_final %>%
  distinct(methodology, scenario_key) %>%
  mutate(methodology = factor(methodology, levels = c("original", "revised"))) %>%
  arrange(methodology, scenario_key)

plots14_final_WAplusplusFixed <- list()
plot_diagnostics <- list()

for (i in seq_len(nrow(scenario_order))) {
  method <- as.character(scenario_order$methodology[i])
  scen <- scenario_order$scenario_key[i]
  key <- paste(method, scen, sep = "__")

  p <- plot_one(method, scen)
  plots14_final_WAplusplusFixed[[key]] <- p
  plot_diagnostics[[key]] <- tibble(
    plot_key = key,
    methodology = method,
    scenario_key = scen,
    n_anchor_points = attr(p, "n_anchor_points"),
    n_ribbon_rows = attr(p, "n_ribbon_rows"),
    n_ribbon_nonzero_width = attr(p, "n_ribbon_nonzero_width")
  )
  print(p)
}

plot_diagnostics <- bind_rows(plot_diagnostics)

message("Done. Plots stored in object: plots14_final_WAplusplusFixed")
message("Diagnostics stored in object: plot_diagnostics")
print(plot_diagnostics)

message("Final matrices used:")
message(" - ", orig_matrix_path)
message(" - ", rev_matrix_path)
