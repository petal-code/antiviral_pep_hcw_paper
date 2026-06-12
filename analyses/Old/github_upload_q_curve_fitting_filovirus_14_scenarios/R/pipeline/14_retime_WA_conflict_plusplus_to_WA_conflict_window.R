# ============================================================
# 14_retime_WA_conflict_plusplus_to_WA_conflict_window.R
#
# Purpose:
#   Fix only the two West Africa conflict++ scenarios.
#
# Charlie's requested logic:
#   - Standard West Africa-with-conflict remains unchanged.
#   - West Africa conflict++ should NOT add an extra DRC-timed collapse.
#   - Instead, it should replace the old mild WA conflict dip with a
#     severe DRC++-style collapse at the same timing as the existing
#     WA-with-conflict dip.
#
# Inputs:
#   final_*_ipcQscaled_WAplusplusFixed.csv if present, otherwise
#   final_*_ipcQscaled.csv.
#
# Outputs:
#   final_original_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusRetimed.csv
#   final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusRetimed.csv
#
# The script:
#   1. Detects the old WA conflict timing by comparing
#      worst_west_africa_conflict vs worst_west_africa.
#   2. Rebuilds worst_west_africa_conflict_plusplus from ordinary
#      worst_west_africa, with a severe collapse imposed over the
#      detected WA conflict window.
#   3. Writes updated final original/revised CSVs.
#   4. Prints only the two updated WA conflict++ plots for checking.
#
# It does NOT modify DRC scenarios or standard WA-with-conflict.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
base_dir <- getwd()
if (!dir.exists(file.path(base_dir, "v6_run_outputs", "final_730day_outputs"))) {
  base_dir <- "C:/Users/jnstapley/Documents/Efficacy_curves/filovirus_14scenario_v8"
}

out_dir <- file.path(base_dir, "v6_run_outputs", "final_730day_outputs")
audit_path <- file.path(out_dir, "scenario_manifest_output_resolution_audit.csv")

orig_in_candidates <- c(
  file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv"),
  file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled.csv")
)

rev_in_candidates <- c(
  file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv"),
  file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled.csv")
)

first_existing <- function(x) {
  y <- x[file.exists(x)]
  if (length(y) == 0) return(NA_character_)
  y[1]
}

orig_in <- first_existing(orig_in_candidates)
rev_in  <- first_existing(rev_in_candidates)

if (is.na(orig_in) || !file.exists(orig_in)) stop("No suitable original matrix found.")
if (is.na(rev_in)  || !file.exists(rev_in))  stop("No suitable revised matrix found.")

orig_out <- file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusRetimed.csv")
rev_out  <- file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusRetimed.csv")

message("Using input matrices:")
message(" - original: ", orig_in)
message(" - revised:  ", rev_in)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
plot_parameters <- c(
  "prob_hosp",
  "delay_hosp",
  "prob_unsafe_funeral_comm",
  "prob_unsafe_funeral_hosp",
  "prop_etu",
  "ipc_helper"
)

good_response_params <- c("prob_hosp", "prop_etu", "ipc_helper")
adverse_params <- c("delay_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp")

panel_lookup <- c(
  "prob_hosp" = "Probability of hospitalisation",
  "delay_hosp" = "Delay to hospitalisation",
  "prob_unsafe_funeral_comm" = "Unsafe funeral after community death",
  "prob_unsafe_funeral_hosp" = "Unsafe funeral after hospital death",
  "prop_etu" = "Proportion in ETU / ETC",
  "ipc_helper" = "Latent IPC / PPE index"
)
panel_levels <- unname(panel_lookup[plot_parameters])

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

# ------------------------------------------------------------
# Detect old WA conflict timing
# ------------------------------------------------------------
run_lengths <- function(flag) {
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  tibble(value = r$values, start_i = starts, end_i = ends)
}

moving_average <- function(x, k = 21) {
  if (length(x) < k) return(x)
  y <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  y[is.na(y)] <- x[is.na(y)]
  y
}

detect_wa_conflict_window <- function(df, methodology_label) {
  base <- df %>%
    filter(scenario_key == "worst_west_africa") %>%
    arrange(relative_day)

  conflict <- df %>%
    filter(scenario_key == "worst_west_africa_conflict") %>%
    arrange(relative_day)

  if (nrow(base) == 0 || nrow(conflict) == 0) {
    stop("Missing base or conflict West Africa scenario for ", methodology_label)
  }

  joined <- base %>%
    select(relative_day, q_base = q_value, all_of(plot_parameters)) %>%
    rename_with(~ paste0(.x, "_base"), all_of(plot_parameters)) %>%
    inner_join(
      conflict %>%
        select(relative_day, q_conflict = q_value, all_of(plot_parameters)) %>%
        rename_with(~ paste0(.x, "_conflict"), all_of(plot_parameters)),
      by = "relative_day"
    )

  joined <- joined %>% mutate(q_drop = q_base - q_conflict)

  search <- joined %>% filter(relative_day >= 250, relative_day <= 650)

  if (nrow(search) == 0) {
    warning("Could not search WA conflict window; using fallback 350-500 for ", methodology_label)
    return(c(350, 500))
  }

  signal <- search$q_drop

  if (all(!is.finite(signal)) || max(signal, na.rm = TRUE) <= 0) {
    comp <- rep(0, nrow(search))
    for (p in good_response_params) {
      b <- search[[paste0(p, "_base")]]
      c <- search[[paste0(p, "_conflict")]]
      rng <- max(b, c, na.rm = TRUE) - min(b, c, na.rm = TRUE)
      if (is.finite(rng) && rng > 0) comp <- comp + pmax(0, (b - c) / rng)
    }
    for (p in adverse_params) {
      b <- search[[paste0(p, "_base")]]
      c <- search[[paste0(p, "_conflict")]]
      rng <- max(b, c, na.rm = TRUE) - min(b, c, na.rm = TRUE)
      if (is.finite(rng) && rng > 0) comp <- comp + pmax(0, (c - b) / rng)
    }
    signal <- comp
  }

  signal_smooth <- moving_average(signal, k = 21)
  peak <- max(signal_smooth, na.rm = TRUE)

  if (!is.finite(peak) || peak <= 0) {
    warning("Could not detect old WA conflict window; using fallback 350-500 for ", methodology_label)
    return(c(350, 500))
  }

  threshold <- max(0.02, 0.30 * peak)
  flag <- signal_smooth >= threshold

  runs <- run_lengths(flag) %>% filter(value)

  if (nrow(runs) == 0) {
    warning("No contiguous WA conflict window detected; using fallback 350-500 for ", methodology_label)
    return(c(350, 500))
  }

  runs <- runs %>%
    rowwise() %>%
    mutate(area = sum(signal_smooth[start_i:end_i], na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(area))

  chosen <- runs[1, ]

  start_day <- search$relative_day[chosen$start_i]
  end_day   <- search$relative_day[chosen$end_i]

  if ((end_day - start_day) < 60) {
    mid <- round((start_day + end_day) / 2)
    start_day <- mid - 50
    end_day <- mid + 50
  }

  start_day <- max(0, round(start_day))
  end_day <- min(730, round(end_day))

  c(start_day, end_day)
}

# ------------------------------------------------------------
# Patch WA conflict++
# ------------------------------------------------------------
poor_endpoint_values <- function(base_df) {
  out <- list()

  for (p in good_response_params) {
    if (p %in% names(base_df)) {
      out[[p]] <- min(base_df[[p]], na.rm = TRUE)
    }
  }

  for (p in adverse_params) {
    if (p %in% names(base_df)) {
      out[[p]] <- max(base_df[[p]], na.rm = TRUE)
    }
  }

  out
}

patch_one_methodology <- function(df, methodology_label) {
  window <- detect_wa_conflict_window(df, methodology_label)

  base <- df %>%
    filter(scenario_key == "worst_west_africa") %>%
    arrange(relative_day)

  existing_target <- df %>%
    filter(scenario_key == "worst_west_africa_conflict_plusplus") %>%
    arrange(relative_day)

  if (nrow(base) == 0 || nrow(existing_target) == 0) {
    stop("Missing base or existing WA++ scenario for ", methodology_label)
  }

  if (!all(base$relative_day == existing_target$relative_day)) {
    stop("Base and existing WA++ day grids differ for ", methodology_label)
  }

  # Start from ordinary WA, not WA-with-conflict. This removes the old
  # mild conflict dip outside the severe collapse window.
  patched <- base

  patched$scenario_key <- existing_target$scenario_key
  if ("scenario" %in% names(patched) && "scenario" %in% names(existing_target)) {
    patched$scenario <- existing_target$scenario
  }
  if ("methodology" %in% names(patched)) patched$methodology <- methodology_label

  collapse_idx <- patched$relative_day >= window[1] & patched$relative_day <= window[2]

  poor <- poor_endpoint_values(base)

  for (p in names(poor)) {
    patched[[p]][collapse_idx] <- poor[[p]]
  }

  if ("q_value" %in% names(patched)) {
    patched$q_value[collapse_idx] <- 0
  }

  out <- df %>%
    filter(scenario_key != "worst_west_africa_conflict_plusplus") %>%
    bind_rows(patched) %>%
    arrange(methodology, scenario_key, relative_day)

  diag <- patched %>%
    mutate(in_collapse = collapse_idx) %>%
    group_by(in_collapse) %>%
    summarise(
      min_day = min(relative_day, na.rm = TRUE),
      max_day = max(relative_day, na.rm = TRUE),
      mean_prob_hosp = mean(prob_hosp, na.rm = TRUE),
      mean_delay_hosp = mean(delay_hosp, na.rm = TRUE),
      mean_ufc = mean(prob_unsafe_funeral_comm, na.rm = TRUE),
      mean_ufh = mean(prob_unsafe_funeral_hosp, na.rm = TRUE),
      mean_prop_etu = mean(prop_etu, na.rm = TRUE),
      mean_ipc = mean(ipc_helper, na.rm = TRUE),
      mean_q = mean(q_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      methodology = methodology_label,
      detected_start = window[1],
      detected_end = window[2],
      .before = 1
    )

  endpoints <- tibble(
    methodology = methodology_label,
    detected_start = window[1],
    detected_end = window[2],
    parameter = names(poor),
    poor_endpoint = unlist(poor)
  )

  attr(out, "detected_window") <- tibble(
    methodology = methodology_label,
    detected_start = window[1],
    detected_end = window[2]
  )
  attr(out, "diagnostic") <- diag
  attr(out, "poor_endpoints") <- endpoints

  out
}

orig <- readr::read_csv(orig_in, show_col_types = FALSE)
rev  <- readr::read_csv(rev_in, show_col_types = FALSE)

orig_retimed <- patch_one_methodology(orig, "original")
rev_retimed  <- patch_one_methodology(rev, "revised")

readr::write_csv(orig_retimed, orig_out)
readr::write_csv(rev_retimed, rev_out)

detected_windows <- bind_rows(
  attr(orig_retimed, "detected_window"),
  attr(rev_retimed, "detected_window")
)

poor_endpoints <- bind_rows(
  attr(orig_retimed, "poor_endpoints"),
  attr(rev_retimed, "poor_endpoints")
)

diagnostics <- bind_rows(
  attr(orig_retimed, "diagnostic"),
  attr(rev_retimed, "diagnostic")
)

message("\nDetected WA conflict timing from old WA-with-conflict scenario:")
print(detected_windows)

message("\nPoor-response endpoints imposed in WA conflict++ collapse window:")
print(poor_endpoints)

message("\nWA conflict++ diagnostics after retiming:")
print(diagnostics)

message("\nWrote retimed final matrices:")
message(" - ", orig_out)
message(" - ", rev_out)

combined_check <- bind_rows(orig_retimed, rev_retimed) %>%
  group_by(methodology, scenario_key) %>%
  summarise(
    n_rows = n(),
    n_days = n_distinct(relative_day),
    min_day = min(relative_day, na.rm = TRUE),
    max_day = max(relative_day, na.rm = TRUE),
    min_q = min(q_value, na.rm = TRUE),
    max_q = max(q_value, na.rm = TRUE),
    n_missing_q = sum(is.na(q_value)),
    .groups = "drop"
  )

message("\nFull 14-scenario sanity check:")
print(combined_check)

bad <- combined_check %>%
  filter(n_rows != 731 | n_days != 731 | min_day != 0 | max_day != 730 |
           n_missing_q != 0 | min_q < -1e-8 | max_q > 1 + 1e-8)

if (nrow(bad) > 0) {
  print(bad)
  stop("Sanity check failed.")
}

# ------------------------------------------------------------
# Plot only the two updated WA conflict++ scenarios
# ------------------------------------------------------------
audit <- if (file.exists(audit_path)) readr::read_csv(audit_path, show_col_types = FALSE) else tibble()
all_retimed <- bind_rows(orig_retimed, rev_retimed)

get_audit_row <- function(method, scen_key) {
  audit %>% filter(methodology == method, scenario_key == scen_key) %>% slice(1)
}

resolve_anchor_path <- function(method, scen_key) {
  base_row <- get_audit_row(method, "worst_west_africa")
  p <- if (nrow(base_row) > 0 && "anchor_path" %in% names(base_row)) clean_path(base_row$anchor_path[[1]]) else NA_character_
  if (!is.na(p) && file.exists(p)) return(p)

  row <- get_audit_row(method, scen_key)
  p2 <- if (nrow(row) > 0 && "anchor_path" %in% names(row)) clean_path(row$anchor_path[[1]]) else NA_character_
  if (!is.na(p2) && file.exists(p2)) return(p2)

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

read_curve_ribbon <- function(method, final_long_for_scenario) {
  base_row <- get_audit_row(method, "worst_west_africa")
  path <- if (nrow(base_row) > 0 && "curve_path" %in% names(base_row)) clean_path(base_row$curve_path[[1]]) else NA_character_
  x <- safe_read_csv(path)

  fallback <- final_long_for_scenario %>%
    transmute(parameter, panel_title, relative_day, mean, q5 = mean, q95 = mean)

  if (nrow(x) == 0 || !"relative_day" %in% names(x)) return(fallback)

  out <- tibble()

  if ("parameter" %in% names(x)) {
    mean_col <- dplyr::case_when(
      "mean" %in% names(x) ~ "mean",
      "median" %in% names(x) ~ "median",
      "value" %in% names(x) ~ "value",
      TRUE ~ NA_character_
    )

    if (!is.na(mean_col)) {
      q5_col <- dplyr::case_when("q5" %in% names(x) ~ "q5", "lower" %in% names(x) ~ "lower", TRUE ~ NA_character_)
      q95_col <- dplyr::case_when("q95" %in% names(x) ~ "q95", "upper" %in% names(x) ~ "upper", TRUE ~ NA_character_)

      out <- x %>%
        transmute(
          parameter = map_parameter(parameter),
          relative_day = as.numeric(relative_day),
          mean = as.numeric(.data[[mean_col]]),
          q5 = if (!is.na(q5_col)) as.numeric(.data[[q5_col]]) else as.numeric(.data[[mean_col]]),
          q95 = if (!is.na(q95_col)) as.numeric(.data[[q95_col]]) else as.numeric(.data[[mean_col]])
        )
    }
  }

  if (nrow(out) == 0) return(fallback)

  out <- out %>%
    filter(parameter %in% plot_parameters, is.finite(relative_day),
           is.finite(mean), is.finite(q5), is.finite(q95)) %>%
    mutate(q_low = pmin(q5, q95), q_high = pmax(q5, q95)) %>%
    group_by(parameter, relative_day) %>%
    summarise(
      mean = median(mean, na.rm = TRUE),
      q5 = median(q_low, na.rm = TRUE),
      q95 = median(q_high, na.rm = TRUE),
      .groups = "drop"
    ) %>%
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

  ribbon <- bind_rows(out, end_ext) %>%
    arrange(parameter, relative_day)

  method_window <- detected_windows %>% filter(methodology == method) %>% slice(1)
  if (nrow(method_window) > 0) {
    collapse_start <- method_window$detected_start
    collapse_end <- method_window$detected_end

    final_collapse <- fallback %>%
      filter(relative_day >= collapse_start, relative_day <= collapse_end) %>%
      select(parameter, panel_title, relative_day, mean, q5, q95)

    ribbon <- ribbon %>%
      filter(!(relative_day >= collapse_start & relative_day <= collapse_end)) %>%
      bind_rows(final_collapse) %>%
      arrange(parameter, relative_day)
  }

  ribbon
}

plot_one_wa_plusplus <- function(method) {
  scen_key <- "worst_west_africa_conflict_plusplus"

  scen_df <- all_retimed %>%
    filter(methodology == method, scenario_key == scen_key)

  df_long <- scen_df %>%
    select(methodology, scenario_key, scenario, relative_day, all_of(plot_parameters)) %>%
    pivot_longer(cols = all_of(plot_parameters), names_to = "parameter", values_to = "mean") %>%
    mutate(
      parameter = as.character(parameter),
      relative_day = as.numeric(relative_day),
      mean = as.numeric(mean),
      panel_title = as_panel(parameter)
    ) %>%
    arrange(parameter, relative_day)

  ribbon <- read_curve_ribbon(method, df_long)
  anchors <- read_anchors(resolve_anchor_path(method, scen_key)) %>% filter(relative_day <= 730)

  ipc_axis_support <- tibble(
    panel_title = factor("Latent IPC / PPE index", levels = panel_levels),
    relative_day = 0,
    y = c(0, 1)
  )

  window <- detected_windows %>% filter(methodology == method) %>% slice(1)
  subtitle <- paste0(
    "WA conflict++ retimed to existing WA conflict window: day ",
    window$detected_start, "-", window$detected_end,
    "; black line from retimed final matrix"
  )

  scenario_label <- if ("scenario" %in% names(scen_df)) unique(scen_df$scenario)[1] else scen_key

  ggplot() +
    geom_blank(data = ipc_axis_support, aes(x = relative_day, y = y), inherit.aes = FALSE) +
    geom_ribbon(data = ribbon, aes(x = relative_day, ymin = q5, ymax = q95, group = parameter),
                alpha = 0.20, fill = "lightblue", inherit.aes = FALSE) +
    geom_line(data = df_long, aes(x = relative_day, y = mean, group = parameter),
              linewidth = 0.9, colour = "black", inherit.aes = FALSE) +
    geom_errorbar(data = anchors, aes(x = relative_day, ymin = value_low, ymax = value_high),
                  inherit.aes = FALSE, width = 0, colour = "orange3", linewidth = 0.55) +
    geom_point(data = anchors, aes(x = relative_day, y = value_used),
               inherit.aes = FALSE, colour = "orange3", size = 2.1) +
    geom_text(data = anchors, aes(x = relative_day, y = value_used, label = anchor_id),
              inherit.aes = FALSE, colour = "orange4", size = 2.4,
              hjust = -0.05, vjust = -0.6, check_overlap = TRUE) +
    facet_wrap(~panel_title, scales = "free_y", ncol = 2) +
    labs(
      title = paste0(stringr::str_to_title(method), " methodology: ", scenario_label),
      subtitle = subtitle,
      x = "Relative day",
      y = "Parameter value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p_original_WA_conflict_plusplus_retimed <- plot_one_wa_plusplus("original")
p_revised_WA_conflict_plusplus_retimed <- plot_one_wa_plusplus("revised")

p_original_WA_conflict_plusplus_retimed
p_revised_WA_conflict_plusplus_retimed

message("\nPlots stored in:")
message(" - p_original_WA_conflict_plusplus_retimed")
message(" - p_revised_WA_conflict_plusplus_retimed")
message("\nFinal retimed matrices:")
message(" - ", orig_out)
message(" - ", rev_out)
