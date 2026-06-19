# =============================================================================
# helper_functions_figure4.R
#
# Post-hoc helper functions for Figure 4 panels a, b, c.
# All panels are derived from baseline (no antiviral) simulations.
#
# Policy definitions:
#   Policy B : antiviral given only to HCWs who were actually infected
#              (PPE failed). Doses consumed = n_infected_hcw.
#   Policy A : antiviral also given to HCWs whom PPE protected, i.e. the
#              full exposure pool is estimated as
#              n_infected_hcw / (1 - ppe_efficacy).
#              Doses consumed = n_policy_A = n_B / (1 - ppe_efficacy).
#
# DPC-efficacy lookup:
#   curve_d50_dat  <- readRDS(here("data-processed",
#                                  "DPC_fixed_efficacy_varied_d50.rds"))
#   Columns used   : dpc, efficacy (mid), eighty_efficacy_lo, eighty_efficacy_hi
#
# Panel a : stockpile × HCW deaths averted  (DPC 0 and DPC 5, Policy A and B)
# Panel b : supply/demand × % HCW deaths averted  (DPC 0, Policy A and B)
# Panel c : intrinsic efficacy × doses per death averted
#           (DPC 0 vs 5, Policy A vs B)
# =============================================================================

library(dplyr)
library(here)

# =============================================================================
# get_efficacy_at_dpc
#
# Read mid / lo / hi efficacy off the DPC-efficacy lookup curve.
# =============================================================================
get_efficacy_at_dpc <- function(dpc_val,
                                col = "efficacy",
                                curve_dat = NULL) {
  if (is.null(curve_dat))
    curve_dat <- readRDS(here("data-processed",
                              "DPC_fixed_efficacy_varied_d50.rds"))
  approx(x = curve_dat$dpc, y = curve_dat[[col]],
         xout = dpc_val, rule = 2)$y
}

# =============================================================================
# extract_figure4_posthoc
#
# Read all baseline RDS files for one scenario and return a per-run data frame
# with the information needed for all three panels.
#
# Returns one row per (particle_id, rep) with:
#   scenario, particle_id, rep,
#   ppe_efficacy,               # posterior draw for this particle
#   n_hcw_infected,             # Policy B eligible pool (actual infections)
#   n_hcw_exposed_A,            # Policy A eligible pool = n_B / (1-ppe_eff)
#   hcw_infection_times,        # list-col: sorted infection times (Policy B order)
#   hcw_death_flags,            # list-col: logical, did each HCW die (baseline)?
#   n_hcw_deaths_baseline       # total HCW deaths without antiviral
# =============================================================================
extract_figure4_posthoc <- function(sc_name     = "WestAfrica",
                                    n_workers   = 10L,
                                    base_dir    = here("outputs", "simulation",
                                                       "figure4_baseline")) {
  library(future)
  library(future.apply)
  
  dir_path <- file.path(base_dir, sc_name)
  # Match only simulation output files (WestAfrica_pNNN_rNN.rds pattern)
  files    <- list.files(dir_path, pattern = "_p\\d+_r\\d+\\.rds$",
                         full.names = TRUE)
  if (length(files) == 0)
    stop("No RDS files found in: ", dir_path)
  
  # Load ppe_efficacy lookup (particle_id -> ppe_efficacy)
  lookup_path <- file.path(base_dir, sc_name,
                           sprintf("ppe_efficacy_lookup_%s.rds", sc_name))
  ppe_lookup  <- readRDS(lookup_path)
  
  message(sprintf("Extracting post-hoc data from %d baseline files (%s)...",
                  length(files), sc_name))
  
  plan(multisession, workers = min(n_workers, future::availableCores()))
  on.exit(plan(sequential), add = TRUE)
  
  results <- future_lapply(seq_along(files), function(i) {
    f <- files[[i]]
    tryCatch({
      x <- readRDS(f)
      
      # Parse particle_id and rep from filename
      fname <- tools::file_path_sans_ext(basename(f))
      parts <- regmatches(fname,
                          regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
      if (length(parts) != 4) return(NULL)
      sc  <- parts[2]
      pid <- as.integer(parts[3])
      rep <- as.integer(parts[4])
      
      tdf <- x$tdf
      
      # Filter to actual HCW infections only
      is_hcw      <- !is.na(tdf$class) & tdf$class == "HCW"
      is_infected <- !is.na(tdf$time_infection_absolute)
      hcw_inf     <- tdf[is_hcw & is_infected, ]
      
      # Sort by infection time for first-come-first-serve stockpile allocation
      hcw_inf <- hcw_inf[order(hcw_inf$time_infection_absolute), ]
      
      n_hcw_infected       <- nrow(hcw_inf)
      n_hcw_deaths_baseline <- sum(!is.na(hcw_inf$outcome) & hcw_inf$outcome)
      
      # Retrieve particle-level ppe_efficacy
      ppe_eff <- ppe_lookup$ppe_efficacy[ppe_lookup$particle_id == pid]
      if (length(ppe_eff) == 0) ppe_eff <- NA_real_
      
      # Policy A eligible pool: back-calculate the full exposure pool
      # including HCWs whom PPE protected from infection
      n_hcw_exposed_A <- if (!is.na(ppe_eff) && ppe_eff < 1)
        n_hcw_infected / (1 - ppe_eff)
      else
        n_hcw_infected   # fallback if ppe_efficacy == 1
      
      data.frame(
        scenario              = sc_name,
        particle_id           = pid,
        rep                   = rep,
        ppe_efficacy          = ppe_eff,
        n_hcw_infected        = n_hcw_infected,
        n_hcw_exposed_A       = n_hcw_exposed_A,
        n_hcw_deaths_baseline = n_hcw_deaths_baseline,
        # Store as comma-separated strings to avoid list-column I/O issues
        hcw_infection_times   = paste(hcw_inf$time_infection_absolute,
                                      collapse = ","),
        hcw_died              = paste(as.integer(!is.na(hcw_inf$outcome) &
                                                   hcw_inf$outcome),
                                      collapse = ","),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message(sprintf("ERROR in %s: %s", basename(f), conditionMessage(e)))
      NULL
    })
  }, future.packages = c("here"), future.seed = TRUE)
  
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) stop("No valid results extracted.")
  do.call(rbind, results)
}

# =============================================================================
# apply_stockpile_posthoc
#
# Given the per-run data frame from extract_figure4_posthoc(), sweep over
# stockpile sizes and compute HCW deaths averted under Policy A and Policy B.
#
# Arguments:
#   run_df        : output of extract_figure4_posthoc()
#   stockpile_seq : integer vector of stockpile sizes to evaluate
#   efficacy      : scalar antiviral efficacy (Bernoulli probability of
#                   preventing death given drug administered)
#   seed          : RNG seed for Bernoulli sampling
#
# Returns a data frame with columns:
#   scenario, particle_id, rep, policy, stockpile,
#   doses_used, deaths_averted, deaths_baseline
# =============================================================================
apply_stockpile_posthoc <- function(run_df,
                                    stockpile_seq    = seq(1000, 100000, by = 1000),
                                    efficacy         = NULL,
                                    dpc              = 0,
                                    curve_dat        = NULL,
                                    seed             = 42L,
                                    DOSES_PER_COURSE = 20L) {
  # stockpile_seq is in DOSES (not courses).
  # Internally, n_inf and n_exp_A are in courses (persons), so we convert S
  # to courses before comparing: S_courses = S / DOSES_PER_COURSE.
  if (is.null(curve_dat))
    curve_dat <- readRDS(here("data-processed",
                              "DPC_fixed_efficacy_varied_d50.rds"))
  if (is.null(efficacy))
    efficacy <- get_efficacy_at_dpc(dpc, col = "efficacy",
                                    curve_dat = curve_dat)
  
  set.seed(seed)
  
  do.call(rbind, lapply(seq_len(nrow(run_df)), function(i) {
    row <- run_df[i, ]
    
    # Reconstruct per-HCW vectors from stored strings
    died_vec <- as.integer(strsplit(row$hcw_died, ",")[[1]])
    n_inf    <- row$n_hcw_infected       # persons (courses)
    n_exp_A  <- row$n_hcw_exposed_A      # persons (courses), may be fractional
    
    do.call(rbind, lapply(stockpile_seq, function(S_doses) {
      
      # Convert stockpile from doses to courses for person-level comparison
      S_courses <- S_doses / DOSES_PER_COURSE
      
      # --- Policy B ---
      # Courses go to actually infected HCWs in arrival order.
      n_treated_B    <- min(S_courses, n_inf)
      died_treated_B <- died_vec[seq_len(round(n_treated_B))]
      averted_B      <- sum(died_treated_B) * efficacy  # expected value
      
      # --- Policy A ---
      # Courses distributed across full exposure pool (infections + PPE-protected).
      # Coverage fraction for actual infections = min(S, n_exp_A) / n_exp_A.
      n_exp_A_round  <- max(round(n_exp_A), 1L)
      coverage_A     <- min(S_courses, n_exp_A_round) / n_exp_A_round
      n_treated_A    <- round(n_inf * coverage_A)
      died_treated_A <- died_vec[seq_len(min(n_treated_A, n_inf))]
      averted_A      <- sum(died_treated_A) * efficacy  # expected value
      
      rbind(
        data.frame(scenario         = row$scenario,
                   particle_id      = row$particle_id,
                   rep              = row$rep,
                   policy           = "B",
                   dpc              = dpc,
                   stockpile_doses  = S_doses,
                   doses_used       = n_treated_B * DOSES_PER_COURSE,
                   deaths_averted   = averted_B,
                   deaths_baseline  = row$n_hcw_deaths_baseline,
                   stringsAsFactors = FALSE),
        data.frame(scenario         = row$scenario,
                   particle_id      = row$particle_id,
                   rep              = row$rep,
                   policy           = "A",
                   dpc              = dpc,
                   stockpile_doses  = S_doses,
                   doses_used       = min(S_courses, n_exp_A_round) * DOSES_PER_COURSE,
                   deaths_averted   = averted_A,
                   deaths_baseline  = row$n_hcw_deaths_baseline,
                   stringsAsFactors = FALSE)
      )
    }))
  }))
}

# =============================================================================
# summarise_stockpile_panel
#
# Aggregate apply_stockpile_posthoc() output to particle-level medians and IQR
# for plotting panel a and panel b.
#
# Panel a : y = deaths_averted (absolute),   x = stockpile
# Panel b : y = pct_deaths_averted,          x = supply_demand_ratio
#             supply_demand_ratio = stockpile / n_demand
#             where n_demand = median n_hcw_infected (Policy B)
#                            = median n_hcw_exposed_A (Policy A)
# =============================================================================
summarise_stockpile_panel <- function(stockpile_df, run_df,
                                      DOSES_PER_COURSE = 20L) {
  
  # Compute demand per policy per run (for panel b x-axis)
  demand_df <- run_df %>%
    group_by(scenario, particle_id) %>%
    summarise(
      demand_B = mean(n_hcw_infected),     # avg across reps
      demand_A = mean(n_hcw_exposed_A),
      .groups = "drop"
    )
  
  # Average across reps within each particle x policy x stockpile
  per_particle <- stockpile_df %>%
    group_by(scenario, particle_id, policy, dpc, stockpile_doses) %>%
    summarise(
      deaths_averted  = mean(deaths_averted),
      deaths_baseline = mean(deaths_baseline),
      doses_used      = mean(doses_used),
      .groups = "drop"
    ) %>%
    mutate(pct_averted = 100 * deaths_averted / pmax(deaths_baseline, 1))
  
  # Join demand; x-axis for panel b is always Policy B demand in doses
  # so that both policies share the same x scale
  per_particle <- per_particle %>%
    left_join(demand_df, by = c("scenario", "particle_id")) %>%
    mutate(
      demand_B_doses = demand_B * DOSES_PER_COURSE,
      supply_ratio   = stockpile_doses / pmax(demand_B_doses, 1)
    )
  
  # Summarise across particles: median + IQR
  per_particle %>%
    group_by(scenario, policy, dpc, stockpile_doses) %>%
    summarise(
      # Panel a quantities
      deaths_averted_med = median(deaths_averted),
      deaths_averted_lo  = quantile(deaths_averted, 0.25),
      deaths_averted_hi  = quantile(deaths_averted, 0.75),
      # Panel b quantities
      pct_averted_med    = median(pct_averted),
      pct_averted_lo     = quantile(pct_averted, 0.25),
      pct_averted_hi     = quantile(pct_averted, 0.75),
      # Supply/demand ratio (median across particles for x-axis placement)
      supply_ratio_med   = median(supply_ratio),
      .groups = "drop"
    )
}

# =============================================================================
# compute_doses_per_death
#
# Panel c: for each combination of (intrinsic_efficacy_scale, dpc, policy),
# compute doses per death averted assuming unlimited stockpile
# (all eligible HCWs receive the drug).
#
# Arguments:
#   run_df         : output of extract_figure4_posthoc()
#   efficacy_scales: numeric vector in [0, 1]; the max efficacy at DPC=0 is
#                    scaled by this factor before being applied
#   dpc_vals       : numeric vector of DPC values (e.g. c(0, 5))
#   curve_dat      : DPC-efficacy lookup table
#   seed           : RNG seed
#
# Returns a data frame with columns:
#   scenario, policy, dpc, efficacy_scale, intrinsic_efficacy,
#   doses_per_death_med, doses_per_death_lo, doses_per_death_hi
# =============================================================================
compute_doses_per_death <- function(run_df,
                                    efficacy_scales  = seq(0.2, 0.9, by = 0.1),
                                    dpc_vals         = c(0, 5),
                                    curve_dat        = NULL,
                                    seed             = 42L,
                                    DOSES_PER_COURSE = 20L) {
  if (is.null(curve_dat))
    curve_dat <- readRDS(here("data-processed",
                              "DPC_fixed_efficacy_varied_d50.rds"))
  
  # Max efficacy at DPC = 0 (the reference peak from the curve)
  max_eff_dpc0 <- get_efficacy_at_dpc(0, col = "efficacy",
                                      curve_dat = curve_dat)
  
  set.seed(seed)
  
  results <- do.call(rbind, lapply(dpc_vals, function(dpc) {
    # Base efficacy at this DPC (mid)
    eff_base_dpc <- get_efficacy_at_dpc(dpc, col = "efficacy",
                                        curve_dat = curve_dat)
    
    do.call(rbind, lapply(efficacy_scales, function(scale) {
      # Scaled intrinsic efficacy: adjust the DPC-0 peak by scale,
      # then apply the same relative DPC decay
      intrinsic_eff <- max_eff_dpc0 * scale
      # Efficacy at this DPC preserves the ratio eff_base_dpc / max_eff_dpc0
      eff_at_dpc    <- intrinsic_eff * (eff_base_dpc / max_eff_dpc0)
      
      do.call(rbind, lapply(seq_len(nrow(run_df)), function(i) {
        row <- run_df[i, ]
        
        died_vec <- as.integer(strsplit(row$hcw_died, ",")[[1]])
        n_inf    <- row$n_hcw_infected
        n_exp_A  <- round(row$n_hcw_exposed_A)
        
        # Policy B: treat all actually infected HCWs (unlimited stockpile)
        # averted = expected value: n_deaths_baseline * efficacy
        # (consistent with figure5: averted = n_prevented_hcw * efficacy)
        n_deaths_B <- sum(died_vec)
        averted_B  <- n_deaths_B * eff_at_dpc
        doses_B    <- n_inf * DOSES_PER_COURSE    # convert courses to doses
        
        # Policy A: treat the full exposure pool (including PPE-protected HCWs)
        # averted count is identical (only actual infections contribute to deaths)
        # but doses consumed are higher
        n_extra_A  <- max(n_exp_A - n_inf, 0L)
        averted_A  <- averted_B
        doses_A    <- (n_inf + n_extra_A) * DOSES_PER_COURSE
        
        rbind(
          data.frame(scenario       = row$scenario,
                     particle_id    = row$particle_id,
                     rep            = row$rep,
                     policy         = "B",
                     dpc            = dpc,
                     efficacy_scale = scale,
                     intrinsic_efficacy = intrinsic_eff,
                     doses          = doses_B,
                     deaths_averted = averted_B,
                     stringsAsFactors = FALSE),
          data.frame(scenario       = row$scenario,
                     particle_id    = row$particle_id,
                     rep            = row$rep,
                     policy         = "A",
                     dpc            = dpc,
                     efficacy_scale = scale,
                     intrinsic_efficacy = intrinsic_eff,
                     doses          = doses_A,
                     deaths_averted = averted_A,
                     stringsAsFactors = FALSE)
        )
      }))
    }))
  }))
  
  # Average across reps within particle, then summarise across particles
  results %>%
    group_by(scenario, particle_id, policy, dpc, efficacy_scale,
             intrinsic_efficacy) %>%
    summarise(
      doses          = mean(doses),
      deaths_averted = mean(deaths_averted),
      .groups = "drop"
    ) %>%
    mutate(
      doses_per_death = ifelse(deaths_averted > 0,
                               doses / deaths_averted,
                               NA_real_)
    ) %>%
    group_by(scenario, policy, dpc, efficacy_scale, intrinsic_efficacy) %>%
    summarise(
      doses_per_death_med = median(doses_per_death, na.rm = TRUE),
      doses_per_death_lo  = quantile(doses_per_death, 0.25, na.rm = TRUE),
      doses_per_death_hi  = quantile(doses_per_death, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}