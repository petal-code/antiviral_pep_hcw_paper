# =============================================================================
# helper_functions_figure4.R
# =============================================================================

library(dplyr)
library(here)

DOSES_PER_COURSE_DEFAULT <- 20L

# =============================================================================
# get_efficacy_at_dpc
# =============================================================================
get_efficacy_at_dpc <- function(dpc_val, col = "efficacy", curve_dat = NULL) {
  if (is.null(curve_dat))
    curve_dat <- readRDS(here("data-processed",
                              "DPC_fixed_efficacy_varied_d50.rds"))
  approx(x = curve_dat$dpc, y = curve_dat[[col]],
         xout = dpc_val, rule = 2)$y
}

# =============================================================================
# extract_figure4_posthoc
# =============================================================================
extract_figure4_posthoc <- function(sc_name   = "WestAfrica",
                                    n_workers = 10L,
                                    base_dir  = here("outputs", "simulation",
                                                     "figure4_baseline")) {
  library(future)
  library(future.apply)
  
  dir_path <- file.path(base_dir, sc_name)
  files    <- list.files(dir_path, pattern = "_p\\d+_r\\d+\\.rds$",
                         full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", dir_path)
  
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
      x     <- readRDS(f)
      fname <- tools::file_path_sans_ext(basename(f))
      parts <- regmatches(fname, regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
      if (length(parts) != 4) return(NULL)
      pid <- as.integer(parts[3])
      rep <- as.integer(parts[4])
      
      tdf         <- x$tdf
      is_hcw      <- !is.na(tdf$class) & tdf$class == "HCW"
      is_infected <- !is.na(tdf$time_infection_absolute)
      hcw_inf     <- tdf[is_hcw & is_infected, ]
      hcw_inf     <- hcw_inf[order(hcw_inf$time_infection_absolute), ]
      
      n_hcw_infected        <- nrow(hcw_inf)
      n_hcw_deaths_baseline <- sum(!is.na(hcw_inf$outcome) & hcw_inf$outcome)
      
      ppe_eff <- ppe_lookup$ppe_efficacy[ppe_lookup$particle_id == pid]
      if (length(ppe_eff) == 0) ppe_eff <- NA_real_
      
      n_hcw_exposed_A <- if (!is.na(ppe_eff) && ppe_eff < 1)
        n_hcw_infected / (1 - ppe_eff)
      else
        n_hcw_infected
      
      data.frame(
        scenario              = sc_name,
        particle_id           = pid,
        rep                   = rep,
        ppe_efficacy          = ppe_eff,
        n_hcw_infected        = n_hcw_infected,
        n_hcw_exposed_A       = n_hcw_exposed_A,
        n_hcw_deaths_baseline = n_hcw_deaths_baseline,
        hcw_died = paste(as.integer(!is.na(hcw_inf$outcome) & hcw_inf$outcome),
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
# stockpile_seq in DOSES. Vectorised per-run to avoid stack overflow.
# averted = cumulative deaths among first k treated HCWs * efficacy.
# =============================================================================
apply_stockpile_posthoc <- function(run_df,
                                    stockpile_seq    = seq(1000, 100000, by = 1000),
                                    efficacy         = NULL,
                                    dpc              = 0,
                                    curve_dat        = NULL,
                                    DOSES_PER_COURSE = DOSES_PER_COURSE_DEFAULT) {
  if (is.null(curve_dat))
    curve_dat <- readRDS(here("data-processed",
                              "DPC_fixed_efficacy_varied_d50.rds"))
  if (is.null(efficacy))
    efficacy <- get_efficacy_at_dpc(dpc, col = "efficacy", curve_dat = curve_dat)
  
  run_results <- lapply(seq_len(nrow(run_df)), function(i) {
    row           <- run_df[i, ]
    died_vec      <- as.integer(strsplit(row$hcw_died, ",")[[1]])
    n_inf         <- row$n_hcw_infected
    n_exp_A       <- row$n_hcw_exposed_A
    n_exp_A_round <- max(round(n_exp_A), 1L)
    
    cum_deaths    <- cumsum(died_vec)
    S_courses_vec <- stockpile_seq / DOSES_PER_COURSE
    
    # Policy B
    n_treated_B_vec <- pmin(round(S_courses_vec), n_inf)
    averted_B_vec   <- ifelse(
      n_treated_B_vec == 0L, 0,
      cum_deaths[n_treated_B_vec] * efficacy
    )
    
    # Policy A: proportional coverage across full exposure pool
    coverage_A_vec  <- pmin(S_courses_vec, n_exp_A_round) / n_exp_A_round
    n_treated_A_vec <- pmin(round(n_inf * coverage_A_vec), n_inf)
    averted_A_vec   <- ifelse(
      n_treated_A_vec == 0L, 0,
      cum_deaths[n_treated_A_vec] * efficacy
    )
    
    rbind(
      data.frame(scenario        = row$scenario,
                 particle_id     = row$particle_id,
                 rep             = row$rep,
                 policy          = "B",
                 dpc             = dpc,
                 stockpile_doses = stockpile_seq,
                 doses_used      = n_treated_B_vec * DOSES_PER_COURSE,
                 deaths_averted  = averted_B_vec,
                 deaths_baseline = row$n_hcw_deaths_baseline,
                 n_hcw_infected  = n_inf,
                 n_hcw_exposed_A = n_exp_A,
                 stringsAsFactors = FALSE),
      data.frame(scenario        = row$scenario,
                 particle_id     = row$particle_id,
                 rep             = row$rep,
                 policy          = "A",
                 dpc             = dpc,
                 stockpile_doses = stockpile_seq,
                 doses_used      = pmin(S_courses_vec, n_exp_A_round) * DOSES_PER_COURSE,
                 deaths_averted  = averted_A_vec,
                 deaths_baseline = row$n_hcw_deaths_baseline,
                 n_hcw_infected  = n_inf,
                 n_hcw_exposed_A = n_exp_A,
                 stringsAsFactors = FALSE)
    )
  })
  do.call(rbind, run_results)
}

# =============================================================================
# summarise_stockpile_panel
#
# Returns a list with two data frames:
#   $panel_a : deaths_averted vs stockpile_doses
#   $panel_b : pct_averted vs supply_ratio (interpolated onto common grid)
#
# Panel b key insight: each particle's supply_ratio = stockpile / its own
# demand. We interpolate pct_averted onto a common ratio grid BEFORE taking
# the median. This guarantees Policy B plateaus at exactly supply_ratio = 1
# because at that point every particle has stockpile == its own full demand.
# =============================================================================
summarise_stockpile_panel <- function(stockpile_df, run_df,
                                      DOSES_PER_COURSE = DOSES_PER_COURSE_DEFAULT,
                                      ratio_grid = seq(0, 3.5, by = 0.05)) {
  
  # supply_ratio_p uses each rep's own n_hcw_infected so that at ratio=1
  # every single run has its full demand covered (Policy B fully saturated)
  per_particle <- stockpile_df %>%
    mutate(
      demand_B_doses  = n_hcw_infected * DOSES_PER_COURSE,
      supply_ratio_p  = stockpile_doses / pmax(demand_B_doses, 1),
      pct_averted     = 100 * deaths_averted / pmax(deaths_baseline, 1)
    )
  
  # ----- Panel a: summarise at raw stockpile levels -----
  panel_a <- per_particle %>%
    group_by(scenario, policy, dpc, stockpile_doses) %>%
    summarise(
      deaths_averted_med = median(deaths_averted),
      deaths_averted_lo  = quantile(deaths_averted, 0.25),
      deaths_averted_hi  = quantile(deaths_averted, 0.75),
      .groups = "drop"
    )
  
  # ----- Panel b: bin runs by supply_ratio_p, then take median -----
  # Each run already has supply_ratio_p = stockpile / its own demand.
  # At ratio_p = 1, every run has exactly full demand covered -> Policy B
  # is always saturated there by definition. We simply bin all runs by
  # their ratio_p and summarise — no interpolation needed.
  panel_b <- per_particle %>%
    mutate(supply_ratio = round(supply_ratio_p / 0.05) * 0.05) %>%
    group_by(scenario, policy, dpc, supply_ratio) %>%
    summarise(
      pct_averted_med = median(pct_averted),
      pct_averted_lo  = quantile(pct_averted, 0.25),
      pct_averted_hi  = quantile(pct_averted, 0.75),
      .groups = "drop"
    )
  
  list(panel_a = panel_a, panel_b = panel_b)
}

# =============================================================================
# compute_doses_per_death  (vectorised)
# =============================================================================
compute_doses_per_death <- function(run_df,
                                    efficacy_scales  = seq(0.2, 0.9, by = 0.1),
                                    dpc_vals         = c(0, 5),
                                    curve_dat        = NULL,
                                    DOSES_PER_COURSE = DOSES_PER_COURSE_DEFAULT) {
  if (is.null(curve_dat))
    curve_dat <- readRDS(here("data-processed",
                              "DPC_fixed_efficacy_varied_d50.rds"))
  
  eff_dpc0 <- get_efficacy_at_dpc(0, col = "efficacy", curve_dat = curve_dat)
  
  run_summary <- run_df %>%
    mutate(
      n_extra_A = pmax(round(n_hcw_exposed_A) - n_hcw_infected, 0L),
      doses_B   = n_hcw_infected * DOSES_PER_COURSE,
      doses_A   = (n_hcw_infected + n_extra_A) * DOSES_PER_COURSE
    ) %>%
    select(scenario, particle_id, rep,
           n_deaths_baseline = n_hcw_deaths_baseline, doses_B, doses_A)
  
  grid <- expand.grid(
    row_idx        = seq_len(nrow(run_summary)),
    dpc            = dpc_vals,
    efficacy_scale = efficacy_scales,
    stringsAsFactors = FALSE
  )
  
  eff_dpc_vals        <- sapply(dpc_vals, function(d)
    get_efficacy_at_dpc(d, col = "efficacy", curve_dat = curve_dat))
  names(eff_dpc_vals) <- as.character(dpc_vals)
  
  grid$dpc_decay   <- eff_dpc_vals[as.character(grid$dpc)] / eff_dpc0
  grid$eff_at_dpc  <- grid$efficacy_scale * grid$dpc_decay
  grid$scenario    <- run_summary$scenario[grid$row_idx]
  grid$particle_id <- run_summary$particle_id[grid$row_idx]
  grid$rep         <- run_summary$rep[grid$row_idx]
  grid$n_deaths    <- run_summary$n_deaths_baseline[grid$row_idx]
  grid$doses_B     <- run_summary$doses_B[grid$row_idx]
  grid$doses_A     <- run_summary$doses_A[grid$row_idx]
  grid$averted     <- grid$n_deaths * grid$eff_at_dpc
  
  results <- bind_rows(
    grid %>% transmute(scenario, particle_id, rep,
                       policy = "B", dpc, efficacy_scale,
                       intrinsic_efficacy = efficacy_scale,
                       doses = doses_B, deaths_averted = averted),
    grid %>% transmute(scenario, particle_id, rep,
                       policy = "A", dpc, efficacy_scale,
                       intrinsic_efficacy = efficacy_scale,
                       doses = doses_A, deaths_averted = averted)
  )
  
  results %>%
    group_by(scenario, particle_id, policy, dpc,
             efficacy_scale, intrinsic_efficacy) %>%
    summarise(doses          = mean(doses),
              deaths_averted = mean(deaths_averted),
              .groups = "drop") %>%
    mutate(doses_per_death = ifelse(deaths_averted > 0,
                                    doses / deaths_averted, NA_real_)) %>%
    group_by(scenario, policy, dpc, efficacy_scale, intrinsic_efficacy) %>%
    summarise(
      doses_per_death_med = median(doses_per_death, na.rm = TRUE),
      doses_per_death_lo  = quantile(doses_per_death, 0.25, na.rm = TRUE),
      doses_per_death_hi  = quantile(doses_per_death, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}