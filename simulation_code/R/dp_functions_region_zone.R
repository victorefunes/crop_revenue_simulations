# ============================================================================
# DYNAMIC PROGRAMMING FUNCTIONS FOR REGION-ZONE FRAMEWORK
# ============================================================================
# Solves for the optimal 6-year crop sequence (and optionally insurance
# coverage level) using backward induction over the 64-state C/S history
# space, then plugs the resulting plan6 into the existing credit/cashflow
# simulation pipeline.
#
# Structural model:
#   State      : 6-char rotation window (e.g. "CCCSCS") -> 64 possible states
#   Action     : crop in {Corn, Soy} x theta in Thetas
#   Payoff     : E[revenue] + insurance_net - E[cost] - RCI_adjustment_cost
#   Transition : drop oldest char, append new crop char -> next state
#
# Integration points:
#   rollout_dp_plan()                -> plan6 string for use in any simulator
#   simulate_dp_credit_region_zone() -> feeds plan6 into credit simulator
#   compare_dp_vs_fixed_credit()     -> DP-optimal vs CCCCCC/CSCSCS/CCSCCS
# ============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, tidyverse)

# ============================================================
# SECTION 1: STATE-SPACE HELPERS
# ============================================================

# All 64 possible 6-char C/S rotation histories
all_histories6_dp <- function() {
  v <- c("C", "S")
  as.vector(outer(outer(outer(outer(outer(
    v, v, paste0), v, paste0), v, paste0), v, paste0), v, paste0))
}

# RCI = 1 + number of crop switches in the 6-char window
rci_from_string <- function(str6) {
  v <- strsplit(str6, "")[[1]]
  1L + sum(v[-1L] != v[-length(v)])
}

# Deterministic state transition: slide window, append new crop char
dp_next_state <- function(window6, crop_char) {
  paste0(substr(window6, 2L, 6L), toupper(substr(crop_char, 1L, 1L)))
}

# ============================================================
# SECTION 2: SHORTFALL COMPUTATION FROM RESIDUAL POOL
# ============================================================
# S_theta = E[max(0, theta - R)]  where R = Y / mu
# This is the per-unit expected shortfall feeding into:
#   Expected Indemnity = Price * mu * S_theta

compute_S_theta_gen <- function(gen, state, theta) {
  mu_row <- gen$mu[rot6 == state]
  if (nrow(mu_row) == 0L) return(NA_real_)
  mu_val <- mu_row$mu[1L]
  if (!is.finite(mu_val) || mu_val <= 0) return(NA_real_)

  pool_row <- gen$pools[rot6 == state]
  if (nrow(pool_row) == 0L) return(NA_real_)
  pool <- pool_row$resid_pool[[1L]]
  if (length(pool) == 0L) return(NA_real_)

  R <- (mu_val + pool) / mu_val     # relative yield: 1 + resid/mu
  mean(pmax(0, theta - R), na.rm = TRUE)
}

# ============================================================
# SECTION 3: BUILD PRECOMPUTED MU / S MAPS
# ============================================================
# Corn generator has entries only for states ending in "C" (32 of 64).
# Soy generator has entries only for states ending in "S" (32 of 64).
# This matches the simulation convention: the generator is called with the
# POST-planting window6 (which ends in the current year's crop char).
# All DP lookups therefore use next6 = dp_next_state(window6, crop).

build_dp_maps_region_zone <- function(geo_region,
                                       prod_zone,
                                       gens_by_region_zone,
                                       thetas  = c(0.65, 0.75, 0.85),
                                       verbose = FALSE) {
  corn_gen <- gens_by_region_zone$corn[[geo_region]][[prod_zone]]
  soy_gen  <- gens_by_region_zone$soy[[geo_region]][[prod_zone]]
  if (is.null(corn_gen) || is.null(soy_gen))
    stop("Generators not available for: ", geo_region, "-", prod_zone)

  all_states  <- all_histories6_dp()
  corn_states <- all_states[endsWith(all_states, "C")]
  soy_states  <- all_states[endsWith(all_states, "S")]

  if (verbose) cat("Building DP maps for", geo_region, "-", prod_zone, "...\n")

  # --- mu vectors (indexed by post-planting state) ---
  extract_mu <- function(gen_obj, st_vec) {
    setNames(vapply(st_vec, function(st) {
      v <- gen_obj$mu[rot6 == st, mu]
      if (length(v) == 0L || !is.finite(v[1L])) NA_real_ else v[1L]
    }, numeric(1L)), st_vec)
  }

  # Pad to all 64 states (non-matching entries stay NA)
  mu_corn <- setNames(rep(NA_real_, length(all_states)), all_states)
  mu_corn[corn_states] <- extract_mu(corn_gen, corn_states)
  mu_soy  <- setNames(rep(NA_real_, length(all_states)), all_states)
  mu_soy[soy_states]   <- extract_mu(soy_gen,  soy_states)

  # --- S(theta) vectors ---
  make_S_vec <- function(gen_obj, st_vec, all64, th) {
    S_raw <- setNames(vapply(st_vec, function(st)
      compute_S_theta_gen(gen_obj, st, th), numeric(1L)), st_vec)
    out <- setNames(rep(NA_real_, length(all64)), all64)
    out[st_vec] <- S_raw
    out
  }

  S_corn <- setNames(
    lapply(thetas, function(th) make_S_vec(corn_gen, corn_states, all_states, th)),
    as.character(thetas))
  S_soy  <- setNames(
    lapply(thetas, function(th) make_S_vec(soy_gen,  soy_states,  all_states, th)),
    as.character(thetas))

  if (verbose) {
    cat(sprintf("  Corn: %d/32 states with finite mu | Soy: %d/32\n",
                sum(is.finite(mu_corn)), sum(is.finite(mu_soy))))
  }

  list(mu = list(Corn = mu_corn, Soy = mu_soy),
       S  = list(Corn = S_corn,  Soy = S_soy),
       geo_region = geo_region, prod_zone = prod_zone)
}

# ============================================================
# SECTION 4: ACTUARIALLY FAIR K ESTIMATION
# ============================================================
# When external K_table (from pairs_rci$K_i) is unavailable:
#   Fair: EI = Premium  =>  P * mu * S = K * mu^(1+beta)
#   =>  K = P * S / mu^beta   (median across states per crop x theta)

estimate_K_from_maps <- function(dp_maps, price_pars,
                                  thetas = c(0.65, 0.75, 0.85),
                                  beta   = -0.5) {
  K_out <- list()
  for (th in thetas) {
    th_key <- as.character(th)
    for (crop in c("Corn", "Soy")) {
      P_pars <- price_pars[[if (crop == "Corn") "C" else "S"]]
      P_exp  <- exp(P_pars$mu + 0.5 * P_pars$sigma^2)
      mu_vec <- dp_maps$mu[[crop]]
      S_vec  <- dp_maps$S[[crop]][[th_key]]
      ok     <- is.finite(mu_vec) & is.finite(S_vec) & mu_vec > 0
      K_vals <- P_exp * S_vec[ok] / (mu_vec[ok]^beta)
      K_out[[paste(crop, th_key, sep = "_")]] <- median(K_vals, na.rm = TRUE)
    }
  }
  K_out
}

# ============================================================
# SECTION 5: PER-PERIOD EXPECTED PAYOFF
# ============================================================
# Stochastic quantities replaced by their means for the DP objective:
#   Price  ~ LogNormal(mu_P, sigma_P)  =>  E[P] = exp(mu_P + sigma_P^2/2)
#   Cost   ~ LogNormal(mu_C, sigma_C)  =>  E[C] = exp(mu_C + sigma_C^2/2)
#   Yield  =>  mu from generator, looked up at POST-planting state (next6)
# RCI adjustment cost is deterministic: kappa * max(0, a_now - a_prev).

dp_payoff <- function(window6,
                       crop,             # "Corn" or "Soy"
                       theta,            # coverage level (used only when use_insurance=TRUE)
                       dp_maps,
                       price_pars,       # list(C = list(mu, sigma), S = list(mu, sigma))
                       cost_pars,        # list(C = list(mu, sigma), S = list(mu, sigma))
                       K_table,
                       subsidy_val   = 0,
                       beta          = -0.5,
                       kappa         = 5,
                       use_insurance = TRUE) {

  crop_char <- substr(crop, 1L, 1L)
  next6     <- dp_next_state(window6, crop_char)

  a_prev <- rci_from_string(window6)
  a_now  <- rci_from_string(next6)

  # Yield: look up at POST-planting state (next6) to match simulation convention
  mu_val <- dp_maps$mu[[crop]][[next6]]
  if (is.null(mu_val) || !is.finite(mu_val) || mu_val <= 0)
    return(list(payoff = NA_real_, next6 = next6))

  P_pars <- price_pars[[if (crop == "Corn") "C" else "S"]]
  P_val  <- exp(P_pars$mu + 0.5 * P_pars$sigma^2)
  C_pars <- cost_pars[[if (crop == "Corn") "C" else "S"]]
  C_val  <- exp(C_pars$mu + 0.5 * C_pars$sigma^2)

  C_rci <- kappa * max(0L, a_now - a_prev)
  rev   <- P_val * mu_val

  ins_net <- 0
  if (use_insurance) {
    th_key <- as.character(theta)
    S_val  <- dp_maps$S[[crop]][[th_key]][[next6]]
    if (is.null(S_val) || !is.finite(S_val)) S_val <- 0

    EI    <- P_val * mu_val * S_val
    K_key <- paste(crop, th_key, sep = "_")
    K_val <- if (!is.null(K_table[[K_key]])) as.numeric(K_table[[K_key]])
              else as.numeric(median(unlist(K_table), na.rm = TRUE))
    ins_net <- EI - (1 - subsidy_val) * K_val * mu_val^(1 + beta)
  }

  list(payoff = rev + ins_net - C_val - C_rci, next6 = next6)
}

# ============================================================
# SECTION 6: BACKWARD INDUCTION (BELLMAN DP)
# ============================================================

solve_dp_region_zone <- function(geo_region,
                                  prod_zone,
                                  gens_by_region_zone  = yield_gens_by_region_zone,
                                  price_pars_regional  = price_params_regional,
                                  base_cost_pars       = base_cost_params,
                                  K_table              = NULL,   # NULL -> estimate from data
                                  subsidy_tbl          = NULL,   # NULL -> zero subsidy
                                  thetas               = c(0.65, 0.75, 0.85),
                                  T_horizon            = 6L,
                                  r_disc               = 0.05,
                                  beta                 = -0.5,
                                  kappa                = 5,
                                  unit_structure       = "basic_optional",
                                  use_insurance        = TRUE,
                                  verbose              = FALSE) {

  delta      <- 1 / (1 + r_disc)
  states     <- all_histories6_dp()
  price_pars <- price_pars_regional[[geo_region]]
  cost_pars  <- base_cost_pars[[geo_region]][[prod_zone]]
  if (is.null(price_pars)) stop("price_pars not found for: ", geo_region)
  if (is.null(cost_pars))  stop("cost_pars not found for: ", geo_region, "-", prod_zone)

  dp_maps <- build_dp_maps_region_zone(
    geo_region, prod_zone, gens_by_region_zone, thetas, verbose)

  if (is.null(K_table))
    K_table <- estimate_K_from_maps(dp_maps, price_pars, thetas, beta)

  subsidy_for <- function(theta) {
    if (is.null(subsidy_tbl)) return(0)
    row <- subsidy_tbl[abs(subsidy_tbl$theta - theta) < 1e-8, ]
    if (nrow(row) == 0L || !unit_structure %in% names(row)) return(0)
    val <- as.numeric(row[[unit_structure]][1L])
    if (is.finite(val)) val else 0
  }

  crops      <- c("Corn", "Soy")
  th_choices <- if (use_insurance) thetas else thetas[1L]  # theta irrelevant w/o insurance

  V   <- vector("list", T_horizon + 1L)
  Pol <- vector("list", T_horizon)
  V[[T_horizon + 1L]] <- setNames(rep(0, length(states)), states)

  if (verbose) cat("Backward induction:", geo_region, "-", prod_zone, "\n")

  for (t in T_horizon:1L) {
    Vt    <- setNames(rep(NA_real_, length(states)), states)
    Pol_t <- setNames(vector("list", length(states)), states)

    for (st in states) {
      best_val <- -Inf
      best_act <- NULL

      for (crop in crops) {
        for (th in th_choices) {
          step <- dp_payoff(
            window6 = st, crop = crop, theta = th,
            dp_maps = dp_maps, price_pars = price_pars, cost_pars = cost_pars,
            K_table = K_table,
            subsidy_val   = if (use_insurance) subsidy_for(th) else 0,
            beta = beta, kappa = kappa, use_insurance = use_insurance)
          if (!is.finite(step$payoff)) next

          val <- step$payoff + delta * V[[t + 1L]][[step$next6]]
          if (is.finite(val) && val > best_val) {
            best_val <- val
            best_act <- list(crop = crop, theta = th, next6 = step$next6)
          }
        }
      }
      Vt[[st]]    <- if (is.finite(best_val)) best_val else NA_real_
      Pol_t[[st]] <- best_act
    }

    V[[t]]   <- Vt
    Pol[[t]] <- Pol_t
    if (verbose)
      cat(sprintf("  t=%d done  V range [%.1f, %.1f]\n",
                  t, min(Vt, na.rm = TRUE), max(Vt, na.rm = TRUE)))
  }

  structure(
    list(V = V, Pol = Pol,
         dp_maps = dp_maps, K_table = K_table,
         price_pars = price_pars, cost_pars = cost_pars,
         geo_region = geo_region, prod_zone = prod_zone,
         thetas = thetas, T_horizon = T_horizon,
         r_disc = r_disc, beta = beta, kappa = kappa,
         use_insurance = use_insurance),
    class = "dp_region_zone")
}

# ============================================================
# SECTION 7: POLICY ROLLOUT -> PLAN6 STRING
# ============================================================
# Follows the optimal policy from start6 forward, returns:
#   plan6  : 6-char crop sequence (e.g. "CSCSCS") for use in any simulator
#   npv_ev : certainty-equivalent NPV (expected values only, no MC)
#   path   : data.table with one row per year

rollout_dp_plan <- function(dp_result, start6 = "CCCCCC") {
  Pol       <- dp_result$Pol
  T_horizon <- dp_result$T_horizon
  delta     <- 1 / (1 + dp_result$r_disc)

  st         <- start6
  plan_chars <- character(T_horizon)
  path_rows  <- vector("list", T_horizon)
  npv_ev     <- 0

  for (t in seq_len(T_horizon)) {
    act <- Pol[[t]][[st]]
    if (is.null(act)) {
      warning("No optimal action at t=", t, " state=", st, "; defaulting to Corn")
      act <- list(crop = "Corn", theta = dp_result$thetas[1L],
                  next6 = dp_next_state(st, "C"))
    }
    plan_chars[[t]] <- substr(act$crop, 1L, 1L)

    step <- dp_payoff(
      window6 = st, crop = act$crop, theta = act$theta,
      dp_maps    = dp_result$dp_maps,
      price_pars = dp_result$price_pars,
      cost_pars  = dp_result$cost_pars,
      K_table    = dp_result$K_table,
      beta = dp_result$beta, kappa = dp_result$kappa,
      use_insurance = dp_result$use_insurance)

    npv_ev <- npv_ev + delta^(t - 1L) * step$payoff
    path_rows[[t]] <- data.table(
      t = t, state = st, crop = act$crop,
      theta = act$theta, next6 = act$next6, payoff_ev = step$payoff)
    st <- act$next6
  }

  list(plan6  = paste(plan_chars, collapse = ""),
       npv_ev = npv_ev,
       path   = rbindlist(path_rows))
}

# ============================================================
# SECTION 8: CREDIT SIMULATION UNDER DP-OPTIMAL PLAN
# ============================================================
# Extracts the optimal plan6 string, then delegates entirely to
# simulate_plan_credit_region_zone() so all stochastic mechanics
# (price draws, yield draws, credit rollover) are unchanged.

simulate_dp_credit_region_zone <- function(dp_result,
                                            start6           = "CCCCCC",
                                            B                = 20000L,
                                            i_op_annual      = 0.09,
                                            credit_limit     = 2000,
                                            corr_method      = "none",
                                            correlation_corn = 0,
                                            correlation_soy  = 0,
                                            seed             = 1L,
                                            verbose          = FALSE) {

  rollout <- rollout_dp_plan(dp_result, start6 = start6)
  plan6   <- rollout$plan6

  if (verbose) cat("DP-optimal plan (from", start6, "):", plan6, "\n")

  credit_sim <- simulate_plan_credit_region_zone(
    history6         = start6,
    plan6            = plan6,
    geo_region       = dp_result$geo_region,
    prod_zone        = dp_result$prod_zone,
    B                = B,
    i_op_annual      = i_op_annual,
    credit_limit     = credit_limit,
    delinquency_rule = "end_balance_positive",
    corr_method      = corr_method,
    correlation_corn = correlation_corn,
    correlation_soy  = correlation_soy,
    seed             = seed,
    verbose          = verbose)

  list(plan6      = plan6,
       npv_ev     = rollout$npv_ev,
       path       = rollout$path,
       credit_sim = credit_sim)
}

# ============================================================
# SECTION 9: SOLVE DP ACROSS FULL 3x3 GRID
# ============================================================

run_dp_grid <- function(gens_by_region_zone  = yield_gens_by_region_zone,
                         price_pars_regional  = price_params_regional,
                         base_cost_pars       = base_cost_params,
                         K_table              = NULL,
                         subsidy_tbl          = NULL,
                         thetas               = c(0.65, 0.75, 0.85),
                         T_horizon            = 6L,
                         r_disc               = 0.05,
                         use_insurance        = TRUE,
                         verbose              = FALSE) {

  regions <- c("northern", "central", "southern")
  zones   <- c("high", "medium", "low")
  cat("\n=== Solving DP across 3x3 region-zone grid ===\n")

  dp_grid <- list()
  for (reg in regions) {
    dp_grid[[reg]] <- list()
    for (zone in zones) {
      cat(sprintf("  %s - %s ... ", reg, zone))
      dp_grid[[reg]][[zone]] <- tryCatch(
        solve_dp_region_zone(
          geo_region          = reg,     prod_zone = zone,
          gens_by_region_zone = gens_by_region_zone,
          price_pars_regional = price_pars_regional,
          base_cost_pars      = base_cost_pars,
          K_table             = K_table, subsidy_tbl = subsidy_tbl,
          thetas = thetas,    T_horizon = T_horizon,
          r_disc = r_disc,    use_insurance = use_insurance,
          verbose = verbose),
        error = function(e) { cat("ERROR:", e$message, "\n"); NULL })
      if (!is.null(dp_grid[[reg]][[zone]])) cat("done\n")
    }
  }

  starts       <- c("CCCCCC", "CSCSCS")
  plan_summary <- rbindlist(lapply(starts, function(start) {
    rbindlist(lapply(regions, function(reg) {
      rbindlist(lapply(zones, function(zone) {
        dp <- dp_grid[[reg]][[zone]]
        if (is.null(dp)) return(NULL)
        ro <- rollout_dp_plan(dp, start6 = start)
        data.table(region = reg, prod_zone = zone,
                   start_history = start, optimal_plan = ro$plan6,
                   npv_ev = ro$npv_ev)
      }))
    }))
  }))

  cat("\n=== Optimal plans by region-zone ===\n")
  print(plan_summary[order(start_history, region, prod_zone)])
  list(dp_grid = dp_grid, plan_summary = plan_summary)
}

# ============================================================
# SECTION 10: COMPARE DP-OPTIMAL vs FIXED PLANS IN CREDIT SPACE
# ============================================================

compare_dp_vs_fixed_credit <- function(dp_grid_result,
                                        fixed_plans      = c("CCCCCC", "CSCSCS", "CCSCCS"),
                                        start6           = "CCCCCC",
                                        B                = 20000L,
                                        corr_method      = "decomposed",
                                        correlation_corn = -0.40,
                                        correlation_soy  = -0.50,
                                        seed_base        = 501L,
                                        verbose          = FALSE) {

  regions <- c("northern", "central", "southern")
  zones   <- c("high", "medium", "low")
  dp_grid <- dp_grid_result$dp_grid
  cat("\n=== Comparing DP-optimal vs fixed plans in credit space ===\n")

  all_cell_results <- list()
  cell_idx <- 0L

  for (reg in regions) {
    all_cell_results[[reg]] <- list()
    for (zone in zones) {
      cell_idx <- cell_idx + 1L
      dp <- dp_grid[[reg]][[zone]]
      if (is.null(dp)) next

      dp_plan6  <- rollout_dp_plan(dp, start6 = start6)$plan6
      plans_run <- unique(c(fixed_plans, dp_plan6))
      cell_sims <- list()

      for (plan in plans_run) {
        cell_seed <- seed_base * 1000L + cell_idx
        cat(sprintf("  %s-%s | %-8s (seed=%d)\n", reg, zone, plan, cell_seed))
        cell_sims[[plan]] <- tryCatch(
          simulate_plan_credit_region_zone(
            history6         = start6,  plan6 = plan,
            geo_region       = reg,     prod_zone = zone,
            B                = B,
            i_op_annual      = 0.09,    credit_limit = 2000,
            delinquency_rule = "end_balance_positive",
            corr_method      = corr_method,
            correlation_corn = correlation_corn,
            correlation_soy  = correlation_soy,
            seed             = cell_seed,
            verbose          = FALSE),
          error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
      }
      all_cell_results[[reg]][[zone]] <- list(sims = cell_sims, dp_plan6 = dp_plan6)
    }
  }

  summary_dt <- rbindlist(lapply(regions, function(reg) {
    rbindlist(lapply(zones, function(zone) {
      cell <- all_cell_results[[reg]][[zone]]
      if (is.null(cell)) return(NULL)
      rbindlist(lapply(names(cell$sims), function(plan) {
        r <- cell$sims[[plan]]
        if (is.null(r)) return(NULL)
        data.table(
          region         = reg,   prod_zone = zone,
          plan           = plan,
          is_dp_optimal  = (plan == cell$dp_plan6),
          pr_delinq      = r$pr_delinquent_any,
          exp_delinq_yrs = r$exp_delinquent_years,
          mean_peak_debt = r$mean_peak_debt,
          pv_net_cf_mean = mean(r$pv_net_cf, na.rm = TRUE))
      }))
    }))
  }))

  cat("\n=== Delinquency + NPV: DP-optimal vs fixed ===\n")
  print(summary_dt[order(region, prod_zone, -is_dp_optimal, plan)])
  list(results = all_cell_results, summary = summary_dt)
}

# ============================================================
# EXAMPLE USAGE (set FALSE -> TRUE to run)
# ============================================================
if (FALSE) {

  source("dp_functions_region_zone.R")

  # --- Single cell: solve and inspect ---
  dp_cm <- solve_dp_region_zone(
    geo_region    = "central",
    prod_zone     = "medium",
    use_insurance = FALSE,     # pure crop-rotation DP, no insurance in objective
    verbose       = TRUE
  )

  ro <- rollout_dp_plan(dp_cm, start6 = "CCCCCC")
  cat("Optimal plan:", ro$plan6, "\n")
  cat("EV-NPV ($/ac):", round(ro$npv_ev, 2), "\n")
  print(ro$path)

  # Contrast: same cell, different starting history
  rollout_dp_plan(dp_cm, start6 = "CCCCCC")$plan6
  rollout_dp_plan(dp_cm, start6 = "CSCSCS")$plan6

  # --- Simulate credit model under optimal plan ---
  dp_credit <- simulate_dp_credit_region_zone(
    dp_cm, start6 = "CCCCCC", B = 20000, corr_method = "decomposed")
  cat("Pr(delinquent):", dp_credit$credit_sim$pr_delinquent_any, "\n")
  cat("Mean peak debt:", dp_credit$credit_sim$mean_peak_debt, "\n")

  # --- Full 3x3 grid (no insurance in objective) ---
  dp_grid <- run_dp_grid(use_insurance = FALSE, verbose = FALSE)
  print(dp_grid$plan_summary)

  # --- Compare DP-optimal vs benchmark plans in credit space ---
  comp <- compare_dp_vs_fixed_credit(
    dp_grid_result   = dp_grid,
    fixed_plans      = c("CCCCCC", "CSCSCS", "CCSCCS"),
    start6           = "CCCCCC",
    B                = 20000,
    corr_method      = "decomposed"
  )
  print(comp$summary)

  # Which cells chose something different from the three benchmarks?
  comp$summary[is_dp_optimal == TRUE & !plan %in% c("CCCCCC","CSCSCS","CCSCCS")]

  # --- With insurance in the DP objective ---
  # Supply subsidy_tbl from price_adj.R's subsidy_tbl object
  dp_grid_ins <- run_dp_grid(
    use_insurance = TRUE,
    subsidy_tbl   = subsidy_tbl,    # bring in from price_adj.R environment
    thetas        = c(0.65, 0.75, 0.85)
  )
  # Optimal coverage level by year for central/medium
  rollout_dp_plan(dp_grid_ins$dp_grid$central$medium, "CCCCCC")$path[, .(t, crop, theta)]
}
