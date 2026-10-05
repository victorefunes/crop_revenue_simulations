# ============================================================================
# ADAPTED CREDIT/OPERATING LOAN FUNCTIONS FOR REGION-ZONE FRAMEWORK
# ============================================================================
# Intra-year cashflow timing with operating credit
# ============================================================================

# ============================================================
# CREDIT MODEL WITH REGION-ZONE PARAMETERS
# ============================================================

simulate_plan_credit_region_zone <- function(history6,
                                             plan6,
                                             geo_region = "central",
                                             prod_zone = "medium",
                                             B = 20000L,
                                             r_disc = 0.05,
                                             gens_by_region_zone = yield_gens_by_region_zone,
                                             price_pars_regional = price_params_regional,
                                             base_cost_pars = base_cost_params,
                                             seed = 1L,
                                             # --- Credit terms ---
                                             i_op_annual = 0.09,            # Operating note APR
                                             months_borrowed = 9,           # Borrow ~Jan, repay ~Oct
                                             credit_limit = 1200,           # $/ac maximum operating debt
                                             repay_fraction_required = 1.0, # 1.0 = must fully repay
                                             penalty_rate = 0.00,           # Extra rate on carried balance
                                             delinquency_rule = c("end_balance_positive", "over_limit"),
                                             # --- Price-yield correlation ---
                                             correlation_corn = 0,
                                             correlation_soy  = 0,
                                             aggregate_fraction = 0.40,
                                             corr_method = c("none", "decomposed", "simple", "t_copula"),
                                             copula_df   = 5L,
                                             # --- Price shock (optional) ---
                                             # Named list with elements $S and/or $C: length-6 numeric
                                             # vectors of multiplicative factors applied to drawn prices
                                             # after all price draws. E.g. list(S = c(1,1.25,1,1,1,1))
                                             # raises soy prices 25% in year 2 only.
                                             price_shock = NULL,
                                             # --- Yield shock (optional) ---
                                             # Same structure as price_shock but applied to Y post-draw.
                                             # Factors are applied only to years whose crop matches the
                                             # list element ($C for corn years, $S for soy years).
                                             # E.g. list(C = c(1,1,0.9,1,1,1)) cuts corn yield 10% in yr 3.
                                             yield_shock = NULL,
                                             # --- Uncertainty freezing (Shapley decomposition) ---
                                             # fix_price = TRUE: replace all price draws with E[P] = Pbar
                                             # fix_yield = TRUE: replace all yield draws with E[Y] per year
                                             fix_price = FALSE,
                                             fix_yield = FALSE,
                                             verbose = FALSE) {
 # Simulate operating credit model with intra-year timing
  
 #  Mechanics:
 #    - Borrow at planting to finance costs
 #    - Interest accrues for months_borrowed
 #    - Repay at harvest from revenue
 #     - Carry shortfall (rollover) into next year
 #    - Track delinquency events
  
 # Returns:
 #    Year-by-year debt, repayment, net cash, and delinquency flags
  
  delinquency_rule <- match.arg(delinquency_rule)
  corr_method      <- match.arg(corr_method)
  use_corr <- corr_method != "none" &&
              (abs(correlation_corn) > 0 || abs(correlation_soy) > 0)

  stopifnot(nchar(history6) == 6, nchar(plan6) == 6)
  set.seed(seed)

  # Build transition path
  path <- make_transition_path(history6, plan6)
  disc <- disc_vec(r_disc)

  # Select region-specific prices
  if (geo_region %in% names(price_pars_regional)) {
    price_pars <- price_pars_regional[[geo_region]]
  } else {
    stop("Region not found: ", geo_region)
  }

  # Select region-zone specific costs
  if (geo_region %in% names(base_cost_pars) &&
      prod_zone %in% names(base_cost_pars[[geo_region]])) {
    cost_pars <- base_cost_pars[[geo_region]][[prod_zone]]
  } else {
    stop("Region-zone combination not found: ", geo_region, "-", prod_zone)
  }

  # Select region-zone specific yield generators
  corn_gen <- gens_by_region_zone$corn[[geo_region]][[prod_zone]]
  soy_gen  <- gens_by_region_zone$soy[[geo_region]][[prod_zone]]

  if (is.null(corn_gen) || is.null(soy_gen)) {
    stop("Yield generators not available: ", geo_region, "-", prod_zone)
  }

  # Crop indicators (needed by both draw paths)
  chars    <- strsplit(plan6, "")[[1]]
  corn_yrs <- which(chars == "C")
  soy_yrs  <- which(chars == "S")
  IC <- matrix(as.numeric(chars == "C"), nrow = B, ncol = 6, byrow = TRUE)
  IS <- matrix(as.numeric(chars == "S"), nrow = B, ncol = 6, byrow = TRUE)

  # Costs drawn independently regardless of correlation setting
  Cc <- matrix(rlnorm(B * 6, cost_pars$C$mu, cost_pars$C$sigma), nrow = B, ncol = 6)
  Cs <- matrix(rlnorm(B * 6, cost_pars$S$mu, cost_pars$S$sigma), nrow = B, ncol = 6)

  Pc <- matrix(0, nrow = B, ncol = 6)
  Ps <- matrix(0, nrow = B, ncol = 6)
  Y  <- matrix(NA_real_, nrow = B, ncol = 6)

  if (!use_corr) {
    # Independent draws (original behaviour)
    for (t in 1:6) {
      crop_t <- path$crop[t]
      st_t   <- path$window6[t]
      if (crop_t == "C") {
        Pc[, t] <- rlnorm(B, price_pars$C$mu, price_pars$C$sigma)
        Y[, t]  <- corn_gen$simulate(st_t, B = B)
      } else if (crop_t == "S") {
        Ps[, t] <- rlnorm(B, price_pars$S$mu, price_pars$S$sigma)
        Y[, t]  <- soy_gen$simulate(st_t, B = B)
      }
    }
  } else {
    # Correlated price-yield draws
    for (i in seq_along(corn_yrs)) {
      t  <- corn_yrs[i]
      st <- path$window6[t]
      mu <- corn_gen$mu[rot6 == st, mu]
      sg <- sd(corn_gen$pools[rot6 == st, resid_pool][[1]])
      if (corr_method == "decomposed") {
        res <- draw_correlated_price_yield_decomposed(
                 B = B, T = 1, price_params = price_pars$C,
                 yield_mean = mu, yield_sd_total = sg,
                 price_yield_correlation = correlation_corn,
                 aggregate_fraction = aggregate_fraction)
      } else if (corr_method == "t_copula") {
        res <- draw_correlated_price_yield_tcopula(
                 B = B, T = 1, price_params = price_pars$C,
                 yield_mean = mu, yield_sd_total = sg,
                 price_yield_correlation = correlation_corn,
                 aggregate_fraction = aggregate_fraction,
                 df = copula_df)
      } else {
        res <- draw_correlated_price_yield_simple(
                 B = B, T = 1,
                 price_mean = price_pars$C$mu, price_sd = price_pars$C$sigma,
                 yield_mean = mu, yield_sd = sg,
                 correlation = correlation_corn)
      }
      Pc[, t] <- as.vector(res$prices)
      Y[, t]  <- as.vector(res$yields)
    }

    for (i in seq_along(soy_yrs)) {
      t  <- soy_yrs[i]
      st <- path$window6[t]
      mu <- soy_gen$mu[rot6 == st, mu]
      sg <- sd(soy_gen$pools[rot6 == st, resid_pool][[1]])
      if (corr_method == "decomposed") {
        res <- draw_correlated_price_yield_decomposed(
                 B = B, T = 1, price_params = price_pars$S,
                 yield_mean = mu, yield_sd_total = sg,
                 price_yield_correlation = correlation_soy,
                 aggregate_fraction = aggregate_fraction)
      } else if (corr_method == "t_copula") {
        res <- draw_correlated_price_yield_tcopula(
                 B = B, T = 1, price_params = price_pars$S,
                 yield_mean = mu, yield_sd_total = sg,
                 price_yield_correlation = correlation_soy,
                 aggregate_fraction = aggregate_fraction,
                 df = copula_df)
      } else {
        res <- draw_correlated_price_yield_simple(
                 B = B, T = 1,
                 price_mean = price_pars$S$mu, price_sd = price_pars$S$sigma,
                 yield_mean = mu, yield_sd = sg,
                 correlation = correlation_soy)
      }
      Ps[, t] <- as.vector(res$prices)
      Y[, t]  <- as.vector(res$yields)
    }

    if (verbose) {
      rc <- if (length(corn_yrs) > 0)
        cor(as.vector(Pc[, corn_yrs, drop = FALSE]),
            as.vector(Y[,  corn_yrs, drop = FALSE])) else NA
      rs <- if (length(soy_yrs) > 0)
        cor(as.vector(Ps[, soy_yrs, drop = FALSE]),
            as.vector(Y[,  soy_yrs, drop = FALSE])) else NA
      cat("Realized price-yield correlations — corn:", round(rc, 3),
          " soy:", round(rs, 3), "\n")
    }
  }

  if (!is.null(price_shock)) {
    if (!is.null(price_shock$S))
      Ps <- Ps * matrix(price_shock$S, nrow = B, ncol = 6, byrow = TRUE)
    if (!is.null(price_shock$C))
      Pc <- Pc * matrix(price_shock$C, nrow = B, ncol = 6, byrow = TRUE)
  }

  if (!is.null(yield_shock)) {
    if (!is.null(yield_shock$C))
      for (t in corn_yrs) Y[, t] <- Y[, t] * yield_shock$C[t]
    if (!is.null(yield_shock$S))
      for (t in soy_yrs)  Y[, t] <- Y[, t] * yield_shock$S[t]
  }

  if (fix_price) {
    Pc[] <- exp(price_pars$C$mu + 0.5 * price_pars$C$sigma^2)
    Ps[] <- exp(price_pars$S$mu + 0.5 * price_pars$S$sigma^2)
  }
  if (fix_yield) {
    for (t in 1:6) Y[, t] <- mean(Y[, t], na.rm = TRUE)
  }

  P <- Pc * IC + Ps * IS
  K <- Cc * IC + Cs * IS
  
  rev  <- P * Y   # Harvest revenue ($/ac)
  cost <- K       # Operating cost financed up front ($/ac)
  
  # --- Operating note mechanics ---
  # Borrow cost at planting. Interest accrues for months_borrowed.
  i_month <- i_op_annual / 12
  growth_factor <- (1 + i_month)^months_borrowed
  
  debt_start <- numeric(B)      # Debt carried into the year (from prior shortfall)
  debt_end   <- matrix(NA_real_, nrow = B, ncol = 6)
  
  repay_amt  <- matrix(NA_real_, nrow = B, ncol = 6)
  cf_net     <- matrix(NA_real_, nrow = B, ncol = 6)  # Net cash after repayment
  
  delinquent <- matrix(FALSE, nrow = B, ncol = 6)
  
  for (t in 1:6) {
    # Amount that must be financed this year = new costs + carried debt
    principal_needed <- cost[, t] + debt_start
    
    # Debt at harvest before repayment (interest on financed principal)
    debt_at_harvest <- principal_needed * growth_factor * (1 + penalty_rate)
    
    # Repayment from harvest revenue
    repayment_capacity <- rev[, t]
    
    # How much lender requires you to repay
    required_payment <- repay_fraction_required * debt_at_harvest
    
    actual_payment <- pmin(repayment_capacity, required_payment)
    
    # Remaining balance after harvest
    remaining_balance <- debt_at_harvest - actual_payment
    
    # Net cash after making the payment (if revenue exceeds payment)
    net_cash <- repayment_capacity - actual_payment
    
    # Delinquency flags
    if (delinquency_rule == "end_balance_positive") {
      delin_flag <- remaining_balance > 0
    } else { # "over_limit"
      delin_flag <- remaining_balance > credit_limit
    }
    
    debt_end[, t] <- remaining_balance
    repay_amt[, t] <- actual_payment
    cf_net[, t] <- net_cash
    delinquent[, t] <- delin_flag
    
    # Rollover to next year
    debt_start <- remaining_balance
  }
  
  # Discounted metrics
  pv_net_cf <- as.numeric(cf_net %*% disc)
  
  list(
    history6 = history6,
    plan6 = plan6,
    region = geo_region,
    productivity = prod_zone,
    path = path,
    rev = rev,
    cost = cost,
    net_cash = cf_net,
    debt_end = debt_end,
    repay_amt = repay_amt,
    delinquent = delinquent,
    pv_net_cf = pv_net_cf,
    # Summary metrics
    pr_delinquent_any = mean(rowSums(delinquent) > 0),
    exp_delinquent_years = mean(rowSums(delinquent)),
    pr_over_limit_any = mean(rowSums(debt_end > credit_limit) > 0),
    peak_debt = apply(debt_end, 1, max, na.rm = TRUE),
    mean_peak_debt = mean(apply(debt_end, 1, max, na.rm = TRUE)),
    p90_peak_debt = as.numeric(quantile(apply(debt_end, 1, max, na.rm = TRUE), 0.9))
  )
}

# ============================================================
# COMPARE TWO PLANS IN CREDIT SPACE
# ============================================================

compare_plans_credit_region_zone <- function(planA, 
                                             planB,
                                             geo_region = "central",
                                             prod_zone = "medium",
                                             B = 20000L,
                                             i_op_annual = 0.09,
                                             credit_limit = 1200,
                                             months_borrowed = 9,
                                             seedA = 11,
                                             seedB = 22) {
  # Compare two plans in credit/delinquency space for a region-zone
  
  outA <- simulate_plan_credit_region_zone(
    history6 = planA, plan6 = planA,
    geo_region = geo_region, prod_zone = prod_zone,
    B = B,
    i_op_annual = i_op_annual,
    credit_limit = credit_limit,
    months_borrowed = months_borrowed,
    seed = seedA,
    delinquency_rule = "end_balance_positive",
    verbose = FALSE
  )
  
  outB <- simulate_plan_credit_region_zone(
    history6 = planB, plan6 = planB,
    geo_region = geo_region, prod_zone = prod_zone,
    B = B,
    i_op_annual = i_op_annual,
    credit_limit = credit_limit,
    months_borrowed = months_borrowed,
    seed = seedB,
    delinquency_rule = "end_balance_positive",
    verbose = FALSE
  )
  
  list(
    region = geo_region,
    productivity = prod_zone,
    A = outA,
    B = outB,
    delta_pv_net_cf = outB$pv_net_cf - outA$pv_net_cf,
    delta_peak_debt = outB$peak_debt - outA$peak_debt,
    comparison = data.table(
      region = geo_region,
      productivity = prod_zone,
      plan = c(planA, planB),
      pr_delinquent_any = c(outA$pr_delinquent_any, outB$pr_delinquent_any),
      exp_delinquent_years = c(outA$exp_delinquent_years, outB$exp_delinquent_years),
      mean_peak_debt = c(outA$mean_peak_debt, outB$mean_peak_debt),
      p90_peak_debt = c(outA$p90_peak_debt, outB$p90_peak_debt)
    )
  )
}

# ============================================================
# GRID CREDIT ANALYSIS
# ============================================================

compare_credit_across_grid <- function(planA,
                                      planB,
                                      B = 20000L,
                                      i_op_annual = 0.09,
                                      credit_limit = 1200) {
  # Run credit analysis across all 9 region-zone combinations
  
  cat("\n=== CREDIT ANALYSIS ACROSS GRID ===\n")
  cat("Comparing:", planA, "vs", planB, "\n\n")
  
  results <- list()
  
  for (reg in c("northern", "central", "southern")) {
    results[[reg]] <- list()
    
    for (zone in c("high", "medium", "low")) {
      cat(sprintf("Running %s-%s... ", reg, zone))
      
      tryCatch({
        results[[reg]][[zone]] <- compare_plans_credit_region_zone(
          planA = planA,
          planB = planB,
          geo_region = reg,
          prod_zone = zone,
          B = B,
          i_op_annual = i_op_annual,
          credit_limit = credit_limit
        )
        cat("✓\n")
      }, error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        results[[reg]][[zone]] <- NULL
      })
    }
  }
  
  # Create summary table
  summary_list <- list()
  for (reg in names(results)) {
    for (zone in names(results[[reg]])) {
      if (!is.null(results[[reg]][[zone]])) {
        comp <- results[[reg]][[zone]]$comparison
        comp[, region := reg]
        comp[, productivity := zone]
        summary_list[[paste0(reg, "_", zone)]] <- comp
      }
    }
  }
  
  full_summary <- rbindlist(summary_list)
  
  cat("\n=== DELINQUENCY RISK SUMMARY ===\n")
  print(full_summary[, .(region, productivity, plan, pr_delinquent_any, mean_peak_debt)])
  
  list(
    results = results,
    summary = full_summary
  )
}

# ============================================================
# VISUALIZATION FUNCTIONS
# ============================================================

plot_debt_paths <- function(credit_result, n_paths = 100) {
  # Plot sample debt paths over time
  
  debt_long <- to_long_mat(
    credit_result$debt_end[1:n_paths, , drop = FALSE],
    plan = credit_result$plan6,
    value_name = "debt"
  )
  
  ggplot(debt_long, aes(x = year, y = debt, group = sim)) +
    geom_line(alpha = 0.1, color = "steelblue") +
    stat_summary(
      aes(group = 1), 
      fun = median, 
      geom = "line", 
      color = "darkblue", 
      linewidth = 1.5
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen") +
    labs(
      title = sprintf("Debt Paths Over Time\n%s IL - %s Productivity",
                      tools::toTitleCase(credit_result$region),
                      tools::toTitleCase(credit_result$productivity)),
      subtitle = paste("Plan:", credit_result$plan6, 
                      sprintf("(Pr(delinquent) = %.2f)", credit_result$pr_delinquent_any)),
      x = "Year (1-6)",
      y = "End-of-Year Debt ($/acre)"
    ) +
    theme_minimal()
}

plot_delinquency_comparison <- function(credit_comparison) {
  # Compare delinquency risk between two plans
  
  dt <- copy(credit_comparison$comparison)
  
  # Add region and productivity if not present
  if (!"region" %in% names(dt)) {
    dt[, region := credit_comparison$region]
  }
  if (!"productivity" %in% names(dt)) {
    dt[, productivity := credit_comparison$productivity]
  }
  
  dt_long <- melt(dt, 
                  id.vars = c("plan"),
                  measure.vars = c("pr_delinquent_any", "mean_peak_debt"),
                  variable.name = "metric",
                  value.name = "value")
  
  ggplot(dt_long, aes(x = plan, y = value, fill = plan)) +
    geom_col() +
    facet_wrap(~ metric, scales = "free_y", 
               labeller = as_labeller(c(
                 pr_delinquent_any = "Pr(Any Delinquency)",
                 mean_peak_debt = "Mean Peak Debt ($/ac)"
               ))) +
    labs(
      title = sprintf("Credit Risk Comparison\n%s IL - %s Productivity",
                      tools::toTitleCase(credit_comparison$region),
                      tools::toTitleCase(credit_comparison$productivity)),
      x = "Plan",
      y = "Value",
      fill = "Plan"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

# ============================================================
# EXAMPLE USAGE
# ============================================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Example 1: Single region-zone credit analysis
  credit_result <- simulate_plan_credit_region_zone(
    history6 = "CCCCCC",
    plan6 = "CCCCCC",
    geo_region = "central",
    prod_zone = "medium",
    B = 20000,
    i_op_annual = 0.09,
    credit_limit = 1200
  )
  
  # View delinquency metrics
  cat("Pr(any delinquency):", credit_result$pr_delinquent_any, "\n")
  cat("Expected delinquent years:", credit_result$exp_delinquent_years, "\n")
  cat("Mean peak debt:", credit_result$mean_peak_debt, "\n")
  
  # Visualize debt paths
  plot_debt_paths(credit_result, n_paths = 100)
  
  # Example 2: Compare two plans
  credit_comp <- compare_plans_credit_region_zone(
    planA = "CCCCCC",
    planB = "CSCSCS",
    geo_region = "central",
    prod_zone = "high",
    B = 20000
  )
  
  # View comparison
  print(credit_comp$comparison)
  
  # Visualize
  plot_delinquency_comparison(credit_comp)
  
  # Example 3: Full grid credit analysis
  grid_credit <- compare_credit_across_grid(
    planA = "CCCCCC",
    planB = "CSCSCS",
    B = 20000,
    i_op_annual = 0.09,
    credit_limit = 1200
  )
  
  # Which region-zones have highest delinquency risk for continuous corn?
  grid_credit$summary[plan == "CCCCCC"][order(-pr_delinquent_any)]
}
