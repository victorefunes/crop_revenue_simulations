# ============================================================================
# ADAPTED CASHFLOW ANALYSIS FUNCTIONS FOR REGION-ZONE FRAMEWORK
# ============================================================================
# These functions provide year-by-year cashflow analysis across the 3×3 grid
# ============================================================================

# ---- Required packages ----------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, tidyverse, ggplot2)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

disc_vec <- function(r) 1 / (1 + r)^(0:5)

draw_lognorm_mat <- function(B, T, mu_log, sd_log) {
  matrix(rlnorm(B * T, meanlog = mu_log, sdlog = sd_log), nrow = B, ncol = T)
}

quantile_vec <- function(x, probs = c(0.1, 0.5, 0.9), na.rm = TRUE) {
  as.numeric(stats::quantile(x, probs = probs, na.rm = na.rm))
}

to_long_mat <- function(mat, plan, value_name = "value") {
  stopifnot(is.matrix(mat))
  data.table(
    plan = plan,
    sim  = rep(seq_len(nrow(mat)), times = ncol(mat)),
    year = rep(seq_len(ncol(mat)), each  = nrow(mat)),
    value = as.vector(mat)
  )[, (value_name) := value][, value := NULL]
}

# ============================================================
# CORE: SIMULATE FULL CASHFLOWS FOR A REGION-ZONE
# ============================================================

simulate_plan_cashflows_region_zone <- function(history6,
                                                plan6,
                                                geo_region = "central",
                                                prod_zone = "medium",
                                                B = 20000L,
                                                r_disc = 0.05,
                                                gens_by_region_zone = yield_gens_by_region_zone,
                                                price_pars_regional = price_params_regional,
                                                base_cost_pars = base_cost_params,
                                                use_rotation_costs = TRUE,
                                                seed = 123,
                                                verbose = FALSE) {

  # Simulate year-by-year cashflows for a specific region-zone combination
  
  # Returns:
  #  Matrices of prices, costs, yields, revenue, profit by year
  #  Plus cumulative profit and NPV metrics
  
  stopifnot(nchar(history6) == 6, nchar(plan6) == 6)
  set.seed(seed)
  
  # Build year-by-year transition path
  path <- make_transition_path(history6, plan6)
  T <- nrow(path)
  stopifnot(T == 6)
  
  if (verbose) {
    cat("\n--- Cashflow path ---\n")
    cat("history6: ", history6, "\n", sep="")
    cat("plan6   : ", plan6, "\n", sep="")
    cat("region  : ", geo_region, "\n", sep="")
    cat("zone    : ", prod_zone, "\n", sep="")
    print(path)
  }
  
  # Select region-specific prices
  if (geo_region %in% names(price_pars_regional)) {
    price_pars <- price_pars_regional[[geo_region]]
  } else {
    stop("Region not found in price parameters: ", geo_region)
  }
  
  # Select region-zone specific costs
  if (geo_region %in% names(base_cost_pars) && 
      prod_zone %in% names(base_cost_pars[[geo_region]])) {
    cost_pars_base <- base_cost_pars[[geo_region]][[prod_zone]]
  } else {
    stop("Region-zone combination not found: ", geo_region, "-", prod_zone)
  }
  
  # Select region-zone specific yield generators
  corn_gen <- gens_by_region_zone$corn[[geo_region]][[prod_zone]]
  soy_gen <- gens_by_region_zone$soy[[geo_region]][[prod_zone]]
  
  if (is.null(corn_gen) || is.null(soy_gen)) {
    stop("Yield generators not available for: ", geo_region, "-", prod_zone)
  }
  
  # Draw year-specific prices (region-specific)
  Pc <- draw_lognorm_mat(B, T, price_pars$C$mu, price_pars$C$sigma)
  Ps <- draw_lognorm_mat(B, T, price_pars$S$mu, price_pars$S$sigma)
  
  # Draw base costs (region-zone specific)
  Cc_base <- draw_lognorm_mat(B, T, cost_pars_base$C$mu, cost_pars_base$C$sigma)
  Cs_base <- draw_lognorm_mat(B, T, cost_pars_base$S$mu, cost_pars_base$S$sigma)
  
  # Apply rotation cost adjustments if requested
  if (use_rotation_costs) {
    # This would call the rotation adjustment function
    # For now, using base costs (rotation adjustments handled in main simulation)
    Cc <- Cc_base
    Cs <- Cs_base
  } else {
    Cc <- Cc_base
    Cs <- Cs_base
  }
  
  # Draw yields: one per year, conditioned on that year's window6 and crop
  Y <- matrix(NA_real_, nrow = B, ncol = T)
  for (t in 1:T) {
    crop_t <- path$crop[t]
    st_t   <- path$window6[t]
    
    if (crop_t == "C") {
      Y[, t] <- corn_gen$simulate(st_t, B = B)
    } else if (crop_t == "S") {
      Y[, t] <- soy_gen$simulate(st_t, B = B)
    } else {
      stop("Unknown crop char in plan: ", crop_t)
    }
  }
  
  # Crop indicators for each year
  chars <- strsplit(plan6, "")[[1]]
  isC <- as.numeric(chars == "C")
  isS <- as.numeric(chars == "S")
  IC <- matrix(isC, nrow = B, ncol = T, byrow = TRUE)
  IS <- matrix(isS, nrow = B, ncol = T, byrow = TRUE)
  
  # Actual prices and costs faced each year
  P <- (Pc * IC) + (Ps * IS)
  K <- (Cc * IC) + (Cs * IS)
  
  # Revenue, profit matrices
  revenue <- P * Y
  profit  <- revenue - K
  
  # Discount factors
  disc <- disc_vec(r_disc)[1:T]
  
  # Present values
  pv_rev   <- as.numeric(revenue %*% disc)
  pv_cost  <- as.numeric(K %*% disc)
  npv_prof <- as.numeric(profit %*% disc)
  
  # Cumulative (undiscounted) cash flow paths
  cum_profit <- t(apply(profit, 1, cumsum))
  
  list(
    history6 = history6,
    plan6 = plan6,
    region = geo_region,
    productivity = prod_zone,
    path = path,
    prices = list(Pc = Pc, Ps = Ps),
    costs  = list(Cc = Cc, Cs = Cs),
    yields = Y,
    revenue = revenue,
    cost = K,
    profit = profit,
    cum_profit = cum_profit,
    pv_revenue = pv_rev,
    pv_cost = pv_cost,
    npv_profit = npv_prof
  )
}

# ============================================================
# SUMMARIZE CASHFLOWS
# ============================================================

summarize_cashflows_region_zone <- function(cf_obj, probs = c(0.1, 0.5, 0.9)) {
  # Calculate annual and cumulative profit summaries plus liquidity metrics
  
  profit <- cf_obj$profit
  cumP   <- cf_obj$cum_profit
  T <- ncol(profit)
  
  # Annual profit summaries
  ann <- rbindlist(lapply(1:T, function(t) {
    x <- profit[, t]
    q <- quantile_vec(x, probs = probs, na.rm = TRUE)
    data.table(
      year = t,
      mean = mean(x, na.rm = TRUE),
      sd   = stats::sd(x, na.rm = TRUE),
      p_lo = q[1], 
      p_mid = q[2], 
      p_hi = q[3],
      pr_loss = mean(x < 0, na.rm = TRUE)
    )
  }))
  
  # Cumulative profit summaries (undiscounted)
  cum <- rbindlist(lapply(1:T, function(t) {
    x <- cumP[, t]
    q <- quantile_vec(x, probs = probs, na.rm = TRUE)
    data.table(
      year = t,
      mean = mean(x, na.rm = TRUE),
      sd   = stats::sd(x, na.rm = TRUE),
      p_lo = q[1], 
      p_mid = q[2], 
      p_hi = q[3],
      pr_cum_negative = mean(x < 0, na.rm = TRUE)
    )
  }))
  
  # Liquidity / risk metrics (pathwise)
  # - Ever negative cumulative CF
  ever_neg <- mean(apply(cumP, 1, function(v) any(v < 0)), na.rm = TRUE)
  
  # - Max drawdown: max(peak - current) over time
  max_dd <- apply(cumP, 1, function(v) {
    peak <- cummax(v)
    max(peak - v, na.rm = TRUE)
  })
  dd_summary <- data.table(
    max_drawdown_mean = mean(max_dd, na.rm = TRUE),
    max_drawdown_p50  = as.numeric(stats::quantile(max_dd, 0.5, na.rm = TRUE)),
    max_drawdown_p90  = as.numeric(stats::quantile(max_dd, 0.9, na.rm = TRUE))
  )
  
  # - First year cumulative becomes positive (payback), NA if never
  payback <- apply(cumP, 1, function(v) {
    idx <- which(v >= 0)
    if (length(idx) == 0) NA_integer_ else idx[1]
  })
  payback_summary <- data.table(
    pr_payback_by6 = mean(!is.na(payback) & payback <= T, na.rm = TRUE),
    payback_p50 = as.numeric(stats::quantile(payback, 0.5, na.rm = TRUE)),
    payback_p90 = as.numeric(stats::quantile(payback, 0.9, na.rm = TRUE))
  )
  
  list(
    plan = cf_obj$plan6,
    history = cf_obj$history6,
    region = cf_obj$region,
    productivity = cf_obj$productivity,
    ann_profit = ann,
    cum_profit = cum,
    liquidity = data.table(
      pr_ever_cum_negative = ever_neg,
      npv_mean = mean(cf_obj$npv_profit, na.rm = TRUE),
      npv_p10  = as.numeric(stats::quantile(cf_obj$npv_profit, 0.1, na.rm = TRUE)),
      npv_p50  = as.numeric(stats::quantile(cf_obj$npv_profit, 0.5, na.rm = TRUE)),
      npv_p90  = as.numeric(stats::quantile(cf_obj$npv_profit, 0.9, na.rm = TRUE))
    ),
    drawdown = dd_summary,
    payback = payback_summary
  )
}

# ============================================================
# COMPARE TWO PLANS IN CASHFLOW SPACE (REGION-ZONE SPECIFIC)
# ============================================================

compare_two_plans_cashflow_region_zone <- function(history6,
                                                   planA,
                                                   planB,
                                                   geo_region = "central",
                                                   prod_zone = "medium",
                                                   B = 20000L,
                                                   r_disc = 0.05,
                                                   gens_by_region_zone = yield_gens_by_region_zone,
                                                   price_pars_regional = price_params_regional,
                                                   base_cost_pars = base_cost_params,
                                                   use_rotation_costs = TRUE,
                                                   seedA = 101,
                                                   seedB = 202,
                                                   verbose = FALSE) {
  
  # Compare two plans from same starting history in cashflow space
  # for a specific region-zone
  
  cfA <- simulate_plan_cashflows_region_zone(
    history6, planA, geo_region, prod_zone, B, r_disc,
    gens_by_region_zone, price_pars_regional, base_cost_pars,
    use_rotation_costs, seedA, verbose
  )
  
  cfB <- simulate_plan_cashflows_region_zone(
    history6, planB, geo_region, prod_zone, B, r_disc,
    gens_by_region_zone, price_pars_regional, base_cost_pars,
    use_rotation_costs, seedB, verbose
  )
  
  sumA <- summarize_cashflows_region_zone(cfA)
  sumB <- summarize_cashflows_region_zone(cfB)
  
  # Delta NPV distribution (B - A)
  delta_npv <- cfB$npv_profit - cfA$npv_profit
  delta_summary <- data.table(
    mean = mean(delta_npv, na.rm = TRUE),
    sd   = sd(delta_npv, na.rm = TRUE),
    p10  = as.numeric(quantile(delta_npv, 0.1, na.rm = TRUE)),
    p50  = as.numeric(quantile(delta_npv, 0.5, na.rm = TRUE)),
    p90  = as.numeric(quantile(delta_npv, 0.9, na.rm = TRUE)),
    pr_positive = mean(delta_npv > 0, na.rm = TRUE)
  )
  
  list(
    region = geo_region,
    productivity = prod_zone,
    A = list(cf = cfA, summary = sumA),
    B = list(cf = cfB, summary = sumB),
    delta = list(delta_npv = delta_npv, summary = delta_summary)
  )
}

# ============================================================
# VISUALIZATION FUNCTIONS
# ============================================================

plot_annual_profit_bands_region_zone <- function(summary_obj, 
                                                 title = NULL) {
  # Plot annual profit distribution with uncertainty bands
  
  if (is.null(title)) {
    title <- sprintf("Annual Profit Distribution\n%s IL - %s Productivity",
                     tools::toTitleCase(summary_obj$region),
                     tools::toTitleCase(summary_obj$productivity))
  }
  
  ggplot(summary_obj$ann_profit, aes(x = year)) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.25, fill = "steelblue") +
    geom_line(aes(y = p_mid), linewidth = 1, color = "steelblue") +
    geom_line(aes(y = mean), linetype = "dashed", linewidth = 0.8, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
    labs(
      title = title,
      subtitle = paste("Plan:", summary_obj$plan),
      x = "Year (1-6)",
      y = "Profit ($/acre)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
}

plot_cum_profit_bands_region_zone <- function(summary_obj, 
                                              title = NULL) {
  # Plot cumulative profit distribution with uncertainty bands
  
  if (is.null(title)) {
    title <- sprintf("Cumulative Profit Distribution\n%s IL - %s Productivity",
                     tools::toTitleCase(summary_obj$region),
                     tools::toTitleCase(summary_obj$productivity))
  }
  
  ggplot(summary_obj$cum_profit, aes(x = year)) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.25, fill = "darkgreen") +
    geom_line(aes(y = p_mid), linewidth = 1, color = "darkgreen") +
    geom_line(aes(y = mean), linetype = "dashed", linewidth = 0.8, color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
    labs(
      title = title,
      subtitle = paste("Plan:", summary_obj$plan),
      x = "Year (1-6)",
      y = "Cumulative Profit ($/acre)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
}

plot_cashflow_comparison <- function(comp_result) {
  # Compare annual profits between two plans side by side
  
  # Combine data from both plans
  dt <- rbind(
    comp_result$A$summary$ann_profit[, plan := comp_result$A$summary$plan],
    comp_result$B$summary$ann_profit[, plan := comp_result$B$summary$plan]
  )
  
  ggplot(dt, aes(x = year, color = plan, fill = plan)) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.2, color = NA) +
    geom_line(aes(y = mean), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    labs(
      title = sprintf("Annual Profit Comparison\n%s IL - %s Productivity",
                      tools::toTitleCase(comp_result$region),
                      tools::toTitleCase(comp_result$productivity)),
      x = "Year (1-6)",
      y = "Mean Annual Profit ($/acre)",
      color = "Plan",
      fill = "Plan"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
}

# ============================================================
# GRID CASHFLOW ANALYSIS
# ============================================================

compare_cashflows_across_grid <- function(history6,
                                         planA,
                                         planB,
                                         B = 20000L,
                                         verbose = FALSE) {
  # Run cashflow comparison across all 9 region-zone combinations
  
  cat("\n=== CASHFLOW COMPARISON ACROSS GRID ===\n")
  cat("From:", planA, "To:", planB, "\n")
  cat("Starting history:", history6, "\n\n")
  
  results <- list()
  
  for (reg in c("northern", "central", "southern")) {
    results[[reg]] <- list()
    
    for (zone in c("high", "medium", "low")) {
      cat(sprintf("Running %s-%s... ", reg, zone))
      
      tryCatch({
        results[[reg]][[zone]] <- compare_two_plans_cashflow_region_zone(
          history6 = history6,
          planA = planA,
          planB = planB,
          geo_region = reg,
          prod_zone = zone,
          B = B,
          verbose = FALSE
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
        comp <- results[[reg]][[zone]]
        
        summary_list[[paste0(reg, "_", zone)]] <- data.table(
          region = reg,
          productivity = zone,
          plan = c(planA, planB),
          rbind(
            comp$A$summary$liquidity,
            comp$B$summary$liquidity
          )
        )
      }
    }
  }
  
  full_summary <- rbindlist(summary_list)
  
  cat("\n=== LIQUIDITY METRICS SUMMARY ===\n")
  print(full_summary)
  
  list(
    results = results,
    liquidity_summary = full_summary
  )
}

# ============================================================
# EXAMPLE USAGE
# ============================================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Example 1: Single region-zone cashflow analysis
  cf_result <- simulate_plan_cashflows_region_zone(
    history6 = "CCCCCC",
    plan6 = "CSCSCS",
    geo_region = "central",
    prod_zone = "medium",
    B = 20000
  )
  
  # Summarize
  cf_summary <- summarize_cashflows_region_zone(cf_result)
  
  # View annual profits
  print(cf_summary$ann_profit)
  
  # View liquidity metrics
  print(cf_summary$liquidity)
  
  # Visualize
  plot_annual_profit_bands_region_zone(cf_summary)
  plot_cum_profit_bands_region_zone(cf_summary)
  
  # Example 2: Compare two plans
  comparison <- compare_two_plans_cashflow_region_zone(
    history6 = "CCCCCC",
    planA = "CCCCCC",
    planB = "CSCSCS",
    geo_region = "central",
    prod_zone = "high",
    B = 20000
  )
  
  # View delta NPV
  print(comparison$delta$summary)
  
  # Compare liquidity
  print(comparison$A$summary$liquidity)
  print(comparison$B$summary$liquidity)
  
  # Visualize comparison
  plot_cashflow_comparison(comparison)
  
  # Example 3: Full grid cashflow analysis
  grid_cf <- compare_cashflows_across_grid(
    history6 = "CCCCCC",
    planA = "CCCCCC",
    planB = "CSCSCS",
    B = 20000
  )
  
  # Access results by region-zone
  grid_cf$results$central$medium$delta$summary
}
