# ============================================================================
# ADAPTED PRICE-YIELD CORRELATION FUNCTIONS FOR REGION-ZONE FRAMEWORK
# ============================================================================
# Three methods for imposing negative price-yield correlation
# ============================================================================

# ---- Required packages ----------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, mvtnorm, fitdistrplus)

# ============================================================================
# METHOD 1: SIMPLE CORRELATED NORMAL APPROACH
# ============================================================================

draw_correlated_price_yield_simple <- function(B, 
                                               T = 6,
                                               price_mean,
                                               price_sd,
                                               yield_mean,
                                               yield_sd,
                                               correlation = -0.4) {
 
  # Draw correlated price and yield using bivariate normal
  
  # Simple but effective approach for imposing correlation
  
  # Correlation matrix
  cor_matrix <- matrix(c(1, correlation,
                         correlation, 1), 
                       nrow = 2, ncol = 2)
  
  # Covariance matrix
  sds <- c(price_sd, yield_sd)
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
  
  # Draw correlated normal random variables
  draws <- matrix(NA_real_, nrow = B * T, ncol = 2)
  
  for (i in 1:(B * T)) {
    draws[i, ] <- mvtnorm::rmvnorm(1, 
                                   mean = c(0, 0),
                                   sigma = cov_matrix)
  }
  
  # Transform to price and yield
  # Price: lognormal
  log_price <- price_mean + draws[, 1]
  price_mat <- matrix(exp(log_price), nrow = B, ncol = T)
  
  # Yield: normal (truncate at zero)
  yield_mat <- matrix(pmax(0, yield_mean + draws[, 2]), nrow = B, ncol = T)
  
  list(
    prices = price_mat,
    yields = yield_mat,
    correlation = cor(draws[, 1], draws[, 2])
  )
}

# ============================================================================
# METHOD 2: AGGREGATE + IDIOSYNCRATIC DECOMPOSITION (More Realistic)
# ============================================================================

draw_correlated_price_yield_decomposed <- function(B,
                                                   T = 6,
                                                   price_params,
                                                   yield_mean,
                                                   yield_sd_total,
                                                   price_yield_correlation = -0.45,
                                                   aggregate_fraction = 0.4) {
  
  # More realistic model decomposing yield variance
  
  # Yield = mean + aggregate_shock + idiosyncratic_shock
  # Price = f(aggregate_shock)  [correlated only with aggregate]
  
  # Captures the economic reality:
  # - Regional/national weather shocks affect price (everyone has high/low yields)
  # - Farm-specific shocks don't affect market price (too small)
  
  # Decompose yield variance
  var_total <- yield_sd_total^2
  var_aggregate <- var_total * aggregate_fraction
  var_idiosyncratic <- var_total * (1 - aggregate_fraction)
  
  sd_aggregate <- sqrt(var_aggregate)
  sd_idiosyncratic <- sqrt(var_idiosyncratic)
  
  # Correlation between log(price) and aggregate yield shock
  cor_matrix <- matrix(c(1, price_yield_correlation,
                         price_yield_correlation, 1),
                       nrow = 2, ncol = 2)
  
  # Covariance matrix
  sds <- c(price_params$sigma, sd_aggregate)
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
  
  # Draw correlated shocks (efficient vectorized approach)
  all_shocks <- mvtnorm::rmvnorm(B * T, mean = c(0, 0), sigma = cov_matrix)
  
  # Reshape into matrices
  log_price_shocks      <- matrix(all_shocks[, 1], nrow = B, ncol = T)
  aggregate_yield_shock <- matrix(all_shocks[, 2], nrow = B, ncol = T)
  
  # Price: lognormal with correlated shocks
  price_mat <- exp(price_params$mu + log_price_shocks)
  
  # Yield: aggregate shock (correlated) + idiosyncratic (uncorrelated)
  idiosyncratic_shock <- matrix(rnorm(B * T, 0, sd_idiosyncratic),
                                nrow = B, ncol = T)
  yield_mat <- yield_mean + aggregate_yield_shock + idiosyncratic_shock
  
  # Ensure non-negative yields
  yield_mat <- pmax(yield_mat, 0)
  
  list(
    prices = price_mat,
    yields = yield_mat,
    aggregate_shocks = aggregate_yield_shock,
    correlation_realized = cor(as.vector(log_price_shocks), 
                              as.vector(aggregate_yield_shock))
  )
}

# ============================================================================
# METHOD 3: COPULA APPROACH (Most Flexible)
# ============================================================================

draw_correlated_copula <- function(B,
                                   T = 6,
                                   price_params,
                                   yield_gen,
                                   state,
                                   correlation = -0.45) {
  # Draw price and yield using Gaussian copula
  
  # Preserves marginal distributions while imposing correlation structure
  # Most flexible but computationally intensive
  
  price_mat <- matrix(NA_real_, nrow = B, ncol = T)
  yield_mat <- matrix(NA_real_, nrow = B, ncol = T)
  
  # Correlation matrix
  cor_matrix <- matrix(c(1, correlation,
                         correlation, 1),
                       nrow = 2, ncol = 2)
  
  for (t in 1:T) {
    # Draw correlated standard normal pairs
    z <- mvtnorm::rmvnorm(B, mean = c(0, 0), sigma = cor_matrix)
    
    # Transform to uniform [0,1] using normal CDF (copula step)
    u <- pnorm(z)
    
    # Price marginal: lognormal
    price_mat[, t] <- qlnorm(u[, 1],
                             meanlog = price_params$mu,
                             sdlog = price_params$sigma)
    
    # Yield marginal: empirical from generator
    raw_yields <- yield_gen$simulate(state, B)
    yield_mat[, t] <- qlnorm(u[, 2], 
                             mean = mean(raw_yields), 
                             sd = sd(raw_yields))
  }
  
  list(
    prices = price_mat,
    yields = yield_mat,
    correlation_check = cor(as.vector(price_mat), as.vector(yield_mat))
  )
}

# ============================================================================
# METHOD 4: t-COPULA (Tail Dependence)
# ============================================================================

draw_correlated_price_yield_tcopula <- function(B,
                                                T = 6,
                                                price_params,
                                                yield_mean,
                                                yield_sd_total,
                                                price_yield_correlation = -0.45,
                                                aggregate_fraction = 0.4,
                                                df = 5) {
  # Same decomposed structure as Method 2 but with a Student-t copula.
  # Marginal distributions are identical to the decomposed Gaussian method;
  # the t-copula adds positive tail dependence — joint extremes (e.g. low
  # price AND low aggregate yield) occur more often than under Gaussian.
  # Effect is strongest at low df (df=5 is a common robustness check).

  sd_aggregate      <- sqrt(yield_sd_total^2 * aggregate_fraction)
  sd_idiosyncratic  <- sqrt(yield_sd_total^2 * (1 - aggregate_fraction))

  cor_matrix <- matrix(c(1, price_yield_correlation,
                         price_yield_correlation, 1), nrow = 2)

  # Draw correlated t variates; sigma = correlation matrix gives standard-t marginals
  t_shocks <- mvtnorm::rmvt(B * T, sigma = cor_matrix, df = df, delta = c(0, 0))

  # Copula step: t-CDF maps each marginal to Uniform(0,1)
  u <- pt(t_shocks, df = df)

  # Apply target marginals via normal quantile function then rescale
  log_price_shocks      <- matrix(qnorm(u[, 1]) * price_params$sigma,
                                  nrow = B, ncol = T)
  aggregate_yield_shock <- matrix(qnorm(u[, 2]) * sd_aggregate,
                                  nrow = B, ncol = T)

  price_mat <- exp(price_params$mu + log_price_shocks)

  idiosyncratic_shock <- matrix(rnorm(B * T, 0, sd_idiosyncratic),
                                nrow = B, ncol = T)
  yield_mat <- pmax(yield_mean + aggregate_yield_shock + idiosyncratic_shock, 0)

  list(
    prices               = price_mat,
    yields               = yield_mat,
    correlation_realized = cor(as.vector(log_price_shocks),
                               as.vector(aggregate_yield_shock))
  )
}

# ============================================================================
# MAIN FUNCTION: SIMULATE NPV WITH PRICE-YIELD CORRELATION (REGION-ZONE)
# ============================================================================

simulate_npv_correlated_region_zone <- function(history6,
                                                plan6,
                                                geo_region = "central",
                                                prod_zone = "medium",
                                                B = 20000L,
                                                r_disc = 0.05,
                                                gens_by_region_zone = yield_gens_by_region_zone,
                                                price_pars_regional = price_params_regional,
                                                base_cost_pars = base_cost_params,
                                                correlation_corn = -0.40,
                                                correlation_soy = -0.50,
                                                aggregate_fraction = 0.40,
                                                method = c("decomposed", "simple", "copula"),
                                                use_rotation_costs = TRUE,
                                                seed = 123,
                                                verbose = TRUE) {
 
  # Simulate NPV with negative price-yield correlation for specific region-zone
  
  # Args:
  #  method: "decomposed" (recommended), "simple", or "copula"
  #  correlation_corn: Negative correlation for corn (-0.3 to -0.5 typical)
  #  correlation_soy: Negative correlation for soybeans
  #  aggregate_fraction: Fraction of yield variance from aggregate shocks (0.3-0.5)
  
  # Returns:
  #  Full simulation results including correlation diagnostics
  
  method <- match.arg(method)
  set.seed(seed)
  stopifnot(nchar(history6) == 6, nchar(plan6) == 6)
  
  # Select region-specific prices
  if (geo_region %in% names(price_pars_regional)) {
    price_pars <- price_pars_regional[[geo_region]]
  } else {
    stop("Region not found: ", geo_region)
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
    stop("Yield generators not available: ", geo_region, "-", prod_zone)
  }
  
  path <- make_transition_path(history6, plan6)
  
  if (verbose) {
    cat("\n=== Correlated Price-Yield Simulation ===\n")
    cat("Region:", geo_region, "  Productivity:", prod_zone, "\n")
    cat("Method:", method, "\n")
    cat("Corn correlation:", correlation_corn, "\n")
    cat("Soy correlation:", correlation_soy, "\n")
    if (method == "decomposed") {
      cat("Aggregate fraction:", aggregate_fraction, "\n")
    }
    cat("\n")
  }
  
  disc <- disc_vec(r_disc)
  
  # Identify which years carry each crop
  chars     <- strsplit(plan6, "")[[1]]
  corn_yrs  <- which(chars == "C")
  soy_yrs   <- which(chars == "S")
  n_corn    <- length(corn_yrs)
  n_soy     <- length(soy_yrs)
  
  # Initialize output matrices
  Pc <- matrix(0, nrow = B, ncol = 6)
  Ps <- matrix(0, nrow = B, ncol = 6)
  Yc <- matrix(0, nrow = B, ncol = 6)
  Ys <- matrix(0, nrow = B, ncol = 6)
  
  # ---- CORN years ----
  if (n_corn > 0) {
    
    # Collect yield parameters
    mu_corn <- numeric(n_corn)
    sd_corn <- numeric(n_corn)
    for (i in seq_along(corn_yrs)) {
      st <- path$window6[corn_yrs[i]]
      mu_corn[i] <- corn_gen$mu[rot6 == st, mu]
      pool       <- corn_gen$pools[rot6 == st, resid_pool][[1]]
      sd_corn[i] <- sd(pool)
    }
    
    if (method == "decomposed") {
      for (i in seq_along(corn_yrs)) {
        t <- corn_yrs[i]
        res <- draw_correlated_price_yield_decomposed(
          B = B, T = 1,
          price_params            = price_pars$C,
          yield_mean              = mu_corn[i],
          yield_sd_total          = sd_corn[i],
          price_yield_correlation = correlation_corn,
          aggregate_fraction      = aggregate_fraction
        )
        Pc[, t] <- as.vector(res$prices)
        Yc[, t] <- as.vector(res$yields)
      }
      
    } else if (method == "simple") {
      for (i in seq_along(corn_yrs)) {
        t <- corn_yrs[i]
        res <- draw_correlated_price_yield_simple(
          B = B, T = 1,
          price_mean = price_pars$C$mu,
          price_sd   = price_pars$C$sigma,
          yield_mean = mu_corn[i],
          yield_sd   = sd_corn[i],
          correlation = correlation_corn
        )
        Pc[, t] <- as.vector(res$prices)
        Yc[, t] <- as.vector(res$yields)
      }
      
    } else if (method == "copula") {
      for (i in seq_along(corn_yrs)) {
        t <- corn_yrs[i]
        st  <- path$window6[t]
        res <- draw_correlated_copula(
          B = B, T = 1,
          price_params = price_pars$C,
          yield_gen    = corn_gen,
          state        = st,
          correlation  = correlation_corn
        )
        Pc[, t] <- as.vector(res$prices)
        Yc[, t] <- as.vector(res$yields)
      }
    }
  }
  
  # ---- SOY years ----
  if (n_soy > 0) {
    
    # Collect yield parameters
    mu_soy <- numeric(n_soy)
    sd_soy <- numeric(n_soy)
    for (i in seq_along(soy_yrs)) {
      st <- path$window6[soy_yrs[i]]
      mu_soy[i] <- soy_gen$mu[rot6 == st, mu]
      pool       <- soy_gen$pools[rot6 == st, resid_pool][[1]]
      sd_soy[i]  <- sd(pool)
    }
    
    if (method == "decomposed") {
      for (i in seq_along(soy_yrs)) {
        t <- soy_yrs[i]
        res <- draw_correlated_price_yield_decomposed(
          B = B, T = 1,
          price_params            = price_pars$S,
          yield_mean              = mu_soy[i],
          yield_sd_total          = sd_soy[i],
          price_yield_correlation = correlation_soy,
          aggregate_fraction      = aggregate_fraction
        )
        Ps[, t] <- as.vector(res$prices)
        Ys[, t] <- as.vector(res$yields)
      }
      
    } else if (method == "simple") {
      for (i in seq_along(soy_yrs)) {
        t <- soy_yrs[i]
        res <- draw_correlated_price_yield_simple(
          B = B, T = 1,
          price_mean = price_pars$S$mu,
          price_sd   = price_pars$S$sigma,
          yield_mean = mu_soy[i],
          yield_sd   = sd_soy[i],
          correlation = correlation_soy
        )
        Ps[, t] <- as.vector(res$prices)
        Ys[, t] <- as.vector(res$yields)
      }
      
    } else if (method == "copula") {
      for (i in seq_along(soy_yrs)) {
        t <- soy_yrs[i]
        st  <- path$window6[t]
        res <- draw_correlated_copula(
          B = B, T = 1,
          price_params = price_pars$S,
          yield_gen    = soy_gen,
          state        = st,
          correlation  = correlation_soy
        )
        Ps[, t] <- as.vector(res$prices)
        Ys[, t] <- as.vector(res$yields)
      }
    }
  }
  
  # ---- Costs: drawn independently ----
  Cc <- matrix(rlnorm(B * 6, cost_pars_base$C$mu, cost_pars_base$C$sigma), 
               nrow = B, ncol = 6)
  Cs <- matrix(rlnorm(B * 6, cost_pars_base$S$mu, cost_pars_base$S$sigma), 
               nrow = B, ncol = 6)
  
  # ---- Combine into revenue/cost/profit ----
  IC <- matrix(as.numeric(chars == "C"), nrow = B, ncol = 6, byrow = TRUE)
  IS <- matrix(as.numeric(chars == "S"), nrow = B, ncol = 6, byrow = TRUE)
  
  Pmat <- (Pc * IC) + (Ps * IS)
  Ymat <- (Yc * IC) + (Ys * IS)
  Kmat <- (Cc * IC) + (Cs * IS)
  
  rev_mat  <- Pmat * Ymat
  cost_mat <- Kmat
  prof_mat <- rev_mat - cost_mat
  
  pv_rev  <- as.numeric(rev_mat  %*% disc)
  pv_cost <- as.numeric(cost_mat %*% disc)
  npv     <- as.numeric(prof_mat %*% disc)
  
  # ---- Realized correlations (diagnostic) ----
  realized_cor_corn <- if (n_corn > 0)
    cor(as.vector(Pc[, corn_yrs, drop = FALSE]),
        as.vector(Yc[, corn_yrs, drop = FALSE])) else NA
  
  realized_cor_soy <- if (n_soy > 0)
    cor(as.vector(Ps[, soy_yrs, drop = FALSE]),
        as.vector(Ys[, soy_yrs, drop = FALSE])) else NA
  
  if (verbose) {
    cat("Realized price-yield correlations:\n")
    if (!is.na(realized_cor_corn)) cat("  Corn:", round(realized_cor_corn, 3), "\n")
    if (!is.na(realized_cor_soy))  cat("  Soy: ", round(realized_cor_soy,  3), "\n")
    cat("\nPrice variation check (should all be > 0):\n")
    for (t in 1:6) {
      cat("  Year", t, "(", chars[t], "): ",
          "n_unique =", length(unique(Pmat[,t])),
          " SD =", round(sd(Pmat[,t]), 4), "\n")
    }
  }
  
  list(
    history6 = history6,
    plan6 = plan6,
    region = geo_region,
    productivity = prod_zone,
    path = path,
    prices_corn = Pc,
    prices_soy  = Ps,
    yields_corn = Yc,
    yields_soy  = Ys,
    revenues  = rev_mat,
    costs     = cost_mat,
    profits   = prof_mat,
    pv_revenue = pv_rev,
    pv_cost    = pv_cost,
    npv_profit = npv,
    realized_cor_corn = realized_cor_corn,
    realized_cor_soy  = realized_cor_soy,
    method = method,
    target_cor_corn = correlation_corn,
    target_cor_soy = correlation_soy
  )
}

# ============================================================================
# COMPARISON: INDEPENDENT vs CORRELATED
# ============================================================================

compare_correlation_impact_region_zone <- function(history6,
                                                   plan6,
                                                   geo_region = "central",
                                                   prod_zone = "medium",
                                                   B = 20000L,
                                                   correlation = -0.45,
                                                   method = "decomposed") {
  
  # Compare NPV results with and without price-yield correlation
  
  cat("\n")
  cat("========================================================\n")
  cat("COMPARING INDEPENDENT vs CORRELATED PRICE-YIELD\n")
  cat("Region:", geo_region, "  Productivity:", prod_zone, "\n")
  cat("========================================================\n\n")
  
  # Independent (original method)
  cat("Running INDEPENDENT simulation...\n")
  result_indep <- simulate_plan_pv_npv_region_zone(
    history6 = history6,
    plan6 = plan6,
    geo_region = geo_region,
    prod_zone = prod_zone,
    B = B,
    seed = 111,
    verbose = FALSE
  )
  
  # Correlated
  cat("Running CORRELATED simulation (ρ =", correlation, ")...\n")
  result_corr <- simulate_npv_correlated_region_zone(
    history6 = history6,
    plan6 = plan6,
    geo_region = geo_region,
    prod_zone = prod_zone,
    B = B,
    correlation_corn = correlation,
    correlation_soy = correlation,
    method = method,
    seed = 111,
    verbose = FALSE
  )
  
  # Compare NPV distributions
  cat("\n--- NPV Comparison ---\n")
  
  indep_summary <- data.table(
    model = "Independent",
    region = geo_region,
    productivity = prod_zone,
    mean = mean(result_indep$npv_profit, na.rm = TRUE),
    sd = sd(result_indep$npv_profit, na.rm = TRUE),
    p05 = quantile(result_indep$npv_profit, 0.05, na.rm = TRUE),
    p50 = quantile(result_indep$npv_profit, 0.50, na.rm = TRUE),
    p95 = quantile(result_indep$npv_profit, 0.95, na.rm = TRUE),
    cv = sd(result_indep$npv_profit, na.rm = TRUE) / 
         mean(result_indep$npv_profit, na.rm = TRUE)
  )
  
  corr_summary <- data.table(
    model = "Correlated",
    region = geo_region,
    productivity = prod_zone,
    mean = mean(result_corr$npv_profit, na.rm = TRUE),
    sd = sd(result_corr$npv_profit, na.rm = TRUE),
    p05 = quantile(result_corr$npv_profit, 0.05, na.rm = TRUE),
    p50 = quantile(result_corr$npv_profit, 0.50, na.rm = TRUE),
    p95 = quantile(result_corr$npv_profit, 0.95, na.rm = TRUE),
    cv = sd(result_corr$npv_profit, na.rm = TRUE) / 
         mean(result_corr$npv_profit, na.rm = TRUE)
  )
  
  comparison <- rbind(indep_summary, corr_summary)
  print(comparison)
  
  cat("\n--- Key Insights ---\n")
  cat("Mean NPV change:", 
      sprintf("%.1f", corr_summary$mean - indep_summary$mean), "$/ac\n")
  cat("SD change:", 
      sprintf("%.1f%%", (corr_summary$sd / indep_summary$sd - 1) * 100), "\n")
  cat("Coefficient of variation change:",
      sprintf("%.1f%%", (corr_summary$cv / indep_summary$cv - 1) * 100), "\n\n")
  
  cat("Interpretation:\n")
  if (corr_summary$sd < indep_summary$sd) {
    cat("  ✓ Negative price-yield correlation REDUCES revenue risk\n")
    cat("    (High yields offset by lower prices; low yields by higher prices)\n")
  } else {
    cat("  ✗ Positive correlation would INCREASE revenue risk\n")
  }
  
  list(
    independent = result_indep,
    correlated = result_corr,
    comparison = comparison
  )
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Example 1: Basic correlated simulation
  result <- simulate_npv_correlated_region_zone(
    history6 = "CCCCCC",
    plan6 = "CSCSCS",
    geo_region = "central",
    prod_zone = "high",
    B = 20000,
    correlation_corn = -0.40,
    correlation_soy = -0.50,
    method = "decomposed"
  )
  
  cat("\nMean NPV:", mean(result$npv_profit), "\n")
  cat("SD NPV:", sd(result$npv_profit), "\n")
  cat("Realized corn correlation:", result$realized_cor_corn, "\n")
  cat("Realized soy correlation:", result$realized_cor_soy, "\n")
  
  # Example 2: Compare independent vs correlated
  comparison <- compare_correlation_impact_region_zone(
    history6 = "CCCCCC",
    plan6 = "CSCSCS",
    geo_region = "central",
    prod_zone = "high",
    B = 20000,
    correlation = -0.45,
    method = "decomposed"
  )
  
  # Check risk reduction
  cat("\nRisk reduction from correlation:",
      sprintf("%.1f%%", 
              (1 - comparison$comparison[model == "Correlated", sd] /
                   comparison$comparison[model == "Independent", sd]) * 100),
      "\n")
}
