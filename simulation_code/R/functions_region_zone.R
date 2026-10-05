# ============================================================================
# ADAPTED FUNCTIONS FOR REGION-ZONE FRAMEWORK
# ============================================================================
# These functions are adapted from functions.R to work with the comprehensive
# region-zone framework (3 regions × 3 productivity zones)
# ============================================================================

# ---- packages ------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table, tidyverse, fixest, arrow, progressr, statar, MASS, mvtnorm
)

# ============================================================
# HELPER FUNCTIONS (unchanged from original)
# ============================================================

## Discount vector
disc_vec <- function(r) 1 / (1 + r)^(0:5)

## Draws from a lognormal distribution
draw_lognorm_6y <- function(B, mu_log, sd_log) {
  matrix(rlnorm(B * 6, meanlog = mu_log, sdlog = sd_log), nrow = B, ncol = 6)
}

## Attach simulate function to generator (if needed)
attach_simulate <- function(gen) {
  data.table::setDT(gen$mu)
  data.table::setDT(gen$pools)
  
  simulate <- function(state, B = 5000L) {
    mu_s <- gen$mu[rot6 == state, mu]
    if (length(mu_s) != 1 || is.na(mu_s)) stop("State not found in mu: ", state)
    
    pool <- gen$pools[rot6 == state, resid_pool][[1]]
    if (length(pool) == 0) stop("Empty residual pool for state: ", state)
    
    mu_s + sample(pool, B, replace = TRUE)
  }
  
  gen$simulate <- simulate
  gen
}

# ============================================================
# TRANSITION PATH BUILDER (unchanged)
# ============================================================
make_transition_path <- function(history6, plan6) {
  stopifnot(nchar(history6) == 6, nchar(plan6) == 6)
  h <- strsplit(history6, "")[[1]]
  p <- strsplit(plan6, "")[[1]]
  
  out <- data.table(t = 1:6, crop = p, window6 = NA_character_)
  for (t in 1:6) {
    seq_now <- c(h, p[1:t])
    win <- tail(seq_now, 6)
    out[t, window6 := paste0(win, collapse="")]
  }
  out
}

# ============================================================
# DIAGNOSTIC: MEAN TERMINAL YIELDS FOR KEY STATES
# Adapted for region-zone generators
# ============================================================
mean_terminal_yield_region_zone <- function(gen, states) {

  # Get mean yields for specific rotation states
  
  #Args:
  #  gen: A region-zone-specific generator (e.g., corn_generators$central$high)
  #  states: Vector of rotation state strings (e.g., c("CCCCCC", "SCSCSC"))
  
  #Returns:
  #  Data.table with states and their mean yields
  
  data.table(
    state = states,
    region = gen$geo_region,
    productivity = gen$prod_zone,
    mu = sapply(states, function(st) {
      mu_val <- gen$mu[rot6 == st, mu]
      if (length(mu_val) == 0) NA_real_ else mu_val[1]
    })
  )
}

# ============================================================
# STATES USED IN A PLAN
# ============================================================
states_used_in_plan <- function(history6, plan6) {
  # Identify which rotation states are used in a transition plan
  
  path <- make_transition_path(history6, plan6)
  list(
    path = path,
    corn_states = unique(path[crop == "C", window6]),
    soy_states  = unique(path[crop == "S", window6])
  )
}

# ============================================================
# SUMMARY STATISTICS FUNCTION
# ============================================================
summarize_dist <- function(x, label = NULL) {
  # Calculate summary statistics for a distribution
  
  out <- tibble::tibble(
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    cv   = sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE)),
    p05  = as.numeric(quantile(x, 0.05, na.rm = TRUE)),
    p25  = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
    p50  = as.numeric(quantile(x, 0.50, na.rm = TRUE)),
    p75  = as.numeric(quantile(x, 0.75, na.rm = TRUE)),
    p95  = as.numeric(quantile(x, 0.95, na.rm = TRUE)),
    pr_positive = mean(x > 0, na.rm = TRUE)
  )
  
  if (!is.null(label)) {
    out <- cbind(label = label, out)
  }
  
  out
}

# ============================================================
# STABLE (STEADY-STATE) ROTATION COMPARISON
# Adapted for region-zone framework
# ============================================================
stable_compare_region_zone <- function(planA,
                                       planB,
                                       geo_region = "central",
                                       prod_zone = "medium",
                                       B = 20000L,
                                       r_disc = 0.05,
                                       gens_by_region_zone = yield_gens_by_region_zone,
                                       price_pars_regional = price_params_regional,
                                       base_cost_pars = base_cost_params,
                                       use_rotation_costs = TRUE,
                                       seedA = 501,
                                       seedB = 502,
                                       verbose = FALSE) {
  
  # Compare two stable rotation strategies for a specific region-zone
  
  # Assumption: history6 = plan6 (steady-state rotation already established)
  
  # Args:
  #  planA, planB: 6-character rotation strings
  #  geo_region: "northern", "central", or "southern"
  #  prod_zone: "high", "medium", or "low"
  
  stopifnot(nchar(planA) == 6, nchar(planB) == 6)
  
  outA <- simulate_plan_pv_npv_region_zone(
    history6   = planA,
    plan6      = planA,
    geo_region = geo_region,
    prod_zone  = prod_zone,
    B          = B,
    r_disc     = r_disc,
    gens_by_region_zone = gens_by_region_zone,
    price_pars_regional = price_pars_regional,
    base_cost_pars = base_cost_pars,
    use_rotation_costs = use_rotation_costs,
    seed       = seedA,
    verbose    = verbose
  )
  
  outB <- simulate_plan_pv_npv_region_zone(
    history6   = planB,
    plan6      = planB,
    geo_region = geo_region,
    prod_zone  = prod_zone,
    B          = B,
    r_disc     = r_disc,
    gens_by_region_zone = gens_by_region_zone,
    price_pars_regional = price_pars_regional,
    base_cost_pars = base_cost_pars,
    use_rotation_costs = use_rotation_costs,
    seed       = seedB,
    verbose    = verbose
  )
  
  list(
    A = outA,
    B = outB,
    region = geo_region,
    productivity = prod_zone,
    delta_pv  = outB$pv_revenue - outA$pv_revenue,
    delta_npv = outB$npv_profit - outA$npv_profit
  )
}

# ============================================================
# TRANSITION DELTA COMPARISON
# Adapted for region-zone framework
# ============================================================
simulate_transition_delta_region_zone <- function(history6,
                                                  plan_from,
                                                  plan_to,
                                                  geo_region = "central",
                                                  prod_zone = "medium",
                                                  B = 20000L,
                                                  r_disc = 0.05,
                                                  gens_by_region_zone = yield_gens_by_region_zone,
                                                  price_pars_regional = price_params_regional,
                                                  base_cost_pars = base_cost_params,
                                                  use_rotation_costs = TRUE,
                                                  seed_from = 123,
                                                  seed_to   = 124,
                                                  verbose = TRUE) {

  # Compare two rotation plans from the same starting history
  # for a specific region-zone combination
  
  out_from <- simulate_plan_pv_npv_region_zone(
    history6, plan_from, B, r_disc,
    gens_by_region_zone, price_pars_regional, base_cost_pars,
    geo_region, prod_zone, use_rotation_costs,
    seed = seed_from, verbose = verbose
  )
  
  out_to <- simulate_plan_pv_npv_region_zone(
    history6, plan_to, B, r_disc,
    gens_by_region_zone, price_pars_regional, base_cost_pars,
    geo_region, prod_zone, use_rotation_costs,
    seed = seed_to, verbose = verbose
  )
  
  list(
    from = out_from,
    to = out_to,
    region = geo_region,
    productivity = prod_zone,
    delta_pv  = out_to$pv_revenue - out_from$pv_revenue,
    delta_npv = out_to$npv_profit - out_from$npv_profit
  )
}

# ============================================================
# GRID COMPARISON: ALL REGION-ZONE COMBINATIONS
# ============================================================
compare_across_grid <- function(history6,
                                plan_from,
                                plan_to,
                                B = 20000L,
                                r_disc = 0.05,
                                gens_by_region_zone = yield_gens_by_region_zone,
                                price_pars_regional = price_params_regional,
                                base_cost_pars = base_cost_params,
                                use_rotation_costs = TRUE,
                                verbose = FALSE) {
  
  # Run comparison across all 9 region-zone combinations
  
  # Returns:
  #  List of results by region and zone, plus summary tables
  
  cat("\n=== COMPARING ACROSS ALL REGION × ZONE COMBINATIONS ===\n")
  cat("From:", plan_from, "To:", plan_to, "\n")
  cat("Starting history:", history6, "\n\n")
  
  results <- list()
  
  for (reg in c("northern", "central", "southern")) {
    results[[reg]] <- list()
    
    for (zone in c("high", "medium", "low")) {
      cat(sprintf("Running %s-%s... ", reg, zone))
      
      tryCatch({
        results[[reg]][[zone]] <- simulate_transition_delta_region_zone(
          history6 = history6,
          plan_from = plan_from,
          plan_to = plan_to,
          geo_region = reg,
          prod_zone = zone,
          B = B,
          r_disc = r_disc,
          gens_by_region_zone = gens_by_region_zone,
          price_pars_regional = price_pars_regional,
          base_cost_pars = base_cost_pars,
          use_rotation_costs = use_rotation_costs,
          seed_from = 123,
          seed_to = 124,
          verbose = FALSE
        )
        cat("✓\n")
      }, error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        results[[reg]][[zone]] <- NULL
      })
    }
  }
  
  # Create summary tables
  summary_list <- list()
  for (reg in names(results)) {
    for (zone in names(results[[reg]])) {
      if (!is.null(results[[reg]][[zone]])) {
        summary_list[[paste0(reg, "_", zone)]] <- data.table(
          region = reg,
          productivity = zone,
          plan = c(plan_from, plan_to),
          rbind(
            summarize_dist(results[[reg]][[zone]]$from$npv_profit),
            summarize_dist(results[[reg]][[zone]]$to$npv_profit)
          )
        )
      }
    }
  }
  
  full_summary <- rbindlist(summary_list)
  
  # Delta summary
  delta_list <- list()
  for (reg in names(results)) {
    for (zone in names(results[[reg]])) {
      if (!is.null(results[[reg]][[zone]])) {
        delta_list[[paste0(reg, "_", zone)]] <- data.table(
          region = reg,
          productivity = zone,
          summarize_dist(results[[reg]][[zone]]$delta_npv)
        )
      }
    }
  }
  
  delta_summary <- rbindlist(delta_list)
  
  cat("\n=== SUMMARY: NPV BY REGION × PRODUCTIVITY × PLAN ===\n")
  print(full_summary)
  
  cat("\n=== SUMMARY: DELTA NPV (", plan_to, " - ", plan_from, ") ===\n", sep="")
  print(delta_summary)
  
  list(
    results = results,
    npv_summary = full_summary,
    delta_summary = delta_summary
  )
}

# ============================================================
# VISUALIZATION FUNCTIONS
# ============================================================

plot_npv_grid <- function(grid_results, metric = "mean") {
  
  # Plot NPV or delta NPV across the 3×3 grid
  
  # Args:
  # grid_results: Output from compare_across_grid()
  # metric: "mean", "median", "sd", or "pr_positive"
  
  # Get delta summary
  dt <- grid_results$delta_summary
  
  # Reorder factors for plotting
  dt[, region := factor(region, levels = c("northern", "central", "southern"))]
  dt[, productivity := factor(productivity, levels = c("low", "medium", "high"))]
  
  # Select metric
  y_var <- switch(metric,
                  "mean" = "mean",
                  "median" = "p50",
                  "sd" = "sd",
                  "pr_positive" = "pr_positive")
  
  y_label <- switch(metric,
                    "mean" = "Mean Delta NPV ($/acre)",
                    "median" = "Median Delta NPV ($/acre)",
                    "sd" = "Std Dev of Delta NPV ($/acre)",
                    "pr_positive" = "Pr(Delta NPV > 0)")
  
  ggplot(dt, aes(x = productivity, y = get(y_var), fill = region)) +
    geom_col(position = position_dodge(width = 0.8), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ region, ncol = 3) +
    labs(
      title = "Rotation Benefits Across Illinois",
      subtitle = paste("Metric:", metric),
      x = "Productivity Zone",
      y = y_label,
      fill = "Region"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_npv_heatmap <- function(grid_results, metric = "mean") {
  
  # Create heatmap of NPV metric across region × productivity grid
  
  dt <- grid_results$delta_summary
  
  # Reorder factors
  dt[, region := factor(region, levels = c("southern", "central", "northern"))]
  dt[, productivity := factor(productivity, levels = c("low", "medium", "high"))]
  
  # Select metric
  y_var <- switch(metric,
                  "mean" = "mean",
                  "median" = "p50",
                  "sd" = "sd",
                  "pr_positive" = "pr_positive")
  
  title_text <- switch(metric,
                       "mean" = "Mean Delta NPV ($/acre)",
                       "median" = "Median Delta NPV ($/acre)",
                       "sd" = "Risk (Std Dev)",
                       "pr_positive" = "Probability of Profit")
  
  ggplot(dt, aes(x = productivity, y = region, fill = get(y_var))) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = sprintf("%.0f", get(y_var))), 
              color = "white", size = 5, fontface = "bold") +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      title = title_text,
      subtitle = "Across Illinois Region × Productivity Grid",
      x = "Productivity Zone",
      y = "Region",
      fill = metric
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold")
    )
}

plot_npv_distributions_by_region <- function(grid_results, 
                                             regions = c("northern", "central", "southern"),
                                             prod_zone = "medium") {
  
  # Compare NPV distributions across regions for a given productivity zone
  
  # Extract distributions
  dist_data <- rbindlist(lapply(regions, function(reg) {
    res <- grid_results$results[[reg]][[prod_zone]]
    if (!is.null(res)) {
      data.table(
        region = reg,
        plan = rep(c("from", "to"), each = length(res$from$npv_profit)),
        npv = c(res$from$npv_profit, res$to$npv_profit)
      )
    }
  }))
  
  ggplot(dist_data, aes(x = npv, fill = plan)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ region, ncol = 1, scales = "free_y") +
    labs(
      title = paste("NPV Distributions by Region -", prod_zone, "Productivity"),
      x = "NPV ($/acre)",
      y = "Density",
      fill = "Plan"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# ============================================================
# BATCH ANALYSIS HELPER
# ============================================================

analyze_multiple_scenarios <- function(scenarios, 
                                      B = 20000,
                                      verbose = FALSE) {
  # Run multiple scenario comparisons
  
  # Args:
    # scenarios: List of lists, each containing:
    #  - history6
    #  - plan_from
    #  - plan_to
    #  - label (optional)
  
  # Returns:
  #  Combined results from all scenarios
  
  all_results <- list()
  
  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    label <- if (!is.null(scenario$label)) scenario$label else paste("Scenario", i)
    
    cat("\n")
    cat("="*60, "\n", sep="")
    cat("Running:", label, "\n")
    cat("="*60, "\n", sep="")
    
    result <- compare_across_grid(
      history6 = scenario$history6,
      plan_from = scenario$plan_from,
      plan_to = scenario$plan_to,
      B = B,
      verbose = verbose
    )
    
    result$scenario_label <- label
    all_results[[label]] <- result
  }
  
  all_results
}

# ============================================================
# EXAMPLE USAGE
# ============================================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Example 1: Single region-zone comparison
  result_central_med <- simulate_transition_delta_region_zone(
    history6 = "CCCCCC",
    plan_from = "CCCCCC",
    plan_to = "CSCSCS",
    geo_region = "central",
    prod_zone = "medium",
    B = 20000
  )
  
  # Summary statistics
  cat("\nFrom (Continuous Corn):\n")
  print(summarize_dist(result_central_med$from$npv_profit))
  
  cat("\nTo (Corn-Soy Rotation):\n")
  print(summarize_dist(result_central_med$to$npv_profit))
  
  cat("\nDelta (To - From):\n")
  print(summarize_dist(result_central_med$delta_npv))
  
  # Example 2: Full grid comparison
  grid_results <- compare_across_grid(
    history6 = "CCCCCC",
    plan_from = "CCCCCC",
    plan_to = "CSCSCS",
    B = 20000
  )
  
  # Visualizations
  plot_npv_grid(grid_results, metric = "mean")
  plot_npv_heatmap(grid_results, metric = "mean")
  plot_npv_distributions_by_region(grid_results, prod_zone = "medium")
  
  # Example 3: Multiple scenarios
  scenarios <- list(
    list(
      history6 = "CCCCCC",
      plan_from = "CCCCCC",
      plan_to = "CSCSCS",
      label = "Transition from Monoculture"
    ),
    list(
      history6 = "CSCSCS",
      plan_from = "CSCSCS",
      plan_to = "CCCCCC",
      label = "Transition to Monoculture"
    )
  )
  
  multi_results <- analyze_multiple_scenarios(scenarios, B = 20000)
}
