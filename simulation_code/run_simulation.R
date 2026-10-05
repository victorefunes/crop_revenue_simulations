library(tidyverse)
library(data.table)
library(patchwork)

source("config.R")   # paths, seed, B; run from the simulation_code/ folder

source("R/functions_region_zone.R")
source("R/correlation_functions_region_zone.R")
source("R/cashflow_functions_region_zone.R")
source("R/credit_functions_region_zone.R")


if (file.exists(YGEN_CACHE)) {
  yield_gens_by_region_zone <- readRDS(YGEN_CACHE)
  base_cost_params          <- readRDS(COST_PARAMS_RDS)
  price_params_regional     <- readRDS(PRICE_PARAMS_RDS)
  cat("Yield generators loaded from cache.\n")
} else {
  message("Cache not found — building from source data (slow).")
  message("Tip: run build_yield_generators.R once to create the cache.")
  source("R/yield_generator.R")
}

if (!dir.exists("results/credit")) dir.create("results/credit", recursive = TRUE)

if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive = TRUE)

# ---- Shared color palettes --------------------------------------------------
plan_colors <- setNames(RColorBrewer::brewer.pal(5, "RdYlBu")[c(1, 2, 5)],
                        c("CCCCCC", "CSCSCS", "CCSCCS"))
zone_colors   <- setNames(RColorBrewer::brewer.pal(3, "BrBG"),
                          c("high", "medium", "low"))
region_colors <- setNames(RColorBrewer::brewer.pal(3, "Spectral"),
                          c("northern", "central", "southern"))

# ---- County-level price-yield correlations ----------------------------------
# Source: corn_correlation.csv / soy_correlation.csv (FIPS x corr).
# Aggregated to region level via the IL county→region map in counties_fips.csv.
# Only IL counties (FIPS 17xxx) are mapped; other states in the CSVs are
# ignored. Correlations are median across counties within each region.
# Fallback to DEFAULT_CORR_* wherever a region has no usable county data.
DEFAULT_CORR_CORN <- -0.40
DEFAULT_CORR_SOY  <- -0.50

.corn_raw <- fread(CORN_CORR_CSV,  select = c("fips", "corr")) |>
  filter(fips > 17000 & fips < 18000)
.soy_raw  <- fread(SOY_CORR_CSV,   select = c("fips", "corr")) |>
  filter(fips > 17000 & fips < 18000)
setnames(.corn_raw, "corr", "corr_corn")
setnames(.soy_raw,  "corr", "corr_soy")

.fips_map <- fread(FBFM_COUNTIES_CSV,
                   select = c("fips", "region"))
.fips_map[, region := tolower(region)]   # "Northern" -> "northern"

.county_corr <- Reduce(function(a, b) merge(a, b, by = "fips", all = TRUE),
                       list(.fips_map, .corn_raw, .soy_raw))

region_corr <- .county_corr[
  !is.na(region),
  .(corr_corn = if (all(is.na(corr_corn))) DEFAULT_CORR_CORN
                else median(corr_corn, na.rm = TRUE),
    corr_soy  = if (all(is.na(corr_soy)))  DEFAULT_CORR_SOY
                else median(corr_soy,  na.rm = TRUE)),
  by = region
]
# Ensure all three regions are present; fill with defaults if any are absent
#region_corr <- merge(
#  data.table(region = regions_vec),
#  region_corr, by = "region", all.x = TRUE
#)
region_corr[is.na(corr_corn), corr_corn := DEFAULT_CORR_CORN]
region_corr[is.na(corr_soy),  corr_soy  := DEFAULT_CORR_SOY]
setkey(region_corr, region)

cat("\n=== County-level price-yield correlations aggregated to region ===\n")
print(region_corr)
rm(.corn_raw, .soy_raw, .fips_map, .county_corr)

# Reproducibility: fix the global seed so the entire grid is replicable.
# If simulate_plan_credit_region_zone() accepts a seed argument, pass it
# through inside run_credit_grid() as well.
# NOTE on seeds: simulate_plan_credit_region_zone() calls set.seed() internally
# from its `seed` argument (default 1L). That means a global set.seed() here is
# overridden by the simulator. We instead thread distinct seeds through
# run_credit_grid() so that (a) cells are reproducible, and (b) Monte Carlo
# errors across cells are independent rather than perfectly correlated.
# Within a single cell, plans share a seed deliberately — common random numbers
# reduces variance of the plan-to-plan difference. Across cells we want
# independent draws so heatmap variation is interpretable.
GLOBAL_SEED <- 21212
set.seed(GLOBAL_SEED)

B <- SIM_B   # was 50000; for Pr ~ 0.05-0.20 the MC SE at B=20k
                          # is ~0.003, well below substantive effect sizes.
                          # Bump back up for final paper tables if needed.
rotation_plans <- c("CCCCCC", "CSCSCS", "CCSCCS")
regions_vec    <- c("northern", "central", "southern")
zones_vec      <- c("high", "medium", "low")
# Credit limit sensitivity: applied post-hoc to peak-debt distributions.
# Delinquency (cash-flow shortfall) is estimated once with a non-binding limit;
# credit rationing is then computed as Pr(peak debt > threshold).
# NOTE: $700-$1200/ac may be low for IL corn operating loans; revisit after
# inspecting the empirical peak-debt distribution in Block H.
credit_limits  <- c(700, 900, 1200)

# History conditioning: starting every plan from "CCCCCC" loads the dice
# against continuous corn (no rotation reset in year 1) and inflates the
# early-year gap between CCCCCC and CSCSCS. Run the full grid from both
# a continuous-corn history (the "transition out of corn-on-corn" framing)
# and an alternating history (a neutral baseline) and report both.
histories_vec  <- c("CCCCCC", "CSCSCS")

# ============================================================================
# PAPER ANALYSIS: DELINQUENCY GRID ACROSS ALL ROTATION STRATEGIES
# ============================================================================
# Three objectives from the paper:
#   1. Compare financial outcomes across rotation strategies. NOTE: CCSCCS
#      and CSCSCS have the same long-run corn share and the same number of
#      unique crops; they differ in the *transition pattern* (consecutive
#      corn years vs strict alternation), not in rotation complexity in the
#      RCI sense. Frame the comparison accordingly. A genuinely-more-complex
#      plan (e.g. corn-soy-wheat) would require extending the yield generator
#      and is left for follow-up work.
#   2. Probability of at least one delinquency over 6-year horizon,
#      conditioning yields on the full rotation history.
#   3. Compare delinquency probabilities with and without the negative
#      price-yield correlation (natural hedging effect).
#      NOTE: corr_method = "decomposed" splits yield variance into a
#      systematic component (correlated with price, fraction =
#      aggregate_fraction = 0.40) and an idiosyncratic component
#      (uncorrelated). This addresses the concern that an aggregate
#      price-yield correlation overstates the natural hedge for an
#      individual field, since field-level idiosyncratic shocks get
#      no price cushion.
#
# Credit limit sensitivity: Blocks A/B use "end_balance_positive" with a
# non-binding ceiling (2000) to measure cash-flow delinquency independently
# of the limit. Block H then sweeps credit_limits as post-hoc thresholds on
# the peak-debt distribution, answering: "at what rate would the lender cut
# off this farm under each limit?" The two questions are kept separate
# because they measure different risks.
# ============================================================================

# Helper: long-format a B x T matrix
to_long_mat <- function(mat, plan_label, var_name) {
  dt <- as.data.table(mat)
  setnames(dt, paste0("V", seq_len(ncol(mat))))
  dt[, sim := .I]
  dt_long <- melt(dt, id.vars = "sim",
                  variable.name = "year_chr", value.name = var_name)
  dt_long[, year     := as.integer(sub("V", "", year_chr))]
  dt_long[, plan     := plan_label]
  dt_long[, year_chr := NULL]
  dt_long
}

# Helper: time to first delinquency (NA = never)
time_to_first <- function(delinq_mat) {
  apply(delinq_mat, 1, function(row) {
    w <- which(row)
    if (length(w) == 0) NA_integer_ else w[1L]
  })
}

# Helper: band summary (p10 / p50 / p90) by plan and year
summ_band <- function(dt, var) {
  dt[, .(
    p10 = quantile(get(var), 0.10, na.rm = TRUE),
    p50 = quantile(get(var), 0.50, na.rm = TRUE),
    p90 = quantile(get(var), 0.90, na.rm = TRUE)
  ), by = .(plan, year)]
}

# Helper: run the full region x zone x plan grid for a given history and
# correlation method. Wrapping the grid loop avoids the duplicated block
# structure of the original and makes adding the second history cheap.
#
# Seeding strategy: each (region, zone) gets a distinct seed derived from
# `seed_base`; all plans within that cell share the seed (common random
# numbers — reduces variance of plan-to-plan differences, which are the
# objects of interest for objective 1). Different `seed_base` values for
# different grids (indep / corr / highcost / histories) give independent
# Monte Carlo error across grids.
run_credit_grid <- function(history6,
                            corr_method = c("none", "decomposed"),
                            cost_pars   = NULL,
                            price_shock = NULL,
                            yield_shock = NULL,
                            fix_price   = FALSE,
                            fix_yield   = FALSE,
                            label       = "",
                            seed_base   = 1L) {
  corr_method <- match.arg(corr_method)
  cat(sprintf("\n=== Running credit grid: history=%s, corr=%s %s===\n",
              history6, corr_method,
              if (nzchar(label)) paste0("(", label, ") ") else ""))
  out <- list()
  cell_idx <- 0L
  for (reg in regions_vec) {
    out[[reg]] <- list()
    for (zone in zones_vec) {
      cell_idx <- cell_idx + 1L
      cell_seed <- seed_base * 1000L + cell_idx
      out[[reg]][[zone]] <- list()
      for (plan in rotation_plans) {
        cat(sprintf("  %s | %s | %s  (seed=%d)\n", reg, zone, plan, cell_seed))
        args <- list(
          history6         = history6,
          plan6            = plan,
          geo_region       = reg,
          prod_zone        = zone,
          B                = B,
          i_op_annual      = 0.09,
          credit_limit     = 2000,   # non-binding; sensitivity in Block H
          delinquency_rule = "end_balance_positive",
          corr_method      = corr_method,
          seed             = cell_seed   # CRN within cell, independent across cells
        )
        if (corr_method == "decomposed") {
          # Look up county-aggregated correlations for this region;
          # fall back to defaults if the region is absent or NA.
          .rc <- region_corr[reg]
          args$correlation_corn <- if (nrow(.rc) && !is.na(.rc$corr_corn))
            .rc$corr_corn else DEFAULT_CORR_CORN
          args$correlation_soy  <- if (nrow(.rc) && !is.na(.rc$corr_soy))
            .rc$corr_soy  else DEFAULT_CORR_SOY
        }
        if (!is.null(cost_pars))   args$base_cost_pars <- cost_pars
        if (!is.null(price_shock)) args$price_shock    <- price_shock
        if (!is.null(yield_shock)) args$yield_shock    <- yield_shock
        if (fix_price)             args$fix_price      <- TRUE
        if (fix_yield)             args$fix_yield      <- TRUE
        out[[reg]][[zone]][[plan]] <- tryCatch(
          do.call(simulate_plan_credit_region_zone, args),
          error = function(e) {
            cat(sprintf("  ERROR: %s\n", e$message)); NULL
          }
        )
      }
    }
  }
  out
}

# Helper: collapse a result list into a per-cell summary table
summarise_grid <- function(res) {
  rbindlist(lapply(regions_vec, function(reg) {
    rbindlist(lapply(zones_vec, function(zone) {
      rbindlist(lapply(rotation_plans, function(plan) {
        r <- res[[reg]][[zone]][[plan]]
        if (is.null(r)) return(NULL)
        data.table(
          region               = reg,
          prod_zone            = zone,
          plan                 = plan,
          pr_delinq            = r$pr_delinquent_any,
          exp_delinquent_years = r$exp_delinquent_years,
          mean_peak_debt       = r$mean_peak_debt,
          pv_net_cf            = mean(r$pv_net_cf)
        )
      }))
    }))
  }))
}

# ---- A/B. Credit grids: independent and correlated, both histories ---------
# Four grids total: {indep, corr} x {history CCCCCC, history CSCSCS}.
# Distinct seed_base per grid -> independent MC error across grids.
seed_bases <- list(
  CCCCCC = list(indep = 101L, corr = 102L),
  CSCSCS = list(indep = 201L, corr = 202L)
)
credit_results <- list()
for (h in histories_vec) {
  credit_results[[h]] <- list(
    indep = run_credit_grid(history6 = h, corr_method = "none",
                            label = "indep prices",
                            seed_base = seed_bases[[h]]$indep),
    corr  = run_credit_grid(history6 = h, corr_method = "decomposed",
                            label = "corr prices",
                            seed_base = seed_bases[[h]]$corr)
  )
}
saveRDS(credit_results, "results/credit/credit_results_all_histories.rds")

# Aliases for downstream blocks written against the original two-object
# structure (history = CCCCCC).
credit_indep <- credit_results[["CCCCCC"]]$indep
credit_corr  <- credit_results[["CCCCCC"]]$corr

# ---- C. Delinquency summary table ------------------------------------------
# Long format across history x correlation x cell, then pivot for the
# main paper table (correlated prices, both histories side-by-side).
delinq_long <- rbindlist(lapply(histories_vec, function(h) {
  rbindlist(lapply(c("indep", "corr"), function(cm) {
    s <- summarise_grid(credit_results[[h]][[cm]])
    s[, history    := h]
    s[, corr_label := cm]
    s
  }))
}))

# Wide table: indep vs corr Pr(delinq), for the CCCCCC starting history.
# (Kept in the original column-name shape so downstream code in F/I works
# without changes.)
delinq_summary <- dcast(
  delinq_long[history == "CCCCCC"],
  region + prod_zone + plan ~ corr_label,
  value.var = c("pr_delinq", "exp_delinquent_years",
                "mean_peak_debt", "pv_net_cf")
)
setnames(delinq_summary,
         old = c("pr_delinq_indep", "pr_delinq_corr",
                 "exp_delinquent_years_indep", "exp_delinquent_years_corr",
                 "mean_peak_debt_indep", "mean_peak_debt_corr",
                 "pv_net_cf_indep", "pv_net_cf_corr"),
         new = c("pr_delinq_indep", "pr_delinq_corr",
                 "exp_years_indep", "exp_years_corr",
                 "mean_peak_debt_indep", "mean_peak_debt_corr",
                 "pv_net_cf_indep", "pv_net_cf_corr"),
         skip_absent = TRUE)

cat("\n=== Delinquency summary across all region-zones (history=CCCCCC) ===\n")
print(delinq_summary)
saveRDS(delinq_summary, "results/credit/delinquency_summary.rds")
fwrite(delinq_summary,  "results/credit/delinquency_summary.csv")

# History sensitivity: how much of the gap between plans is driven by the
# corn-on-corn starting state vs the plan itself?
history_sensitivity <- delinq_long[corr_label == "corr",
  .(region, prod_zone, plan, history, pr_delinq)]
history_sens_wide <- dcast(history_sensitivity,
                           region + prod_zone + plan ~ history,
                           value.var = "pr_delinq")
history_sens_wide[, gap_due_to_history := CCCCCC - CSCSCS]
cat("\n=== History sensitivity: Pr(delinq) starting from CCCCCC vs CSCSCS ===\n")
print(history_sens_wide[order(plan, region, prod_zone)])
fwrite(history_sens_wide, "results/credit/history_sensitivity.csv")

# ---- D. Year-by-year delinquency (all regions × zones) ----------------------
dt_yby_all <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    outs <- credit_corr[[reg]][[zone]]
    dt_delq <- rbindlist(lapply(rotation_plans, \(p)
      to_long_mat(outs[[p]]$delinquent * 1L, p, "delinq")))
    dt_delq[, region    := reg]
    dt_delq[, prod_zone := zone]
    dt_delq[, .(pr_delinq = mean(delinq)), by = .(region, prod_zone, plan, year)]
  }))
}))

dt_any_all <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    outs <- credit_corr[[reg]][[zone]]
    rbindlist(lapply(rotation_plans, function(p) {
      d <- outs[[p]]$delinquent
      data.table(region = reg, prod_zone = zone, plan = p,
                 pr_delinq_any = mean(rowSums(d) > 0))
    }))
  }))
}))

dt_yby_all[, region    := factor(region,    levels = regions_vec)]
dt_yby_all[, prod_zone := factor(prod_zone, levels = zones_vec)]
dt_any_all[, region    := factor(region,    levels = regions_vec)]
dt_any_all[, prod_zone := factor(prod_zone, levels = zones_vec)]

cat("\n=== Pr(any delinquency) — all regions x zones ===\n")
print(dt_any_all[order(region, prod_zone, plan)])

p_bar_all <- dt_any_all |>
  ggplot(aes(x = plan, y = pr_delinq_any, fill = plan)) +
  geom_col() +
  facet_grid(region ~ prod_zone, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = plan_colors, guide = "none") +
  labs(title = "Pr(at least one delinquency over 6 years) — all regions x zones",
       x = NULL, y = "Probability") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_year_all <- dt_yby_all |>
  ggplot(aes(x = year, y = pr_delinq, color = plan)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_grid(region ~ prod_zone, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = plan_colors) +
  labs(title = "Delinquency probability by year — all regions x zones",
       x = "Year", y = "Pr(delinquent)") +
  theme_bw(base_size = 9)

delinq_any_long <- melt(
  delinq_summary[, .(region, prod_zone, plan,
                     Independent = pr_delinq_indep,
                     Correlated  = pr_delinq_corr)],
  id.vars       = c("region", "prod_zone", "plan"),
  variable.name = "price_model",
  value.name    = "pr_delinq_any"
)
delinq_any_long[, price_model := factor(price_model,
                                        levels = c("Independent", "Correlated"))]
delinq_any_long[, region    := factor(region,    levels = regions_vec)]
delinq_any_long[, prod_zone := factor(prod_zone, levels = zones_vec)]

p_bar_corr <- delinq_any_long |>
  ggplot(aes(x = plan, y = pr_delinq_any,
             fill = plan, alpha = price_model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  facet_grid(region ~ prod_zone, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = plan_colors, guide = "none") +
  scale_alpha_manual(values = c(Independent = 0.40, Correlated = 1.00),
                     name = "Price model") +
  labs(title = "Pr(at least one delinquency): independent vs. correlated prices",
       subtitle = "Faded = independent | Solid = correlated | Gap = natural hedge benefit",
       x = NULL, y = "Probability") +
  theme_bw(base_size = 9) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

print(p_bar_all)
print(p_bar_corr)
print(p_year_all)

ggsave(file.path(FIGS_DIR, "fig_p_bar_all.png"), p_bar_all,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIGS_DIR, "fig_p_bar_corr.png"), p_bar_corr,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIGS_DIR, "fig_p_bar_corr.pdf"), p_bar_corr,
       width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_p_year_all.png"), p_year_all,
       width = 9, height = 7, dpi = 300)

fwrite(dt_yby_all, "results/credit/delinquency_by_year_all.csv")
fwrite(dt_any_all, "results/credit/delinquency_any_all.csv")

# ---- D2. Publication figure: delinquency probability by year, zone, region --
plan_labels <- c(
  CCCCCC  = "Continuous corn",
  CSCSCS  = "Strict alternation (CS)",
  CCSCCS  = "2-yr corn / 1-yr soy"
)
zone_labels   <- c(high = "High zone", medium = "Medium zone", low = "Low zone")
region_labels <- c(northern = "Northern", central = "Central", southern = "Southern")

fig_delinq_year <- dt_yby_all |>
  ggplot(aes(x = year, y = pr_delinq,
             color = plan, group = plan)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(
    region ~ prod_zone,
    labeller = labeller(region = region_labels, prod_zone = zone_labels)
  ) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = plan_colors, name = "Rotation",
                     labels = plan_labels) +
  labs(
    title    = "Delinquency probability by year, zone, and region",
    subtitle = "Correlated prices | history = CCCCCC",
    x        = "Year",
    y        = "Pr(delinquent in year t)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold"),
    strip.background  = element_rect(fill = "grey92"),
    strip.text        = element_text(face = "bold"),
    panel.grid.minor  = element_blank()
  )

print(fig_delinq_year)
ggsave(file.path(FIGS_DIR, "fig_delinquency_by_year.pdf"), fig_delinq_year,
       width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_delinquency_by_year.png"), fig_delinq_year,
       width = 9, height = 7, dpi = 300)

# ---- E. Survival curves (all regions × zones) ------------------------------
# survival(t) = fraction of paths whose first delinquency is strictly after t
# (or never). Reformulated with the complement of "ever delinquent by t" for
# transparency; numerically equivalent to the original formulation.
dt_surv_all <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    outs <- credit_corr[[reg]][[zone]]
    rbindlist(lapply(rotation_plans, function(p) {
      tfirst <- time_to_first(outs[[p]]$delinquent)
      data.table(
        region    = reg,
        prod_zone = zone,
        plan      = p,
        year      = 0:6,
        surv      = sapply(0:6, function(t)
          1 - mean(!is.na(tfirst) & tfirst <= t))
      )
    }))
  }))
}))

dt_surv_all[, region    := factor(region,    levels = regions_vec)]
dt_surv_all[, prod_zone := factor(prod_zone, levels = zones_vec)]

p_surv_all <- dt_surv_all |>
  ggplot(aes(x = year, y = surv, color = plan)) +
  geom_step(linewidth = 0.8) +
  facet_grid(region ~ prod_zone, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = 0:6) +
  scale_color_manual(values = plan_colors) +
  labs(title = "Survival: Pr(no delinquency yet) — all regions x zones",
       x = "Year", y = "Pr(not yet delinquent)") +
  theme_bw(base_size = 9)
print(p_surv_all)

ggsave(file.path(FIGS_DIR, "fig_p_surv_all.pdf"), p_surv_all,
       width = 9, height = 7, device = grDevices::cairo_pdf)

ggsave(file.path(FIGS_DIR, "fig_p_surv_all.png"), p_surv_all,
       width = 9, height = 7, dpi = 300)

fwrite(dt_surv_all, "results/credit/survival_curves_all.csv")

# ---- F. Correlation impact on delinquency (objective 3) --------------------
corr_impact <- delinq_summary[, .(
  region, prod_zone, plan,
  pr_delinq_indep = pr_delinq_indep,
  pr_delinq_corr  = pr_delinq_corr,
  reduction_pp    = pr_delinq_indep - pr_delinq_corr,
  reduction_pct   = (pr_delinq_indep - pr_delinq_corr) / pr_delinq_indep * 100
)]

corr_impact[, region    := factor(region,    levels = regions_vec)]
corr_impact[, prod_zone := factor(prod_zone, levels = zones_vec)]
cat("\n=== Correlation impact on delinquency probability ===\n")
print(corr_impact[order(plan, region, prod_zone)])
fwrite(corr_impact, "results/credit/correlation_impact_delinquency.csv")

# Heterogeneity in the natural-hedge effect: the abstract claims correlation
# reduces delinquency, but magnitude should differ by productivity zone and
# plan. Pull out the heterogeneity explicitly.
corr_hetero <- corr_impact[, .(
  mean_reduction_pp   = mean(reduction_pp),
  median_reduction_pp = median(reduction_pp),
  min_reduction_pp    = min(reduction_pp),
  max_reduction_pp    = max(reduction_pp)
), by = .(plan, prod_zone)]
cat("\n=== Heterogeneity in natural-hedge effect by plan x zone ===\n")
print(corr_hetero[order(plan, prod_zone)])
fwrite(corr_hetero, "results/credit/correlation_impact_heterogeneity.csv")

p_heat <- corr_impact |>
  ggplot(aes(x = prod_zone, y = region, fill = reduction_pp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f pp", reduction_pp)), size = 3) +
  scale_fill_gradient2(low = "white", high = "steelblue", name = "Reduction\n(pp)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "Delinquency reduction from price-yield correlation",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_heat)

ggsave(file.path(FIGS_DIR, "fig_p_heat.pdf"), p_heat,
       width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_p_heat.png"), p_heat,
       width = 9, height = 7, dpi = 300)

# ---- F2. NPV distribution: full correlation impact by plan and history ------
# Extracts the raw pv_net_cf vector (length B) for every region x zone cell
# under both price assumptions and overlays density curves. Dashed lines mark
# the median of each distribution so the shift is immediately readable.
extract_pv_long <- function(history6, plan) {
  rbindlist(lapply(c("indep", "corr"), function(cm) {
    rbindlist(lapply(regions_vec, function(reg) {
      rbindlist(lapply(zones_vec, function(zone) {
        rc <- credit_results[[history6]][[cm]][[reg]][[zone]][[plan]]
        if (is.null(rc)) return(NULL)
        data.table(
          history    = history6,
          corr_label = if (cm == "indep") "Independent prices" else "Correlated prices",
          region     = reg,
          prod_zone  = zone,
          plan       = plan,
          pv_net_cf  = rc$pv_net_cf
        )
      }))
    }))
  }))
}

corr_colors <- c("Independent prices" = "#d73027", "Correlated prices" = "#4575b4")

for (focal_hist in histories_vec) {
  for (focal_plan in rotation_plans) {
    dt_pv <- extract_pv_long(focal_hist, focal_plan)
    dt_pv[, region     := factor(region,     levels = regions_vec)]
    dt_pv[, prod_zone  := factor(prod_zone,  levels = zones_vec)]
    dt_pv[, corr_label := factor(corr_label,
                                 levels = c("Independent prices", "Correlated prices"))]

    dt_med <- dt_pv[, .(med = median(pv_net_cf)), by = .(region, prod_zone, corr_label)]

    fig_npv <- dt_pv |>
      ggplot(aes(x = pv_net_cf, color = corr_label, fill = corr_label)) +
      geom_density(alpha = 0.20, linewidth = 0.8) +
      geom_vline(data = dt_med,
                 aes(xintercept = med, color = corr_label),
                 linetype = "dashed", linewidth = 0.6) +
      facet_grid(
        region ~ prod_zone,
        labeller = labeller(region = region_labels, prod_zone = zone_labels)
      ) +
      scale_x_continuous(labels = scales::dollar_format(suffix = "/ac")) +
      scale_color_manual(values = corr_colors, name = NULL) +
      scale_fill_manual(values  = corr_colors, name = NULL) +
      labs(
        title    = sprintf("NPV distribution: correlation impact — plan %s", focal_plan),
        subtitle = sprintf("History = %s | dashed lines = medians", focal_hist),
        x        = "PV net cash flow ($/ac)",
        y        = "Density"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    print(fig_npv)
    out_stem <- sprintf(file.path(FIGS_DIR, "fig_npv_corr_impact_%s_%s"), focal_plan, focal_hist)
    ggsave(paste0(out_stem, ".pdf"), fig_npv, width = 9, height = 7, device = grDevices::cairo_pdf)
    ggsave(paste0(out_stem, ".png"), fig_npv, width = 9, height = 7, dpi = 300)
  }
}

# ---- G. Cash flow and carried debt bands (all regions × zones) -------------
# net_cash is a flow: cumsum across years gives cumulative position.
# debt_end is a stock (rollover balance): plot directly by year, not cumsum.
dt_cf_all <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    outs <- credit_corr[[reg]][[zone]]
    dt <- rbindlist(lapply(rotation_plans, \(p)
      to_long_mat(outs[[p]]$net_cash, p, "net_cash")))
    dt[, region    := reg]
    dt[, prod_zone := zone]
    dt
  }))
}))
dt_cf_all[, cum_cf := ave(net_cash, region, prod_zone, plan, sim, FUN = cumsum)]

dt_debt_all <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    outs <- credit_corr[[reg]][[zone]]
    dt <- rbindlist(lapply(rotation_plans, \(p)
      to_long_mat(outs[[p]]$debt_end, p, "debt_end")))
    dt[, region    := reg]
    dt[, prod_zone := zone]
    dt
  }))
}))

dt_cf_all[,   region    := factor(region,    levels = regions_vec)]
dt_cf_all[,   prod_zone := factor(prod_zone, levels = zones_vec)]
dt_debt_all[, region    := factor(region,    levels = regions_vec)]
dt_debt_all[, prod_zone := factor(prod_zone, levels = zones_vec)]

band_cf_all <- dt_cf_all[, .(
  p10 = quantile(cum_cf, 0.10),
  p50 = quantile(cum_cf, 0.50),
  p90 = quantile(cum_cf, 0.90)
), by = .(region, prod_zone, plan, year)]

band_debt_all <- dt_debt_all[, .(
  mean       = mean(debt_end, na.rm = TRUE),
  p90        = quantile(debt_end, 0.90, na.rm = TRUE),
  p99        = quantile(debt_end, 0.99, na.rm = TRUE),
  pr_pos     = mean(debt_end > 0, na.rm = TRUE)   # fraction of paths with unpaid debt
), by = .(region, prod_zone, plan, year)]

p_cf_all <- ggplot(band_cf_all, aes(x = year, color = plan)) +
  geom_ribbon(aes(ymin = p10, ymax = p90, fill = plan), alpha = 0.15, color = NA) +
  geom_line(aes(y = p50), linewidth = 0.8) +
  facet_grid(region ~ prod_zone, labeller = label_both) +
  scale_color_manual(values = plan_colors) +
  scale_fill_manual(values = plan_colors, guide = "none") +
  labs(title = "Cumulative net cash (median + 10-90 band) — all regions x zones",
       x = "Year", y = "Cumulative net cash ($/ac)") +
  theme_bw(base_size = 9)

# Ribbon uses p90-p99 to capture the delinquent tail.
# For Southern IL, >90% of paths repay in full each year, so p75/p90 = 0;
# p90-p99 shows the regime where debt actually accumulates.
#p_debt_all <- ggplot(band_debt_all, aes(x = year, color = plan)) +
#  geom_ribbon(aes(ymin = p90, ymax = p99, fill = plan), alpha = 0.20, color = NA) +
#  geom_line(aes(y = mean), linewidth = 0.8) +
#  facet_grid(prod_zone ~ region, labeller = label_both, scales = "free") +
#  labs(title = "Year-end carried debt (mean + 90-99 tail band) — all regions x zones",
#       subtitle = "Band shows top 9% of paths by debt level",
#       x = "Year", y = "Rollover debt ($/ac)") +
#  guides(fill = "none") +
#  theme_bw(base_size = 9)

# Build the same debt summary for independent-price draws to show hedge effect
dt_debt_indep <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    outs <- credit_indep[[reg]][[zone]]
    dt <- rbindlist(lapply(rotation_plans, \(p)
      to_long_mat(outs[[p]]$debt_end, p, "debt_end")))
    dt[, region    := reg]
    dt[, prod_zone := zone]
    dt
  }))
}))
dt_debt_indep[, region    := factor(region,    levels = regions_vec)]
dt_debt_indep[, prod_zone := factor(prod_zone, levels = zones_vec)]

band_debt_indep <- dt_debt_indep[, .(mean = mean(debt_end, na.rm = TRUE)),
                                 by = .(region, prod_zone, plan, year)]

band_debt_hedge <- rbind(
  band_debt_all[,   .(region, prod_zone, plan, year, mean, model = "Correlated (hedge)")],
  band_debt_indep[, .(region, prod_zone, plan, year, mean, model = "Independent")]
)
band_debt_hedge[, model := factor(model,
  levels = c("Independent", "Correlated (hedge)"))]

p_debt_all <- ggplot(band_debt_hedge,
                     aes(x = year, y = mean, color = plan, linetype = model)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = plan_colors) +
  scale_linetype_manual(values = c("Independent" = "dashed",
                                   "Correlated (hedge)" = "solid"),
                        name = NULL) +
  facet_grid(region ~ prod_zone, labeller = label_both, scales = "free_y") +
  labs(title = "Natural hedge from crop rotation: mean year-end debt",
       subtitle = paste("Solid = correlated prices/yields (hedge active)",
                        "| Dashed = independent (no hedge)",
                        "\nGap = hedge benefit; wider for CSCSCS (stronger soy correlation)"),
       x = "Year", y = "Mean rollover debt ($/ac)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

print(p_cf_all)
print(p_debt_all)

ggsave(file.path(FIGS_DIR, "fig_p_cf_all.png"), p_cf_all,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIGS_DIR, "fig_p_debt_all.png"), p_debt_all,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIGS_DIR, "fig_p_cf_all.pdf"), p_cf_all,
       width = 9, height = 7, dpi = 300, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_p_debt_all.pdf"), p_debt_all,
       width = 9, height = 7, dpi = 300, device = grDevices::cairo_pdf)       

fwrite(band_cf_all,   "results/credit/cashflow_bands_all.csv")
fwrite(band_debt_all, "results/credit/debt_bands_all.csv")

# ---- H. Credit-rationing sensitivity (post-hoc on peak debt) ---------------
# For each simulation path, peak debt = max year-end balance over 6 years.
# Pr(peak debt > threshold) answers: how often would this farm exceed each
# lender's credit ceiling, triggering a credit-rationing event?
cat("\n=== Credit-rationing sensitivity ===\n")

# Inspect the empirical peak-debt distribution before deciding whether the
# chosen credit_limits are informative. If even the 99.9th percentile is
# below the smallest limit, the rationing block is a null finding and the
# limits should be revisited.
# peak_gross_loan(): gross spring loan = operating costs + prior-year rollover.
# This is the quantity a lender evaluates against a credit ceiling, not
# debt_end (which is the post-harvest residual and is 0 for solvent paths).
peak_gross_loan <- function(rc) {
  debt_carry <- cbind(0, rc$debt_end[, -6])   # prior year residual; year 1 = 0
  apply(rc$cost + debt_carry, 1, max)
}

cat("\n--- Peak gross loan diagnostics (corr prices, history=CCCCCC) ---\n")
peak_diag <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    rbindlist(lapply(rotation_plans, function(plan) {
      rc <- credit_corr[[reg]][[zone]][[plan]]
      if (is.null(rc)) return(NULL)
      pd <- peak_gross_loan(rc)
      data.table(region = reg, prod_zone = zone, plan = plan,
                 p95  = quantile(pd, 0.95),
                 p99  = quantile(pd, 0.99),
                 p999 = quantile(pd, 0.999),
                 max  = max(pd))
    }))
  }))
}))
print(peak_diag[order(-max)][1:10])  # 10 worst cells
fwrite(peak_diag, "results/credit/peak_debt_diagnostics.csv")

rationing_summary <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    rbindlist(lapply(rotation_plans, function(plan) {
      rc <- credit_corr[[reg]][[zone]][[plan]]
      if (is.null(rc)) return(NULL)
      peak_loan <- peak_gross_loan(rc)
      rbindlist(lapply(credit_limits, function(lim) {
        data.table(
          region       = reg,
          prod_zone    = zone,
          plan         = plan,
          credit_limit = lim,
          pr_rationed  = mean(peak_loan > lim),
          p50_peak     = median(peak_loan),
          p90_peak     = quantile(peak_loan, 0.90)
        )
      }))
    }))
  }))
}))

print(rationing_summary[order(credit_limit, plan, region, prod_zone)])
saveRDS(rationing_summary, "results/credit/credit_rationing_sensitivity.rds")
fwrite(rationing_summary,  "results/credit/credit_rationing_sensitivity.csv")

# Peak gross loan distribution — correlated and independent models together
# to expose the natural hedge: negative price-yield correlation compresses
# the upper tail of spring borrowing demand.
peak_dt_all <- rbindlist(lapply(
  list(list(label = "Correlated",  sim = credit_corr),
       list(label = "Independent", sim = credit_indep)),
  function(item) {
    rbindlist(lapply(regions_vec, function(reg) {
      rbindlist(lapply(zones_vec, function(zone) {
        rbindlist(lapply(rotation_plans, function(p) {
          rc <- item$sim[[reg]][[zone]][[p]]
          if (is.null(rc)) return(NULL)
          data.table(region = reg, prod_zone = zone, plan = p,
                     model = item$label,
                     peak_loan = peak_gross_loan(rc))
        }))
      }))
    }))
  }
))
peak_dt_all[, region    := factor(region,    levels = regions_vec)]
peak_dt_all[, prod_zone := factor(prod_zone, levels = zones_vec)]
peak_dt_all[, model     := factor(model, levels = c("Correlated", "Independent"))]

limit_labels <- data.table(
  credit_limit = credit_limits,
  label        = paste0("$", credit_limits, "/ac"),
  color        = c("#d73027", "#fc8d59", "#4393c3")
)

# Blue (correlated) violin narrower in the upper tail than red (independent)
# = hedge value from negative price-yield co-movement.
p_ration_all <- peak_dt_all |>
  ggplot(aes(x = plan, y = peak_loan, fill = model)) +
  geom_violin(alpha = 0.45, position = position_dodge(width = 0.9),
              draw_quantiles = c(0.50, 0.90), scale = "width") +
  geom_hline(data = limit_labels,
             aes(yintercept = credit_limit, color = label),
             linetype = "dashed", linewidth = 0.7) +
  facet_grid(region ~ prod_zone, labeller = label_both) +
  scale_fill_manual(values = c("Correlated" = "#2166ac", "Independent" = "#d6604d"),
                    name = "Model") +
  scale_color_manual(values = setNames(limit_labels$color, limit_labels$label),
                     name = "Credit limit") +
  labs(title = "Natural hedge: peak gross spring loan — correlated vs. independent",
       subtitle = "Blue (correlated) upper tail narrower than red (independent) = hedge value",
       x = NULL, y = "Peak loan demanded ($/ac)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_ration_all)

ggsave(file.path(FIGS_DIR, "fig_peak_loan_hedge.png"), p_ration_all,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIGS_DIR, "fig_peak_loan_all.png"),   p_ration_all,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIGS_DIR, "fig_peak_loan_hedge.pdf"), p_ration_all,
       width = 9, height = 7, dpi = 300, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_peak_loan_all.pdf"),   p_ration_all,
       width = 9, height = 7, dpi = 300, device = grDevices::cairo_pdf)       

fwrite(peak_dt_all, "results/credit/peak_gross_loan_all.csv")

# ---- I. Cost sensitivity: +25% direct costs -------------------------------
# Mechanism: higher costs raise principal_needed at planting → larger
# debt_at_harvest → more paths with remaining_balance > 0 → delinquency.
# Baseline for comparison is credit_corr (correlated, non-binding limit).
cat("\n=== Cost sensitivity sweep: +15% and +25% direct costs ===\n")

# Assumes Cbar = E[C] = exp(mu + 0.5*sigma^2). Verify against base_cost_params
# construction if results look off (median convention would require mu update
# without the -0.5*sigma^2 term).
.cb_check <- base_cost_params[[ regions_vec[1] ]][[ zones_vec[1] ]]$C
.cb_mean  <- exp(.cb_check$mu + 0.5 * .cb_check$sigma^2)
.cb_med   <- exp(.cb_check$mu)
if (!isTRUE(all.equal(.cb_check$Cbar, .cb_mean, tolerance = 1e-6))) {
  if (isTRUE(all.equal(.cb_check$Cbar, .cb_med, tolerance = 1e-6))) {
    warning("base_cost_params$Cbar appears to be the MEDIAN (exp(mu)), ",
            "not the mean. The mu update below assumes Cbar = mean. ",
            "Adjust before running.")
  } else {
    warning("base_cost_params$Cbar matches neither exp(mu+0.5*sigma^2) ",
            "nor exp(mu); verify the cost-distribution parameterization.")
  }
}
rm(.cb_check, .cb_mean, .cb_med)

cost_impact_list <- list()

for (cost_mult in c(1.15, 1.25)) {
  mult_label <- sprintf("%+.0f%%", (cost_mult - 1) * 100)   # "+15%", "+25%"
  tag        <- sprintf("%.0f",    (cost_mult - 1) * 100)   # "15",   "25"

  shocked_cost <- base_cost_params
  for (.reg in names(shocked_cost)) {
    for (.zone in names(shocked_cost[[.reg]])) {
      shocked_cost[[.reg]][[.zone]]$C$Cbar <-
        base_cost_params[[.reg]][[.zone]]$C$Cbar * cost_mult
      shocked_cost[[.reg]][[.zone]]$S$Cbar <-
        base_cost_params[[.reg]][[.zone]]$S$Cbar * cost_mult
      shocked_cost[[.reg]][[.zone]]$C$mu <-
        log(shocked_cost[[.reg]][[.zone]]$C$Cbar) -
        0.5 * base_cost_params[[.reg]][[.zone]]$C$sigma^2
      shocked_cost[[.reg]][[.zone]]$S$mu <-
        log(shocked_cost[[.reg]][[.zone]]$S$Cbar) -
        0.5 * base_cost_params[[.reg]][[.zone]]$S$sigma^2
    }
  }
  rm(.reg, .zone)

  grid_out <- run_credit_grid(
    history6    = "CCCCCC",
    corr_method = "decomposed",
    cost_pars   = shocked_cost,
    label       = paste0(mult_label, " costs"),
    seed_base   = 300L + as.integer(tag)   # 315L for +15%, 325L for +25%
  )
  saveRDS(grid_out,
          sprintf("results/credit/credit_corr_highcost_%spct.rds", tag))

  impact <- merge(
    delinq_summary[, .(region, prod_zone, plan,
                       pr_delinq_base = pr_delinq_corr,
                       pv_net_cf_base = pv_net_cf_corr)],
    summarise_grid(grid_out)[, .(region, prod_zone, plan,
                                  pr_delinq_shock = pr_delinq,
                                  pv_net_cf_shock = pv_net_cf)],
    by = c("region", "prod_zone", "plan")
  )
  impact[, delta_pr_delinq_pp := (pr_delinq_shock - pr_delinq_base) * 100]
  impact[, delta_pv_net_cf    :=  pv_net_cf_shock - pv_net_cf_base]
  impact[, cost_shock         := mult_label]

  saveRDS(impact, sprintf("results/credit/cost_sensitivity_%spct.rds", tag))
  fwrite(impact,  sprintf("results/credit/cost_sensitivity_%spct.csv", tag))
  cost_impact_list[[mult_label]] <- impact
  cat(sprintf("  %s cost shock complete.\n", mult_label))
}

cost_impact_all <- rbindlist(cost_impact_list)
cost_impact_all[, region    := factor(region,    levels = regions_vec)]
cost_impact_all[, prod_zone := factor(prod_zone, levels = zones_vec)]
cat("\n=== Cost sensitivity sweep results ===\n")
print(cost_impact_all[order(cost_shock, plan, region, prod_zone)])
saveRDS(cost_impact_all, "results/credit/cost_sensitivity_sweep.rds")
fwrite(cost_impact_all,  "results/credit/cost_sensitivity_sweep.csv")

p_cost_heat <- cost_impact_all |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pr_delinq_pp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("+%.1f pp", delta_pr_delinq_pp)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick",
                      name = "Chg. Pr(delinq)\n(pp)") +
  scale_y_discrete(limits = rev) +
  facet_grid(cost_shock ~ plan) +
  labs(title = "Delinquency increase from direct cost shock (vs baseline)",
       subtitle = "Correlated prices | history = CCCCCC",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_cost_heat)
ggsave(file.path(FIGS_DIR, "fig_cost_sensitivity_heat.pdf"),
       p_cost_heat, width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_cost_sensitivity_heat.png"),
       p_cost_heat, width = 9, height = 7, dpi = 300)

p_cost_pv <- cost_impact_all |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pv_net_cf)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("$%.0f", delta_pv_net_cf)), size = 3) +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                       midpoint = 0, name = "Chg. PV net\ncash ($/ac)") +
  scale_y_discrete(limits = rev) +
  facet_grid(cost_shock ~ plan) +
  labs(title = "PV net cash change from direct cost shock (vs baseline)",
       subtitle = "Correlated prices | history = CCCCCC",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_cost_pv)
ggsave(file.path(FIGS_DIR, "fig_cost_sensitivity_pv.pdf"),
       p_cost_pv, width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_cost_sensitivity_pv.png"),
       p_cost_pv, width = 9, height = 7, dpi = 300)

# ---- J. Tail-dependence robustness: t-copula vs Gaussian -------------------
# Gaussian copula has zero tail dependence by construction. For multi-year
# ruin probabilities, joint extremes in (low yield, high price) and
# (high yield, low price) drive the result, and tail dependence can shift
# Pr(delinq) materially even at the same linear correlation. Re-run a
# representative cell (central / medium) with a t-copula at low df and
# compare. Block J fails gracefully if the simulator has not yet been
# extended to support t-copula draws.
cat("\n=== Tail-dependence robustness: t-copula at central/medium ===\n")

run_tcopula_cell <- function(plan, df = 5) {
  args <- list(
    history6         = "CCCCCC",
    plan6            = plan,
    geo_region       = "central",
    prod_zone        = "medium",
    B                = B,
    i_op_annual      = 0.09,
    credit_limit     = 2000,
    delinquency_rule = "end_balance_positive",
    corr_method      = "t_copula",   # requires simulator support
    correlation_corn = -0.40,
    correlation_soy  = -0.50,
    copula_df        = df,
    seed             = 401L
  )
  tryCatch(do.call(simulate_plan_credit_region_zone, args),
           error = function(e) {
             message("t-copula path not available: ", e$message,
                     ". Skipping robustness block.")
             NULL
           })
}

tcop_results <- setNames(lapply(rotation_plans, run_tcopula_cell),
                         rotation_plans)

if (!all(sapply(tcop_results, is.null))) {
  tcop_summary <- rbindlist(lapply(rotation_plans, function(p) {
    rc <- tcop_results[[p]]
    if (is.null(rc)) return(NULL)
    data.table(
      plan      = p,
      copula    = "t (df=5)",
      pr_delinq = rc$pr_delinquent_any,
      mean_peak = rc$mean_peak_debt,
      p99_peak  = quantile(peak_gross_loan(rc), 0.99)
    )
  }))
  gauss_summary <- rbindlist(lapply(rotation_plans, function(p) {
    rc <- credit_corr$central$medium[[p]]
    data.table(
      plan      = p,
      copula    = "Gaussian",
      pr_delinq = rc$pr_delinquent_any,
      mean_peak = rc$mean_peak_debt,
      p99_peak  = quantile(peak_gross_loan(rc), 0.99)
    )
  }))
  copula_compare <- rbind(gauss_summary, tcop_summary)
  cat("\n=== Gaussian vs t-copula comparison (central / medium) ===\n")
  print(copula_compare[order(plan, copula)])
  fwrite(copula_compare, "results/credit/copula_robustness.csv")

  p_copula <- copula_compare |>
    melt(id.vars = c("plan", "copula"),
         measure.vars = c("pr_delinq", "p99_peak"),
         variable.name = "metric", value.name = "value") |>
    ggplot(aes(x = plan, y = value, fill = copula)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Gaussian" = "steelblue", "t (df=5)" = "firebrick"),
                      name = "Copula") +
    facet_wrap(~metric, scales = "free_y",
               labeller = as_labeller(c(
                 pr_delinq = "Pr(any delinquency)",
                 p99_peak  = "P99 peak debt ($/ac)"
               ))) +
    labs(title = "Tail-dependence robustness: Gaussian vs t-copula (df=5)",
         subtitle = "Central IL, medium productivity | history = CCCCCC",
         x = "Rotation plan", y = NULL) +
    theme_bw()
  print(p_copula)
  ggsave(file.path(FIGS_DIR, "fig_copula_robustness.pdf"),
         p_copula, width = 8, height = 4, device = grDevices::cairo_pdf)
  ggsave(file.path(FIGS_DIR, "fig_copula_robustness.png"),
         p_copula, width = 8, height = 4, dpi = 300)
} else {
  cat("\nt-copula corr_method not implemented in simulator; skipping Block J.\n")
  cat("To enable: add a 't_copula' branch to correlation_functions_region_zone.R\n")
  cat("that draws standard t variates with the same correlation matrix and\n")
  cat("transforms via pt() to uniforms before applying the marginals.\n")
}

# ---- K. Yield sensitivity: corn yield -12% in year 3 ----------------------
# Structural roles by plan:
#   CCCCCC (C,C,C,C,C,C): corn in yr 3 → direct loss in revenue.
#   CSCSCS (C,S,C,S,C,S): corn in yr 3 → direct loss in revenue.
#   CCSCCS (C,C,S,C,C,S): soy  in yr 3 → null effect (clean control).
# The shock is applied post-draw (Y multiplied after all yield draws), so it
# is additive to the price-yield correlation structure. Only years where the
# plan actually plants corn are affected; soy-year positions in yield_shock$C
# are ignored because the application loop iterates over corn_yrs only.
cat("\n=== Yield sensitivity: -12% corn yield in year 3 ===\n")

corn_yr3_yield_shock <- list(C = c(1, 1, 0.88, 1, 1, 1))

credit_corr_cornyield_yr3 <- run_credit_grid(
  history6    = "CCCCCC",
  corr_method = "decomposed",
  yield_shock = corn_yr3_yield_shock,
  label       = "-12% corn yield yr3",
  seed_base   = 501L
)
saveRDS(credit_corr_cornyield_yr3, "results/credit/credit_corr_cornyield_yr3.rds")

delinq_summary_cornyield_yr3 <- summarise_grid(credit_corr_cornyield_yr3)

yield_impact_yr3 <- merge(
  delinq_summary[, .(region, prod_zone, plan,
                     pr_delinq_base  = pr_delinq_corr,
                     pv_net_cf_base  = pv_net_cf_corr)],
  delinq_summary_cornyield_yr3[, .(region, prod_zone, plan,
                                    pr_delinq_shock = pr_delinq,
                                    pv_net_cf_shock = pv_net_cf)],
  by = c("region", "prod_zone", "plan")
)
yield_impact_yr3[, delta_pr_delinq_pp := (pr_delinq_shock - pr_delinq_base) * 100]
yield_impact_yr3[, delta_pv_net_cf    :=  pv_net_cf_shock - pv_net_cf_base]
yield_impact_yr3[, region    := factor(region,    levels = regions_vec)]
yield_impact_yr3[, prod_zone := factor(prod_zone, levels = zones_vec)]

cat("\n=== Yield impact: -12% corn yr3 vs baseline ===\n")
print(yield_impact_yr3[order(plan, region, prod_zone)])
saveRDS(yield_impact_yr3, "results/credit/yield_sensitivity_corn_yr3.rds")
fwrite(yield_impact_yr3,  "results/credit/yield_sensitivity_corn_yr3.csv")

p_yield_heat_yr3 <- yield_impact_yr3 |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pr_delinq_pp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f pp", delta_pr_delinq_pp)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Chg. Pr(delinq)\n(pp)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "-12% corn yield (year 3 only): delinquency change vs baseline",
       subtitle = "Positive = worsening | correlated prices | history = CCCCCC\nCCSCCS: null (soy in year 3)",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_yield_heat_yr3)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_corn_yr3.pdf"),
       p_yield_heat_yr3, width = 9, height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_corn_yr3.png"),
       p_yield_heat_yr3, width = 9, height = 5, dpi = 300)

p_yield_pv_yr3 <- yield_impact_yr3 |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pv_net_cf)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("$%.0f", delta_pv_net_cf)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Chg. PV net\ncash ($/ac)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "-12% corn yield (year 3 only): PV net cash change vs baseline",
       subtitle = "Correlated prices | history = CCCCCC",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_yield_pv_yr3)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_pv_corn_yr3.pdf"),
       p_yield_pv_yr3, width = 9, height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_pv_corn_yr3.png"),
       p_yield_pv_yr3, width = 9, height = 5, dpi = 300)

# ---- K2. Yield sensitivity (independent prices): corn yield -12% in year 3 --
# Same shock as K but with corr_method = "none" (independent prices).
# Baseline is the independent-price grid (pr_delinq_indep / pv_net_cf_indep).
# Structural roles are identical: CCCCCC and CSCSCS have corn in yr 3;
# CCSCCS has soy → null control.
cat("\n=== Yield sensitivity (independent prices): -12% corn yield in year 3 ===\n")

credit_indep_cornyield_yr3 <- run_credit_grid(
  history6    = "CCCCCC",
  corr_method = "none",
  yield_shock = corn_yr3_yield_shock,
  label       = "-12% corn yield yr3 (indep)",
  seed_base   = 502L
)
saveRDS(credit_indep_cornyield_yr3, "results/credit/credit_indep_cornyield_yr3.rds")

delinq_summary_cornyield_yr3_indep <- summarise_grid(credit_indep_cornyield_yr3)

yield_impact_yr3_indep <- merge(
  delinq_summary[, .(region, prod_zone, plan,
                     pr_delinq_base  = pr_delinq_indep,
                     pv_net_cf_base  = pv_net_cf_indep)],
  delinq_summary_cornyield_yr3_indep[, .(region, prod_zone, plan,
                                          pr_delinq_shock = pr_delinq,
                                          pv_net_cf_shock = pv_net_cf)],
  by = c("region", "prod_zone", "plan")
)
yield_impact_yr3_indep[, delta_pr_delinq_pp := (pr_delinq_shock - pr_delinq_base) * 100]
yield_impact_yr3_indep[, delta_pv_net_cf    :=  pv_net_cf_shock - pv_net_cf_base]
yield_impact_yr3_indep[, region    := factor(region,    levels = regions_vec)]
yield_impact_yr3_indep[, prod_zone := factor(prod_zone, levels = zones_vec)]

cat("\n=== Yield impact (indep): -12% corn yr3 vs baseline ===\n")
print(yield_impact_yr3_indep[order(plan, region, prod_zone)])
saveRDS(yield_impact_yr3_indep, "results/credit/yield_sensitivity_corn_yr3_indep.rds")
fwrite(yield_impact_yr3_indep,  "results/credit/yield_sensitivity_corn_yr3_indep.csv")

p_yield_heat_yr3_indep <- yield_impact_yr3_indep |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pr_delinq_pp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f pp", delta_pr_delinq_pp)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Chg. Pr(delinq)\n(pp)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "-12% corn yield (year 3 only): delinquency change vs baseline",
       subtitle = "Positive = worsening | independent prices | history = CCCCCC\nCCSCCS: null (soy in year 3)",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_yield_heat_yr3_indep)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_corn_yr3_indep.pdf"),
       p_yield_heat_yr3_indep, width = 9, height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_corn_yr3_indep.png"),
       p_yield_heat_yr3_indep, width = 9, height = 5, dpi = 300)

p_yield_pv_yr3_indep <- yield_impact_yr3_indep |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pv_net_cf)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("$%.0f", delta_pv_net_cf)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Chg. PV net\ncash ($/ac)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "-10% corn yield (year 3 only): PV net cash change vs baseline",
       subtitle = "Independent prices | history = CCCCCC",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_yield_pv_yr3_indep)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_pv_corn_yr3_indep.pdf"),
       p_yield_pv_yr3_indep, width = 9, height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_yield_sensitivity_pv_corn_yr3_indep.png"),
       p_yield_pv_yr3_indep, width = 9, height = 5, dpi = 300)

cat("\n=== Paper analysis complete. Results saved to results/credit/ ===\n")

# ---- K3. Comparison plot: correlated vs. independent prices -----------------
yield_impact_combined <- rbind(
  copy(yield_impact_yr3)[,       price_model := "Correlated"],
  copy(yield_impact_yr3_indep)[, price_model := "Independent"]
)
yield_impact_combined[, price_model := factor(price_model,
                                              levels = c("Correlated", "Independent"))]

p_compare_delinq <- yield_impact_combined |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pr_delinq_pp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f pp", delta_pr_delinq_pp)), size = 2.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Chg. Pr(delinq)\n(pp)") +
  scale_y_discrete(limits = rev) +
  facet_grid(price_model ~ plan) +
  labs(title = "Delinquency change",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"))

p_compare_pv <- yield_impact_combined |>
  ggplot(aes(x = prod_zone, y = region, fill = delta_pv_net_cf)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("$%.0f", delta_pv_net_cf)), size = 2.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Chg. PV net\ncash ($/ac)") +
  scale_y_discrete(limits = rev) +
  facet_grid(price_model ~ plan) +
  labs(title = "PV net cash change",
       x = "Productivity zone", y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"))

p_yield_compare <- p_compare_delinq / p_compare_pv +
  plot_annotation(
    title    = "-12% corn yield shock (year 3 only): correlated vs. independent prices",
    subtitle = "Positive = worsening | history = CCCCCC | CCSCCS: null control (soy in year 3)",
    theme    = theme(plot.title    = element_text(face = "bold"),
                     plot.subtitle = element_text(color = "grey40"))
  )

print(p_yield_compare)
ggsave(file.path(FIGS_DIR, "fig_yield_compare_corr_vs_indep.pdf"),
       p_yield_compare, width = 10, height = 8, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_yield_compare_corr_vs_indep.png"),
       p_yield_compare, width = 10, height = 8, dpi = 300)


# ---- L. Shapley decomposition: price vs. yield uncertainty -----------------
# Four variants of the correlated grid:
#   both        : baseline (reuse credit_corr from Block A)
#   price_only  : fix_yield=TRUE  — price varies, yield fixed at E[Y|window]
#   yield_only  : fix_price=TRUE  — yield varies, price fixed at Pbar
#   none        : both fixed      — residual structural delinquency
#
# Shapley values for two factors are additive and sum to (both - none):
#   shapley_price = 0.5*(price_only - none) + 0.5*(both - yield_only)
#   shapley_yield = 0.5*(yield_only - none) + 0.5*(both - price_only)
cat("\n=== Shapley decomposition: price vs. yield uncertainty ===\n")

credit_shapley_price <- run_credit_grid(
  history6    = "CCCCCC",
  corr_method = "decomposed",
  fix_yield   = TRUE,
  label       = "price-only",
  seed_base   = 601L
)
saveRDS(credit_shapley_price, "results/credit/credit_shapley_price.rds")

credit_shapley_yield <- run_credit_grid(
  history6    = "CCCCCC",
  corr_method = "decomposed",
  fix_price   = TRUE,
  label       = "yield-only",
  seed_base   = 602L
)
saveRDS(credit_shapley_yield, "results/credit/credit_shapley_yield.rds")

credit_shapley_none <- run_credit_grid(
  history6    = "CCCCCC",
  corr_method = "decomposed",
  fix_price   = TRUE,
  fix_yield   = TRUE,
  label       = "deterministic",
  seed_base   = 603L
)
saveRDS(credit_shapley_none, "results/credit/credit_shapley_none.rds")

credit_shapley_price_2 <- run_credit_grid(
  history6    = "CSCSCS",
  corr_method = "decomposed",
  fix_yield   = TRUE,
  label       = "price-only",
  seed_base   = 601L
)
saveRDS(credit_shapley_price, "results/credit/credit_shapley_price_2.rds")

credit_shapley_yield_2 <- run_credit_grid(
  history6    = "CSCSCS",
  corr_method = "decomposed",
  fix_price   = TRUE,
  label       = "yield-only",
  seed_base   = 602L
)
saveRDS(credit_shapley_yield, "results/credit/credit_shapley_yield_2.rds")

credit_shapley_none_2 <- run_credit_grid(
  history6    = "CSCSCS",
  corr_method = "decomposed",
  fix_price   = TRUE,
  fix_yield   = TRUE,
  label       = "deterministic",
  seed_base   = 603L
)
saveRDS(credit_shapley_none, "results/credit/credit_shapley_none_2.rds")

# Summarise all four grids and compute Shapley values
delinq_both  <- summarise_grid(credit_corr)[,  variant := "both"]
delinq_price <- summarise_grid(credit_shapley_price)[, variant := "price_only"]
delinq_yield <- summarise_grid(credit_shapley_yield)[, variant := "yield_only"]
delinq_none  <- summarise_grid(credit_shapley_none)[,  variant := "none"]

shapley_long <- rbindlist(list(delinq_both, delinq_price, delinq_yield, delinq_none))

shapley_wide <- dcast(
  shapley_long[, .(region, prod_zone, plan, variant, pr_delinq)],
  region + prod_zone + plan ~ variant,
  value.var = "pr_delinq"
)

shapley_wide[, shapley_price := 0.5 * (price_only - none) + 0.5 * (both - yield_only)]
shapley_wide[, shapley_yield := 0.5 * (yield_only - none) + 0.5 * (both - price_only)]
shapley_wide[, share_price   := shapley_price / (shapley_price + shapley_yield)]
shapley_wide[, share_yield   := shapley_yield / (shapley_price + shapley_yield)]
shapley_wide[, region    := factor(region,    levels = regions_vec)]
shapley_wide[, prod_zone := factor(prod_zone, levels = zones_vec)]

cat("\n=== Shapley decomposition results ===\n")
print(shapley_wide[order(plan, region, prod_zone),
  .(region, prod_zone, plan,
    shapley_price_pp = round(shapley_price * 100, 2),
    shapley_yield_pp = round(shapley_yield * 100, 2),
    share_price_pct  = round(share_price  * 100, 1))])
saveRDS(shapley_wide, "results/credit/shapley_decomposition.rds")
fwrite(shapley_wide,  "results/credit/shapley_decomposition.csv")

# Figure 1: price share heatmap across region x zone x plan
p_shapley_share <- shapley_wide |>
  ggplot(aes(x = prod_zone, y = region, fill = share_price * 100)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", share_price * 100)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 50, limits = c(0, 100),
                       name = "Price share\n(%)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "Shapley decomposition: share of delinquency risk from price uncertainty",
       subtitle = "Red = price-dominated | Blue = yield-dominated | Correlated prices | history = CCCCCC",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_shapley_share)
ggsave(file.path(FIGS_DIR, "fig_shapley_price_share.pdf"),
       p_shapley_share, width = 9, height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_shapley_price_share.png"),
       p_shapley_share, width = 9, height = 5, dpi = 300)

# Figure 2: stacked bar of absolute Shapley values (pp) by region x zone x plan
shapley_melt <- melt(shapley_wide,
  id.vars      = c("region", "prod_zone", "plan"),
  measure.vars = c("shapley_price", "shapley_yield"),
  variable.name = "source", value.name = "shapley_pp")
shapley_melt[, source := ifelse(source == "shapley_price", "Price", "Yield")]
shapley_melt[, region    := factor(region,    levels = regions_vec)]
shapley_melt[, prod_zone := factor(prod_zone, levels = zones_vec)]

p_shapley_abs <- shapley_melt |>
  ggplot(aes(x = prod_zone, y = shapley_pp * 100, fill = source)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(Price = "firebrick", Yield = "steelblue"),
                    name = "Uncertainty\nsource") +
  facet_grid(region ~ plan) +
  labs(title = "Shapley decomposition: absolute delinquency risk (pp) by uncertainty source",
       subtitle = "Correlated prices | history = CCCCCC",
       x = "Productivity zone", y = "Shapley value (pp)") +
  theme_bw()
print(p_shapley_abs)
ggsave(file.path(FIGS_DIR, "fig_shapley_abs.pdf"),
       p_shapley_abs, width = 9, height = 6, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_shapley_abs.png"),
       p_shapley_abs, width = 9, height = 6, dpi = 300)



# ---- M. Revenue variance decomposition: rotation as natural hedge -----------
# Structural model:
#   log R_t = log P_t + log Y_t
#   Var(log R_t) = sigma_P^2 + sigma_Y^2 + 2*rho*sqrt(f)*sigma_P*sigma_Y
#   Covariance term is negative when rho < 0 (price-yield hedge).
#   Plans with more soy years accumulate more negative covariance across the
#   6-year horizon => lower effective revenue variance => lower delinquency.
#   hedge_efficiency = -cov_total / (price_total + yield_total): fraction of
#   gross variance offset by the natural hedge.
#   Computed analytically from parameters — no new simulations.
aggregate_fraction <- 0.40   # matches simulate_plan_credit_region_zone default

decompose_revenue_variance <- function(history6, plan6, geo_region, prod_zone,
                                       r_disc = 0.05) {
  price_pars <- price_params_regional[[geo_region]]
  corn_gen   <- yield_gens_by_region_zone$corn[[geo_region]][[prod_zone]]
  soy_gen    <- yield_gens_by_region_zone$soy[[geo_region]][[prod_zone]]
  rc         <- region_corr[geo_region]
  rho_corn   <- rc$corr_corn
  rho_soy    <- rc$corr_soy
  disc       <- disc_vec(r_disc)

  path <- make_transition_path(history6, plan6)

  rows <- lapply(1:6, function(t) {
    crop_t <- path$crop[t]
    st_t   <- path$window6[t]
    if (crop_t == "C") {
      sigma_P <- price_pars$C$sigma
      mu_Y    <- corn_gen$mu[rot6 == st_t, mu]
      sigma_Y <- sd(corn_gen$pools[rot6 == st_t, resid_pool][[1]]) / mu_Y  # CV
      rho     <- rho_corn
    } else {
      sigma_P <- price_pars$S$sigma
      mu_Y    <- soy_gen$mu[rot6 == st_t, mu]
      sigma_Y <- sd(soy_gen$pools[rot6 == st_t, resid_pool][[1]]) / mu_Y   # CV
      rho     <- rho_soy
    }
    # Decomposition in CV² units: yield residuals are in level space (bu/ac),
    # so dividing by mean yield converts to CV, making price and yield terms
    # comparable. Var(R)/E[R]² ≈ sigma_P² + CV_Y² + 2·rho·sqrt(f)·sigma_P·CV_Y
    var_P   <- sigma_P^2
    var_Y   <- sigma_Y^2   # CV_Y^2
    var_cov <- 2 * rho * sqrt(aggregate_fraction) * sigma_P * sigma_Y
    data.table(year = t, crop = crop_t, sigma_P, sigma_Y, rho,
               var_P, var_Y, var_cov,
               var_total = var_P + var_Y + var_cov,
               disc2 = disc[t]^2)   # variance of discounted sum uses disc^2
  })

  dt <- rbindlist(rows)

  # Plan-level: discount-weighted sum of per-year variances
  totals <- dt[, .(
    var_price = sum(var_P   * disc2),
    var_yield = sum(var_Y   * disc2),
    var_cov   = sum(var_cov * disc2)
  )]
  totals[, var_gross       := var_price + var_yield]
  totals[, var_net         := var_gross + var_cov]
  totals[, hedge_efficiency := -var_cov / var_gross]   # >0 when rho<0
  totals[, share_price     := var_price / var_gross]
  totals[, share_yield     := var_yield / var_gross]
  totals[, region := geo_region][, prod_zone := prod_zone][, plan := plan6]
  totals
}

rev_decomp <- rbindlist(lapply(regions_vec, function(reg) {
  rbindlist(lapply(zones_vec, function(zone) {
    rbindlist(lapply(rotation_plans, function(plan) {
      tryCatch(
        decompose_revenue_variance("CCCCCC", plan, reg, zone),
        error = function(e) {
          cat(sprintf("  ERROR %s/%s/%s: %s\n", reg, zone, plan, e$message)); NULL
        }
      )
    }))
  }))
}))

cat("\n=== Revenue variance decomposition by rotation plan ===\n")
print(rev_decomp[order(plan, region, prod_zone),
  .(region, prod_zone, plan,
    var_price  = round(var_price,  4),
    var_yield  = round(var_yield,  4),
    var_cov    = round(var_cov,    4),
    hedge_eff  = round(hedge_efficiency, 3),
    share_P    = round(share_price, 3))])
rev_decomp[, region    := factor(region,    levels = regions_vec)]
rev_decomp[, prod_zone := factor(prod_zone, levels = zones_vec)]
saveRDS(rev_decomp, "results/credit/revenue_variance_decomp.rds")
fwrite(rev_decomp,  "results/credit/revenue_variance_decomp.csv")

# Figure 1a: normalized shares — what fraction of gross variance is price vs. yield?
# Covariance shown as hedge efficiency (separate, avoids mixing positive/negative bars)
rev_decomp[, share_price_pct := share_price * 100]
rev_decomp[, share_yield_pct := (1 - share_price) * 100]

rev_shares <- melt(rev_decomp,
  id.vars      = c("region", "prod_zone", "plan"),
  measure.vars = c("share_price_pct", "share_yield_pct"),
  variable.name = "component", value.name = "share")
rev_shares[, component := factor(component,
  levels = c("share_yield_pct", "share_price_pct"),
  labels = c("Yield", "Price"))]
rev_shares[, region    := factor(region,    levels = regions_vec)]
rev_shares[, prod_zone := factor(prod_zone, levels = zones_vec)]

p_rev_decomp <- rev_shares |>
  ggplot(aes(x = plan, y = share, fill = component)) +
  geom_col() +
  scale_fill_manual(values = c("Price" = "firebrick", "Yield" = "steelblue"),
                    name = "Component") +
  facet_grid(region ~ prod_zone) +
  labs(title = "Revenue risk composition: price vs. yield share of gross variance",
       subtitle = "Excludes covariance offset — see hedge efficiency figure",
       x = "Rotation plan", y = "Share of gross variance (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_rev_decomp)
ggsave(file.path(FIGS_DIR, "fig_revenue_variance_decomp.pdf"),
       p_rev_decomp, width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_revenue_variance_decomp.png"),
       p_rev_decomp, width = 9, height = 7, dpi = 300)

# Figure 1b: total CV of discounted revenue by plan — absolute risk scale
# CV(R) = sqrt(exp(var_net) - 1) for log-normal; var_net already discounted
rev_decomp[, cv_revenue := sqrt(pmax(var_net, 0)) * 100]   # CV in %; var_net = Var(R)/E[R]^2

p_rev_cv <- rev_decomp |>
  ggplot(aes(x = plan, y = cv_revenue, fill = plan)) +
  geom_col() +
  scale_fill_manual(values = plan_colors, guide = "none") +
  facet_grid(region ~ prod_zone) +
  labs(title = "Total revenue risk: CV of discounted revenue by rotation plan",
       subtitle = "Includes price-yield hedge offset (var_net = gross variance + covariance)",
       x = "Rotation plan", y = "CV of discounted revenue (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_rev_cv)
ggsave(file.path(FIGS_DIR, "fig_revenue_cv.pdf"),
       p_rev_cv, width = 9, height = 7, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_revenue_cv.png"),
       p_rev_cv, width = 9, height = 7, dpi = 300)

# Figure 2: hedge efficiency heatmap — fraction of gross variance offset per plan
p_hedge_eff <- rev_decomp |>
  ggplot(aes(x = prod_zone, y = region, fill = hedge_efficiency * 100)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", hedge_efficiency * 100)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Hedge\nefficiency (%)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~plan) +
  labs(title = "Natural hedge efficiency: % of gross revenue variance offset by price-yield covariance",
       subtitle = "Higher = more protection from rotation | history = CCCCCC",
       x = "Productivity zone", y = "Region") +
  theme_bw()
print(p_hedge_eff)
ggsave(file.path(FIGS_DIR, "fig_hedge_efficiency.pdf"),
       p_hedge_eff, width = 9, height = 5, device = grDevices::cairo_pdf)
ggsave(file.path(FIGS_DIR, "fig_hedge_efficiency.png"),
       p_hedge_eff, width = 9, height = 5, dpi = 300)
gc()

       

## Dynamic programming for seeding strategy optimization: run on a single cell for inspection, then run across the full grid and compare to fixed plans in credit space. The DP
if (isTRUE(RUN_DP_EXTENSION)) {   # extension not reported in the paper; see README
source("R/dp_functions_region_zone.R")

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

subsidy_tbl <- tribble(
  ~theta, ~basic_optional, ~enterprise,
  0.50,   0.67,            0.80,
  0.55,   0.64,            0.80,
  0.60,   0.64,            0.80,
  0.65,   0.59,            0.80,
  0.70,   0.59,            0.80,
  0.75,   0.55,            0.77,
  0.80,   0.48,            0.68,
  0.85,   0.38,            0.53
)

# --- With insurance in the DP objective ---
# Supply subsidy_tbl from price_adj.R's subsidy_tbl object
dp_grid_ins <- run_dp_grid(
    use_insurance = TRUE,
    subsidy_tbl   = subsidy_tbl,    # bring in from price_adj.R environment
    thetas        = c(0.65, 0.75, 0.85)
  )
# Optimal coverage level by year for central/medium
rollout_dp_plan(dp_grid_ins$dp_grid$central$medium, "CCCCCC")$path[, .(t, crop, theta)]
}  # end RUN_DP_EXTENSION
