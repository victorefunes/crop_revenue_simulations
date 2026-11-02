# ============================================================================
# STAGE 2 (breakeven paper): probability that revenue falls below cost, state by state.
#
# Structural model
#   inputs    : state x NCCPI-zone yield generator conditioned on the 6-year rotation history,
#               state price (lognormal, mean + log-sd), state total non-land cost (lognormal,
#               sd 0.10) plus a same-crop penalty (corn after corn +$15, soy after soy +$5),
#               price-yield dependence (state median county correlation; decomposed method)
#   mechanism : R_t = P_t * Y_t ;  K_t = cost_t ;  year t is a loss if R_t < K_t
#   outcome   : Pr(loss) per year, Pr(at least one loss in 6 years), expected loss years,
#               compared across rotations CCCCCC, CCSCCS, CSCSCS, within each state x zone.
# Each state is simulated on its own (never pooled).
#
# Run from simulation_code/ after build_state_generators.R.
# Env: SIM_CACHE_DIR (generators), SIM_B (paths per cell, default 50,000), SIM_STATES, SIM_OUT_DIR.
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
source("config.R")
source("R/functions_region_zone.R")
source("R/correlation_functions_region_zone.R")
source("R/credit_functions_region_zone.R")

OUT_DIR <- .env("SIM_OUT_DIR", "results/breakeven")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
STATES <- strsplit(.env("SIM_STATES", "IA,IL,IN,MN,MO,OH,SD,WI"), ",")[[1]]
ZONES  <- c("high", "medium", "low")
PLANS  <- c("CCCCCC", "CCSCCS", "CSCSCS")
B      <- SIM_B
COST_SD <- 0.10
SAME_CROP_PENALTY <- c(C = 15, S = 5)    # $/ac extra direct cost when the crop follows itself (farmdoc 2025)
GLOBAL_SEED <- 21212

gens   <- readRDS(file.path(CACHE_DIR, "state_generators.rds"))
costs  <- fread("data/state_costs.csv")
prices <- fread("data/state_prices.csv")

# ---- parameter lists in the layout the simulator expects --------------------
price_pars <- lapply(setNames(STATES, STATES), function(st) {
  mk <- function(cr) { p <- prices[state == st & crop == cr]
    list(Pbar = p$price_per_bu, sigma = p$log_sd, mu = log(p$price_per_bu) - 0.5 * p$log_sd^2) }
  list(C = mk("C"), S = mk("S"))
})
cost_pars <- lapply(setNames(STATES, STATES), function(st) {
  mk <- function(cr) { k <- costs[state == st & crop == cr]$total_nonland_cost
    list(Cbar = k, sigma = COST_SD, mu = log(k) - 0.5 * COST_SD^2) }
  z <- list(C = mk("C"), S = mk("S"))
  setNames(list(z, z, z), ZONES)
})

# ---- state price-yield correlations (median of county correlations) ---------
fips_state <- c(IL = 17, IN = 18, IA = 19, MN = 27, MO = 29, OH = 39, SD = 46, WI = 55)
med_corr <- function(f, def) { d <- fread(f, select = c("fips", "corr"))
  sapply(STATES, function(st) { x <- d[floor(fips / 1000) == fips_state[[st]], corr]
    if (length(x)) median(x, na.rm = TRUE) else def }) }
corr_corn <- med_corr(CORN_CORR_CSV, -0.40)
corr_soy  <- med_corr(SOY_CORR_CSV,  -0.50)
print(data.table(state = STATES, corr_corn, corr_soy))

# ---- one cell ---------------------------------------------------------------
run_cell <- function(st, zone, plan, corr_method, seed) {
  o <- simulate_plan_credit_region_zone(
    history6 = plan, plan6 = plan, geo_region = st, prod_zone = zone, B = B,
    gens_by_region_zone = gens, price_pars_regional = price_pars, base_cost_pars = cost_pars,
    seed = seed, credit_limit = 2000, delinquency_rule = "end_balance_positive",
    correlation_corn = corr_corn[[st]], correlation_soy = corr_soy[[st]],
    corr_method = corr_method, aggregate_fraction = 0.40)
  ch <- strsplit(plan, "")[[1]]
  prev <- c(ch[6], ch[1:5])                         # steady state: history equals the plan
  pen <- ifelse(ch == prev, SAME_CROP_PENALTY[ch], 0)
  K <- sweep(o$cost, 2, pen, "+")
  loss <- o$rev < K
  net  <- o$rev - K
  data.table(state = st, zone = zone, plan = plan, corr = corr_method,
             p_loss_any = mean(rowSums(loss) > 0), exp_loss_years = mean(rowSums(loss)),
             p_loss_year = mean(loss),
             p_loss_corn_yr = if (any(ch == "C")) mean(loss[, ch == "C"]) else NA_real_,
             p_loss_soy_yr  = if (any(ch == "S")) mean(loss[, ch == "S"]) else NA_real_,
             mean_rev = mean(o$rev), mean_cost = mean(K), mean_net = mean(net),
             p10_net = quantile(net, 0.10), mean_cum_net = mean(rowSums(net)),
             p_loss_cum = mean(rowSums(net) < 0),
             t(setNames(colMeans(loss), paste0("p_loss_t", 1:6))))
}

res <- list(); k <- 0L
for (si in seq_along(STATES)) for (zi in seq_along(ZONES)) {
  seed <- si * 1000L + zi                           # shared across plans in a cell (common random numbers)
  for (cm in c("decomposed", "none")) for (pl in PLANS) {
    k <- k + 1L
    res[[k]] <- run_cell(STATES[si], ZONES[zi], pl, cm, seed)
  }
  cat(sprintf("done %s-%s\n", STATES[si], ZONES[zi]))
}
res <- rbindlist(res)
fwrite(res, file.path(OUT_DIR, "breakeven_cells.csv"))

# ---- analytical breakeven thresholds ----------------------------------------
be <- merge(costs, prices[, .(state, crop, price_per_bu)], by = c("state", "crop"))
be[, breakeven_yield := total_nonland_cost / price_per_bu]
be[, breakeven_price_at_assumed_yield := total_nonland_cost / assumed_yield_bu]
fwrite(be, file.path(OUT_DIR, "breakeven_thresholds.csv"))
print(res[corr == "decomposed", .(state, zone, plan, p_loss_any = round(p_loss_any, 3),
                                  p_loss_year = round(p_loss_year, 3), mean_net = round(mean_net))])
