# ============================================================================
# Direct (operating) cost variants. Cost of a drawn year = nondirect_state + direct_state * r_y (non-direct kept at the state budget).
#   baseline       : state 2025 budget, direct and non-direct fixed
#   direct_avg     : direct = ERS regional real operating cost averaged over the state's panel years (constant); non-direct fixed
#   direct_raw     : direct varies by year with the ERS real operating cost (mean-one index), state direct level
#   direct_detrend : as direct_raw after removing a linear time trend from the ERS real operating cost (mean kept)
# Medium-zone cell of each state, seeds of run_breakeven_years.R. Run from simulation_code/ (SIM_CACHE_DIR, SIM_B).
# Needs data/ers_costs_real.csv. Writes results/breakeven_years/cost_direct.csv and cost_direct_trend.csv.
# ============================================================================
suppressPackageStartupMessages(library(data.table))
source("config.R")
source("R/functions_region_zone.R")
source("R/breakeven_sim.R")
OUT_DIR <- "results/breakeven_years"
STATES <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
PLANS <- c("CCCCCC", "CCSCCS", "CSCSCS"); B <- SIM_B; COST_SD <- 0.10
rmap <- c(IA = "Heartland", IL = "Heartland", IN = "Heartland", MO = "Heartland", OH = "Heartland",
          WI = "Northern Crescent", MN = "Northern Great Plains", SD = "Northern Great Plains")
gens <- readRDS(file.path(CACHE_DIR, "state_generators.rds"))
costs <- fread("data/state_costs.csv"); prices <- fread("data/state_prices.csv"); ers <- fread("data/ers_costs_real.csv")
fips_state <- c(IL = 17, IN = 18, IA = 19, MN = 27, MO = 29, OH = 39, SD = 46, WI = 55)
med_corr <- function(f, def) { d <- fread(f, select = c("fips", "corr"))
  sapply(STATES, function(st) { x <- d[floor(fips / 1000) == fips_state[[st]], corr]
    if (length(x)) median(x, na.rm = TRUE) else def }) }
rho_c <- med_corr(CORN_CORR_CSV, -0.40); rho_s <- med_corr(SOY_CORR_CSV, -0.50)
scen <- c("baseline", "direct_avg", "direct_raw", "direct_detrend")

out <- list(); trend <- list()
for (si in seq_along(STATES)) {
  st <- STATES[si]; seed <- si * 1000L + 2L
  cg <- gens$corn[[st]][["medium"]]; sg <- gens$soy[[st]][["medium"]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  yrs <- wt$years
  pp <- lapply(c(C = "C", S = "S"), function(cr) { p <- prices[state == st & crop == cr]
    list(Pbar = p$price_per_bu, sigma = p$log_sd, mu = log(p$price_per_bu) - 0.5 * p$log_sd^2) })
  for (cr in c("C", "S")) {
    x <- ers[region == rmap[[st]] & crop == cr & year %in% yrs]
    fit <- lm(operating_real ~ year, data = x)
    trend[[length(trend) + 1]] <- data.table(state = st, crop = cr, mean_real = mean(x$operating_real),
      slope = coef(fit)[["year"]], r2 = summary(fit)$r.squared, cv = sd(x$operating_real) / mean(x$operating_real),
      cv_detrended = sd(resid(fit)) / mean(x$operating_real))
  }
  for (sn in scen) {
    cp <- lapply(c(C = "C", S = "S"), function(cr) {
      b <- costs[state == st & crop == cr]; K <- b$total_nonland_cost; d <- b$direct_cost; nd <- K - d
      x <- ers[region == rmap[[st]] & crop == cr]; x <- x[match(yrs, year)]
      op <- setNames(x$operating_real, as.character(yrs))
      op_dt <- setNames(as.numeric(resid(lm(op ~ yrs))) + mean(op), as.character(yrs))
      switch(sn, baseline = list(Cbar = K, sigma = COST_SD),
        direct_avg = list(Cbar = nd + mean(op), sigma = COST_SD),
        direct_raw = list(Cbar = K, sigma = COST_SD, idx = (nd + d * op / mean(op)) / K),
        direct_detrend = list(Cbar = K, sigma = COST_SD, idx = (nd + d * op_dt / mean(op_dt)) / K))
    })
    for (pl in PLANS) {
      o <- simulate_breakeven(cg, sg, pl, pp, cp, rho_c[[st]], rho_s[[st]], B = B, seed = seed, wt = wt)
      loss <- o$rev < o$cost; net <- o$rev - o$cost; ty <- matrix(wt$type[o$j], nrow(o$j))
      out[[length(out) + 1]] <- data.table(state = st, scenario = sn, plan = pl, p_year = mean(loss), p_cum = mean(rowSums(net) < 0),
        p_bad = mean(loss[ty == "bad"]), p_normal = mean(loss[ty == "normal"]), p_good = mean(loss[ty == "good"]), mean_cost = mean(o$cost))
    }
  }
  cat("done", st, "\n")
}
res <- rbindlist(out); fwrite(res, file.path(OUT_DIR, "cost_direct.csv"))
tr <- rbindlist(trend); fwrite(tr, file.path(OUT_DIR, "cost_direct_trend.csv"))
print(tr[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
w <- dcast(res, state + scenario ~ plan, value.var = c("p_year", "p_bad", "p_good", "p_cum"))
w[, gap := p_year_CSCSCS - p_year_CCCCCC]
sh <- function(v) print(dcast(w, state ~ scenario, value.var = v)[, lapply(.SD, function(x) if (is.numeric(x)) round(100 * x, 1) else x)])
sh("p_year_CCCCCC"); sh("gap"); sh("p_bad_CCCCCC"); sh("p_good_CCCCCC"); sh("p_cum_CCCCCC")
print(dcast(res[plan == "CCCCCC"], state ~ scenario, value.var = "mean_cost")[, lapply(.SD, function(x) if (is.numeric(x)) round(x) else x)])
