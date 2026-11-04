# ============================================================================
# Does using year-specific costs (ERS regional real non-land cost, 2008-2022) change the results?
# Medium-zone cell of each state, seeds of run_breakeven_years.R. Cost of a drawn calendar year j is Cbar * idx[j] * lognormal noise.
#   baseline     : one 2025 budget cost per state (idx = 1)
#   ers_timing   : state 2025 budget level, ERS real-cost timing (idx = ERS_j / mean over the state's panel years)
#   ers_level    : state 2025 budget scaled back by ERS real cost of year j relative to 2025 (historical real level)
#   ers_region   : ERS regional real cost itself (replaces the state budget; assumed non-direct costs not used)
#   ers_heartland: as ers_timing, with the Heartland series for every state
# State -> ERS region: IA/IL/IN/MO/OH Heartland, WI Northern Crescent, MN/SD Northern Great Plains.
# Run from simulation_code/. Env: SIM_CACHE_DIR, SIM_B. Needs data/ers_costs_real.csv (build_cost_index.R).
# Writes results/breakeven_years/cost_timing.csv.
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
scen <- c("baseline", "ers_timing", "ers_level", "ers_region", "ers_heartland")

out <- list()
for (si in seq_along(STATES)) {
  st <- STATES[si]; seed <- si * 1000L + 2L
  cg <- gens$corn[[st]][["medium"]]; sg <- gens$soy[[st]][["medium"]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  yrs <- as.character(wt$years)
  pp <- lapply(c(C = "C", S = "S"), function(cr) { p <- prices[state == st & crop == cr]
    list(Pbar = p$price_per_bu, sigma = p$log_sd, mu = log(p$price_per_bu) - 0.5 * p$log_sd^2) })
  for (sn in scen) {
    cp <- lapply(c(C = "C", S = "S"), function(cr) {
      budget <- costs[state == st & crop == cr]$total_nonland_cost
      reg <- if (sn == "ers_heartland") "Heartland" else rmap[[st]]
      x <- ers[ers$region == reg & ers$crop == cr]
      R <- setNames(x$nonland_real, as.character(x$year))[yrs]
      R25 <- x[year == 2025]$nonland_real
      switch(sn, baseline = list(Cbar = budget, sigma = COST_SD),
        ers_timing = , ers_heartland = list(Cbar = budget, sigma = COST_SD, idx = R / mean(R)),
        ers_level = list(Cbar = budget, sigma = COST_SD, idx = R / R25),
        ers_region = list(Cbar = mean(R), sigma = COST_SD, idx = R / mean(R)))
    })
    for (pl in PLANS) {
      o <- simulate_breakeven(cg, sg, pl, pp, cp, rho_c[[st]], rho_s[[st]], B = B, seed = seed, wt = wt)
      loss <- o$rev < o$cost; net <- o$rev - o$cost; ty <- matrix(wt$type[o$j], nrow(o$j))
      out[[length(out) + 1]] <- data.table(state = st, scenario = sn, plan = pl, p_year = mean(loss), p_cum = mean(rowSums(net) < 0),
        p_bad = mean(loss[ty == "bad"]), p_normal = mean(loss[ty == "normal"]), p_good = mean(loss[ty == "good"]),
        mean_cost = mean(o$cost), mean_net = mean(net))
    }
  }
  cat("done", st, "\n")
}
res <- rbindlist(out); fwrite(res, file.path(OUT_DIR, "cost_timing.csv"))
w <- dcast(res, state + scenario ~ plan, value.var = c("p_year", "p_bad", "p_good"))
w[, `:=`(gap = p_year_CSCSCS - p_year_CCCCCC)]
sh <- function(v) print(dcast(w, state ~ scenario, value.var = v)[, lapply(.SD, function(x) if (is.numeric(x)) round(100 * x, 1) else x)])
sh("p_year_CCCCCC"); sh("gap"); sh("p_bad_CCCCCC"); sh("p_good_CCCCCC")
