# ============================================================================
# Appendix robustness (breakeven paper): re-run the medium-zone cell of each state with one input changed.
#   baseline | costs -10% | costs +10% | costs +20% | price log-sd x0.75 | price log-sd x1.25 | no same-crop penalty
# Each scenario changes one parameter that enters R_t < K_t directly (Cbar, sigma_P, or the penalty added to K_t),
# keeps the seeds of run_breakeven_years.R (so "baseline" reproduces the main results), and reports the annual and
# six-year probabilities by plan. Also writes the state inputs (rho) for the parameter table.
# Run from simulation_code/. Env: SIM_CACHE_DIR, SIM_B. Writes results/breakeven_years/robustness_inputs.csv and
# results/breakeven_years/state_inputs.csv.
# ============================================================================
suppressPackageStartupMessages(library(data.table))
source("config.R")
source("R/functions_region_zone.R")
source("R/breakeven_sim.R")
OUT_DIR <- "results/breakeven_years"
STATES <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI"); ZONES <- c("high", "medium", "low")
PLANS <- c("CCCCCC", "CCSCCS", "CSCSCS"); B <- SIM_B; COST_SD <- 0.10
gens <- readRDS(file.path(CACHE_DIR, "state_generators.rds"))
costs <- fread("data/state_costs.csv"); prices <- fread("data/state_prices.csv")
fips_state <- c(IL = 17, IN = 18, IA = 19, MN = 27, MO = 29, OH = 39, SD = 46, WI = 55)
med_corr <- function(f, def) { d <- fread(f, select = c("fips", "corr"))
  sapply(STATES, function(st) { x <- d[floor(fips / 1000) == fips_state[[st]], corr]
    if (length(x)) median(x, na.rm = TRUE) else def }) }
rho_c <- med_corr(CORN_CORR_CSV, -0.40); rho_s <- med_corr(SOY_CORR_CSV, -0.50)
fwrite(data.table(state = STATES, rho_corn = rho_c, rho_soy = rho_s), file.path(OUT_DIR, "state_inputs.csv"))

scen <- list(baseline = list(), cost_m10 = list(cm = 0.9), cost_p10 = list(cm = 1.1), cost_p20 = list(cm = 1.2),
             vol_x075 = list(vm = 0.75), vol_x125 = list(vm = 1.25), no_penalty = list(pen = c(C = 0, S = 0)))
out <- list()
for (si in seq_along(STATES)) {
  st <- STATES[si]; zi <- match("medium", ZONES); seed <- si * 1000L + zi
  cg <- gens$corn[[st]][["medium"]]; sg <- gens$soy[[st]][["medium"]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  for (sn in names(scen)) {
    s <- scen[[sn]]; cm <- if (is.null(s$cm)) 1 else s$cm; vm <- if (is.null(s$vm)) 1 else s$vm
    pen <- if (is.null(s$pen)) c(C = 15, S = 5) else s$pen
    pp <- lapply(c(C = "C", S = "S"), function(cr) { p <- prices[state == st & crop == cr]; sg_ <- p$log_sd * vm
      list(Pbar = p$price_per_bu, sigma = sg_, mu = log(p$price_per_bu) - 0.5 * sg_^2) })
    cp <- lapply(c(C = "C", S = "S"), function(cr) list(Cbar = costs[state == st & crop == cr]$total_nonland_cost * cm, sigma = COST_SD))
    for (pl in PLANS) {
      o <- simulate_breakeven(cg, sg, pl, pp, cp, rho_c[[st]], rho_s[[st]], B = B, seed = seed, wt = wt, same_crop_penalty = pen)
      loss <- o$rev < o$cost; net <- o$rev - o$cost
      out[[length(out) + 1]] <- data.table(state = st, scenario = sn, plan = pl, p_year = mean(loss), p_cum = mean(rowSums(net) < 0))
    }
  }
  cat("done", st, "\n")
}
res <- rbindlist(out); fwrite(res, file.path(OUT_DIR, "robustness_inputs.csv"))
w <- dcast(res, state + scenario ~ plan, value.var = "p_year")[, gap := CCCCCC - CSCSCS]
print(dcast(w, state ~ scenario, value.var = "gap")[, lapply(.SD, function(x) if (is.numeric(x)) round(100 * x, 1) else x)])
