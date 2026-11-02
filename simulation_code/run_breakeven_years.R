# ============================================================================
# STAGE 2 (breakeven paper, year-shock version): Pr(revenue < cost) in good vs bad years, and a Shapley
# decomposition of price, yield and cost risk. State by state, never pooled.
#
# Structural model
#   inputs    : state x NCCPI-zone generator (field-effect mean by rotation state + year effects + residual
#               pool), state price and total non-land cost, price-yield dependence through the year shock
#   mechanism : R_t = P_t * Y_t ; year t is a shortfall if R_t < K_t ; a year is good/normal/bad by the tercile
#               of its weather index (standardized corn + soybean year effects)
#   outcome   : Pr(R_t < K_t | year type), by calendar year, annual, six-year cumulative; Shapley shares of
#               the variance of net return and of the shortfall probability
# Run from simulation_code/ after build_state_generators.R (generators must contain year_fe).
# Env: SIM_CACHE_DIR, SIM_B (default 50,000), SIM_STATES, SIM_OUT_DIR (default results/breakeven_years).
# ============================================================================
suppressPackageStartupMessages(library(data.table))
source("config.R")
source("R/functions_region_zone.R")      # make_transition_path()
source("R/breakeven_sim.R")

OUT_DIR <- .env("SIM_OUT_DIR", "results/breakeven_years")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
STATES <- strsplit(.env("SIM_STATES", "IA,IL,IN,MN,MO,OH,SD,WI"), ",")[[1]]
ZONES  <- c("high", "medium", "low")
PLANS  <- c("CCCCCC", "CCSCCS", "CSCSCS")
B <- SIM_B; COST_SD <- 0.10

gens   <- readRDS(file.path(CACHE_DIR, "state_generators.rds"))
costs  <- fread("data/state_costs.csv"); prices <- fread("data/state_prices.csv")
price_pars <- lapply(setNames(STATES, STATES), function(st) {
  mk <- function(cr) { p <- prices[state == st & crop == cr]
    list(Pbar = p$price_per_bu, sigma = p$log_sd, mu = log(p$price_per_bu) - 0.5 * p$log_sd^2) }
  list(C = mk("C"), S = mk("S")) })
cost_pars <- lapply(setNames(STATES, STATES), function(st) {
  mk <- function(cr) list(Cbar = costs[state == st & crop == cr]$total_nonland_cost, sigma = COST_SD)
  list(C = mk("C"), S = mk("S")) })
fips_state <- c(IL = 17, IN = 18, IA = 19, MN = 27, MO = 29, OH = 39, SD = 46, WI = 55)
med_corr <- function(f, def) { d <- fread(f, select = c("fips", "corr"))
  sapply(STATES, function(st) { x <- d[floor(fips / 1000) == fips_state[[st]], corr]
    if (length(x)) median(x, na.rm = TRUE) else def }) }
rho_c <- med_corr(CORN_CORR_CSV, -0.40); rho_s <- med_corr(SOY_CORR_CSV, -0.50)
print(data.table(state = STATES, rho_corn = rho_c, rho_soy = rho_s))

main <- list(); byyr <- list(); shap <- list(); dep <- list()
for (si in seq_along(STATES)) for (zi in seq_along(ZONES)) {
  st <- STATES[si]; z <- ZONES[zi]; seed <- si * 1000L + zi
  cg <- gens$corn[[st]][[z]]; sg <- gens$soy[[st]][[z]]
  wt <- weather_types(cg$year_fe, sg$year_fe)
  sim <- function(plan, rc = rho_c[[st]], rs = rho_s[[st]], ...)
    simulate_breakeven(cg, sg, plan, price_pars[[st]], cost_pars[[st]], rc, rs, B = B, seed = seed, wt = wt, ...)
  for (pl in PLANS) {
    o <- sim(pl); loss <- o$rev < o$cost; net <- o$rev - o$cost; ty <- matrix(wt$type[o$j], nrow(o$j))
    row <- data.table(state = st, zone = z, plan = pl, p_loss_year = mean(loss),
                      p_loss_cum = mean(rowSums(net) < 0), p_loss_any = mean(rowSums(loss) > 0), mean_net = mean(net))
    for (tp in c("bad", "normal", "good"))
      row[, (paste0("p_", tp)) := mean(loss[ty == tp])][, (paste0("net_", tp)) := mean(net[ty == tp])]
    main[[length(main) + 1]] <- row
    byyr[[length(byyr) + 1]] <- data.table(state = st, zone = z, plan = pl, year = wt$years, type = wt$type,
      p_loss = sapply(seq_along(wt$years), function(k) mean(loss[o$j == k])),
      mean_net = sapply(seq_along(wt$years), function(k) mean(net[o$j == k])))
    o0 <- sim(pl, 0, 0)                                  # no price-yield dependence
    dep[[length(dep) + 1]] <- data.table(state = st, zone = z, plan = pl, p_loss_year_nodep = mean(o0$rev < o0$cost))
    # Shapley over price / yield / cost uncertainty, two functionals
    v_p <- function(random) { o <- sim(pl, fix_price = !"price" %in% random, fix_yield = !"yield" %in% random,
                                       fix_cost = !"cost" %in% random); mean(o$rev < o$cost) }
    v_v <- function(random) { o <- sim(pl, fix_price = !"price" %in% random, fix_yield = !"yield" %in% random,
                                       fix_cost = !"cost" %in% random); var(as.vector(o$rev - o$cost)) }
    sp <- shapley3(v_p); sv <- shapley3(v_v)
    shap[[length(shap) + 1]] <- data.table(state = st, zone = z, plan = pl,
      p_empty = sp$v_empty, p_all = sp$v_all, phi_p_price = sp$phi[["price"]], phi_p_yield = sp$phi[["yield"]], phi_p_cost = sp$phi[["cost"]],
      var_all = sv$v_all, phi_v_price = sv$phi[["price"]], phi_v_yield = sv$phi[["yield"]], phi_v_cost = sv$phi[["cost"]])
  }
  cat(sprintf("done %s-%s\n", st, z))
}
fwrite(rbindlist(main),  file.path(OUT_DIR, "breakeven_year_types.csv"))
fwrite(rbindlist(byyr),  file.path(OUT_DIR, "breakeven_by_calendar_year.csv"))
fwrite(rbindlist(shap),  file.path(OUT_DIR, "shapley.csv"))
fwrite(rbindlist(dep),   file.path(OUT_DIR, "dependence_check.csv"))
print(rbindlist(main)[zone == "medium", .(state, plan, annual = round(p_loss_year, 2), bad = round(p_bad, 2),
                                          normal = round(p_normal, 2), good = round(p_good, 2))])
