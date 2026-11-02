# Breakeven probability by VPD regime (low / mid / high July-August VPD), state by state.
# Two variants per cell, on common random numbers for the year draw:
#   reclass : generator unchanged (window pools); calendar years are labelled by their VPD regime
#   pools   : the field residual is drawn from the pool for (six-year window, VPD regime of the drawn year), so VPD enters the
#             generator through both the year effects drawn and the regime-specific residual distribution
# Reuses the setup of run_breakeven_years.R. Needs generators built with VPD regimes (build_state_generators.R, v3 cache).
# Env: SIM_CACHE_DIR (v3 cache), SIM_B, SIM_OUT_DIR (default results/breakeven_vpd). Run from simulation_code/.
suppressPackageStartupMessages(library(data.table))
source("config.R")
Sys.setenv(SIM_OUT_DIR = .env("SIM_OUT_DIR", "results/breakeven_vpd"))
src <- readLines("run_breakeven_years.R")
eval(parse(text = src[1:(grep("^main <- list", src)[1] - 1)]))

regs <- fread("results/vpd_regimes_full.csv")
out <- list(); ann <- list()
for (si in seq_along(STATES)) for (zi in seq_along(ZONES)) {
  st <- STATES[si]; z <- ZONES[zi]; seed <- si * 1000L + zi
  cg <- gens$corn[[st]][[z]]; sg <- gens$soy[[st]][[z]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  ry <- regs[state == st]$regime[match(wt$years, regs[state == st]$year)]
  for (pl in c("CCCCCC", "CCSCCS", "CSCSCS")) for (vp in c(FALSE, TRUE)) {
    o <- simulate_breakeven(cg, sg, pl, price_pars[[st]], cost_pars[[st]], rho_c[[st]], rho_s[[st]], B = B, seed = seed, wt = wt,
                            regime_years = ry, vpd_pools = vp)
    loss <- o$rev < o$cost; net <- o$rev - o$cost; rg <- o$regime
    variant <- if (vp) "pools" else "reclass"
    ann[[length(ann) + 1]] <- data.table(state = st, zone = z, plan = pl, variant = variant, p_loss_year = mean(loss),
                                         p_loss_cum = mean(rowSums(net) < 0), mean_net = mean(net))
    for (r in c("low", "mid", "high"))
      out[[length(out) + 1]] <- data.table(state = st, zone = z, plan = pl, variant = variant, regime = r,
        share = mean(rg == r, na.rm = TRUE), p_loss = mean(loss[rg == r], na.rm = TRUE), mean_net = mean(net[rg == r], na.rm = TRUE),
        sd_net = sd(net[rg == r], na.rm = TRUE))
  }
  cat(sprintf("done %s-%s\n", st, z))
}
fwrite(rbindlist(out), file.path(OUT_DIR, "breakeven_by_vpd.csv"))
fwrite(rbindlist(ann), file.path(OUT_DIR, "breakeven_vpd_annual.csv"))
cat("wrote", OUT_DIR, "\n")
