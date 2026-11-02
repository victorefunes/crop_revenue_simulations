# Breakeven probability by type of year with and without the price-yield dependence (all zones, all plans).
# Reuses the setup of run_breakeven_years.R (same seeds, so the "with" columns reproduce breakeven_year_types.csv).
# Env: SIM_CACHE_DIR (v2 cache with year effects), SIM_B, SIM_OUT_DIR (default results/breakeven_years). Run from simulation_code/.
suppressPackageStartupMessages(library(data.table))
src <- readLines("run_breakeven_years.R")
eval(parse(text = src[1:(grep("^main <- list", src)[1] - 1)]))
out <- list()
for (si in seq_along(STATES)) for (zi in seq_along(ZONES)) {
  st <- STATES[si]; z <- ZONES[zi]; seed <- si * 1000L + zi
  cg <- gens$corn[[st]][[z]]; sg <- gens$soy[[st]][[z]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  for (pl in PLANS) for (dp in c(TRUE, FALSE)) {
    o <- simulate_breakeven(cg, sg, pl, price_pars[[st]], cost_pars[[st]], if (dp) rho_c[[st]] else 0, if (dp) rho_s[[st]] else 0,
                            B = B, seed = seed, wt = wt)
    loss <- o$rev < o$cost; ty <- matrix(wt$type[o$j], nrow(o$j))
    out[[length(out) + 1]] <- data.table(state = st, zone = z, plan = pl, dependence = dp, p_loss_year = mean(loss),
      p_bad = mean(loss[ty == "bad"]), p_normal = mean(loss[ty == "normal"]), p_good = mean(loss[ty == "good"]))
  }
  cat(sprintf("done %s-%s\n", st, z))
}
fwrite(rbindlist(out), file.path(OUT_DIR, "dependence_year_types.csv"))
