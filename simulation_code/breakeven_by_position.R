# Breakeven probability by position in the six-year path (year 1..6), for each state, zone and plan.
# Under strict alternation this shows the corn/soybean zig-zag; under continuous corn it is flat.
# Reuses the setup of run_breakeven_years.R (everything before its main loop). Run from simulation_code/.
suppressPackageStartupMessages(library(data.table))
source("config.R")
src <- readLines("run_breakeven_years.R")
eval(parse(text = src[1:(grep("^main <- list", src)[1] - 1)]))

out <- list()
for (si in seq_along(STATES)) for (zi in seq_along(ZONES)) {
  st <- STATES[si]; z <- ZONES[zi]; seed <- si * 1000L + zi
  cg <- gens$corn[[st]][[z]]; sg <- gens$soy[[st]][[z]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  for (pl in PLANS) {
    o <- simulate_breakeven(cg, sg, pl, price_pars[[st]], cost_pars[[st]], rho_c[[st]], rho_s[[st]],
                            B = B, seed = seed, wt = wt)
    out[[length(out) + 1]] <- data.table(state = st, zone = z, plan = pl, t = 1:6,
                                         crop = strsplit(pl, "")[[1]], p = colMeans(o$rev < o$cost))
  }
}
fwrite(rbindlist(out), file.path(OUT_DIR, "breakeven_by_position.csv"))
cat("wrote", file.path(OUT_DIR, "breakeven_by_position.csv"), "\n")
