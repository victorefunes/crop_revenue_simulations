# Shapley decomposition of the variance of net return, computed separately within bad / normal / good years.
# V(S) = variance of net return among simulated path-years of the given type when only the sources in S are random
# (price, yield, cost); the other sources are fixed at their means. The calendar year j is drawn in every coalition,
# so the type of each path-year is the same across coalitions and the coalitions are compared on common random numbers.
# Reuses the setup of run_breakeven_years.R. Run from simulation_code/.
suppressPackageStartupMessages(library(data.table))
source("config.R")
src <- readLines("run_breakeven_years.R")
eval(parse(text = src[1:(grep("^main <- list", src)[1] - 1)]))

players <- c("price", "yield", "cost")
subsets <- expand.grid(price = c(FALSE, TRUE), yield = c(FALSE, TRUE), cost = c(FALSE, TRUE))
key <- function(r) paste(as.integer(r), collapse = "")
perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))

out <- list()
for (si in seq_along(STATES)) for (zi in seq_along(ZONES)) for (pl in c("CCCCCC", "CSCSCS")) {
  st <- STATES[si]; z <- ZONES[zi]; seed <- si * 1000L + zi
  cg <- gens$corn[[st]][[z]]; sg <- gens$soy[[st]][[z]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  val <- list()                                   # coalition -> named vector of conditional variances
  for (i in seq_len(nrow(subsets))) {
    s <- as.logical(subsets[i, ])
    o <- simulate_breakeven(cg, sg, pl, price_pars[[st]], cost_pars[[st]], rho_c[[st]], rho_s[[st]], B = B, seed = seed, wt = wt,
                            fix_price = !s[1], fix_yield = !s[2], fix_cost = !s[3])
    net <- as.vector(o$rev - o$cost); ty <- wt$type[as.vector(o$j)]
    val[[key(s)]] <- sapply(c(bad = "bad", normal = "normal", good = "good"), function(tp) var(net[ty == tp]))
  }
  for (tp in c("bad", "normal", "good")) {
    phi <- setNames(numeric(3), players)
    for (p in perms) { cur <- c(FALSE, FALSE, FALSE)
      for (i in p) { before <- val[[key(cur)]][[tp]]; cur[i] <- TRUE; phi[i] <- phi[i] + (val[[key(cur)]][[tp]] - before) / 6 } }
    out[[length(out) + 1]] <- data.table(state = st, zone = z, plan = pl, type = tp, var_all = val[[key(c(TRUE, TRUE, TRUE))]][[tp]],
      phi_price = phi[["price"]], phi_yield = phi[["yield"]], phi_cost = phi[["cost"]])
  }
}
res <- rbindlist(out)
res[, tot := phi_price + phi_yield + phi_cost]
res[, `:=`(sh_price = phi_price / tot, sh_yield = phi_yield / tot, sh_cost = phi_cost / tot)]
fwrite(res, file.path(OUT_DIR, "shapley_by_year_type.csv"))
cat("wrote", file.path(OUT_DIR, "shapley_by_year_type.csv"), "\n")
