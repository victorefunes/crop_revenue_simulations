# ============================================================================
# Year-shock breakeven simulator.
#
#   Y_t = mu[window_t] + (yfe[j_t] - mean(yfe)) + e_t     j_t ~ U{calendar years}, e_t from the residual pool
#   P_t = Pbar-lognormal, log-shock z_p = rho * z_a(j_t) + sqrt(1 - rho^2) * eps,
#         z_a = Phi^-1 of the rank of crop-specific year effect yfe[j_t]  (price moves against the year shock)
#   K_t = Kbar * lognormal(unit mean, sd 0.10) + same-crop penalty
#   loss_t = 1{P_t * Y_t < K_t}
# Year effects (yfe) are the fitted year fixed effects of the cell's two-way FE model; they are the
# common weather/market shock that the residual pool does not contain. A calendar year's "type" (bad /
# normal / good) is its tercile of a weather index: the sum of the standardized corn and soybean year effects.
# All random numbers are drawn every call (even when a factor is switched off) so coalitions that differ in
# which factors are random are compared on common random numbers.
# ============================================================================

weather_types <- function(corn_fe, soy_fe) {
  yrs <- intersect(corn_fe$year, soy_fe$year)
  ce <- corn_fe[match(yrs, year), fe]; se <- soy_fe[match(yrs, year), fe]
  w  <- as.numeric(scale(ce)) + as.numeric(scale(se))
  r  <- rank(w, ties.method = "first"); J <- length(yrs); k <- round(J / 3)
  list(years = yrs, ce = ce, se = se, w = w,
       type = ifelse(r <= k, "bad", ifelse(r > J - k, "good", "normal")),
       rank_c = rank(ce, ties.method = "first"), rank_s = rank(se, ties.method = "first"))
}

simulate_breakeven <- function(corn_gen, soy_gen, plan, price_par, cost_par, rho_c, rho_s,
                               B = 50000L, seed = 1L, same_crop_penalty = c(C = 15, S = 5),
                               fix_price = FALSE, fix_yield = FALSE, fix_cost = FALSE, wt = NULL,
                               regime_years = NULL, vpd_pools = FALSE) {
  # regime_years: character vector aligned with wt$years giving each calendar year's VPD regime (optional).
  # vpd_pools: draw the field residual from the pool for (window, VPD regime of the drawn year) instead of the
  # window's pooled residuals; needs generators built with regime_tab (gen$pools_vpd).
  set.seed(seed)
  if (is.null(wt)) wt <- weather_types(corn_gen$year_fe, soy_gen$year_fe)
  J  <- length(wt$years)
  ch <- strsplit(plan, "")[[1]]
  path <- make_transition_path(plan, plan)               # steady state: history equals the plan
  prev <- c(ch[6], ch[1:5]); pen <- ifelse(ch == prev, same_crop_penalty[ch], 0)
  rev <- K <- Jm <- matrix(NA_real_, B, 6)
  regm <- if (!is.null(regime_years)) matrix(NA_character_, B, 6) else NULL
  for (t in 1:6) {
    c_t <- ch[t]; gen <- if (c_t == "C") corn_gen else soy_gen
    pp <- price_par[[c_t]]; cp <- cost_par[[c_t]]; rho <- if (c_t == "C") rho_c else rho_s
    efe <- if (c_t == "C") wt$ce else wt$se; rk <- if (c_t == "C") wt$rank_c else wt$rank_s
    j  <- sample.int(J, B, replace = TRUE)
    pool <- gen$pools[rot6 == path$window6[t], resid_pool][[1]]
    mu_s <- gen$mu[rot6 == path$window6[t], mu]
    reg_j <- if (!is.null(regime_years)) regime_years[j] else NULL
    if (vpd_pools && !is.null(gen$pools_vpd) && !is.null(reg_j)) {
      e <- numeric(B)
      for (r in unique(reg_j[!is.na(reg_j)])) {
        ii <- which(reg_j == r)
        pr <- gen$pools_vpd[rot6 == path$window6[t] & regime == r, resid_pool]
        e[ii] <- sample(if (length(pr)) pr[[1]] else pool, length(ii), replace = TRUE)
      }
      e[is.na(reg_j)] <- sample(pool, sum(is.na(reg_j)), replace = TRUE)
    } else e <- sample(pool, B, replace = TRUE)
    eps <- rnorm(B); ec <- rnorm(B)
    Y <- if (fix_yield) rep(mu_s, B) else pmax(mu_s + (efe[j] - mean(efe)) + e, 0)
    za <- qnorm((rk[j] - 0.5) / J)
    zp <- rho * za + sqrt(1 - rho^2) * eps
    P  <- if (fix_price) rep(pp$Pbar, B) else exp(pp$mu + pp$sigma * zp)
    # cp$idx (optional): named numeric by calendar year, multiplies the cost level of the drawn year (cost timing)
    ki <- if (is.null(cp$idx)) 1 else cp$idx[as.character(wt$years[j])]
    Kt <- if (fix_cost) rep(cp$Cbar, B) else cp$Cbar * ki * exp(cp$sigma * ec - 0.5 * cp$sigma^2)
    rev[, t] <- P * Y; K[, t] <- Kt + pen[t]; Jm[, t] <- j
    if (!is.null(regm)) regm[, t] <- reg_j
  }
  list(rev = rev, cost = K, j = Jm, wt = wt, plan = plan, regime = regm)
}

# Shapley values of three uncertainty sources (price, yield, cost) for a scalar functional v(random = <names>)
shapley3 <- function(v_of) {
  players <- c("price", "yield", "cost")
  subsets <- expand.grid(price = c(FALSE, TRUE), yield = c(FALSE, TRUE), cost = c(FALSE, TRUE))
  key <- function(r) paste(as.integer(r), collapse = "")
  val <- list()
  for (i in seq_len(nrow(subsets))) { s <- as.logical(subsets[i, ])
    val[[key(s)]] <- v_of(random = players[s]) }
  perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  phi <- setNames(numeric(3), players)
  for (p in perms) { cur <- c(FALSE, FALSE, FALSE)
    for (i in p) { before <- val[[key(cur)]]; cur[i] <- TRUE; phi[i] <- phi[i] + (val[[key(cur)]] - before) / 6 } }
  list(phi = phi, v_empty = val[[key(c(FALSE, FALSE, FALSE))]], v_all = val[[key(c(TRUE, TRUE, TRUE))]])
}
