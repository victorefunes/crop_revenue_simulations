# Empirical checks of the conditions in the theoretical framework section.
#  (1) Proposition 1: dominance of corn yield after soybeans over continuous corn at the breakeven yield K/P and the mean-yield
#      condition, from the estimated generators (year shock + residual).
#  (2) Proposition 2: per-year breakeven probability under rotation <= under continuous corn, in every year of the path.
#  (4) Proposition 4: yield elasticity of the year shock (sd of year effects / mean yield) against the price hedge |rho| * sigma_P.
# Reuses the setup of run_breakeven_years.R (v2 generators with year effects). Run from simulation_code/.
suppressPackageStartupMessages(library(data.table))
source("config.R")
Sys.setenv(SIM_OUT_DIR = .env("SIM_OUT_DIR", "results/theory"))
src <- readLines("run_breakeven_years.R")
eval(parse(text = src[1:(grep("^main <- list", src)[1] - 1)]))
TAB_DIR <- .env("SIM_TABLES_DIR", "../tables")
f1 <- function(x) formatC(x, format = "f", digits = 1); f3 <- function(x) formatC(x, format = "f", digits = 3)
f0 <- function(x) formatC(x, format = "f", digits = 0)
write_tab <- function(rows, spec, head, file) writeLines(c(sprintf("\\begin{tabular}{%s}", spec), "\\toprule", head,
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file)
Bq <- 400000L; set.seed(123)

draw_yield <- function(gen, w, wt, efe, B) {
  j <- sample.int(length(wt$years), B, replace = TRUE)
  mu <- gen$mu[rot6 == w, mu]; pool <- gen$pools[rot6 == w, resid_pool][[1]]
  pmax(mu + (efe[j] - mean(efe)) + sample(pool, B, replace = TRUE), 0)
}

res1 <- list(); res4 <- list()
for (st in STATES) for (z in ZONES) {
  cg <- gens$corn[[st]][[z]]; sg <- gens$soy[[st]][[z]]; wt <- weather_types(cg$year_fe, sg$year_fe)
  Kc <- cost_pars[[st]]$C$Cbar; Pc <- price_pars[[st]]$C$Pbar
  ycc <- draw_yield(cg, "CCCCCC", wt, wt$ce, Bq); ycs <- draw_yield(cg, "SCSCSC", wt, wt$ce, Bq)
  ybar <- Kc / Pc; ybar_cc <- (Kc + 15) / Pc
  res1[[length(res1) + 1]] <- data.table(state = st, zone = z, ybar = ybar, ybar_cc = ybar_cc,
    F_cc = mean(ycc < ybar), F_cs = mean(ycs < ybar), F_cc_pen = mean(ycc < ybar_cc),
    mean_cc = mean(ycc), mean_cs = mean(ycs))
  for (cr in c("C", "S")) {
    gen <- if (cr == "C") cg else sg; efe <- if (cr == "C") wt$ce else wt$se
    w <- if (cr == "C") "SCSCSC" else "CSCSCS"
    mu <- gen$mu[rot6 == w, mu]; sdA <- sd(efe)
    rho <- if (cr == "C") rho_c[[st]] else rho_s[[st]]; sP <- price_pars[[st]][[cr]]$sigma
    res4[[length(res4) + 1]] <- data.table(state = st, zone = z, crop = cr, mu = mu, sd_year = sdA, elasticity = sdA / mu,
      rho = rho, sigma_P = sP, hedge = -rho * sP)
  }
}
r1 <- rbindlist(res1); r4 <- rbindlist(res4)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(r1, file.path(OUT_DIR, "prop1_check.csv")); fwrite(r4, file.path(OUT_DIR, "prop4_check.csv"))

# (2) per-year dominance of the breakeven probability (rotation vs continuous corn), from the by-position output
pos <- fread("results/breakeven_years/breakeven_by_position.csv")
cc <- pos[plan == "CCCCCC", .(state, zone, t, p_cc = p)]
cmp <- merge(pos[plan != "CCCCCC"], cc, by = c("state", "zone", "t"))
cmp[, excess := p - p_cc]
s2 <- cmp[, .(cells = .N, n_violate = sum(excess > 0), max_excess_pp = 100 * max(excess)), by = plan]
fwrite(s2, file.path(OUT_DIR, "prop2_check.csv"))
cat("Proposition 2 check (rotation minus continuous corn, by year; excess > 0 means rotation worse):\n"); print(s2)
cat("largest violations:\n"); print(cmp[excess > 0][order(-excess)][1:6, .(state, zone, plan, t, crop, excess_pp = round(100 * excess, 2))])

# tables for the medium zone
m1 <- r1[zone == "medium"]
write_tab(sapply(STATES, function(s) { r <- m1[state == s]
  sprintf("%s & %s & %s & %s & %s & %s \\\\", s, f0(r$ybar), f1(100 * r$F_cc), f1(100 * r$F_cs), f1(r$mean_cs - r$mean_cc), f1(100 * r$F_cc_pen)) }),
  "lccccc", c("State & Breakeven yield & $F_{CC}(\\bar y)$ & $F_{CS}(\\bar y)$ & $E[Y_{CS}]-E[Y_{CC}]$ & $F_{CC}$ with penalty\\\\"),
  file.path(TAB_DIR, "th_ltd.tex"))
m4 <- r4[zone == "medium"]
yt <- fread("results/breakeven_years/breakeven_year_types.csv")[zone == "medium" & plan == "CCCCCC"]
write_tab(sapply(STATES, function(s) { a <- m4[state == s & crop == "C"]; b <- m4[state == s & crop == "S"]; y <- yt[state == s]
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", s, f3(a$elasticity), f3(a$hedge), f3(a$elasticity - a$hedge),
          f3(b$elasticity), f3(b$hedge), f1(100 * (y$p_bad - y$p_good))) }),
  "lcccccc", c("& \\multicolumn{3}{c}{Corn} & \\multicolumn{2}{c}{Soybeans} & \\\\", "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}",
  "State & $\\sigma_A/\\mu$ & $|\\rho|\\sigma_P$ & Net & $\\sigma_A/\\mu$ & $|\\rho|\\sigma_P$ & Bad minus good (pp)\\\\"),
  file.path(TAB_DIR, "th_hedge.tex"))
# Proposition 3: mean margin under each plan and the reduction in the annual breakeven probability
yy <- fread("results/breakeven_years/breakeven_year_types.csv")[zone == "medium" & plan %in% c("CCCCCC", "CSCSCS")]
write_tab(sapply(STATES, function(s) { a <- yy[state == s & plan == "CCCCCC"]; b <- yy[state == s & plan == "CSCSCS"]
  sprintf("%s & %s & %s & %s & %s \\\\", s, f0(a$mean_net), f0(b$mean_net), f0(b$mean_net - a$mean_net), f1(100 * (a$p_loss_year - b$p_loss_year))) }),
  "lcccc", c("State & Continuous corn & Alternation & Gain & Fall in Pr (pp)\\\\"), file.path(TAB_DIR, "th_margin.tex"))
cat("\nProposition 3 (medium zone: mean margin, \\$/ac):\n"); print(yy[, .(state, plan, mean_net = round(mean_net), p = round(100 * p_loss_year, 1))][order(plan, state)])
cat("\nProposition 1 check (medium zone):\n"); print(m1[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\nProposition 4 check (medium zone):\n"); print(m4[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
