# ============================================================================
# Lower-tail dominance (LTD) check: corn yields after soybeans in a strict
# alternation (SCSCSC) vs. continuous corn (CCCCCC), by region x productivity zone.
# Reads the slim yield generators (config.R); writes figs/ltd_cdf_comparison.{pdf,png},
# results/ltd_dominance_results.csv and results/sample_sizes.csv.
#   setwd("<path>/simulation_code"); source("ltd_dominance_check.R")
# Delinquency threshold (Section 3): ybar = C (1 + i/12)^9 / P at the median
# cost and median price, i = 9%.
# ============================================================================
suppressMessages({ library(data.table); library(ggplot2) })
source("config.R")

yg <- readRDS(YGEN_CACHE); cp <- readRDS(COST_PARAMS_RDS); pp <- readRDS(PRICE_PARAMS_RDS)
regions <- c("northern", "central", "southern"); zones <- c("high", "medium", "low")
i_rate <- 0.09; B <- 50000L
FIGS_DIR <- Sys.getenv("SIM_FIGS_DIR", "figs")
dir.create(FIGS_DIR, showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
set.seed(42)

res <- list(); cdf <- list(); ns <- list()
for (r in regions) for (z in zones) {
  g   <- yg$corn[[r]][[z]]
  ycc <- g$simulate("CCCCCC", B = B)
  ycs <- g$simulate("SCSCSC", B = B)
  ybar <- exp(cp[[r]][[z]]$C$mu) * (1 + i_rate / 12)^9 / exp(pp[[r]]$C$mu)
  Fcc <- ecdf(ycc); Fcs <- ecdf(ycs)

  grid <- sort(unique(c(ycc, ycs))); lo <- grid[grid <= ybar]
  max_diff_below <- if (length(lo)) max(Fcs(lo) - Fcc(lo)) else NA_real_   # > 0 would violate LTD

  gg <- seq(min(ycc, ycs), max(ycc, ycs), length.out = 2000)
  dd <- Fcs(gg) - Fcc(gg)
  keep <- abs(dd) > 1.5 * sqrt(2 * 0.25 / B)                  # ignore Monte Carlo noise
  s <- sign(dd[keep]); ncross <- sum(diff(s) != 0)

  res[[paste(r, z)]] <- data.table(
    region = r, zone = z, ybar = ybar,
    p10_cc = quantile(ycc, .10), p10_cs = quantile(ycs, .10), F_cs_at_p10cc = mean(ycs <= quantile(ycc, .10)),
    p25_cc = quantile(ycc, .25), p25_cs = quantile(ycs, .25), F_cs_at_p25cc = mean(ycs <= quantile(ycc, .25)),
    F_cc_at_ybar = Fcc(ybar), F_cs_at_ybar = Fcs(ybar),
    max_Fcs_minus_Fcc_below_ybar = max_diff_below,
    mean_gap = mean(ycs) - mean(ycc), n_crossings = ncross)
  cdf[[paste(r, z)]] <- rbind(data.table(region = r, zone = z, plan = "CCCCCC", y = gg, F = Fcc(gg)),
                              data.table(region = r, zone = z, plan = "CSCSCS", y = gg, F = Fcs(gg)),
                              data.table(region = r, zone = z, plan = NA_character_, y = ybar, F = NA_real_))
}
tab <- rbindlist(res)
tab[, `:=`(ltd_p10 = F_cs_at_p10cc <= 0.10, ltd_p25 = F_cs_at_p25cc <= 0.25, ltd_ybar = F_cs_at_ybar <= F_cc_at_ybar)]
fwrite(tab, "results/ltd_dominance_results.csv")

# Sample sizes behind the generators: field-year observations per rotation state
# (terminal crop = the crop being modelled, complete six-year history).
for (crop in c("corn", "soy")) for (r in regions) for (z in zones) {
  p <- yg[[crop]][[r]][[z]]$pools
  cc <- if (crop == "corn") "CCCCCC" else "CCCCCS"; alt <- if (crop == "corn") "SCSCSC" else "CSCSCS"
  ns[[paste(crop, r, z)]] <- data.table(crop = crop, region = r, zone = z, n_field_years = sum(p$n), n_states = nrow(p),
                                        n_continuous = p[rot6 == cc, n], n_alternation = p[rot6 == alt, n],
                                        min_state_n = min(p$n), states_below_100 = sum(p$n < 100),
                                        states_below_1000 = sum(p$n < 1000))
}
fwrite(rbindlist(ns), "results/sample_sizes.csv")

cat(sprintf("LTD at P10: %d/9 | at P25: %d/9 | at ybar: %d/9\n", sum(tab$ltd_p10), sum(tab$ltd_p25), sum(tab$ltd_ybar)))
print(tab[, .(region, zone, ybar = round(ybar, 1), F_cs_at_p10cc = round(F_cs_at_p10cc, 3), F_cs_at_p25cc = round(F_cs_at_p25cc, 3),
              F_cc_at_ybar = round(F_cc_at_ybar, 4), F_cs_at_ybar = round(F_cs_at_ybar, 4),
              max_below = round(max_Fcs_minus_Fcc_below_ybar, 4), mean_gap = round(mean_gap, 2), n_crossings)])

lines <- rbindlist(cdf)[!is.na(plan)]; vl <- rbindlist(cdf)[is.na(plan), .(region, zone, ybar = y)]
for (d in list(lines, vl)) { d[, region := factor(region, regions)]; d[, zone := factor(zone, zones)] }
p <- ggplot(lines, aes(y, F, color = plan, linetype = plan)) + geom_line(linewidth = 0.7) +
  geom_vline(data = vl, aes(xintercept = ybar), linetype = "dotted", color = "grey40") +
  facet_grid(zone ~ region, scales = "free_x", labeller = label_both) +
  scale_color_manual(values = c(CCCCCC = "#D55E00", CSCSCS = "#0072B2"), labels = c("Continuous corn", "Corn-soy rotation")) +
  scale_linetype_manual(values = c(CCCCCC = "solid", CSCSCS = "dashed"), labels = c("Continuous corn", "Corn-soy rotation")) +
  labs(x = "Corn yield (bu/ac)", y = "F(y)", color = NULL, linetype = NULL,
       title = "Empirical CDFs of corn yield: rotation vs. continuous corn",
       subtitle = "Dotted line: delinquency threshold at median cost and price (i = 9%, 9-month accrual).") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(FIGS_DIR, "ltd_cdf_comparison.pdf"), p, width = 10, height = 9)
ggsave(file.path(FIGS_DIR, "ltd_cdf_comparison.png"), p, width = 10, height = 9, dpi = 200)
