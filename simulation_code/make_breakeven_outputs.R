# Figures and LaTeX tables for the breakeven paper, from results/breakeven/*.csv.
# Figures go to FIGS_DIR (set SIM_FIGS_DIR=../figs to write into the manuscript folder); all names start "be_".
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("config.R")
IN_DIR <- .env("SIM_OUT_DIR", "results/breakeven")
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive = TRUE)

d  <- fread(file.path(IN_DIR, "breakeven_cells.csv"))
be <- fread(file.path(IN_DIR, "breakeven_thresholds.csv"))
pr <- fread("data/state_prices.csv")
d[, zone := factor(zone, c("high", "medium", "low"))]
d[, plan := factor(plan, c("CCCCCC", "CCSCCS", "CSCSCS"))]
plan_lab <- c(CCCCCC = "Continuous corn", CCSCCS = "Corn-corn-soy", CSCSCS = "Strict alternation")
cols <- c(CCCCCC = "#B2182B", CCSCCS = "#EF8A62", CSCSCS = "#2166AC")
# states ordered by continuous-corn annual loss probability (medium zone)
ord <- d[corr == "decomposed" & plan == "CCCCCC" & zone == "medium"][order(p_loss_year), as.character(state)]
d[, state := factor(state, ord)]
th <- theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
sv <- function(p, nm, w = 7.5, h = 4.6) ggsave(file.path(FIGS_DIR, paste0(nm, ".pdf")), p, width = w, height = h)

main <- d[corr == "decomposed"]

# 1. annual breakeven probability by state, rotation, zone
p1 <- ggplot(main, aes(state, p_loss_year, colour = plan, shape = plan)) +
  geom_point(size = 2.2, position = position_dodge(0.5)) + facet_wrap(~zone, nrow = 1) +
  scale_colour_manual(values = cols, labels = plan_lab, name = NULL) +
  scale_shape_discrete(labels = plan_lab, name = NULL) +
  scale_y_continuous(limits = c(0, 0.8), labels = scales::percent) +
  labs(x = NULL, y = "Pr(revenue < cost) in a given year") + th
sv(p1, "be_p_loss_year", 8.5, 4)

# 2. six-year cumulative: probability that cumulative net return is negative
p2 <- ggplot(main, aes(state, p_loss_cum, colour = plan, shape = plan)) +
  geom_point(size = 2.2, position = position_dodge(0.5)) + facet_wrap(~zone, nrow = 1) +
  scale_colour_manual(values = cols, labels = plan_lab, name = NULL) +
  scale_shape_discrete(labels = plan_lab, name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = NULL, y = "Pr(six-year cumulative net return < 0)") + th
sv(p2, "be_p_loss_cum", 8.5, 4)

# 3. rotation gap (alternation minus continuous corn), percentage points
g <- dcast(main, state + zone ~ plan, value.var = "p_loss_year")
g[, `:=`(gap = 100 * (CSCSCS - CCCCCC), gap_cc = 100 * (CCSCCS - CCCCCC))]
p3 <- ggplot(g, aes(state, gap, fill = zone)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(high = "#01665E", medium = "#80CDC1", low = "#BF812D"), name = "NCCPI zone") +
  labs(x = NULL, y = "Change in annual Pr(revenue < cost), pp\n(strict alternation vs continuous corn)") + th
sv(p3, "be_rotation_gap", 7, 4)

# ---- LaTeX tables -----------------------------------------------------------
f0 <- function(x) formatC(x, format = "f", digits = 0); f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3); f1 <- function(x) formatC(x, format = "f", digits = 1)
sto <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
# Each table is written as a complete tabular. (\input directly inside a tabular body breaks under the
# 2025 LaTeX kernel.)
write_tab <- function(rows, spec, head, file) {
  writeLines(c(sprintf("\\begin{tabular}{%s}", spec), "\\toprule", head, "\\midrule", rows,
               "\\bottomrule", "\\end{tabular}"), file)
}

# Inputs table
w <- dcast(merge(be, pr[, .(state, crop, log_sd)], by = c("state", "crop")), state ~ crop,
           value.var = c("price_per_bu", "log_sd", "total_nonland_cost", "breakeven_yield"))
write_tab(sapply(sto, function(s) { r <- w[state == s]
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", s,
          f2(r$price_per_bu_C), f3(r$log_sd_C), f0(r$total_nonland_cost_C), f0(r$breakeven_yield_C),
          f2(r$price_per_bu_S), f3(r$log_sd_S), f0(r$total_nonland_cost_S), f1(r$breakeven_yield_S)) }),
  "lcccccccc",
  c("& \\multicolumn{4}{c}{Corn} & \\multicolumn{4}{c}{Soybeans}\\\\",
    "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
    "State & Price & $\\sigma$ & Cost & BE yield & Price & $\\sigma$ & Cost & BE yield\\\\"),
  file.path(IN_DIR, "tab_inputs.tex"))

# Main results table (medium zone)
m <- dcast(main, state + zone ~ plan, value.var = c("p_loss_year", "p_loss_cum", "mean_net"))
write_tab(sapply(sto, function(s) { r <- m[state == s & zone == "medium"]
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", s,
          f1(100 * r$p_loss_year_CCCCCC), f1(100 * r$p_loss_year_CCSCCS), f1(100 * r$p_loss_year_CSCSCS),
          f1(100 * r$p_loss_cum_CCCCCC), f1(100 * r$p_loss_cum_CSCSCS),
          f1(100 * (r$p_loss_year_CSCSCS - r$p_loss_year_CCCCCC))) }),
  "lcccccc",
  c("& \\multicolumn{3}{c}{Annual} & \\multicolumn{2}{c}{Six-year cumulative} & \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}",
    "State & CCCCCC & CCSCCS & CSCSCS & CCCCCC & CSCSCS & Gap (pp)\\\\"),
  file.path(IN_DIR, "tab_main.tex"))

# Zone table: annual loss probability by zone, continuous corn and alternation
write_tab(sapply(sto, function(s) { r <- m[state == s]; r <- r[match(c("high", "medium", "low"), zone)]
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", s,
          f1(100 * r$p_loss_year_CCCCCC[1]), f1(100 * r$p_loss_year_CCCCCC[2]), f1(100 * r$p_loss_year_CCCCCC[3]),
          f1(100 * r$p_loss_year_CSCSCS[1]), f1(100 * r$p_loss_year_CSCSCS[2]), f1(100 * r$p_loss_year_CSCSCS[3])) }),
  "lcccccc",
  c("& \\multicolumn{3}{c}{Continuous corn} & \\multicolumn{3}{c}{Strict alternation}\\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    "State & High & Medium & Low & High & Medium & Low\\\\"),
  file.path(IN_DIR, "tab_zone.tex"))

# Price-yield dependence check
h <- dcast(d[zone == "medium"], state + plan ~ corr, value.var = "p_loss_year")
h[, diff_pp := 100 * (decomposed - none)]
fwrite(h, file.path(IN_DIR, "hedge_check.csv"))
cat(sprintf("max |dependence effect| on annual Pr(loss): %.2f pp\n", max(abs(h$diff_pp))))
fwrite(g, file.path(IN_DIR, "rotation_gap.csv"))
cat("wrote figures to", FIGS_DIR, "and tables to", IN_DIR, "\n")
