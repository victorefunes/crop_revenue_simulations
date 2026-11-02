# Figures and LaTeX tables for the year-shock breakeven paper, from results/breakeven_years/*.csv.
# Figures -> FIGS_DIR (SIM_FIGS_DIR=../figs for the manuscript), tables -> SIM_TABLES_DIR (default ../tables).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("config.R")
IN_DIR <- .env("SIM_OUT_DIR", "results/breakeven_years")
TAB_DIR <- .env("SIM_TABLES_DIR", "../tables")
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive = TRUE)
if (!dir.exists(TAB_DIR)) dir.create(TAB_DIR, recursive = TRUE)

m   <- fread(file.path(IN_DIR, "breakeven_year_types.csv"))
cy  <- fread(file.path(IN_DIR, "breakeven_by_calendar_year.csv"))
sh  <- fread(file.path(IN_DIR, "shapley.csv"))
dep <- fread(file.path(IN_DIR, "dependence_check.csv"))
be  <- fread("results/breakeven/breakeven_thresholds.csv"); pr <- fread("data/state_prices.csv")
plan_lab <- c(CCCCCC = "Continuous corn", CCSCCS = "Corn-corn-soy", CSCSCS = "Strict alternation")
cols <- c(CCCCCC = "#B2182B", CCSCCS = "#EF8A62", CSCSCS = "#2166AC")
sto <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
th <- theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
sv <- function(p, nm, w, h) ggsave(file.path(FIGS_DIR, paste0(nm, ".pdf")), p, width = w, height = h)
f1 <- function(x) formatC(x, format = "f", digits = 1); f0 <- function(x) formatC(x, format = "f", digits = 0)
f2 <- function(x) formatC(x, format = "f", digits = 2); f3 <- function(x) formatC(x, format = "f", digits = 3)
write_tab <- function(rows, spec, head, file) writeLines(c(sprintf("\\begin{tabular}{%s}", spec), "\\toprule", head,
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file)

med <- m[zone == "medium"]
long <- melt(med, id.vars = c("state", "plan"), measure.vars = c("p_bad", "p_normal", "p_good"),
             variable.name = "type", value.name = "p")
long[, type := factor(sub("p_", "", type), c("bad", "normal", "good"))]
long[, state := factor(state, sto)]

# 1. PRIMARY: breakeven probability in bad, normal and good years
p1 <- ggplot(long, aes(type, p, colour = plan, group = plan, shape = plan)) +
  geom_line(linewidth = 0.5) + geom_point(size = 2) + facet_wrap(~state, nrow = 2) +
  scale_colour_manual(values = cols, labels = plan_lab, name = NULL) + scale_shape_discrete(labels = plan_lab, name = NULL) +
  scale_y_continuous(limits = c(0, 0.9), labels = scales::percent) +
  labs(x = "Type of year (tercile of the weather index)", y = "Pr(revenue < cost) in a year of that type") + th
sv(p1, "be_year_types", 8.5, 5)

# 2. rotation gap by year type
g <- dcast(long, state + type ~ plan, value.var = "p")[, gap := 100 * (CSCSCS - CCCCCC)]
p2 <- ggplot(g, aes(state, gap, fill = type)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(bad = "#8C510A", normal = "#C7C7C7", good = "#01665E"), name = "Year type") +
  labs(x = NULL, y = "Change in Pr(revenue < cost), pp\n(strict alternation vs continuous corn)") + th
sv(p2, "be_gap_by_year_type", 7, 4)

# 3. calendar years (medium zone)
cym <- cy[zone == "medium" & plan %in% c("CCCCCC", "CSCSCS")]
cym[, state := factor(state, sto)]
p3 <- ggplot(cym, aes(year, p_loss, colour = plan, group = plan)) + geom_line(linewidth = 0.5) +
  geom_point(aes(shape = type), size = 1.8) + facet_wrap(~state, nrow = 2) +
  scale_colour_manual(values = cols, labels = plan_lab, name = NULL) +
  scale_shape_manual(values = c(bad = 25, normal = 1, good = 17), name = "Year type") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(breaks = seq(2008, 2022, 4)) +
  labs(x = NULL, y = "Pr(revenue < cost) in that calendar year's conditions") + th
sv(p3, "be_calendar_years", 9, 5)

# 4. Shapley: shares of net-return variance
# Shares are taken over the Shapley total phi_price + phi_yield + phi_cost = Var(all) - Var(none). With a rotation
# the deterministic margin differs between corn and soybean years, so Var(none) > 0 and is not risk.
sh[, phi_v_tot := phi_v_price + phi_v_yield + phi_v_cost]
sm <- sh[zone == "medium" & plan %in% c("CCCCCC", "CSCSCS")]
sl <- melt(sm[, .(state, plan, price = phi_v_price / phi_v_tot, yield = phi_v_yield / phi_v_tot, cost = phi_v_cost / phi_v_tot)],
           id.vars = c("state", "plan"), variable.name = "source", value.name = "share")
sl[, state := factor(state, sto)]
p4 <- ggplot(sl, aes(state, 100 * share, fill = source)) + geom_col(width = 0.7) + facet_wrap(~plan, labeller = as_labeller(plan_lab)) +
  scale_fill_manual(values = c(price = "#4575B4", yield = "#D73027", cost = "#BDBDBD"), name = "Source of risk") +
  labs(x = NULL, y = "Shapley share of the variance of net return (%)") + th
sv(p4, "be_shapley", 8, 3.8)

# 5. breakeven probability by position in the six-year path (run breakeven_by_position.R first)
pos_file <- file.path(IN_DIR, "breakeven_by_position.csv")
if (file.exists(pos_file)) {
  pos <- fread(pos_file)[zone == "medium"]; pos[, state := factor(state, sto)]
  p5 <- ggplot(pos, aes(t, 100 * p, colour = plan, group = plan)) + geom_line(linewidth = 0.5) +
    geom_point(aes(shape = crop), size = 1.9) + facet_wrap(~state, nrow = 2) +
    scale_colour_manual(values = cols, labels = plan_lab, name = NULL) +
    scale_shape_manual(values = c(C = 16, S = 17), labels = c(C = "Corn year", S = "Soybean year"), name = NULL) +
    scale_x_continuous(breaks = 1:6) + scale_y_continuous(limits = c(0, 80)) +
    labs(x = "Year in the six-year path", y = "Pr(revenue < cost) in that year (%)") + th
  sv(p5, "be_by_position", 9, 5)
}

# ---- tables ------------------------------------------------------------------
w <- dcast(merge(be, pr[, .(state, crop, log_sd)], by = c("state", "crop")), state ~ crop,
           value.var = c("price_per_bu", "log_sd", "total_nonland_cost", "breakeven_yield"))
write_tab(sapply(sto, function(s) { r <- w[state == s]
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", s, f2(r$price_per_bu_C), f3(r$log_sd_C), f0(r$total_nonland_cost_C),
          f0(r$breakeven_yield_C), f2(r$price_per_bu_S), f3(r$log_sd_S), f0(r$total_nonland_cost_S), f1(r$breakeven_yield_S)) }),
  "lcccccccc", c("& \\multicolumn{4}{c}{Corn} & \\multicolumn{4}{c}{Soybeans}\\\\", "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  "State & Price & $\\sigma$ & Cost & BE yield & Price & $\\sigma$ & Cost & BE yield\\\\"), file.path(TAB_DIR, "by_inputs.tex"))

# primary table: bad / normal / good years, continuous corn and strict alternation, with gap in bad and good years
write_tab(sapply(sto, function(s) { cc <- med[state == s & plan == "CCCCCC"]; cs <- med[state == s & plan == "CSCSCS"]
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", s, f1(100 * cc$p_bad), f1(100 * cc$p_normal), f1(100 * cc$p_good),
          f1(100 * cs$p_bad), f1(100 * cs$p_normal), f1(100 * cs$p_good),
          f1(100 * (cs$p_bad - cc$p_bad)), f1(100 * (cs$p_good - cc$p_good))) }),
  "lcccccccc", c("& \\multicolumn{3}{c}{Continuous corn} & \\multicolumn{3}{c}{Strict alternation} & \\multicolumn{2}{c}{Gap (pp)}\\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-9}", "State & Bad & Normal & Good & Bad & Normal & Good & Bad & Good\\\\"),
  file.path(TAB_DIR, "by_yeartypes.tex"))

# annual and six-year cumulative
write_tab(sapply(sto, function(s) { cc <- med[state == s & plan == "CCCCCC"]; cm <- med[state == s & plan == "CCSCCS"]; cs <- med[state == s & plan == "CSCSCS"]
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", s, f1(100 * cc$p_loss_year), f1(100 * cm$p_loss_year), f1(100 * cs$p_loss_year),
          f1(100 * cc$p_loss_cum), f1(100 * cs$p_loss_cum), f1(100 * (cs$p_loss_year - cc$p_loss_year))) }),
  "lcccccc", c("& \\multicolumn{3}{c}{Annual} & \\multicolumn{2}{c}{Six-year cumulative} & \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}", "State & CCCCCC & CCSCCS & CSCSCS & CCCCCC & CSCSCS & Gap (pp)\\\\"),
  file.path(TAB_DIR, "by_main.tex"))

# zones: annual probability by zone (continuous corn and alternation)
write_tab(sapply(sto, function(s) { r <- function(pl) sapply(c("high", "medium", "low"), function(z) f1(100 * m[state == s & zone == z & plan == pl]$p_loss_year))
  sprintf("%s & %s & %s \\\\", s, paste(r("CCCCCC"), collapse = " & "), paste(r("CSCSCS"), collapse = " & ")) }),
  "lcccccc", c("& \\multicolumn{3}{c}{Continuous corn} & \\multicolumn{3}{c}{Strict alternation}\\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}", "State & High & Medium & Low & High & Medium & Low\\\\"), file.path(TAB_DIR, "by_zone.tex"))

# Shapley table: variance shares (continuous corn and alternation) plus signed probability contributions for alternation
write_tab(sapply(sto, function(s) { a <- sh[state == s & zone == "medium" & plan == "CCCCCC"]; b <- sh[state == s & zone == "medium" & plan == "CSCSCS"]
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", s,
          f0(100 * a$phi_v_price / a$phi_v_tot), f0(100 * a$phi_v_yield / a$phi_v_tot), f0(100 * a$phi_v_cost / a$phi_v_tot),
          f0(100 * b$phi_v_price / b$phi_v_tot), f0(100 * b$phi_v_yield / b$phi_v_tot), f0(100 * b$phi_v_cost / b$phi_v_tot),
          f1(100 * b$phi_p_price), f1(100 * b$phi_p_yield), f1(100 * b$phi_p_cost)) }),
  "lccccccccc", c("& \\multicolumn{3}{c}{Variance share, cont. corn (\\%)} & \\multicolumn{3}{c}{Variance share, alternation (\\%)} & \\multicolumn{3}{c}{Contribution to Pr (pp), alternation}\\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}", "State & Price & Yield & Cost & Price & Yield & Cost & Price & Yield & Cost\\\\"),
  file.path(TAB_DIR, "by_shapley.tex"))

d <- merge(m, dep, by = c("state", "zone", "plan")); d[, diff_pp := 100 * (p_loss_year - p_loss_year_nodep)]
fwrite(d[, .(state, zone, plan, p_loss_year, p_loss_year_nodep, diff_pp)], file.path(IN_DIR, "dependence_summary.csv"))
cat(sprintf("max |dependence effect| %.1f pp\n", max(abs(d$diff_pp))))
cat("figures ->", FIGS_DIR, "; tables ->", TAB_DIR, "\n")
