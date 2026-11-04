# Appendix tables for the breakeven paper: robustness to inputs (medium zone) and the full state-zone grid.
# Reads results/breakeven_years/robustness_inputs.csv and breakeven_year_types.csv. Run from simulation_code/
# (SIM_TABLES_DIR=../tables for the manuscript). Writes rob_inputs.tex and app_grid.tex.
suppressPackageStartupMessages(library(data.table))
TAB <- Sys.getenv("SIM_TABLES_DIR", "../tables"); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
STATES <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
f1 <- function(x) sprintf("%.1f", 100 * x)

r <- fread("results/breakeven_years/robustness_inputs.csv")
scn <- c("baseline", "cost_m10", "cost_p10", "cost_p20", "vol_x075", "vol_x125", "no_penalty")
hdr <- c("Baseline", "Cost $-10\\%$", "Cost $+10\\%$", "Cost $+20\\%$", "Vol. $\\times 0.75$", "Vol. $\\times 1.25$", "No penalty")
w <- dcast(r, state + scenario ~ plan, value.var = "p_year")
w[, gap := CSCSCS - CCCCCC]
row <- function(st, col) paste(c(st, sapply(scn, function(s) f1(w[state == st & scenario == s][[col]]))), collapse = " & ")
lines <- c("\\begin{tabular}{lccccccc}", "\\toprule", paste0("State & ", paste(hdr, collapse = " & "), "\\\\"), "\\midrule",
  "\\multicolumn{8}{l}{\\emph{Panel A. Continuous corn, annual probability (percent)}}\\\\",
  paste0(sapply(STATES, row, col = "CCCCCC"), " \\\\"), "\\midrule",
  "\\multicolumn{8}{l}{\\emph{Panel B. Strict alternation, annual probability (percent)}}\\\\",
  paste0(sapply(STATES, row, col = "CSCSCS"), " \\\\"), "\\midrule",
  "\\multicolumn{8}{l}{\\emph{Panel C. Difference, strict alternation minus continuous corn (points)}}\\\\",
  paste0(sapply(STATES, row, col = "gap"), " \\\\"), "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(TAB, "rob_inputs.tex"))

# Year-specific costs (cost_timing.R): continuous corn, baseline vs ERS cost timing, medium zone
ct <- fread("results/breakeven_years/cost_timing.csv")
cw <- dcast(ct[scenario %in% c("baseline", "ers_timing")], state + scenario ~ plan, value.var = c("p_year", "p_bad", "p_good"))
cw[, gap := p_year_CSCSCS - p_year_CCCCCC]
cget <- function(st, s, col) cw[state == st & scenario == s][[col]]
crow <- function(st) paste(c(st, f1(cget(st, "baseline", "p_year_CCCCCC")), f1(cget(st, "ers_timing", "p_year_CCCCCC")),
  f1(cget(st, "baseline", "gap")), f1(cget(st, "ers_timing", "gap")),
  f1(cget(st, "baseline", "p_bad_CCCCCC")), f1(cget(st, "ers_timing", "p_bad_CCCCCC")),
  f1(cget(st, "baseline", "p_good_CCCCCC")), f1(cget(st, "ers_timing", "p_good_CCCCCC"))), collapse = " & ")
cl <- c("\\begin{tabular}{lcccccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{Annual} & \\multicolumn{2}{c}{Gap (pp)} & \\multicolumn{2}{c}{Bad years} & \\multicolumn{2}{c}{Good years}\\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "State & Fixed & ERS & Fixed & ERS & Fixed & ERS & Fixed & ERS\\\\", "\\midrule",
  paste0(sapply(STATES, crow), " \\\\"), "\\bottomrule", "\\end{tabular}")
writeLines(cl, file.path(TAB, "rob_cost_timing.tex"))

m <- fread("results/breakeven_years/breakeven_year_types.csv")
Z <-c(high = "High", medium = "Medium", low = "Low")
g <- function(st, z) {
  a <- m[state == st & zone == z & plan == "CCCCCC"]; b <- m[state == st & zone == z & plan == "CSCSCS"]
  paste(c(st, Z[[z]], f1(a$p_loss_year), f1(a$p_bad), f1(a$p_good), f1(a$p_loss_cum),
          f1(b$p_loss_year), f1(b$p_bad), f1(b$p_good), f1(b$p_loss_cum)), collapse = " & ")
}
gl <- c("\\begin{tabular}{llcccccccc}", "\\toprule",
  " & & \\multicolumn{4}{c}{Continuous corn} & \\multicolumn{4}{c}{Strict alternation}\\\\",
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
  "State & Zone & Annual & Bad & Good & Six-year & Annual & Bad & Good & Six-year\\\\", "\\midrule")
for (st in STATES) { gl <- c(gl, paste0(sapply(names(Z), function(z) g(st, z)), " \\\\")); if (st != "WI") gl <- c(gl, "\\addlinespace") }
writeLines(c(gl, "\\bottomrule", "\\end{tabular}"), file.path(TAB, "app_grid.tex"))
cat("ok\n")
