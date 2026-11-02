# Field-level NCCPI vs mean yield, by state, with the tercile cut-points that define the high / medium / low zones.
# One point = one field (mean yield over its observed years). Zones are terciles of field NCCPI computed within the
# state from the pooled corn and soybean fields (same rule as build_state_generators.R); cut-points are read from
# results/state_zone_cuts.csv. Run from simulation_code/. Figures -> FIGS_DIR (SIM_FIGS_DIR=../figs).
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(arrow) })
source("config.R")
PANEL_DIR <- .env("SIM_PANEL_DIR", "D:/Crop data")
STATES <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
KG_CORN <- 0.01593; KG_SOY <- 0.01487
cuts <- fread("results/state_zone_cuts.csv")
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive = TRUE)

read_fields <- function(path, ycol, conv) {
  d <- as.data.table(read_parquet(path, col_select = c("tile", "field_id", ycol, "nccpi3all_mean")))
  d[, uid := paste(tile, field_id, sep = ":")]
  d[, .(nccpi = mean(nccpi3all_mean, na.rm = TRUE), yield = mean(get(ycol), na.rm = TRUE) * conv), by = uid]
}
res <- list()
for (st in STATES) {
  cat(st, "")
  cn <- read_fields(file.path(PANEL_DIR, sprintf("d_igis13_%s.parquet", st)), "corn_yield", KG_CORN)
  sy <- read_fields(file.path(PANEL_DIR, sprintf("d_igis13soy_%s.parquet", st)), "soy_yield", KG_SOY)
  res[[length(res) + 1]] <- rbind(cn[, crop := "Corn"], sy[, crop := "Soybeans"])[, state := st]
}
cat("\n")
d <- rbindlist(res)[is.finite(nccpi) & is.finite(yield)]
d <- merge(d, cuts, by = "state")
d[, zone := factor(fifelse(nccpi >= p67, "high", fifelse(nccpi >= p33, "medium", "low")), c("low", "medium", "high"))]
d[, state := factor(state, STATES)]
zone_cols <- c(low = "#BF812D", medium = "#C7C7C7", high = "#01665E")

# zone summary table
zs <- d[, .(n_fields = .N, nccpi_mean = mean(nccpi), yield_mean = mean(yield), yield_sd = sd(yield)), by = .(state, crop, zone)]
fwrite(zs[order(state, crop, zone)], "results/nccpi_zone_summary.csv")
print(zs[crop == "Corn"][order(state, zone)][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

th <- theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
cuts_f <- copy(cuts)[state %in% STATES][, state := factor(state, STATES)]
for (cr in c("Corn", "Soybeans")) {
  x <- d[crop == cr]
  # binned mean yield in 25 equal-width NCCPI bins within each state (>= 30 fields per bin)
  x[, bin := cut(nccpi, 25), by = state]
  bm <- x[, .(nccpi = mean(nccpi), yield = mean(yield), n = .N), by = .(state, bin)][n >= 30]
  ymax <- quantile(x$yield, 0.999)
  p <- ggplot(x, aes(nccpi, yield)) +
    geom_bin2d(bins = 70) + scale_fill_viridis_c(trans = "log10", name = "Fields", option = "mako", direction = -1) +
    geom_vline(data = cuts_f, aes(xintercept = p33), linetype = "dashed", colour = "black", linewidth = 0.35) +
    geom_vline(data = cuts_f, aes(xintercept = p67), linetype = "dashed", colour = "black", linewidth = 0.35) +
    geom_line(data = bm, aes(nccpi, yield), colour = "#D73027", linewidth = 0.7) +
    facet_wrap(~state, nrow = 2, scales = "free_x") + coord_cartesian(ylim = c(0, ymax)) +
    labs(x = "Field NCCPI (dashed lines: tercile cut-points separating low, medium, high zones)",
         y = sprintf("Mean %s yield of the field (bu/ac)", tolower(cr)), subtitle = "Red line: mean yield in 25 equal-width NCCPI bins") + th
  ggsave(file.path(FIGS_DIR, sprintf("be_nccpi_yield_%s.pdf", if (cr == "Corn") "corn" else "soy")), p, width = 10, height = 5.4)
}

# distribution of field NCCPI with the zones shaded
fl <- d[, .SD[!duplicated(uid)], by = state][, .(state, uid, nccpi, p33, p67)]
fl[, zone := factor(fifelse(nccpi >= p67, "high", fifelse(nccpi >= p33, "medium", "low")), c("low", "medium", "high"))]
p3 <- ggplot(fl, aes(nccpi, fill = zone)) + geom_histogram(bins = 60, colour = NA) +
  scale_fill_manual(values = zone_cols, name = "Zone (NCCPI tercile)") + facet_wrap(~state, nrow = 2, scales = "free") +
  labs(x = "Field NCCPI", y = "Fields") + th
ggsave(file.path(FIGS_DIR, "be_nccpi_dist.pdf"), p3, width = 10, height = 5.2)
cat("done\n")
