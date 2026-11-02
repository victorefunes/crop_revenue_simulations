# VPD regimes by state and calendar year, and how well the VPD index explains the estimated yield year effects.
# Index: mean over fields of the July-August maximum vapor-pressure deficit (vpdmax_7, vpdmax_8; suffix = calendar month),
# averaged over the corn and soybean panels. Regimes are terciles of the index across the state's calendar years
# (low / mid / high VPD; k = round(J/3) years in each tail, same rule as weather_types()).
# Reads the v2 generators only for their year effects. Run from simulation_code/.
# Writes results/vpd_regimes.csv and results/vpd_year_effect_fit.csv.
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("config.R")
source("R/functions_region_zone.R")
source("R/breakeven_sim.R")        # weather_types()
PANEL_DIR <- .env("SIM_PANEL_DIR", "D:/Crop data")
STATES <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
gens <- readRDS(file.path(CACHE_DIR, "state_generators.rds"))

state_year_vpd <- function(path) {
  cols <- intersect(c("year", "vpdmax_7", "vpdmax_8", "vpd_7", "vpd_8"), names(open_dataset(path)$schema))
  d <- as.data.table(read_parquet(path, col_select = all_of(cols)))
  d[, idx := rowMeans(cbind(vpdmax_7, vpdmax_8), na.rm = TRUE)]
  d[, idx_mean := if (all(c("vpd_7", "vpd_8") %in% names(d))) rowMeans(cbind(vpd_7, vpd_8), na.rm = TRUE) else NA_real_]
  d[, .(idx = mean(idx, na.rm = TRUE), idx_mean = mean(idx_mean, na.rm = TRUE)), by = year]
}

out <- list(); fit <- list()
for (st in STATES) {
  a <- state_year_vpd(file.path(PANEL_DIR, sprintf("d_igis13_%s.parquet", st)))
  b <- state_year_vpd(file.path(PANEL_DIR, sprintf("d_igis13soy_%s.parquet", st)))
  v <- merge(a, b, by = "year", suffixes = c("_c", "_s"))
  v[, `:=`(idx = (idx_c + idx_s) / 2, idx_mean = (idx_mean_c + idx_mean_s) / 2)]
  wt <- weather_types(gens$corn[[st]][["medium"]]$year_fe, gens$soy[[st]][["medium"]]$year_fe)
  y <- data.table(year = wt$years, ce = wt$ce, se = wt$se, w = wt$w, type = wt$type)
  m <- merge(v[, .(year, idx, idx_mean)], y, by = "year")
  J <- nrow(m); k <- round(J / 3); r <- rank(m$idx, ties.method = "first")
  m[, regime := ifelse(r <= k, "low", ifelse(r > J - k, "high", "mid"))]
  m[, state := st]
  out[[st]] <- m
  fit[[st]] <- data.table(state = st, years = J,
    cor_idx_corn_fe = cor(m$idx, m$ce), cor_idx_soy_fe = cor(m$idx, m$se), cor_idx_weather = cor(m$idx, m$w),
    cor_idxmean_weather = cor(m$idx_mean, m$w, use = "complete.obs"),
    high_vpd_is_bad = mean(m$regime == "high" & m$type == "bad") / pmax(mean(m$regime == "high"), 1e-9),
    low_vpd_is_good = mean(m$regime == "low" & m$type == "good") / pmax(mean(m$regime == "low"), 1e-9))
}
fwrite(rbindlist(out), "results/vpd_regimes.csv"); fit <- rbindlist(fit); fwrite(fit, "results/vpd_year_effect_fit.csv")
print(fit[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
print(rbindlist(out)[state == "IA", .(year, idx = round(idx, 2), regime, type)])
