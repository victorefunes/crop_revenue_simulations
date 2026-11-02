# ============================================================================
# Rotation-conditioned yield generators, built state by state.
# Same estimator as R/yield_generator.R (two-way field + year FE per cell; per-rotation-state
# residual pool; mean = mean field effect + one common cell-wide year effect), but:
#   * the cell is state x NCCPI zone (no sub-state regions),
#   * zones are terciles of field-level NCCPI (nccpi3all_mean), computed within the state,
#   * the cell filter is applied to a pre-filtered table, so no name shadowing is possible.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(fixest); library(stringr) })

KGHA_TO_BUAC_CORN <- 0.01593
KGHA_TO_BUAC_SOY  <- 0.01487

# 6-year crop history ending in the observation year, e.g. "CSCSCS" (NA if incomplete)
make_rot6_full <- function(DT) {
  crop_cols  <- sort(names(DT)[str_detect(names(DT), "^crop_\\d{4}$")])
  crop_years <- as.integer(str_remove(crop_cols, "^crop_"))
  pos_y <- match(DT[["year"]], crop_years)
  vals  <- unlist(DT[, ..crop_cols], use.names = FALSE)
  crop_bin <- ifelse(vals == "Corn", 0L, ifelse(vals == "Soybeans", 1L, NA_integer_))
  crop_bin <- matrix(crop_bin, nrow = nrow(DT), ncol = length(crop_cols))
  seq6 <- matrix(NA_integer_, nrow(DT), 6)
  for (j in 1:6) {
    cj <- pos_y + (j - 6)
    good <- !is.na(cj) & cj >= 1 & cj <= ncol(crop_bin)
    if (any(good)) seq6[good, j] <- crop_bin[cbind(which(good), cj[good])]
  }
  incomplete <- apply(is.na(seq6), 1, any)
  out <- rep(NA_character_, nrow(DT))
  if (any(!incomplete))
    out[!incomplete] <- apply(ifelse(seq6[!incomplete, , drop = FALSE] == 1L, "S", "C"), 1, paste0, collapse = "")
  out
}

# Read one state's panel for one crop, keep the columns the generator needs, convert to bu/ac
read_state_panel <- function(path, crop = c("corn", "soy")) {
  crop <- match.arg(crop)
  y_col <- if (crop == "corn") "corn_yield" else "soy_yield"
  sch <- names(arrow::open_dataset(path)$schema)
  keep <- c("tile", "field_id", "year", y_col, "nccpi3all_mean", intersect(c("vpdmax_7", "vpdmax_8"), sch),
            sch[str_detect(sch, "^crop_\\d{4}$")])
  DT <- as.data.table(arrow::read_parquet(path, col_select = dplyr::all_of(keep)))
  DT[, (y_col) := get(y_col) * if (crop == "corn") KGHA_TO_BUAC_CORN else KGHA_TO_BUAC_SOY]
  DT[, field_uid := paste(tile, field_id, sep = ":")]
  DT
}

# NCCPI terciles from the pooled field-level distribution of both crops' panels.
# A field's zone is fixed by its mean NCCPI, so it never changes across years.
assign_nccpi_zones <- function(corn_dt, soy_dt) {
  f <- rbind(corn_dt[, .(field_uid, nccpi = nccpi3all_mean)],
             soy_dt[,  .(field_uid, nccpi = nccpi3all_mean)])
  f <- f[, .(nccpi = mean(nccpi, na.rm = TRUE)), by = field_uid][is.finite(nccpi)]
  cuts <- quantile(f$nccpi, c(1/3, 2/3))
  f[, zone := fcase(nccpi >= cuts[2], "high", nccpi >= cuts[1], "medium", default = "low")]
  corn_dt[, prod_zone := f$zone[match(field_uid, f$field_uid)]]
  soy_dt[,  prod_zone := f$zone[match(field_uid, f$field_uid)]]
  list(corn = corn_dt, soy = soy_dt, cuts = cuts)
}

# VPD regimes by calendar year for one state, with the fixed bins of Pizzo (2025, Chapter 1, Appendix C), calibrated to
# Illinois from Xu et al. (2021): July maximum VPD (vpdmax_7) below 19 hPa = "normal" (stored as low), 19 to 21 hPa =
# "somewhat dry" (mid), above 21 hPa = "dry" (high). The panel stores VPD in kPa, so the cut-points are 1.9 and 2.1.
# The index is the field-mean vpdmax_7 of the calendar year, averaged over the corn and soybean panels.
# Returns data.table(year, idx, regime) with regime in low / mid / high.
VPD_CUTS_KPA <- c(1.9, 2.1)
vpd_regime_table <- function(corn_dt, soy_dt, cuts = VPD_CUTS_KPA) {
  yr <- function(d) d[, .(idx = mean(vpdmax_7, na.rm = TRUE)), by = year]
  v <- merge(yr(corn_dt), yr(soy_dt), by = "year", suffixes = c("_c", "_s"))
  v[, idx := (idx_c + idx_s) / 2]
  v[, regime := ifelse(idx < cuts[1], "low", ifelse(idx <= cuts[2], "mid", "high"))]
  v[, .(year, idx, regime)]
}

# regime_tab (optional): data.table(year, regime). When given, the generator also stores residual pools by
# (six-year window, VPD regime): pools_vpd. A (window, regime) cell with fewer than 200 observations falls back to the
# window's pooled residuals. The window-by-regime residual mean is the rotation-by-VPD interaction that the common
# year effect does not absorb.
build_cell_generator <- function(DT, yield_col, terminal_crop, label,
                                 seed = 1L, B_per_state = 20000L, regime_tab = NULL) {
  DT <- copy(DT)
  DT[, rot6 := make_rot6_full(DT)]
  DT_sub <- DT[!is.na(rot6) & substr(rot6, 6, 6) == terminal_crop]
  DT_fy <- DT_sub[!is.na(get(yield_col)),
                  .(y = mean(get(yield_col), na.rm = TRUE), n_tiles = .N),
                  by = .(field_uid, year, rot6)]
  if (nrow(DT_fy) == 0) stop("0 field-years in ", label)
  m  <- feols(y ~ 1 | field_uid + year, data = DT_fy, notes = FALSE)
  fe <- fixef(m)
  fe_year <- unname(fe$year[as.character(DT_fy$year)])
  fe_fld  <- unname(fe$field_uid[as.character(DT_fy$field_uid)])
  DT_fy[, resid := y - (fe_fld + fe_year)]
  ok <- is.finite(DT_fy$resid)
  DT_use <- DT_fy[ok]
  set.seed(seed)
  pools <- DT_use[, .(n = .N, resid_pool = list({
    x <- resid[is.finite(resid)]
    if (length(x) >= B_per_state) sample(x, B_per_state) else sample(x, B_per_state, replace = TRUE)
  })), by = rot6]
  DT_use[, fe_f := y - resid - fe_year[ok]]
  ybar_year <- mean(fe_year[ok])
  mu <- DT_use[, .(mu = mean(fe_f) + ybar_year, n = .N), by = rot6]
  cat(sprintf("  %-22s field-years %7d | fields %6d | rotation states %2d\n",
              label, nrow(DT_fy), uniqueN(DT_fy$field_uid), nrow(mu)))
  # Fitted year effects (one per calendar year in the cell): the common weather/market shock that the
  # residual pool does not contain. Used by the year-shock simulator (run_breakeven_years.R).
  year_fe <- unique(data.table(year = DT_fy$year, fe = fe_year))[order(year)]
  pools_vpd <- NULL
  if (!is.null(regime_tab)) {
    DT_use[, regime := regime_tab$regime[match(year, regime_tab$year)]]
    pools_vpd <- DT_use[!is.na(regime), .(n = .N, resid_pool = list({
      x <- resid[is.finite(resid)]
      if (length(x) >= B_per_state) sample(x, B_per_state) else sample(x, B_per_state, replace = TRUE)
    })), by = .(rot6, regime)]
    thin <- which(pools_vpd$n < 200)
    for (i in thin) pools_vpd$resid_pool[[i]] <- pools[rot6 == pools_vpd$rot6[i], resid_pool][[1]]
  }
  list(yield_col = yield_col, terminal_crop = terminal_crop, geo_region = label,
       mu = mu[, .(rot6, mu)], pools = pools[, .(rot6, n, resid_pool)], year_fe = year_fe, pools_vpd = pools_vpd)
}
