# Check on the price volatility in sigma_estimation.r: how much of the variance of elevator-demeaned log
# harvest prices is inter-annual (common regional year effect) versus cross-elevator dispersion?
# Needs the elevator price file (3.4 GB): set SIM_PRICE_CSV; counties_fips.csv comes from ../data.
#   Sys.setenv(SIM_PRICE_CSV = "D:/Crop data/cash_price_with_counties.csv"); source("calibration/sigma_decomposition.R")
# Result (2015-2024): pooled sigma reproduces sigma_regions.csv; the sd of region-year means is 0.213-0.235
# (corn) and 0.176-0.182 (soy); year effects are >99% of the variance.
suppressMessages({ library(data.table); library(fixest) })
regions <- fread(file.path("data", "counties_fips.csv"))
p <- fread(Sys.getenv("SIM_PRICE_CSV", "raw/cash_price_with_counties.csv"),
           select = c("value", "location", "month", "year", "cropLongName", "fips"))
p[, fips := suppressWarnings(as.numeric(fips))]
dt <- merge(p[cropLongName %in% c("Corn", "Soybeans") & fips > 17000 & fips < 17999], regions[, .(fips, region)], by = "fips", all.x = TRUE)
harvest <- dt[(cropLongName == "Corn" & month == 10) | (cropLongName == "Soybeans" & month %in% c(10, 11)),
              .(price = mean(value, na.rm = TRUE)), by = .(location, region, cropLongName, year)]
harvest[, log_price := log(price)]
harvest[, lp_dm := log_price - mean(log_price), by = .(location, cropLongName)]      # as in sigma_estimation.r
cur <- harvest[, .(sigma = sd(lp_dm, na.rm = TRUE), n_obs = .N), by = .(region, cropLongName)]
ym <- harvest[!is.na(lp_dm), .(m = mean(lp_dm)), by = .(region, cropLongName, year)][, .(sd_year_mean = sd(m)), by = .(region, cropLongName)]
res <- merge(cur, ym, by = c("region", "cropLongName"))
fe <- rbindlist(lapply(split(harvest[is.finite(log_price)], by = c("region", "cropLongName")), function(d) {
  m <- feols(log_price ~ 1 | location + year, data = d, notes = FALSE)
  data.table(region = d$region[1], cropLongName = d$cropLongName[1], sd_resid = sqrt(mean(resid(m)^2))) }))
res <- merge(res, fe, by = c("region", "cropLongName")); res[, year_share_of_var := 1 - sd_resid^2 / sd_year_mean^2]
print(res); dir.create("results", showWarnings = FALSE); fwrite(res, "results/sigma_decomposition.csv")
