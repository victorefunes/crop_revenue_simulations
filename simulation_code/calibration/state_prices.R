# State-level price parameters for the breakeven analysis.
#   log_sd      : sd of elevator-demeaned log harvest-season cash prices, by state x crop
#                 (same method as sigma_estimation.r, extended from IL to all 8 states)
#   price_per_bu: common reference level (the 2025 budget prices, $4.20 corn / $10.25 soy, so that
#                 revenue and cost are in the same year) times a state basis relative to IL. The basis is the state
#                 coefficient from log(price) ~ state + year FE on harvest-season elevator prices.
# Run from simulation_code/. Writes data/state_prices.csv.
library(data.table)

PRICES_CSV <- Sys.getenv("SIM_CASH_PRICES", "D:/Crop data/cash_price_with_counties.csv")
ref <- c(Corn = 4.20, Soybeans = 10.25)   # 2025 farmdoc budget benchmark prices, matching the 2025 cost year
st_fips <- c(IL = 17, IN = 18, IA = 19, MN = 27, MO = 29, OH = 39, SD = 46, WI = 55)

p <- fread(PRICES_CSV, select = c("value", "location", "month", "year", "cropLongName", "fips"))
p[, fips := as.numeric(fips)]
p[, st_code := floor(fips / 1000)]
p <- p[cropLongName %in% c("Corn", "Soybeans") & st_code %in% st_fips & is.finite(value) & value > 0]
p <- p[(cropLongName == "Corn" & month == 10) | (cropLongName == "Soybeans" & month %in% c(10, 11))]
p[, state := names(st_fips)[match(st_code, st_fips)]]
h <- p[, .(price = mean(value)), by = .(location, state, crop = cropLongName, year)]
h[, lp := log(price)]
h[, lp_dm := lp - mean(lp), by = .(location, crop)]

out <- rbindlist(lapply(c("Corn", "Soybeans"), function(cr) {
  d <- h[crop == cr]
  d[, state := relevel(factor(state), ref = "IL")]
  m <- lm(lp ~ state + factor(year), data = d)
  b <- c(IL = 0, coef(m)[grep("^state", names(coef(m)))])
  names(b) <- sub("^state", "", names(b))
  sg <- d[, .(log_sd = sd(lp_dm), n_obs = .N, n_elev = uniqueN(location), n_years = uniqueN(year)), by = state]
  sg[, `:=`(state = as.character(state), crop = ifelse(cr == "Corn", "C", "S"))]
  sg[, price_per_bu := round(ref[cr] * exp(b[state]), 2)]
  sg
}))
setcolorder(out, c("state", "crop", "price_per_bu", "log_sd"))
setorder(out, state, crop)
print(out)
fwrite(out, "data/state_prices.csv")
