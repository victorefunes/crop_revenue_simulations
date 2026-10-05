# ============================================================================
# Interest-rate and accrual-window sensitivity of six-year delinquency.
#   Inputs -> mechanism -> outcome:
#   i_op_annual, months_borrowed -> growth = (1 + i/12)^months (credit_functions_region_zone.R)
#   -> harvest debt H = L * growth, B = H - min(R, H) -> delinquency B > 0.
# Correlated price-yield draws, initial history CCCCCC, all 27 cells, same seeds as the
# baseline corr grid (seed_base = 102), so the (9%, 9 months) scenario reproduces the
# baseline exactly. Writes results/credit/interest_rate_sensitivity.csv.
#   setwd("<path>/simulation_code"); source("interest_rate_sensitivity.R")
# ============================================================================
suppressMessages({ library(data.table); library(tidyverse) })
source("config.R")
source("R/functions_region_zone.R"); source("R/correlation_functions_region_zone.R")
source("R/cashflow_functions_region_zone.R"); source("R/credit_functions_region_zone.R")

yield_gens_by_region_zone <- readRDS(YGEN_CACHE)
base_cost_params          <- readRDS(COST_PARAMS_RDS)
price_params_regional     <- readRDS(PRICE_PARAMS_RDS)

# Same regional correlations as run_simulation.R (median of county correlations)
.c <- fread(CORN_CORR_CSV, select = c("fips", "corr"))[fips > 17000 & fips < 18000]; setnames(.c, "corr", "corr_corn")
.s <- fread(SOY_CORR_CSV,  select = c("fips", "corr"))[fips > 17000 & fips < 18000]; setnames(.s, "corr", "corr_soy")
.m <- fread(FBFM_COUNTIES_CSV, select = c("fips", "region")); .m[, region := tolower(region)]
region_corr <- Reduce(function(a, b) merge(a, b, by = "fips", all = TRUE), list(.m, .c, .s))[
  !is.na(region), .(corr_corn = median(corr_corn, na.rm = TRUE), corr_soy = median(corr_soy, na.rm = TRUE)), by = region]
setkey(region_corr, region)

B <- as.integer(Sys.getenv("SIM_B", "50000"))
regions <- c("northern", "central", "southern"); zones <- c("high", "medium", "low")
plans <- c("CCCCCC", "CSCSCS", "CCSCCS")
scenarios <- data.table(scenario = c("i=6%", "i=9% (baseline)", "i=12%", "9% / 6 months", "9% / 12 months"),
                        i = c(.06, .09, .12, .09, .09), months = c(9, 9, 9, 6, 12))

out <- list(); cell_idx <- 0L
for (r in regions) for (z in zones) {
  cell_idx <- cell_idx + 1L; seed <- 102L * 1000L + cell_idx
  for (p in plans) for (k in seq_len(nrow(scenarios))) {
    sc <- scenarios[k]
    res <- simulate_plan_credit_region_zone(
      history6 = "CCCCCC", plan6 = p, geo_region = r, prod_zone = z, B = B,
      i_op_annual = sc$i, months_borrowed = sc$months, credit_limit = 2000,
      delinquency_rule = "end_balance_positive", corr_method = "decomposed", seed = seed,
      correlation_corn = region_corr[r]$corr_corn, correlation_soy = region_corr[r]$corr_soy)
    out[[length(out) + 1]] <- data.table(scenario = sc$scenario, i = sc$i, months = sc$months, region = r,
      prod_zone = z, plan = p, pr_delinq = res$pr_delinquent_any, mean_peak_debt = res$mean_peak_debt)
    rm(res)
  }
  cat(sprintf("done %s-%s\n", r, z))
}
tab <- rbindlist(out)
dir.create("results/credit", recursive = TRUE, showWarnings = FALSE)
fwrite(tab, "results/credit/interest_rate_sensitivity.csv")
