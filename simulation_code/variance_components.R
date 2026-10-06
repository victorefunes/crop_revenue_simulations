# Variance components behind tab:decomp: R2 of field effects, year effects and both, on corn field-year yields.
# Run from simulation_code/ with SIM_CORN_PARQUET / SIM_SOY_PARQUET set (needs the raw panels, ~10 min).
# Sample B (generator sample, pooled) gives field 0.195, year (both - field) 0.620, residual 0.184.
suppressMessages({library(data.table); library(tidyverse); library(fixest)}); options(width=200)
source("config.R"); source("R/functions_region_zone.R"); source("R/slim_generators.R")
L <- readLines("R/yield_generator.R"); cut <- grep("^# BUILD GENERATORS FOR ALL REGION", L)[1] - 2
eval(parse(text = L[1:cut]))
r2s <- function(DT) { c(field=unname(r2(feols(y ~ 1 | field_uid, DT), "r2")), year=unname(r2(feols(y ~ 1 | year, DT), "r2")),
                        both=unname(r2(feols(y ~ 1 | field_uid + year, DT), "r2")), n=nrow(DT)) }
DT <- as.data.table(corn_arranged)[!is.na(geo_region) & !is.na(prod_zone) & !is.na(corn_yield)]
DT[, field_uid := paste(tile, field_id, sep=":")]
# sample A: every corn field-year
A <- DT[, .(y=mean(corn_yield)), by=.(field_uid, year, geo_region, prod_zone)]
# sample B: generator sample (terminal crop corn, complete six-year history), as in build_yield_generator_region_zone
res <- list()
cat("== sample A: all corn field-years, statewide\n"); a <- r2s(A); print(round(a,4))
cat("== sample A: per cell\n"); pa <- A[, as.list(r2s(.SD)), by=.(geo_region, prod_zone)]; print(pa[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)])
B <- list()
for (r in unique(DT$geo_region)) for (z in unique(DT$prod_zone)) { d <- DT[geo_region==r & prod_zone==z]; d[, rot6 := make_rot6_full(d)]
  d <- d[!is.na(rot6) & substr(rot6,6,6)=="C", .(y=mean(corn_yield)), by=.(field_uid, year, rot6)]; d[, `:=`(geo_region=r, prod_zone=z)]; B[[length(B)+1]] <- d }
B <- rbindlist(B)
cat("== sample B: generator sample, pooled across cells\n"); b <- r2s(B[, .(field_uid, year, y)]); print(round(b,4))
cat("== sample B: per cell\n"); pb <- B[, as.list(r2s(.SD)), by=.(geo_region, prod_zone)]; print(pb[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)])
cat("additivity check A: field+year =", round(a["field"]+a["year"],3), "vs both =", round(a["both"],3), "\n")
fwrite(pa, "results/variance_components_allfieldyears.csv"); fwrite(pb, "results/variance_components_generator_sample.csv")
