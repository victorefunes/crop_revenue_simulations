library(data.table)
library(tidyverse)
library(fixest)
library(statar)
library(broom)
library(arrow)
library(progressr)


setwd("C:/Users/vf006/Box/premiums")

prices <- fread("./county_prices/cash_price_with_counties.csv")
regions <- fread("./data-raw/counties_fips.csv")

prices |>
  mutate(fips = as.numeric(fips)) |>
  filter(cropLongName %in% c("Corn", "Soybeans") & month %in% 10:12 &
           fips > 17000 & fips < 17999) |>
  group_by(cropLongName, fips, month, year) |>
  summarise(price = mean(value, na.rm = TRUE)) |>
  ungroup() |>
  filter(!is.nan(price)) ->
  county_prices

county_prices |> 
    left_join(regions, by = "fips")  ->
    region_prices

region_prices <- region_prices |>
    data.table()


# dt columns: elevator_id, region, crop, year, price
# "region" = your Central/Northern/Southern IL assignment

prices |>
  mutate(fips = as.numeric(fips)) |>
  filter(cropLongName %in% c("Corn", "Soybeans") & month %in% 10:12 &
           fips > 17000 & fips < 17999) ->
  dt

# Step 1: de-mean each elevator to remove its permanent basis
dt[, log_price := log(value)]
dt[, log_price_dm := log_price - mean(log_price, na.rm = TRUE), by = .(location, crop)]

dt <- dt |> 
    left_join(regions, by = "fips")

# Step 2: sigma = sd of demeaned log prices, by region × crop
sigma_est <- dt[, .(
  sigma     = sd(log_price_dm, na.rm = TRUE),
  n_obs     = .N,
  n_years   = uniqueN(year),
  n_elevs   = uniqueN(location)
), by = .(region, crop)]

# Step 3: attach a 95% CI (chi-squared, n_obs-1 df)
sigma_est[, sigma_lo := sigma * sqrt((n_obs - 1) / qchisq(0.975, n_obs - 1))]
sigma_est[, sigma_hi := sigma * sqrt((n_obs - 1) / qchisq(0.025, n_obs - 1))]


## Monthly prices

# Step 0: extract the harvest-month price for each elevator-year
# October for corn, October or November for soybeans
dt <- prices |>
    mutate(fips = as.numeric(fips)) |>
    filter(cropLongName %in% c("Corn", "Soybeans") &
           fips > 17000 & fips < 17999) |> 
    left_join(regions, by = "fips") |>
    data.table()


harvest <- dt[
  (cropLongName == "Corn" & month == 10) |
  (cropLongName == "Soybeans" & month %in% c(10, 11)),
  .(price = mean(value, na.rm = TRUE)),   # average Oct-Nov for beans
  by = .(location, region, cropLongName, year)
]

# Then proceed exactly as before
harvest[, log_price := log(price)]
harvest[, log_price_dm := log_price - mean(log_price), by = .(location, cropLongName)]

sigma_est <- harvest[, .(
  sigma   = sd(log_price_dm, na.rm = TRUE),
  n_obs   = .N,
  n_years = uniqueN(year),
  n_elevs = uniqueN(location)
), by = .(region, cropLongName)]

sigma_est |> 
    readr::write_csv("./NPV_sims/yield_regions/sigma_regions.csv")