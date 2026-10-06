# Simulation code

Code behind the simulation in `../crop_rotation_delinquencies.tex`. Everything needed to run the
model is in this folder except two large inputs that are not in the repo (see *Large inputs*).

## The model in three lines

- **Inputs:** rotation-conditioned yield generator (mean + empirical residual pool for each 6-year rotation
  state, by region x productivity zone), regional lognormal prices and costs, price-yield dependence,
  credit terms (9% APR, 9-month accrual).
- **Mechanism:** each year `L = C + B_prev`, `H = L (1+i)^(9/12)`, `R = P*Y`, `M = min(R, H)`,
  `B = H - M`; an unpaid balance rolls into next year's loan. Delinquency in year t is `B_t > 0`.
- **Outcome:** six-year delinquency probability, per-year hazard, survival, peak loan demand, hedge effect
  (correlated vs. independent prices), yield-shock stress test, variance and Shapley decompositions.

## Layout

| Path | What it is |
|---|---|
| `config.R` | All paths and run options (environment-variable overrides). Sourced by every script. |
| `run_simulation.R` | **Stage 2 driver.** Runs the full grid and writes results and figures. |
| `build_yield_generators.R` | **Stage 1** (once): raw parquet panels -> slim `cache/yield_generators.rds`. |
| `slim_yield_generators.R` | Alternative to Stage 1: shrink an existing multi-GB cache, no raw data needed. |
| `R/functions_region_zone.R` | Shared helpers (transition paths, `attach_simulate`, summaries). |
| `R/correlation_functions_region_zone.R` | Price-yield dependence (simple, decomposed, copula, t-copula). |
| `R/cashflow_functions_region_zone.R` | Cash-flow simulation and plotting helpers. |
| `R/credit_functions_region_zone.R` | **Operating-credit recursion** (`simulate_plan_credit_region_zone`). |
| `R/yield_generator.R` | Builds the yield generators, base costs and price parameters (Stage 1). |
| `R/slim_generators.R` | Strips generators to what the simulation reads. |
| `R/dp_functions_region_zone.R` | Dynamic-programming extension. Not used in the paper; off by default. |
| `data/` | Small inputs: price-yield correlations by county, FBFM county-region map, county FIPS map, `cost_params.rds`, `price_params.rds`. |
| `calibration/` | `sigma_estimation.r` and its output `sigma_regions.csv` (source of the price volatilities). |

## Requirements

R (tested with 4.6.1). Packages: tidyverse, data.table, patchwork, RColorBrewer, scales, mvtnorm, fixest,
arrow, progressr, statar, MASS, fitdistrplus, stringr, tibble, pacman. The scripts install `pacman` and
whatever it loads if missing; set a CRAN mirror first when running non-interactively.

## Large inputs (not in the repo)

1. **Yield generators**, a list `$corn` / `$soy` -> `[[region]][[zone]]` (about 4 GB in the original cache,
   far smaller once slimmed). Get them one of two ways:
   - *Existing cache:* `Rscript slim_yield_generators.R "D:/region_sims/yield_generators.rds"`
   - *From raw data:* put the two field-panel parquet files in `raw/` (or set `SIM_CORN_PARQUET`,
     `SIM_SOY_PARQUET`) and run `source("build_yield_generators.R")` (10-20 min).
   Both write `cache/yield_generators.rds` (or `SIM_CACHE_DIR`).
2. Raw parquet panels, only for the "from raw data" route.

## Running

```r
setwd("<path>/simulation_code")        # all scripts assume this is the working directory
Sys.setenv(SIM_B = "500")              # optional smoke test; the paper uses 50,000 (the default)
source("run_simulation.R")
```

Intermediate results go to `results/credit/`, figures to `figs/` (so the paper's `../figs/` is never
overwritten by accident; set `SIM_FIGS_DIR=../figs` to refresh the manuscript).

| Variable | Default | Meaning |
|---|---|---|
| `SIM_CACHE_DIR` | `cache` | folder holding `yield_generators.rds` |
| `SIM_FIGS_DIR` | `figs` | figure output folder |
| `SIM_B` | `50000` | Monte Carlo paths per cell |
| `SIM_RUN_DP` | `FALSE` | run the dynamic-programming extension at the end of the driver |
| `SIM_CORN_PARQUET`, `SIM_SOY_PARQUET` | `raw/...` | raw panels (Stage 1 only) |

Seeds: `GLOBAL_SEED = 21212`; each region x zone cell uses `seed_base * 1000 + cell index`, shared across
plans within a cell (common random numbers).

## Figures produced for the paper

`fig_p_bar_all`, `fig_p_year_all`, `fig_p_surv_all`, `fig_peak_loan_all`, `fig_p_debt_all`, `fig_p_bar_corr`,
`fig_yield_compare_corr_vs_indep`, `fig_revenue_variance_decomp`, `fig_shapley_price_share`, `fig_revenue_cv`,
`fig_p_cf_all`, `fig_shapley_abs`. `ltd_dominance_check.R` (run separately, after Stage 1) produces
`ltd_cdf_comparison.{pdf,png}`, `results/ltd_dominance_results.csv` (lower-tail dominance check) and
`results/sample_sizes.csv` (observations per cell and rotation state). `variance_components.R` (needs the raw
panels) computes the field / year / residual variance shares behind `tab:decomp`. Not produced here: `map_practice.png`
(maps script), which is not part of the simulation model.

## Implementation notes (as coded)

- Credit terms live in `simulate_plan_credit_region_zone()`: `i_op_annual = 0.09`, `months_borrowed = 9`;
  the driver passes a non-binding `credit_limit = 2000` and `delinquency_rule = "end_balance_positive"`.
  Credit-limit sensitivity is applied afterwards to the peak-debt distribution.
- `corr_method = "none"` resamples the empirical residual pool. `corr_method = "decomposed"` draws a normal
  aggregate yield shock (variance share `aggregate_fraction = 0.40`) correlated with the log-price shock, plus a
  normal idiosyncratic shock. Its correlation is the regional median of the county-level correlations in
  `data/corn_correlation.csv` / `soy_correlation.csv` (defaults -0.40 corn, -0.50 soy).
- Yield generators (`R/yield_generator.R`): per region x zone cell, a two-way fixed-effects fit (field and year)
  on field-year yields; each rotation state keeps a pool of residuals (20,000 draws) and a mean equal to the
  mean field effect of the state's observations plus **one common, cell-wide year effect**. Earlier versions
  averaged each state's own year effects, so thin states inherited their drought-year mix (2011-12); that
  produced spurious state means (for example a 12.7 bu/ac penalty for CCCSCC in southern-low).
- Prices are hard-coded regional means with log-sd from `calibration/sigma_regions.csv`; costs are hard-coded
  2025 FBFM direct costs (see `R/yield_generator.R`). `data/cost_params.rds` and `data/price_params.rds` are
  exactly what that code produces.

## Changes from the original scripts

The originals in `yield_regions/` are untouched. This folder differs only in that:

- absolute paths (`setwd`, `D:/...`, `Box/...`) are replaced by `config.R`;
- the Monte Carlo size comes from `SIM_B` (default unchanged, 50,000);
- the dynamic-programming block at the end of the driver runs only if `SIM_RUN_DP=TRUE`;
- `build_yield_generators.R` is rewritten: the original wrote an older list format that the driver does not
  read, and stored generators that embed the whole panel. The new one writes the three files the driver
  loads, with generators slimmed to `mu`, `pools` and `simulate`;
- in `R/yield_generator.R`, a read of `crop_prices_il.csv` / `direct_cost.csv` whose result was discarded is
  removed (it needed a 211 MB file that never entered the model);
- the dead line `rm(ygen)` is removed;
- two stale figure subtitles are corrected (the `fig_p_debt_all` subtitle no longer claims the gap is wider for
  CSCSCS, and the variance-decomposition subtitle no longer points to a figure the paper does not contain).
  Figures must be regenerated to show this.
