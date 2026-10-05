# ============================================================================
# Central configuration for the simulation package.
# All scripts are meant to be run from THIS folder (the one containing config.R):
#   setwd("<path>/simulation_code")   # or open the folder as an RStudio project
# Every path below can be overridden with an environment variable, e.g.
#   Sys.setenv(SIM_CACHE_DIR = "D:/region_sims")
# ============================================================================

if (!file.exists("config.R") || !dir.exists("R")) {
  stop("Run from the simulation_code/ folder: setwd() there first (config.R and R/ must be visible).")
}

.env <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

# ---- Stage 2 inputs (run_simulation.R) -------------------------------------
# Cached yield generators. Too large for the repo; build once with
# build_yield_generators.R (from raw parquet) or slim_yield_generators.R (from an
# existing cache). Expected: a list with $corn / $soy -> [[region]][[zone]] generators.
CACHE_DIR  <- .env("SIM_CACHE_DIR", "cache")
YGEN_CACHE <- file.path(CACHE_DIR, "yield_generators.rds")

# Small parameter files shipped with the package.
COST_PARAMS_RDS   <- file.path("data", "cost_params.rds")
PRICE_PARAMS_RDS  <- file.path("data", "price_params.rds")
CORN_CORR_CSV     <- file.path("data", "corn_correlation.csv")
SOY_CORR_CSV      <- file.path("data", "soy_correlation.csv")
FBFM_COUNTIES_CSV <- file.path("data", "FBFM_counties.csv")

# ---- Stage 1 inputs (build_yield_generators.R only) ------------------------
CORN_PARQUET      <- .env("SIM_CORN_PARQUET", file.path("raw", "d_igis13_12_1_2025.parquet"))
SOY_PARQUET       <- .env("SIM_SOY_PARQUET",  file.path("raw", "d_igis13soy_11_30_2025.parquet"))
COUNTIES_FIPS_CSV <- file.path("data", "counties_fips.csv")

# ---- Outputs ----------------------------------------------------------------
# Figures. Defaults to ./figs so a run never overwrites the paper's figures;
# point SIM_FIGS_DIR at ../figs to refresh the manuscript.
# Intermediate results are written to ./results/credit/ (relative paths in run_simulation.R).
FIGS_DIR <- .env("SIM_FIGS_DIR", "figs")

# ---- Run controls -----------------------------------------------------------
# Monte Carlo paths per cell. The paper uses 50,000; lower it (e.g. 500) for a smoke test.
SIM_B <- as.integer(.env("SIM_B", "50000"))

# The dynamic-programming extension at the end of run_simulation.R is NOT used in the
# paper (the paper holds the rotation plan fixed). Set TRUE to run it.
RUN_DP_EXTENSION <- identical(.env("SIM_RUN_DP", "FALSE"), "TRUE")
