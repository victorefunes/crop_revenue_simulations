# ============================================================================
# STAGE 1 (run once): build the yield generators from the raw field-panel parquet
# files and write the three objects that run_simulation.R loads:
#   cache/yield_generators.rds   rotation-conditioned yield generators (slimmed)
#   data/cost_params.rds         base operating costs by region x zone
#   data/price_params.rds        lognormal price parameters by region
#
# Needs the raw panels (see config.R: SIM_CORN_PARQUET, SIM_SOY_PARQUET).
# Typical run time: 10-20 minutes (parquet I/O + two-way FE fits for 18 cells).
# Run from the simulation_code/ folder.
# ============================================================================

source("config.R")
source("R/functions_region_zone.R")   # attach_simulate()
source("R/slim_generators.R")
source("R/yield_generator.R")         # builds yield_gens_by_region_zone,
                                      # base_cost_params, price_params_regional

yield_gens_slim <- slim_generators(yield_gens_by_region_zone)

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)
saveRDS(yield_gens_slim,       YGEN_CACHE)
saveRDS(base_cost_params,      COST_PARAMS_RDS)
saveRDS(price_params_regional,  PRICE_PARAMS_RDS)

cat(sprintf("\nDone.\n  %s (%.1f MB)\n  %s\n  %s\nNext: source('run_simulation.R')\n",
            YGEN_CACHE, file.size(YGEN_CACHE) / 1e6, COST_PARAMS_RDS, PRICE_PARAMS_RDS))
