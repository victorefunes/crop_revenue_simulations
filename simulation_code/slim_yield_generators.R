# ============================================================================
# Convert an existing (multi-GB) yield-generator cache into the slim cache that
# run_simulation.R loads, without needing the raw parquet panels.
#
# Usage (from the simulation_code/ folder):
#   Rscript slim_yield_generators.R "D:/region_sims/yield_generators.rds"
# or set SRC below. Writes <CACHE_DIR>/yield_generators.rds (see config.R).
# The slim object draws exactly the same yields as the original for a given seed:
# both compute mu[state] + sample(resid_pool[state], B, replace = TRUE).
# ============================================================================

source("config.R")
source("R/functions_region_zone.R")   # attach_simulate()
source("R/slim_generators.R")

args <- commandArgs(trailingOnly = TRUE)
SRC  <- if (length(args) >= 1) args[1] else file.path("D:/region_sims", "yield_generators.rds")
stopifnot(file.exists(SRC))

cat("Reading", SRC, "(this can take a few minutes) ...\n")
gens <- readRDS(SRC)
gens <- slim_generators(gens)

if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)
saveRDS(gens, YGEN_CACHE)
cat(sprintf("Wrote %s (%.1f MB)\n", YGEN_CACHE, file.size(YGEN_CACHE) / 1e6))
