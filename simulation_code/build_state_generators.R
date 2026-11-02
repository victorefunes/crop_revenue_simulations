# ============================================================================
# STAGE 1 (breakeven paper): build rotation-conditioned yield generators, one state at a time.
#   Input : D:/Crop data/d_igis13_<ST>.parquet (corn), d_igis13soy_<ST>.parquet (soy)
#   Output: <CACHE_DIR>/state_generators.rds    list(corn/soy -> [[state]][[zone]])
#           results/state_sample_sizes.csv      observations per state x zone x crop x rotation state
#           results/state_zone_cuts.csv         NCCPI tercile cut-points by state
# Run from simulation_code/. Subset with Sys.setenv(SIM_STATES = "OH,WI").
# ============================================================================
source("config.R")
source("R/functions_region_zone.R")   # attach_simulate()
source("R/state_generator.R")

PANEL_DIR <- .env("SIM_PANEL_DIR", "D:/Crop data")
STATES <- strsplit(.env("SIM_STATES", "IA,IL,IN,MN,MO,OH,SD,WI"), ",")[[1]]
ZONES  <- c("high", "medium", "low")
OUT    <- file.path(CACHE_DIR, "state_generators.rds")
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)
if (!dir.exists("results"))  dir.create("results")

gens <- if (file.exists(OUT)) readRDS(OUT) else list(corn = list(), soy = list())
samp <- list(); cuts_out <- list(); vpd_out <- list()

for (st in STATES) {
  cat(sprintf("\n=== %s ===\n", st))
  corn <- read_state_panel(file.path(PANEL_DIR, sprintf("d_igis13_%s.parquet", st)), "corn")
  soy  <- read_state_panel(file.path(PANEL_DIR, sprintf("d_igis13soy_%s.parquet", st)), "soy")
  z <- assign_nccpi_zones(corn, soy)
  cuts_out[[st]] <- data.table(state = st, p33 = z$cuts[[1]], p67 = z$cuts[[2]])
  regs <- NULL                                    # VPD regimes by calendar year (needs the vpdmax_7 / vpdmax_8 columns)
  if (all(c("vpdmax_7", "vpdmax_8") %in% names(corn)) && all(c("vpdmax_7", "vpdmax_8") %in% names(soy))) {
    regs <- vpd_regime_table(corn, soy); regs[, state := st]; vpd_out[[st]] <- regs
  }
  gens$corn[[st]] <- list(); gens$soy[[st]] <- list()
  for (zone in ZONES) {
    lab <- paste(st, zone, sep = "-")
    gc_ <- build_cell_generator(z$corn[prod_zone == zone & is.finite(corn_yield)], "corn_yield", "C", paste(lab, "corn"), seed = 1L, regime_tab = regs)
    gs_ <- build_cell_generator(z$soy[ prod_zone == zone & is.finite(soy_yield)],  "soy_yield",  "S", paste(lab, "soy"),  seed = 2L, regime_tab = regs)
    samp[[length(samp) + 1]] <- rbind(
      cbind(state = st, zone = zone, crop = "C", gc_$pools[, .(rot6, n)]),
      cbind(state = st, zone = zone, crop = "S", gs_$pools[, .(rot6, n)]))
    gens$corn[[st]][[zone]] <- attach_simulate(gc_)
    gens$soy[[st]][[zone]]  <- attach_simulate(gs_)
  }
  rm(corn, soy, z); gc()
  saveRDS(gens, OUT)    # checkpoint after every state
}

fwrite(rbindlist(samp),     "results/state_sample_sizes.csv")
fwrite(rbindlist(cuts_out), "results/state_zone_cuts.csv")
if (length(vpd_out)) fwrite(rbindlist(vpd_out), "results/vpd_regimes_full.csv")
cat(sprintf("\nDone: %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6))
