# Reduce yield generators to what the simulation actually reads.
#
# Raw generators from build_yield_generator_region_zone() embed the full field-year
# table (DT_fy), the fixest model, and a $simulate closure whose environment holds the
# whole panel, so the serialized object is several GB. The simulation only uses
# $mu, $pools and $simulate (draw mu[state] + resample of the residual pool), plus
# the geo_region / prod_zone labels. Nothing downstream touches DT_fy or the model.
#
# Requires attach_simulate() from R/functions_region_zone.R (sourced before this file).
# $simulate is rebuilt from the slim object, so its environment holds only mu and pools.
slim_generators <- function(gens) {
  keep <- c("yield_col", "terminal_crop", "geo_region", "prod_zone", "mu", "pools")
  lapply(gens, function(by_region)
    lapply(by_region, function(by_zone)
      lapply(by_zone, function(g) attach_simulate(g[intersect(keep, names(g))]))))
}
