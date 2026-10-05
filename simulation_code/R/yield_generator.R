# ============================================================
# CORN + SOY YIELD GENERATORS WITH PRODUCTIVITY ZONES
#   - Three productivity zones: High, Medium, Low
#   - Zone-specific yield distributions (separate FE models)
#   - Zone-specific costs (both base costs AND rotation adjustments)
#   - Regional variation within zones
#
# New Structure:
#   - Productivity zones determine BOTH yield potential AND cost structure
#   - High productivity → higher yields, higher costs
#   - Low productivity → lower yields, lower costs
#   - Rotation effects apply within each zone
# ============================================================

# ---- packages ------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table, stringr, fixest, arrow, tibble, tidyverse, progressr, statar
)

# (paths come from config.R; no setwd here)

# ---- user paths ----------------------------------------------------------
corn_path <- CORN_PARQUET
soy_path  <- SOY_PARQUET
regions <- fread(COUNTIES_FIPS_CSV)

# ---- unit conversions ----------------------------------------------------
KGHA_TO_BUAC_CORN <- 0.01593
KGHA_TO_BUAC_SOY  <- 0.01487

# ============================================================
# 0) GEOGRAPHIC REGION CLASSIFICATION
# ============================================================

classify_geographic_region <- function(DT) {
  #Classify counties into North, Central, or South Illinois regions
  #Based on the 'region' column from the regions file
  
  if (!"region" %in% names(DT)) {
    stop("'region' column not found in data. Check regions file merge.")
  }
  
  # The regions file should have: north, central, south
  # Standardize naming
  DT[, geo_region := tolower(trimws(region))]
  
  # Validate regions
  unique_regions <- unique(DT$geo_region)
  expected <- c("northern", "central", "southern")
  
  if (!all(unique_regions %in% c(expected, NA))) {
    warning("Unexpected region values: ", 
            paste(setdiff(unique_regions, expected), collapse=", "))
  }
  
  # Report distribution
  region_dist <- DT[!is.na(geo_region), .N, by = geo_region]
  cat("\nGeographic region distribution:\n")
  print(region_dist)
  
  DT
}


# Define productivity zone thresholds based on NCCPI or average yields
# Option 1: Use NCCPI (National Commodity Crop Productivity Index) if available
# Option 2: Use average field yields to classify
# Option 3: Use soil series or other field characteristics

classify_productivity_zone_by_region <- function(DT, method = c("nccpi", "yield", "percentile")) {
  # Classify fields into productivity zones WITHIN each geographic region
  # This ensures balanced representation across regions
  
  # Returns: DT with new column 'prod_zone' (high, medium, low) assigned within geo_region
  method <- match.arg(method)
  
  if (!"geo_region" %in% names(DT)) {
    stop("Must run classify_geographic_region() first")
  }
  
  # Determine yield column
  y_col <- if ("corn_yield" %in% names(DT)) "corn_yield" else "soy_yield"
  
  if (method == "nccpi") {
    if ("nccpi" %in% names(DT)) {
      # Classify within each region using NCCPI
      DT[, prod_zone := fcase(
        nccpi >= 0.70, "high",
        nccpi >= 0.50, "medium",
        nccpi < 0.50, "low"
      ), by = geo_region]
    } else {
      warning("NCCPI not found, falling back to percentile method")
      method <- "percentile"
    }
  }
  
  if (method == "yield") {
    # Calculate field-level average yields
    field_avg <- DT[!is.na(get(y_col)), 
                    .(avg_yield = mean(get(y_col), na.rm = TRUE)), 
                    by = .(tile, field_id, geo_region)]
    
    # Define thresholds by crop type
    if (y_col == "corn_yield") {
      field_avg[, prod_zone := fcase(
        avg_yield >= 180, "high",
        avg_yield >= 150, "medium",
        avg_yield < 150, "low"
      ), by = geo_region]
    } else {
      field_avg[, prod_zone := fcase(
        avg_yield >= 55, "high",
        avg_yield >= 45, "medium",
        avg_yield < 45, "low"
      ), by = geo_region]
    }
    
    # Remove prod_zone if it already exists to avoid conflicts
    if ("prod_zone" %in% names(DT)) {
      DT[, prod_zone := NULL]
    }
    
    DT <- merge(DT, field_avg[, .(tile, field_id, prod_zone)], 
                by = c("tile", "field_id"), all.x = TRUE)
  }
  
  if (method == "percentile") {
    # Calculate field-level average yields
    field_avg <- DT[!is.na(get(y_col)), 
                    .(avg_yield = mean(get(y_col), na.rm = TRUE)), 
                    by = .(tile, field_id, geo_region)]
    
    # Calculate percentiles WITHIN each region
    field_avg[, prod_zone := {
      p33 <- quantile(avg_yield, 0.33, na.rm = TRUE)
      p67 <- quantile(avg_yield, 0.67, na.rm = TRUE)
      
      fcase(
        avg_yield >= p67, "high",
        avg_yield >= p33, "medium",
        avg_yield < p33, "low"
      )
    }, by = geo_region]
    
    # Print thresholds by region
    thresholds <- field_avg[, .(
      p33 = quantile(avg_yield, 0.33, na.rm = TRUE),
      p67 = quantile(avg_yield, 0.67, na.rm = TRUE)
    ), by = geo_region]
    
    cat("\nProductivity zone thresholds by region (percentile method):\n")
    for (i in 1:nrow(thresholds)) {
      cat(sprintf("\n%s:\n", toupper(thresholds$geo_region[i])))
      cat(sprintf("  High:   >= %.1f bu/ac\n", thresholds$p67[i]))
      cat(sprintf("  Medium: %.1f - %.1f bu/ac\n", thresholds$p33[i], thresholds$p67[i]))
      cat(sprintf("  Low:    < %.1f bu/ac\n", thresholds$p33[i]))
    }
    
    # Remove prod_zone if it already exists to avoid conflicts
    if ("prod_zone" %in% names(DT)) {
      DT[, prod_zone := NULL]
    }
    
    DT <- merge(DT, field_avg[, .(tile, field_id, prod_zone)], 
                by = c("tile", "field_id"), all.x = TRUE)
  }
  
  # Report distribution by region × productivity
  zone_dist <- DT[!is.na(geo_region) & !is.na(prod_zone), .N, by = .(geo_region, prod_zone)]
  cat("\nRegion × Productivity zone distribution:\n")
  print(dcast(zone_dist, geo_region ~ prod_zone, value.var = "N"))
  
  # Verify prod_zone column exists and check for missing values
  if (!"prod_zone" %in% names(DT)) {
    stop("ERROR: prod_zone column was not created successfully!")
  }
  
  n_missing <- sum(is.na(DT$prod_zone))
  if (n_missing > 0) {
    warning(sprintf("%d observations have missing prod_zone values", n_missing))
  }
  
  DT
}


# ============================================================
# 0) READ + PREP (memory-aware)
# ============================================================
prep_crop_df_dt <- function(raw, crop = c("corn","soy")) {
  crop <- match.arg(crop)
  DT <- as.data.table(raw)
  
  y_col  <- if (crop == "corn") "corn_yield" else "soy_yield"
  ys_col <- if (crop == "corn") "corn_yield_scym" else "soy_yield_scym"
  
  crop_hist_cols <- names(DT)[stringr::str_detect(names(DT), "^crop_\\d{4}$")]
  
  core_keep <- intersect(
    c("X","tile","field_id","STATE_ABBR","COUNTY_FIPS","year","ccyear",
      y_col, ys_col, "CC_probability", "region", "nccpi3all_mean"),  # Added nccpi
    names(DT)
  )
  keep_cols <- unique(c(core_keep, crop_hist_cols))
  DT <- DT[, ..keep_cols]
  
  ord_cols <- intersect(c("STATE_ABBR","COUNTY_FIPS","year","tile","field_id"), names(DT))
  if (length(ord_cols) > 0) setorderv(DT, ord_cols)
  
  conv <- if (crop == "corn") KGHA_TO_BUAC_CORN else KGHA_TO_BUAC_SOY
  if (y_col %in% names(DT))  DT[, (y_col)  := get(y_col)  * conv]
  if (ys_col %in% names(DT)) DT[, (ys_col) := get(ys_col) * conv]
  
  DT
}

cat("\n--- Reading corn parquet ---\n")
corn_raw <- arrow::read_parquet(corn_path)
corn_raw |>
  filter(STATE_ABBR == "IL") |>
  mutate(fips = paste0("17", COUNTY_FIPS),
         fips = as.numeric(fips)) |>
  left_join(regions, by = "fips") ->
  corn_raw
corn_arranged <- prep_crop_df_dt(corn_raw, "corn")
rm(corn_raw); gc()

cat("\n--- Reading soy parquet ---\n")
soy_raw <- arrow::read_parquet(soy_path)
soy_raw |>
  filter(STATE_ABBR == "IL") |>
  mutate(fips = paste0("17", COUNTY_FIPS),
         fips = as.numeric(fips)) |>
  left_join(regions, by = "fips") ->
  soy_raw
soy_arranged <- prep_crop_df_dt(soy_raw, "soy")
rm(soy_raw); gc()

# ============================================================
# CLASSIFY BY REGION AND PRODUCTIVITY ZONE
# ============================================================
cat("\n=== CLASSIFYING FIELDS ===\n")

cat("\n--- Classifying corn fields by geographic region ---\n")
corn_arranged <- classify_geographic_region(corn_arranged)

cat("\n--- Classifying corn fields by productivity zone (within region) ---\n")
corn_arranged <- classify_productivity_zone_by_region(corn_arranged, method = "percentile")

cat("\n--- Classifying soy fields by geographic region ---\n")
soy_arranged <- classify_geographic_region(soy_arranged)

cat("\n--- Classifying soy fields by productivity zone (within region) ---\n")
soy_arranged <- classify_productivity_zone_by_region(soy_arranged, method = "percentile")

stopifnot(all(c("tile","field_id","year","geo_region","prod_zone") %in% names(corn_arranged)))
stopifnot(all(c("tile","field_id","year","geo_region","prod_zone") %in% names(soy_arranged)))

corn_arranged |>
  filter(!is.na(geo_region) & !is.na(prod_zone)) |>
  group_by(geo_region, prod_zone) |>
  summarise(corn_m = mean(corn_yield, na.rm = TRUE),
            corn_sd = sd(corn_yield, na.rm = TRUE))

soy_arranged |>
  filter(!is.na(geo_region) & !is.na(prod_zone)) |>
  group_by(geo_region, prod_zone) |>
  summarise(corn_m = mean(soy_yield, na.rm = TRUE),
            corn_sd = sd(soy_yield, na.rm = TRUE))

# ============================================================
# 1) FULL 6-YEAR HISTORY STATE BUILDER
# ============================================================
make_rot6_full <- function(DT) {
  stopifnot(is.data.table(DT))
  crop_cols <- sort(names(DT)[stringr::str_detect(names(DT), "^crop_\\d{4}$")])
  if (length(crop_cols) == 0) stop("No crop_YYYY columns found.")
  
  crop_years <- as.integer(stringr::str_remove(crop_cols, "^crop_"))
  pos_y <- match(DT[["year"]], crop_years)
  
  vals <- unlist(DT[, ..crop_cols], use.names = FALSE)
  crop_bin <- ifelse(vals == "Corn", 0L,
                     ifelse(vals == "Soybeans", 1L, NA_integer_))
  crop_bin <- matrix(crop_bin, nrow = nrow(DT), ncol = length(crop_cols))
  rm(vals); gc()
  
  idx_mat <- sapply(-5:0, function(o) pos_y + o)
  
  n <- nrow(DT)
  seq6 <- matrix(NA_integer_, n, 6)
  
  for (j in 1:6) {
    cj <- idx_mat[, j]
    good <- !is.na(cj) & cj >= 1 & cj <= ncol(crop_bin)
    if (any(good)) seq6[good, j] <- crop_bin[cbind(which(good), cj[good])]
  }
  
  incomplete <- apply(is.na(seq6), 1, any)
  
  out <- rep(NA_character_, n)
  if (any(!incomplete)) {
    chars <- ifelse(seq6[!incomplete, , drop = FALSE] == 1L, "S", "C")
    out[!incomplete] <- apply(chars, 1, paste0, collapse = "")
  }
  out
}

# ============================================================
# REGION-ZONE-SPECIFIC YIELD GENERATOR BUILDER
# ============================================================
build_yield_generator_region_zone <- function(DT,
                                              yield_col,
                                              terminal_crop = c("C","S"),
                                              geo_region = c("northern", "central", "southern"),
                                              prod_zone = c("high", "medium", "low"),
                                              seed = 1L,
                                              B_per_state = 20000L) {
  # Build yield generator for a specific region × productivity zone combination
  
  terminal_crop <- match.arg(terminal_crop)
  geo_region <- match.arg(geo_region)
  prod_zone <- match.arg(prod_zone)
  
  stopifnot(is.data.table(DT))
  stopifnot(all(c("tile","field_id","year","geo_region","prod_zone") %in% names(DT)))
  stopifnot(yield_col %in% names(DT))
  
  # Filter to specific region × zone.
  # NOTE: the argument names geo_region / prod_zone shadow the identically named
  # columns of DT. Under dplyr data masking the column wins on both sides of `==`,
  # so filter(geo_region == geo_region) is TRUE for every row and the cell filter
  # silently does nothing. Copy the arguments into differently named locals first.
  .reg      <- geo_region
  .zone     <- prod_zone
  .n_before <- nrow(DT)
  DT <- DT[geo_region == .reg & prod_zone == .zone]
  stopifnot(nrow(DT) < .n_before)   # guard: the cell filter must actually bind
  
  if (nrow(DT) == 0) {
    stop("No observations for ", geo_region, "-", prod_zone)
  }
  
  cat(sprintf("  Building for %s-%s: %d observations\n", geo_region, prod_zone, nrow(DT)))
  
  DT[, field_uid := paste(tile, field_id, sep=":")]
  DT[, rot6 := make_rot6_full(DT)]
  DT_sub <- DT[!is.na(rot6) & substr(rot6, 6, 6) == terminal_crop]
  
  DT_fy <- DT_sub[!is.na(get(yield_col)),
                  .(y = mean(get(yield_col), na.rm = TRUE),
                    n_tiles = .N),
                  by = .(field_uid, year, rot6)]
  
  if (nrow(DT_fy) == 0) {
    stop("0 rows after filtering for ", geo_region, "-", prod_zone)
  }
  
  # FE model within this region-zone
  m <- fixest::feols(y ~ 1 | field_uid + year, data = DT_fy)
  
  fe <- fixest::fixef(m)
  fe_field <- unname(fe$field_uid[as.character(DT_fy$field_uid)])
  fe_year  <- unname(fe$year[as.character(DT_fy$year)])
  
  DT_fy[, resid := y - (fe_field + fe_year)]
  DT_fy[!is.finite(resid), resid := NA_real_]
  
  DT_use <- DT_fy[!is.na(resid)]
  if (nrow(DT_use) == 0) stop("All residuals are NA")
  
  set.seed(seed)
  
  pools <- DT_use[, .(
    n = .N,
    resid_pool = list({
      x <- resid
      x <- x[is.finite(x)]
      if (length(x) == 0) numeric(0)
      else if (length(x) >= B_per_state) sample(x, B_per_state, replace = FALSE)
      else sample(x, B_per_state, replace = TRUE)
    })
  ), by = rot6]
  
  # State mean = mean field effect of the state's observations + ONE common year effect
  # (cell-wide, observation-weighted). Averaging each state's own fitted year effects would
  # let a thin state whose observations sit in drought years (2011-12) inherit a spuriously
  # low mean; the year effects are common shocks, not part of the rotation state.
  DT_use[, fe_f := y - resid - fe_year[!is.na(DT_fy$resid)]]
  ybar_year <- mean(fe_year[!is.na(DT_fy$resid)])
  mu <- DT_use[, .(mu = mean(fe_f, na.rm = TRUE) + ybar_year), by = rot6]
  
  simulate <- function(state, B = 5000L) {
    mu_s <- mu[rot6 == state, mu]
    if (length(mu_s) != 1 || is.na(mu_s)) stop("State not found: ", state)
    pool <- pools[rot6 == state, resid_pool][[1]]
    if (length(pool) == 0) stop("Empty residual pool for: ", state)
    mu_s + sample(pool, B, replace = TRUE)
  }
  
  list(
    yield_col = yield_col,
    terminal_crop = terminal_crop,
    geo_region = geo_region,
    prod_zone = prod_zone,
    DT_fy = DT_fy,
    model = m,
    pools = pools,
    mu = mu,
    simulate = simulate
  )
}

# ============================================================
# BUILD GENERATORS FOR ALL REGION × ZONE COMBINATIONS
# ============================================================

cat("\n=== BUILDING REGION × ZONE YIELD GENERATORS ===\n")
corn_arranged |>
  filter(!is.na(prod_zone)) ->
  corn_arranged

soy_arranged |>
  filter(!is.na(prod_zone)) ->
  soy_arranged

# Build corn generators
cat("\n--- CORN generators ---\n")
corn_generators <- list()
for (reg in c("northern", "central", "southern")) {
  corn_generators[[reg]] <- list()
  for (zone in c("high", "medium", "low")) {
    key <- paste0(reg, "_", zone)
    cat(sprintf("\nBuilding corn generator: %s\n", key))
    tryCatch({
      corn_generators[[reg]][[zone]] <- build_yield_generator_region_zone(
        DT = copy(corn_arranged),
        yield_col = "corn_yield",
        terminal_crop = "C",
        geo_region = reg,
        prod_zone = zone,
        seed = 1L,
        B_per_state = 20000L
      )
      cat(sprintf("  ✓ Success: %d states available\n", 
                  length(corn_generators[[reg]][[zone]]$mu$rot6)))
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR: %s\n", e$message))
      corn_generators[[reg]][[zone]] <- NULL
    })
  }
}

# Build soy generators
cat("\n--- SOY generators ---\n")
soy_generators <- list()
for (reg in c("northern", "central", "southern")) {
  soy_generators[[reg]] <- list()
  for (zone in c("high", "medium", "low")) {
    key <- paste0(reg, "_", zone)
    cat(sprintf("\nBuilding soy generator: %s\n", key))
    tryCatch({
      soy_generators[[reg]][[zone]] <- build_yield_generator_region_zone(
        DT = copy(soy_arranged),
        yield_col = "soy_yield",
        terminal_crop = "S",
        geo_region = reg,
        prod_zone = zone,
        seed = 2L,
        B_per_state = 20000L
      )
      cat(sprintf("  ✓ Success: %d states available\n", 
                  length(soy_generators[[reg]][[zone]]$mu$rot6)))
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR: %s\n", e$message))
      soy_generators[[reg]][[zone]] <- NULL
    })
  }
}

yield_gens_by_region_zone <- list(
  corn = corn_generators,
  soy = soy_generators
)
gc()

# ============================================================
# 4) **NEW: ZONE-SPECIFIC COST PARAMETERS**
# ============================================================

# NOTE: price parameters below are hard-coded (regional mean / log-sd; sigma from
# calibration/sigma_estimation.r). The original script also read
# county_prices/crop_prices_il.csv and data-raw/direct_cost.csv here, but the
# result of that block was discarded and never entered the model, so it was removed.

# Region-specific price parameters
# Based on provided data for Illinois regions
price_params_regional <- list(
  northern = list(
    C = list(Pbar = 3.72, sigma = 0.233),   # Northern IL Corn
    S = list(Pbar = 9.57, sigma = 0.172)    # Northern IL Soybeans
  ),
  central = list(
    C = list(Pbar = 3.80, sigma = 0.225),  # Central IL Corn
    S = list(Pbar = 9.62, sigma = 0.173)    # Central IL Soybeans
  ),
  southern = list(
    C = list(Pbar = 3.71, sigma = 0.204),   # Southern IL Corn
    S = list(Pbar = 9.70, sigma = 0.168)    # Southern IL Soybeans
  )
)

# Calculate log-normal parameters for each region
for (reg in names(price_params_regional)) {
  price_params_regional[[reg]]$C$mu <- 
    log(price_params_regional[[reg]]$C$Pbar) - 
    0.5 * price_params_regional[[reg]]$C$sigma^2
  
  price_params_regional[[reg]]$S$mu <- 
    log(price_params_regional[[reg]]$S$Pbar) - 
    0.5 * price_params_regional[[reg]]$S$sigma^2
}

# Print price structure
cat("\n=== REGION-SPECIFIC PRICE PARAMETERS ===\n")
for (reg in names(price_params_regional)) {
  cat("\n", toupper(reg), " ILLINOIS:\n", sep="")
  cat(sprintf("  Corn:     mean=$%.2f/bu, sd=$%.4f\n",
              price_params_regional[[reg]]$C$Pbar,
              price_params_regional[[reg]]$C$sigma))
  cat(sprintf("  Soybeans: mean=$%.2f/bu, sd=$%.4f\n",
              price_params_regional[[reg]]$S$Pbar,
              price_params_regional[[reg]]$S$sigma))
}

# Default to central for backward compatibility
price_params <- price_params_regional$central

#  ============================================================
# REGION × PRODUCTIVITY ZONE COST STRUCTURE
# Based on empirical direct cost data
# ============================================================

# Empirical base costs by region × productivity × crop × previous crop
# This replaces the generic cost structure with actual data

# Create cost lookup table from empirical data
empirical_costs <- data.table(
  region = c(
    rep("northern", 8),
    rep("central", 8),
    rep("southern", 8)
  ),
  productivity = rep(c(
    rep("high", 4), rep("low", 4)
  ), 3),
  crop = rep(c("corn", "corn", "soy", "soy"), 6),
  prev_crop = rep(c("soy", "corn", "corn", "soy"), 6),
  direct_cost = c(
    # Northern high
    468, 483, 207, 212,
    # Northern low
    468, 483, 207, 212,
    # Central high
    486, 501, 225, 230,
    # Central low
    462, 478, 210, 215,
    # Southern high
    452, 468, 221, 226,
    # Southern low
    452, 468, 221, 226
  )
)

# Display the cost structure
cat("\n=== EMPIRICAL COST STRUCTURE ===\n")
cat("\nCosts by Region × Productivity × Rotation:\n")
print(dcast(empirical_costs, region + productivity ~ crop + prev_crop, value.var = "direct_cost"))

# Calculate rotation cost effects from empirical data
# Effect = (after corn - after soy) / after soy
cat("\n=== ROTATION COST EFFECTS (from empirical data) ===\n")
rotation_effects <- empirical_costs[, .(
  region = region[1],
  productivity = productivity[1],
  crop = crop[1],
  cost_after_soy = direct_cost[prev_crop == "soy"],
  cost_after_corn = direct_cost[prev_crop == "corn"],
  cost_increase = direct_cost[prev_crop == "corn"] - direct_cost[prev_crop == "soy"],
  pct_increase = (direct_cost[prev_crop == "corn"] - direct_cost[prev_crop == "soy"]) / 
    direct_cost[prev_crop == "soy"] * 100
), by = .(region, productivity, crop)]

print(rotation_effects)

# Convert to cost parameter structure for simulation
# Use medium productivity as interpolation between high and low
base_cost_params <- list()

for (reg in c("northern", "central", "southern")) {
  base_cost_params[[reg]] <- list()
  
  for (zone in c("high", "medium", "low")) {
    
    # For medium, average high and low
    if (zone == "medium") {
      # Corn after soy (baseline)
      corn_high <- empirical_costs[region == reg & productivity == "high" & 
                                     crop == "corn" & prev_crop == "soy", direct_cost]
      corn_low <- empirical_costs[region == reg & productivity == "low" & 
                                    crop == "corn" & prev_crop == "soy", direct_cost]
      corn_base <- mean(c(corn_high, corn_low))
      
      # Soy after corn (baseline)
      soy_high <- empirical_costs[region == reg & productivity == "high" & 
                                    crop == "soy" & prev_crop == "corn", direct_cost]
      soy_low <- empirical_costs[region == reg & productivity == "low" & 
                                   crop == "soy" & prev_crop == "corn", direct_cost]
      soy_base <- mean(c(soy_high, soy_low))
    } else {
      # Use actual high or low values
      prod_level <- zone
      
      # Corn after soy (this will be our baseline, then adjusted by rotation multipliers)
      corn_base <- empirical_costs[region == reg & productivity == prod_level & 
                                     crop == "corn" & prev_crop == "soy", direct_cost]
      
      # Soy after corn (baseline)
      soy_base <- empirical_costs[region == reg & productivity == prod_level & 
                                    crop == "soy" & prev_crop == "corn", direct_cost]
    }
    
    # Store as cost parameters with assumed volatility
    base_cost_params[[reg]][[zone]] <- list(
      C = list(Cbar = corn_base, sigma = 0.10),  # 10% CV assumed
      S = list(Cbar = soy_base, sigma = 0.10)
    )
  }
}

# Calculate log-normal parameters
for (reg in names(base_cost_params)) {
  for (zone in names(base_cost_params[[reg]])) {
    base_cost_params[[reg]][[zone]]$C$mu <- 
      log(base_cost_params[[reg]][[zone]]$C$Cbar) - 
      0.5 * base_cost_params[[reg]][[zone]]$C$sigma^2
    
    base_cost_params[[reg]][[zone]]$S$mu <- 
      log(base_cost_params[[reg]][[zone]]$S$Cbar) - 
      0.5 * base_cost_params[[reg]][[zone]]$S$sigma^2
  }
}

# Print final cost structure
cat("\n=== BASE COST PARAMETERS (Corn after Soy, Soy after Corn) ===\n")
for (reg in names(base_cost_params)) {
  cat("\n", toupper(reg), " ILLINOIS:\n", sep="")
  for (zone in names(base_cost_params[[reg]])) {
    cat(sprintf("  %s productivity: Corn=$%.0f/ac, Soy=$%.0f/ac\n",
                zone,
                base_cost_params[[reg]][[zone]]$C$Cbar,
                base_cost_params[[reg]][[zone]]$S$Cbar))
  }
}

# ============================================================
# ROTATION COST ADJUSTMENTS (calculated from empirical data)
# ============================================================

# Calculate empirical rotation multipliers
# These are derived from the actual cost differences in the data

# For corn: after corn vs after soy
corn_rotation_effects <- empirical_costs[crop == "corn", .(
  mean_mult = mean(direct_cost[prev_crop == "corn"] / direct_cost[prev_crop == "soy"])
), by = .(region, productivity)]

cat("\n=== EMPIRICAL CORN ROTATION MULTIPLIERS ===\n")
print(corn_rotation_effects)

# For soy: after soy vs after corn  
soy_rotation_effects <- empirical_costs[crop == "soy", .(
  mean_mult = mean(direct_cost[prev_crop == "soy"] / direct_cost[prev_crop == "corn"])
), by = .(region, productivity)]

cat("\n=== EMPIRICAL SOY ROTATION MULTIPLIERS ===\n")
print(soy_rotation_effects)

# Calculate average multipliers across regions/productivity
avg_corn_mult <- empirical_costs[crop == "corn", 
                                 mean(direct_cost[prev_crop == "corn"] / direct_cost[prev_crop == "soy"])]
avg_soy_mult <- empirical_costs[crop == "soy", 
                                mean(direct_cost[prev_crop == "soy"] / direct_cost[prev_crop == "corn"])]

cat("\n=== AVERAGE ROTATION MULTIPLIERS ===\n")
cat(sprintf("Corn after corn vs after soy: %.4f (%.1f%% increase)\n", 
            avg_corn_mult, (avg_corn_mult - 1) * 100))
cat(sprintf("Soy after soy vs after corn: %.4f (%.1f%% increase)\n", 
            avg_soy_mult, (avg_soy_mult - 1) * 100))

# Set up rotation cost adjustments based on empirical data
# Corn after corn = ~3.2% increase
# Soy after soy = ~2.4% increase
# For consecutive years, we'll scale these effects

rotation_cost_adjustments <- list(
  corn = list(
    consecutive_1yr = list(
      pattern = "C$", 
      multiplier = avg_corn_mult,  # ~1.032 from data
      description = "1 year consecutive corn (empirical)"
    ),
    consecutive_2yr = list(
      pattern = "CC$", 
      multiplier = 1 + 2 * (avg_corn_mult - 1),  # ~1.064 (doubled effect)
      description = "2 years consecutive corn (scaled empirical)"
    ),
    consecutive_3plus = list(
      pattern = "CCC", 
      multiplier = 1 + 3 * (avg_corn_mult - 1),  # ~1.096 (tripled effect)
      description = "3+ years consecutive corn (scaled empirical)"
    ),
    after_soy = list(
      pattern = "S$", 
      multiplier = 1.00,  # Baseline from data
      description = "Corn after soy (baseline)"
    )
  ),
  soy = list(
    consecutive_1yr = list(
      pattern = "S$", 
      multiplier = avg_soy_mult,  # ~1.024 from data
      description = "1 year consecutive soy (empirical)"
    ),
    consecutive_2plus = list(
      pattern = "SS", 
      multiplier = 1 + 2 * (avg_soy_mult - 1),  # ~1.048 (doubled effect)
      description = "2+ years consecutive soy (scaled empirical)"
    ),
    after_corn = list(
      pattern = "C$", 
      multiplier = 1.00,  # Baseline from data
      description = "Soy after corn (baseline)"
    )
  )
)

cat("\n=== ROTATION COST ADJUSTMENT RULES (from empirical data) ===\n")
cat("\nCORN:\n")
for (rule_name in names(rotation_cost_adjustments$corn)) {
  rule <- rotation_cost_adjustments$corn[[rule_name]]
  cat(sprintf("  %s: %.3f (%.1f%% change) - %s\n", 
              rule_name, 
              rule$multiplier,
              (rule$multiplier - 1) * 100,
              rule$description))
}

cat("\nSOY:\n")
for (rule_name in names(rotation_cost_adjustments$soy)) {
  rule <- rotation_cost_adjustments$soy[[rule_name]]
  cat(sprintf("  %s: %.3f (%.1f%% change) - %s\n", 
              rule_name, 
              rule$multiplier,
              (rule$multiplier - 1) * 100,
              rule$description))
}

get_rotation_cost_multiplier <- function(window6, current_crop) {
  if (is.null(window6) || is.na(window6) || nchar(window6) != 6) {
    return(1.00)
  }
  if (is.null(current_crop) || is.na(current_crop) || !current_crop %in% c("C", "S")) {
    return(1.00)
  }
  
  history5 <- substr(window6, 1, 5)
  crop_key <- ifelse(current_crop == "C", "corn", "soy")
  
  if (!crop_key %in% names(rotation_cost_adjustments)) return(1.00)
  
  rules <- rotation_cost_adjustments[[crop_key]]
  
  if (current_crop == "C") {
    if (grepl("CCC", history5)) {
      mult <- rules$consecutive_3plus$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
    if (grepl("CC$", history5)) {
      mult <- rules$consecutive_2yr$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
    if (grepl("C$", history5)) {
      mult <- rules$consecutive_1yr$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
    if (grepl("S$", history5)) {
      mult <- rules$after_soy$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
  } else if (current_crop == "S") {
    if (grepl("SS", history5)) {
      mult <- rules$consecutive_2plus$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
    if (grepl("S$", history5)) {
      mult <- rules$consecutive_1yr$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
    if (grepl("C$", history5)) {
      mult <- rules$after_corn$multiplier
      if (!is.null(mult) && !is.na(mult)) return(mult)
    }
  }
  
  return(1.00)
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================

disc_vec <- function(r) 1 / (1 + r)^(0:5)

draw_lognorm_6y <- function(B, mu_log, sd_log) {
  matrix(rlnorm(B * 6, meanlog = mu_log, sdlog = sd_log), nrow = B, ncol = 6)
}

summarize_dist <- function(x) {
  tibble::tibble(
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    p05  = as.numeric(quantile(x, 0.05, na.rm = TRUE)),
    p50  = as.numeric(quantile(x, 0.50, na.rm = TRUE)),
    p95  = as.numeric(quantile(x, 0.95, na.rm = TRUE)),
    pr_positive = mean(x > 0, na.rm = TRUE)
  )
}

make_transition_path <- function(history6, plan6) {
  stopifnot(nchar(history6) == 6, nchar(plan6) == 6)
  h <- strsplit(history6, "")[[1]]
  p <- strsplit(plan6, "")[[1]]
  
  out <- data.table(t = 1:6, crop = p, window6 = NA_character_)
  for (t in 1:6) {
    seq_now <- c(h, p[1:t])
    win <- tail(seq_now, 6)
    out[t, window6 := paste0(win, collapse="")]
  }
  out
}

# ============================================================
# REGION-ZONE SIMULATION FUNCTION
# ============================================================

simulate_plan_pv_npv_region_zone <- function(history6,
                                             plan6,
                                             B = 20000L,
                                             r_disc = 0.05,
                                             gens_by_region_zone = yield_gens_by_region_zone,
                                             price_pars_regional = price_params_regional,
                                             base_cost_pars = base_cost_params,
                                             geo_region = "central",
                                             prod_zone = "medium",
                                             use_rotation_costs = TRUE,
                                             seed = 123,
                                             verbose = TRUE) {
  stopifnot(nchar(history6) == 6, nchar(plan6) == 6)
  set.seed(seed)
  
  # Select region-specific prices
  if (geo_region %in% names(price_pars_regional)) {
    price_pars <- price_pars_regional[[geo_region]]
    if (verbose) cat("Using prices for", geo_region, "Illinois\n")
  } else {
    stop("Region not found in price parameters: ", geo_region)
  }
  
  # Select region-zone specific costs
  if (geo_region %in% names(base_cost_pars) && 
      prod_zone %in% names(base_cost_pars[[geo_region]])) {
    cost_pars_base <- base_cost_pars[[geo_region]][[prod_zone]]
    if (verbose) cat("Using costs for", geo_region, "-", prod_zone, "\n")
  } else {
    stop("Region-zone combination not found: ", geo_region, "-", prod_zone)
  }
  
  # Select region-zone specific yield generators
  corn_gen <- gens_by_region_zone$corn[[geo_region]][[prod_zone]]
  soy_gen <- gens_by_region_zone$soy[[geo_region]][[prod_zone]]
  
  if (is.null(corn_gen) || is.null(soy_gen)) {
    stop("Yield generators not available for: ", geo_region, "-", prod_zone)
  }
  
  path <- make_transition_path(history6, plan6)
  if (verbose) {
    cat("\n--- Transition path ---\n")
    cat("history6:     ", history6, "\n", sep="")
    cat("plan:         ", plan6, "\n", sep="")
    cat("geo_region:   ", geo_region, "\n", sep="")
    cat("prod_zone:    ", prod_zone, "\n", sep="")
    print(path)
  }
  
  disc <- disc_vec(r_disc)
  
  # Draw base prices (region-specific) and costs
  Pc_base <- draw_lognorm_6y(B, price_pars$C$mu, price_pars$C$sigma)
  Ps_base <- draw_lognorm_6y(B, price_pars$S$mu, price_pars$S$sigma)
  Cc_base <- draw_lognorm_6y(B, cost_pars_base$C$mu, cost_pars_base$C$sigma)
  Cs_base <- draw_lognorm_6y(B, cost_pars_base$S$mu, cost_pars_base$S$sigma)
  
  # Apply rotation cost adjustments
  Cc_adjusted <- Cc_base
  Cs_adjusted <- Cs_base
  
  if (use_rotation_costs) {
    for (t in 1:6) {
      crop_t <- path$crop[t]
      window_t <- path$window6[t]
      
      multiplier <- get_rotation_cost_multiplier(window_t, crop_t)
      
      if (is.null(multiplier) || is.na(multiplier) || length(multiplier) == 0) {
        multiplier <- 1.00
      }
      
      multiplier <- as.numeric(multiplier[1])
      
      if (crop_t == "C") {
        Cc_adjusted[, t] <- Cc_base[, t] * multiplier
      } else if (crop_t == "S") {
        Cs_adjusted[, t] <- Cs_base[, t] * multiplier
      }
      
      if (verbose) {
        cat(sprintf("  Year %d (%s): window=%s, multiplier=%.3f\n", 
                    t, crop_t, window_t, multiplier))
      }
    }
  }
  
  # Yield draws from region-zone specific generators
  Ymat <- matrix(NA_real_, nrow = B, ncol = 6)
  
  for (t in 1:6) {
    crop_t <- path$crop[t]
    st_t   <- path$window6[t]
    
    if (crop_t == "C") {
      Ymat[, t] <- corn_gen$simulate(st_t, B = B)
    } else if (crop_t == "S") {
      Ymat[, t] <- soy_gen$simulate(st_t, B = B)
    } else stop("Unknown crop: ", crop_t)
  }
  
  # Revenue and costs
  chars <- strsplit(plan6, "")[[1]]
  isC <- as.numeric(chars == "C")
  isS <- as.numeric(chars == "S")
  IC <- matrix(isC, nrow = B, ncol = 6, byrow = TRUE)
  IS <- matrix(isS, nrow = B, ncol = 6, byrow = TRUE)
  
  Pmat <- (Pc_base * IC) + (Ps_base * IS)
  Kmat <- (Cc_adjusted * IC) + (Cs_adjusted * IS)
  
  rev_mat   <- Pmat * Ymat
  cost_mat  <- Kmat
  prof_mat  <- rev_mat - cost_mat
  
  pv_rev <- as.numeric(rev_mat  %*% disc)
  pv_cost<- as.numeric(cost_mat %*% disc)
  npv    <- as.numeric(prof_mat %*% disc)
  
  list(
    history6 = history6,
    plan6 = plan6,
    geo_region = geo_region,
    prod_zone = prod_zone,
    use_rotation_costs = use_rotation_costs,
    path = path,
    pv_revenue = pv_rev,
    pv_cost = pv_cost,
    npv_profit = npv
  )
}

rm(corn_arranged, soy_arranged)
gc()