# ============================================================================
# Map of rotation practices (paper Figure 2), field level, 2016.
# Rows: corn / soybeans as the terminal crop. Columns: continuous crop, strict alternation, 2-corn/1-soy.
# Each point is a field whose six-year history (2011-2016) matches the panel's practice.
# Projected to NAD83 / Illinois East (EPSG:26971, meters) so the state is not drawn slanted.
# Needs the raw panels (SIM_CORN_PARQUET, SIM_SOY_PARQUET); writes figs/map_practice.png.
#   setwd("<path>/simulation_code"); source("maps_practice.R")
# ============================================================================
suppressMessages({ library(data.table); library(arrow); library(sf); library(usmap); library(ggplot2) })
source("config.R")
dir.create(FIGS_DIR <- Sys.getenv("SIM_FIGS_DIR", "figs"), showWarnings = FALSE)
CRS_MAP <- 26971

crop_cols <- paste0("crop_", 2011:2016)
field_points <- function(path, terminal, patterns) {
  x <- as.data.table(read_parquet(path, col_select = dplyr::any_of(c("STATE_ABBR", "year", "tile", "field_id", "lat", "lon", crop_cols))))
  x <- x[STATE_ABBR == "IL" & year == 2016 & !is.na(lat) & !is.na(lon)]
  x <- unique(x, by = c("tile", "field_id"))
  code <- lapply(crop_cols, function(cc) fifelse(x[[cc]] == "Corn", "C", fifelse(x[[cc]] == "Soybeans", "S", NA_character_)))
  x[, hist := do.call(paste0, code)]
  x <- x[!Reduce(`|`, lapply(code, is.na))]
  n_complete <- nrow(x)
  x[, practice := NA_character_]
  for (p in names(patterns)) x[hist %in% patterns[[p]], practice := p]
  cat(sprintf("%s: %d fields with a complete history; %s\n", terminal, n_complete,
              paste(sprintf("%s %.1f%%", names(patterns), 100 * sapply(patterns, function(pt) mean(x$hist %in% pt))), collapse = ", ")))
  x <- x[!is.na(practice), .(lon, lat, practice)]
  x[, terminal := terminal]
  x
}
corn <- field_points(CORN_PARQUET, "Corn as terminal crop",
  list("Continuous crop" = "CCCCCC", "Strict alternation" = "SCSCSC", "2 corn - 1 soy" = c("SCCSCC", "CSCCSC")))
soy <- field_points(SOY_PARQUET, "Soybeans as terminal crop",
  list("Continuous crop" = "SSSSSS", "Strict alternation" = "CSCSCS", "2 corn - 1 soy" = "CCSCCS"))
pts <- rbind(corn, soy)
pts[, practice := factor(practice, c("Continuous crop", "Strict alternation", "2 corn - 1 soy"))]
pts[, terminal := factor(terminal, c("Corn as terminal crop", "Soybeans as terminal crop"))]
pts <- st_transform(st_as_sf(pts, coords = c("lon", "lat"), crs = 4326), CRS_MAP)

cty <- us_map(regions = "counties")
il <- st_transform(st_as_sf(cty[cty$abbr == "IL", ]), CRS_MAP)

p <- ggplot() +
  geom_sf(data = il, fill = "white", color = "grey70", linewidth = 0.1) +
  geom_sf(data = pts, aes(color = practice), size = 0.45, alpha = 0.45, stroke = 0, show.legend = FALSE) +
  facet_grid(terminal ~ practice) +
  scale_color_manual(values = c("Continuous crop" = "#D73027", "Strict alternation" = "#F4A340", "2 corn - 1 soy" = "#4575B4")) +
  coord_sf(crs = CRS_MAP, datum = NA) +
  labs(caption = "Each point is a field in 2016 with the six-year history (2011-2016) shown in the column title. Other histories are not shown.") +
  theme_void(base_size = 10) + theme(strip.text = element_text(face = "bold"), plot.caption = element_text(size = 7, hjust = 0),
                                     plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(FIGS_DIR, "map_practice.png"), p, width = 9, height = 7.5, dpi = 200, bg = "white")
