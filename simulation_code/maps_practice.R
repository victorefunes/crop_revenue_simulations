# ============================================================================
# Map of rotation practices (paper Figure 2): county-level share of fields in each practice, 2016.
# Rows: corn / soybeans as the terminal crop. Columns: continuous crop, strict alternation, 2-corn/1-soy.
# Share = fields with that six-year history (2011-2016) / fields with a complete corn-soy history.
# Needs the raw panels (SIM_CORN_PARQUET, SIM_SOY_PARQUET); writes figs/map_practice.png and
# results/map_practice_county_shares.csv.
#   setwd("<path>/simulation_code"); source("maps_practice.R")
# ============================================================================
suppressMessages({ library(data.table); library(arrow); library(sf); library(usmap); library(ggplot2) })
source("config.R")
dir.create(FIGS_DIR <- Sys.getenv("SIM_FIGS_DIR", "figs"), showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
MIN_FIELDS <- 30      # counties with fewer complete-history fields are left blank

crop_cols <- paste0("crop_", 2011:2016)
county_shares <- function(path, terminal, patterns) {
  x <- as.data.table(read_parquet(path, col_select = dplyr::any_of(c("STATE_ABBR", "COUNTY_FIPS", "year", "tile", "field_id", crop_cols))))
  x <- x[STATE_ABBR == "IL" & year == 2016]
  x <- unique(x, by = c("tile", "field_id"))
  code <- lapply(crop_cols, function(cc) fifelse(x[[cc]] == "Corn", "C", fifelse(x[[cc]] == "Soybeans", "S", NA_character_)))
  x[, hist := do.call(paste0, code)][, ok := !Reduce(`|`, lapply(code, is.na))]
  x[, fips := sprintf("17%03d", as.integer(COUNTY_FIPS))]
  x <- x[ok == TRUE]
  rbindlist(lapply(names(patterns), function(p) {
    x[, .(n_fields = .N, n_practice = sum(hist %in% patterns[[p]])), by = fips][, `:=`(practice = p, terminal = terminal, share = 100 * n_practice / n_fields)]
  }))
}
corn <- county_shares(CORN_PARQUET, "Corn as terminal crop",
  list("Continuous crop" = "CCCCCC", "Strict alternation" = "SCSCSC", "2 corn - 1 soy" = c("SCCSCC", "CSCCSC")))
soy <- county_shares(SOY_PARQUET, "Soybeans as terminal crop",
  list("Continuous crop" = "SSSSSS", "Strict alternation" = "CSCSCS", "2 corn - 1 soy" = "CCSCCS"))
sh <- rbind(corn, soy)
sh[n_fields < MIN_FIELDS, share := NA_real_]
fwrite(sh, "results/map_practice_county_shares.csv")

cty <- us_map(regions = "counties")
il <- st_as_sf(cty[cty$abbr == "IL", ])
d <- merge(il[, "fips"], sh, by = "fips")
d$practice <- factor(d$practice, c("Continuous crop", "Strict alternation", "2 corn - 1 soy"))
d$terminal <- factor(d$terminal, c("Corn as terminal crop", "Soybeans as terminal crop"))
p <- ggplot(d) + geom_sf(aes(fill = share), color = "white", linewidth = 0.1) +
  facet_grid(terminal ~ practice) +
  scale_fill_viridis_c(name = "% of fields", option = "mako", direction = -1, na.value = "grey90", limits = c(0, NA)) +
  labs(caption = paste0("Share of fields with a complete 2011-2016 corn-soybean history, by county (2016). ",
                        "Counties with fewer than ", MIN_FIELDS, " such fields are blank. Other histories are not shown.")) +
  theme_void(base_size = 10) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom",
                                     plot.caption = element_text(size = 7, hjust = 0), legend.key.width = grid::unit(1.6, "cm"))
ggsave(file.path(FIGS_DIR, "map_practice.png"), p, width = 10, height = 7.5, dpi = 200)
print(sh[!is.na(share), .(counties = .N, mean_share = round(mean(share), 1), max_share = round(max(share), 1)), by = .(terminal, practice)])
