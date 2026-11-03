# Map of rotation practices in the eight states: share of six-year crop histories ending in corn that are CCCCCC,
# SCSCSC (strict alternation), SCCSCC or CSCCSC (two corn, one soybean), on a 0.25-degree grid, from the corn field panels. Run from simulation_code/
# (SIM_FIGS_DIR=../figs for the manuscript). Writes figs/be_rotation_map.pdf and results/rotation_map_cells.csv.
suppressPackageStartupMessages({ library(data.table); library(arrow); library(ggplot2); library(maps) })
source("config.R")
source("R/state_generator.R")   # make_rot6_full
PANEL_DIR <- .env("SIM_PANEL_DIR", "D:/Crop data")
FIGS_DIR <- .env("SIM_FIGS_DIR", "figs"); dir.create(FIGS_DIR, showWarnings = FALSE, recursive = TRUE)
STATES <- c("IA", "IL", "IN", "MN", "MO", "OH", "SD", "WI")
CELL <- 0.25; MIN_N <- 50

out <- list()
for (st in STATES) {
  p <- file.path(PANEL_DIR, sprintf("d_igis13_%s.parquet", st))
  sch <- names(open_dataset(p)$schema)
  d <- as.data.table(read_parquet(p, col_select = c("year", "lat", "lon", sch[grepl("^crop_\\d{4}$", sch)])))
  d[, rot6 := make_rot6_full(d)]
  d <- d[!is.na(rot6) & substr(rot6, 6, 6) == "C" & !is.na(lat) & !is.na(lon)]
  d[, `:=`(gx = floor(lon / CELL) * CELL + CELL / 2, gy = floor(lat / CELL) * CELL + CELL / 2)]
  out[[st]] <- d[, .(n = .N, CCCCCC = mean(rot6 == "CCCCCC"), CCSCCS = mean(rot6 %in% c("SCCSCC", "CSCCSC")), CSCSCS = mean(rot6 == "SCSCSC")),
                   by = .(gx, gy)][, state := st]
  cat(st, nrow(d), "windows\n")
}
cells <- rbindlist(out); fwrite(cells, "results/rotation_map_cells.csv")
long <- melt(cells[n >= MIN_N], id.vars = c("gx", "gy", "n", "state"), variable.name = "rot", value.name = "share")
long[, rot := factor(rot, c("CCCCCC", "CCSCCS", "CSCSCS"),
                     c("Continuous corn", "Two corn, one soybean", "Strict alternation"))]
usa <- map_data("state", region = c("iowa", "illinois", "indiana", "minnesota", "missouri", "ohio", "south dakota", "wisconsin"))
g <- ggplot() +
  geom_tile(data = long, aes(gx, gy, fill = share * 100), width = CELL, height = CELL) +
  geom_polygon(data = usa, aes(long, lat, group = group), fill = NA, colour = "grey30", linewidth = 0.25) +
  facet_wrap(~rot, nrow = 1) + coord_quickmap() +
  scale_fill_viridis_c(option = "C", direction = -1, name = "Share of corn field-years (%)", limits = c(0, NA)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid = element_blank(),
                                        axis.text = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave(file.path(FIGS_DIR, "be_rotation_map.pdf"), g, width = 9, height = 4.2)
print(long[, .(min = min(share), med = median(share), max = max(share), cells = .N), by = rot])
