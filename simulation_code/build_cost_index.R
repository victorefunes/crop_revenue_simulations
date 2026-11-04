# ============================================================================
# Real (2025-dollar) corn and soybean cost per planted acre by ERS Farm Resource Region and year, 1996-2025.
# Source: USDA ERS Commodity Costs and Returns (corn.csv, soybeans.csv; ARMS-based), deflated with the BEA GDP
# price deflator (FRED GDPDEF, annual). Total non-land cost = "Total, costs listed" minus "Opportunity cost of land"
# (operating costs plus allocated overhead other than land). Operating cost is also kept as the direct-cost series.
# Inputs : SIM_ERS_DIR (default D:/sim_run/ers) with corn.csv, soybeans.csv, gdpdef.csv.
# Output : data/ers_costs_real.csv  (region, crop, year, operating_nom, nonland_nom, deflator, operating_real, nonland_real)
# State -> region map used by run_breakeven_years.R: IA/IL/IN/MO/OH Heartland, WI Northern Crescent, MN/SD Northern Great Plains.
# ============================================================================
suppressPackageStartupMessages(library(data.table))
ERS <- Sys.getenv("SIM_ERS_DIR", "D:/sim_run/ers")
rd <- function(cr) { d <- fread(file.path(ERS, paste0(cr, ".csv"))); d[, crop := c(corn = "C", soybeans = "S")[[cr]]]; d }
d <- rbind(rd("corn"), rd("soybeans"))
d <- d[Item %in% c("Total, operating costs", "Total, costs listed", "Opportunity cost of land"),
       .(crop, region = Region, year = Year, item = Item, value = Value)]
w <- dcast(d, crop + region + year ~ item, value.var = "value")
setnames(w, c("Total, operating costs", "Total, costs listed", "Opportunity cost of land"), c("operating_nom", "total_nom", "land_nom"))
w[, nonland_nom := total_nom - land_nom]
g <- fread(file.path(ERS, "gdpdef.csv")); g[, year := as.integer(substr(observation_date, 1, 4))]
base <- g[year == 2025, GDPDEF]
w[, deflator := g$GDPDEF[match(year, g$year)] / base]
w[, `:=`(operating_real = operating_nom / deflator, nonland_real = nonland_nom / deflator)]
w <- w[region != "U.S. total" & !is.na(nonland_real), .(region, crop, year, operating_nom, nonland_nom, deflator, operating_real, nonland_real)][order(region, crop, year)]
fwrite(w, "data/ers_costs_real.csv")
print(dcast(w[region %in% c("Heartland", "Northern Crescent", "Northern Great Plains") & year %in% c(2008, 2012, 2016, 2020, 2022, 2025)],
            crop + region ~ year, value.var = "nonland_real")[, lapply(.SD, function(x) if (is.numeric(x)) round(x) else x)])
