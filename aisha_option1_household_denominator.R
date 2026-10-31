# ============================================================================
# Response to Dr Aisha Adamu (Re-draft_V3, comment 3, Option 1)
# Retrospective wave only (n = 93 deaths). Deaths AND denominator both come
# from the KASSEP household listing (household GPS), so surveyed cells with
# no deaths are true zeros and the numerator/denominator mismatch is removed.
#   A. Cell-level negative binomial / Poisson, log(surveyed live births) offset
#   B. Household-level logistic: maternal death in household vs. not, among
#      households with >= 1 live birth, cluster-robust by EA
# Writes to outputs_n116/aisha_option1/. Does not touch other outputs.
# ============================================================================

setwd("/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR")
HH_XLSX   <- "/Users/babamusa/Downloads/Sharable_Final KASSEP BASELINE DATA.xlsx"
STATE_SHP <- "../NGA_State_Boundaries_V2_-6583428934468491568 (1)/grid3_nga_boundary_vaccstates.shp"
LGA_SHP   <- "../NGA_LGA_Boundaries_2_8280999904383118650/grid3_nga_boundary_vacclgas.shp"
DEPRIV_CSV <- "../dataset-outputs.csv"
OSM_SHP   <- "nga_health_facilities_osm/health_facilities_points.shp"
# BIRTHS=live: offset = babies delivered minus stillbirths (floored at 0), closer to the published live-birth denominator.
# Default offset = Number_babies_delivered, which is TOTAL births (live births + stillbirths).
LIVE <- Sys.getenv("BIRTHS") == "live"
OUT <- if (LIVE) "outputs_n116/aisha_option1_livebirths" else "outputs_n116/aisha_option1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages({library(readxl); library(dplyr); library(sf); library(MASS)
                  library(sandwich); library(lmtest); library(spdep)})
select <- dplyr::select

# ---- 1. Household listing ---------------------------------------------------
hh_raw <- suppressMessages(read_excel(HH_XLSX, sheet = 1, .name_repair = "unique_quiet"))
hh <- hh_raw %>%
  transmute(LGA = trimws(LGA), EA, Settlement_type,
            completed = `survey_type...654` == 1,
            strategy = `survey_type...10`,
            lat = as.numeric(latitude), lon = as.numeric(longitude),
            births = suppressWarnings(as.numeric(Number_babies_delivered)),
            still = suppressWarnings(as.numeric(number_stillbirth)),
            death = suppressWarnings(as.numeric(No_maternal_death)))
hh$LGA[hh$LGA == "Ungoggo"] <- "Ungogo"           # match GRID3 spelling
hh$LGA[hh$LGA == "Garum Mallam"] <- "Garun Mallam"

flow <- data.frame(
  step = c("Listing rows", "Completed households", "...with GPS in Kano box",
           "...with >=1 live birth", "Maternal deaths (completed, GPS)",
           "Live births (completed, GPS)"),
  n = NA_real_)
flow$n[1] <- nrow(hh)
hh <- hh %>% filter(completed %in% TRUE)
flow$n[2] <- nrow(hh)
all_deaths <- sum(hh$death, na.rm = TRUE); all_births <- sum(hh$births, na.rm = TRUE)
hh <- hh %>% filter(!is.na(lat), !is.na(lon), lat > 10, lat < 13, lon > 7, lon < 10)
flow$n[3] <- nrow(hh)
hh$births[is.na(hh$births)] <- 0; hh$death[is.na(hh$death)] <- 0; hh$still[is.na(hh$still)] <- 0
if (LIVE) hh$births <- pmax(hh$births - hh$still, 0)
# Maternal deaths need not follow a live birth (stillbirth, miscarriage, abortion, antepartum death),
# so death households with 0 live births are valid numerator events, not missing data.
cat("Death households with 0 recorded live births:", sum(hh$death > 0 & hh$births == 0), "of", sum(hh$death > 0), "\n")
flow$n[4] <- sum(hh$births > 0)
flow$n[5] <- sum(hh$death); flow$n[6] <- sum(hh$births)
cat("All completed HH: deaths", all_deaths, " births", all_births, "\n")
print(flow); write.csv(flow, file.path(OUT, "flow_households.csv"), row.names = FALSE)
cat("Deaths lost to missing GPS:", all_deaths - sum(hh$death),
    " births lost:", all_births - sum(hh$births), "\n")

# ---- 2. Spatial layers ------------------------------------------------------
st <- st_read(STATE_SHP, quiet = TRUE)
sc <- grep("state|admin", names(st), value = TRUE, ignore.case = TRUE)[1]
kano <- st %>% filter(!!sym(sc) == "Kano") %>% st_transform(32632)
lgas <- st_read(LGA_SHP, quiet = TRUE) %>% filter(!!sym(sc) == "Kano") %>% st_transform(32632)
fac  <- st_read(OSM_SHP, quiet = TRUE) %>% st_transform(32632)
emoc <- read.csv(DEPRIV_CSV) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(32632)
emoc <- emoc[lengths(st_intersects(emoc, kano)) > 0, ]

hh_sf <- st_as_sf(hh, coords = c("lon", "lat"), crs = 4326) %>% st_transform(32632)
hh_sf$dist_km <- as.numeric(st_distance(hh_sf, fac[st_nearest_feature(hh_sf, fac), ],
                                        by_element = TRUE)) / 1000
hh_sf$census <- as.integer(hh_sf$strategy == "Census")
hh_sf$urban <- as.integer(hh_sf$Settlement_type == "Urban")
hh_sf$semiurban <- as.integer(hh_sf$Settlement_type == "Semi-Urban")

# Strategy by LGA
strat <- hh_sf %>% st_drop_geometry() %>% group_by(LGA) %>%
  summarise(strategy = paste(unique(strategy), collapse = "/"), households = n(),
            births = sum(births), deaths = sum(death),
            crude_MMR_per_100k = round(100000 * deaths / births)) %>% arrange(strategy, LGA)
write.csv(strat, file.path(OUT, "strategy_by_lga.csv"), row.names = FALSE)

# ---- 3. A. Cell-level count model ------------------------------------------
build_cells <- function(cell_m) {
  g <- st_sf(geometry = st_make_grid(kano, cellsize = cell_m)) %>%
    st_intersection(kano) %>% mutate(cell_id = row_number()) %>% dplyr::select(cell_id)
  j <- st_join(hh_sf, g)
  d <- j %>% st_drop_geometry() %>% filter(!is.na(cell_id)) %>% group_by(cell_id) %>%
    summarise(hh = n(), births = sum(births), deaths = sum(death),
              census_share = mean(census), urban_share = mean(urban),
              dist_km = mean(dist_km), .groups = "drop")
  g <- inner_join(g, d, by = "cell_id")
  g$depr_n <- lengths(st_intersects(g, emoc))
  dep <- st_join(emoc, g["cell_id"]) %>% st_drop_geometry() %>% filter(!is.na(cell_id)) %>%
    group_by(cell_id) %>% summarise(depr = mean(result == 2), .groups = "drop")  # share high-deprivation points
  g <- left_join(g, dep, by = "cell_id")
  g$lga <- lgas$lganame[st_nearest_feature(st_centroid(g), lgas)]
  g
}

fit_cells <- function(cell_m, label) {
  g <- build_cells(cell_m) %>% filter(births > 0)
  cat("\n=====", label, ": cells with surveyed births =", nrow(g), " deaths =", sum(g$deaths),
      " zero-death cells =", sum(g$deaths == 0), "\n")
  d <- st_drop_geometry(g)
  res <- list()
  f_main <- deaths ~ dist_km + census_share + offset(log(births))
  f_depr <- deaths ~ dist_km + census_share + depr + offset(log(births))
  pois <- glm(f_main, family = poisson, data = d)
  nb   <- tryCatch(glm.nb(f_main, data = d), error = function(e) NULL)
  disp <- sum(residuals(pois, "pearson")^2) / pois$df.residual
  cat("Poisson Pearson dispersion:", round(disp, 2), "\n")
  tidy <- function(m, nm) {
    ct <- coeftest(m, vcov. = vcovCL(m, cluster = d$lga[as.numeric(rownames(model.frame(m)))]))
    data.frame(model = nm, term = rownames(ct), est = ct[, 1], se = ct[, 2], p = ct[, 4],
               RR = exp(ct[, 1]), lo = exp(ct[, 1] - 1.96 * ct[, 2]), hi = exp(ct[, 1] + 1.96 * ct[, 2]),
               row.names = NULL)
  }
  out <- tidy(pois, "Poisson (LGA-cluster robust SE)")
  if (!is.null(nb)) out <- rbind(out, tidy(nb, "Negative binomial (LGA-cluster robust SE)"))
  # Sensitivity: add deprivation, restricted to IDEAMAPS-covered cells
  dd <- d %>% filter(!is.na(depr))
  cat("Cells with IDEAMAPS coverage:", nrow(dd), " deaths in them:", sum(dd$deaths), "\n")
  if (sum(dd$deaths) >= 10) {
    pd <- glm(f_depr, family = poisson, data = dd)
    ct <- coeftest(pd, vcov. = vcovCL(pd, cluster = dd$lga))
    out <- rbind(out, data.frame(model = "Poisson + deprivation (IDEAMAPS-covered cells only)",
      term = rownames(ct), est = ct[, 1], se = ct[, 2], p = ct[, 4], RR = exp(ct[, 1]),
      lo = exp(ct[, 1] - 1.96 * ct[, 2]), hi = exp(ct[, 1] + 1.96 * ct[, 2]), row.names = NULL))
  }
  write.csv(out, file.path(OUT, paste0("cell_model_", label, ".csv")), row.names = FALSE)
  write.csv(d, file.path(OUT, paste0("cell_data_", label, ".csv")), row.names = FALSE)
  print(out %>% mutate(across(where(is.numeric), ~round(.x, 3))))
  # Residual spatial autocorrelation (Moran's I on Pearson residuals, KNN5)
  if (nrow(g) > 10) {
    nbk <- knn2nb(knearneigh(st_coordinates(st_centroid(g)), k = min(5, nrow(g) - 1)))
    mi <- moran.test(residuals(pois, "pearson"), nb2listw(nbk, style = "W"), zero.policy = TRUE)
    cat("Residual Moran's I:", round(mi$estimate[1], 3), " p =", round(mi$p.value, 3), "\n")
  }
  invisible(g)
}
g5  <- fit_cells(5000, "5km")
g10 <- fit_cells(10000, "10km")

# ---- 4. B. Household-level logistic ----------------------------------------
cat("\n===== Household-level logistic (all completed households) =====\n")
hd <- hh_sf %>% st_drop_geometry() %>% mutate(died = as.integer(death > 0)) %>%
  filter(!is.na(EA), !is.na(urban), !is.na(semiurban))
cat("Households:", nrow(hd), " with maternal death:", sum(hd$died), "\n")
lg <- glm(died ~ dist_km + census + urban + semiurban, family = binomial, data = hd)
ct <- coeftest(lg, vcov. = vcovCL(lg, cluster = hd$EA))
hres <- data.frame(term = rownames(ct), OR = exp(ct[, 1]), lo = exp(ct[, 1] - 1.96 * ct[, 2]),
                   hi = exp(ct[, 1] + 1.96 * ct[, 2]), p = ct[, 4], row.names = NULL)
print(hres %>% mutate(across(where(is.numeric), ~round(.x, 3))))
write.csv(hres, file.path(OUT, "household_logistic.csv"), row.names = FALSE)

# Crude summary by strategy
cs <- hh_sf %>% st_drop_geometry() %>% group_by(strategy) %>%
  summarise(households = n(), births = sum(births), deaths = sum(death),
            MMR_per_100k = round(100000 * sum(death) / sum(births)))
print(cs); write.csv(cs, file.path(OUT, "crude_by_strategy.csv"), row.names = FALSE)

# ---- 5. Compare with VA GPS points -----------------------------------------
va <- read.csv("combined_points_n116.csv") %>% filter(wave == "retrospective", !is.na(Lat)) %>%
  st_as_sf(coords = c("Lng", "Lat"), crs = 4326) %>% st_transform(32632)
dh <- hh_sf %>% filter(death > 0)
nn <- st_distance(va, dh)
cat("\nVA retrospective points:", nrow(va), "| household-listing death households:", nrow(dh),
    "| VA points within 100 m of a listing death household:", sum(apply(nn, 1, min) < 100), "\n")

# ---- 6. Save cell geometries for figures ------------------------------------
saveRDS(list(g5 = g5, g10 = g10, lgas = lgas, kano = kano, strat = strat), file.path(OUT, "cells_for_figures.rds"))
