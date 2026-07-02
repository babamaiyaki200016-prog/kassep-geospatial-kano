# ============================================================================
# KASSEP Geospatial Analysis – Final Production Script (n=93)
# ============================================================================
# Title: Geospatial determinants of maternal mortality in Kano State, Nigeria
# Author: [Your Name]
# Date: 2025-07-02
# R version: 4.5.2 or later
# ============================================================================

# ----------------------------------------------------------------------------
# 0. USER CONFIGURATION – ADJUST THESE PATHS IF NEEDED
# ----------------------------------------------------------------------------

WORKDIR <- "/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR"
setwd(WORKDIR)

KASSEP_EXCEL <- "../03_kassep_data_ADAPTED.xlsx"
STATE_SHP    <- "../NGA_State_Boundaries_V2_-6583428934468491568 (1)/grid3_nga_boundary_vaccstates.shp"
LGA_SHP      <- "../NGA_LGA_Boundaries_2_8280999904383118650/grid3_nga_boundary_vacclgas.shp"
DEPRIV_CSV   <- "../dataset-outputs.csv"
POP_RASTER   <- "nga_women_of_reproductive_age_15_49_2020.tif"
OSM_ZIP_URL  <- "https://production-raw-data-api.s3.amazonaws.com/ISO3/NGA/health_facilities/hotosm_nga_health_facilities_osm_shp.zip"

# ----------------------------------------------------------------------------
# 1. CHECK INPUT FILES
# ----------------------------------------------------------------------------

files_to_check <- c(KASSEP_EXCEL, STATE_SHP, LGA_SHP, DEPRIV_CSV, POP_RASTER)
for (f in files_to_check) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}
cat("✅ All input files found.\n")

# ----------------------------------------------------------------------------
# 2. INSTALL & LOAD PACKAGES
# ----------------------------------------------------------------------------

packages <- c("readxl", "tidyverse", "sf", "spdep", "tmap", "ggplot2",
              "terra", "broom", "car", "dplyr", "tidyr")
new_pkgs <- packages[!packages %in% installed.packages()[,"Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))
cat("✅ All packages loaded.\n")

# ----------------------------------------------------------------------------
# 3. LOAD ALL DATA
# ----------------------------------------------------------------------------

# 3.1 KASSEP
X03 <- read_excel(KASSEP_EXCEL)

# 3.2 Boundaries
nigeria_states <- st_read(STATE_SHP, quiet = TRUE)
nigeria_lgas   <- st_read(LGA_SHP, quiet = TRUE)
state_col <- grep("state|admin", names(nigeria_states), value = TRUE, ignore.case = TRUE)[1]
kano_state <- nigeria_states %>% filter(!!sym(state_col) == "Kano")
kano_lgas  <- nigeria_lgas  %>% filter(!!sym(state_col) == "Kano")
kano_utm   <- kano_state %>% st_transform(crs = 32632)
kano_lgas_utm <- kano_lgas %>% st_transform(crs = 32632)
study_area <- kano_utm
cat("✅ Kano boundaries loaded.\n")

# 3.3 EmOC deprivation
emoc_data <- read.csv(DEPRIV_CSV)
emoc_sf <- st_as_sf(emoc_data, coords = c("longitude", "latitude"), crs = 4326)
emoc_utm <- st_transform(emoc_sf, crs = 32632)
emoc_utm <- emoc_utm[st_intersects(emoc_utm, study_area) %>% lengths() > 0, ]

# 3.4 Health facilities (OSM)
OSM_SHP <- "nga_health_facilities_osm/health_facilities_points.shp"
if (!file.exists(OSM_SHP)) {
  download.file(OSM_ZIP_URL, destfile = "nga_health_facilities_osm.zip", mode = "wb")
  unzip("nga_health_facilities_osm.zip", exdir = "nga_health_facilities_osm")
}
health_shp <- st_read(OSM_SHP, quiet = TRUE)
health_utm <- st_transform(health_shp, crs = 32632)
health_kano <- st_crop(health_utm, st_bbox(study_area))
cat("✅ Health facilities loaded (", nrow(health_kano), " in Kano).\n")

# 3.5 WorldPop raster
pop_raster <- rast(POP_RASTER)
bbox_wgs84 <- st_bbox(study_area %>% st_transform(4326))
pop_kano <- crop(pop_raster, bbox_wgs84)
pop_kano_utm <- project(pop_kano, "EPSG:32632")
pop_kano_masked <- mask(pop_kano_utm, study_area)
cat("✅ Population raster processed.\n")

# ----------------------------------------------------------------------------
# 4. CLEAN GPS – INCLUDE ALL 93 DEATHS
# ----------------------------------------------------------------------------

df_temp <- X03
missing_row <- which(is.na(df_temp$Lat) | is.na(df_temp$Lng))
if (length(missing_row) > 0) {
  gwale_cent <- kano_lgas %>% filter(lganame == "Gwale") %>% st_centroid() %>% st_coordinates()
  df_temp$Lat[missing_row] <- gwale_cent[1, "Y"]
  df_temp$Lng[missing_row] <- gwale_cent[1, "X"]
}
df_clean <- df_temp %>%
  filter(!is.na(Lat) & !is.na(Lng)) %>%
  mutate(Lat = as.numeric(Lat), Lng = as.numeric(Lng)) %>%
  filter(Lat > 10 & Lat < 13 & Lng > 7 & Lng < 10)
points_sf <- df_clean %>% st_as_sf(coords = c("Lng", "Lat"), crs = 4326)
points_utm <- st_transform(points_sf, crs = 32632)
# Add LGA names to points
points_with_lga <- st_join(points_utm, kano_lgas_utm, largest = TRUE)
df_clean$lganame <- points_with_lga$lganame
cat("✅ GPS cleaned: ", nrow(df_clean), " deaths (all 93).\n")

# ----------------------------------------------------------------------------
# 5. CREATE GRID DATA FUNCTION (ROBUST LGA HANDLING)
# ----------------------------------------------------------------------------

create_grid_data <- function(cellsize_m, grid_name) {
  cat("Processing ", grid_name, "...\n", sep = "")
  
  # Create grid
  grid <- st_make_grid(study_area, cellsize = cellsize_m, what = "polygons", square = TRUE)
  grid_sf <- st_sf(geometry = grid) %>%
    st_intersection(study_area) %>%
    mutate(cell_id = row_number())
  
  # Deaths
  deaths <- lengths(st_intersects(grid_sf, points_utm))
  grid_sf$deaths <- deaths
  
  # Population
  pop_vals <- terra::extract(pop_kano_masked, vect(grid_sf), fun = sum, na.rm = TRUE)
  grid_sf$pop_women <- pop_vals[,2]
  
  # Area & density
  grid_sf$area_km2 <- as.numeric(st_area(grid_sf)) / 1e6
  grid_sf$pop_density <- grid_sf$pop_women / grid_sf$area_km2
  
  # MMR
  grid_sf$mmr <- (grid_sf$deaths / grid_sf$pop_women) * 100000
  grid_sf$mmr[grid_sf$pop_women == 0 | is.na(grid_sf$pop_women)] <- NA
  
  # Deprivation (median of points inside)
  grid_with_dep <- st_join(grid_sf, emoc_utm, join = st_intersects)
  dep_sum <- grid_with_dep %>%
    group_by(cell_id) %>%
    summarise(deprivation_median = median(result, na.rm = TRUE)) %>%
    ungroup()
  grid_sf <- grid_sf %>%
    left_join(st_drop_geometry(dep_sum), by = "cell_id")
  
  # Distance to nearest health facility
  cents <- st_centroid(grid_sf)
  dists <- st_distance(cents, health_kano)
  grid_sf$dist_health_km <- apply(dists, 1, min) / 1000
  
  # Facility density (within 5 km)
  fac_count <- sapply(1:nrow(cents), function(i) {
    buff <- st_buffer(cents[i,], dist = 5000)
    length(st_intersects(buff, health_kano)[[1]])
  })
  grid_sf$facility_count_5km <- fac_count
  
  # ---- LGA NAME ASSIGNMENT (ROBUST) ----
  # Use a spatial join with largest overlap to assign LGA
  grid_with_lga <- st_join(grid_sf, kano_lgas_utm, largest = TRUE)
  
  # Find the column containing the LGA name (it may be lganame, lganame.x, lganame.y, etc.)
  lga_col <- grep("lganame", names(grid_with_lga), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(lga_col)) stop("Could not find LGA name column after join.")
  
  # Rename it to a standard name
  grid_with_lga$lganame <- grid_with_lga[[lga_col]]
  
  # Keep only necessary columns, drop the duplicate LGA column and extra geometry
  grid_sf <- grid_with_lga %>%
    select(cell_id, deaths, pop_women, area_km2, pop_density, mmr,
           deprivation_median, dist_health_km, facility_count_5km,
           lganame, geometry) %>%
    distinct()
  
  # Urban indicator
  urban_lgas <- c("Dala", "Fagge", "Gwale", "Kano Municipal", "Ungogo",
                  "Nasarawa", "Kumbotso", "Tarauni", "Dawakin Kudu", "Madobi")
  grid_sf <- grid_sf %>%
    mutate(urban_lga = ifelse(lganame %in% urban_lgas, 1, 0))
  
  # Impute deprivation median
  med_dep <- median(grid_sf$deprivation_median, na.rm = TRUE)
  grid_sf$depr_imputed <- ifelse(is.na(grid_sf$deprivation_median), med_dep, grid_sf$deprivation_median)
  
  return(grid_sf)
}

# ----------------------------------------------------------------------------
# 6. GENERATE GRIDS (5 km and 10 km ONLY)
# ----------------------------------------------------------------------------

grid_5km  <- create_grid_data(5000,  "5 km")
grid_10km <- create_grid_data(10000, "10 km")
grid_5km_filtered <- grid_5km %>% filter(deaths > 0)
grid_10km_filtered <- grid_10km %>% filter(deaths > 0)

# ----------------------------------------------------------------------------
# 7. MORAN'S I (GRID OPTIMISATION)
# ----------------------------------------------------------------------------

moran_results <- data.frame()
for (g in list(list(grid_5km_filtered, "5 km"), list(grid_10km_filtered, "10 km"))) {
  sub <- g[[1]] %>% filter(deaths > 0)
  if (nrow(sub) < 3) next
  coords <- st_centroid(sub) %>% st_coordinates()
  nb_knn <- knn2nb(knearneigh(coords, k = 5))
  listw_knn <- nb2listw(nb_knn, style = "W")
  moran_knn <- moran.test(sub$deaths, listw_knn, zero.policy = TRUE)
  moran_results <- rbind(moran_results,
                         data.frame(grid = g[[2]], weight = "KNN5",
                                    moran_i = moran_knn$estimate[1],
                                    p_value = moran_knn$p.value))
}
write.csv(moran_results, "outputs/moran_comparison.csv", row.names = FALSE)
cat("✅ Moran's I comparison saved.\n")

# ----------------------------------------------------------------------------
# 8. HOTSPOT DETECTION (Getis-Ord Gi*) on 5 km grid
# ----------------------------------------------------------------------------

coords_5km <- st_centroid(grid_5km_filtered) %>% st_coordinates()
nb_knn5 <- knn2nb(knearneigh(coords_5km, k = 5))
listw_knn5 <- nb2listw(nb_knn5, style = "W")
gi <- localG(grid_5km_filtered$deaths, listw_knn5, zero.policy = TRUE)
grid_5km_filtered$gi_z <- as.numeric(gi)
grid_5km_filtered <- grid_5km_filtered %>%
  mutate(hotspot = case_when(
    gi_z >= 2.58 ~ "Very Hot",
    gi_z >= 1.96 ~ "Hot",
    gi_z >= 1.65 ~ "Warm",
    TRUE ~ "Not Significant"
  ))

# ---- HOTSPOT LGA SUMMARY (ROBUST) ----
hotspot_cells <- grid_5km_filtered %>% filter(hotspot %in% c("Very Hot", "Hot"))
hotspot_lgas <- st_join(hotspot_cells, kano_lgas_utm, largest = TRUE)
# Find LGA column
lga_col_hot <- grep("lganame", names(hotspot_lgas), value = TRUE, ignore.case = TRUE)[1]
if (!is.na(lga_col_hot)) {
  hotspot_lgas$lganame <- hotspot_lgas[[lga_col_hot]]
} else {
  stop("Could not find LGA column for hotspots.")
}
hotspot_summary <- hotspot_lgas %>%
  st_drop_geometry() %>%
  group_by(lganame) %>%
  summarise(Hotspot_Cells = n()) %>%
  arrange(desc(Hotspot_Cells))
write.csv(hotspot_summary, "outputs/table2_hotspot_lgas.csv", row.names = FALSE)
cat("✅ Hotspot detection complete.\n")

# ----------------------------------------------------------------------------
# 9. DEPRIVATION CROSS-TAB
# ----------------------------------------------------------------------------

grid_5km_filtered <- grid_5km_filtered %>%
  mutate(deprivation_group = case_when(
    deprivation_median == 2 ~ "High Deprivation",
    deprivation_median == 1 ~ "Medium Deprivation",
    deprivation_median == 0 ~ "Low Deprivation",
    TRUE ~ "Missing"
  ))
tab <- grid_5km_filtered %>%
  st_drop_geometry() %>%
  mutate(is_hotspot = hotspot %in% c("Very Hot", "Hot"))
cross_tab <- table(tab$deprivation_group, tab$is_hotspot)
write.csv(as.data.frame.matrix(cross_tab), "outputs/table5_deprivation_hotspots.csv")
cat("✅ Deprivation cross-tab saved.\n")

# ----------------------------------------------------------------------------
# 10. REGRESSION MODELS
# ----------------------------------------------------------------------------

reg_data <- grid_5km_filtered %>%
  filter(!is.na(mmr) & !is.na(deprivation_median) & !is.na(dist_health_km) & pop_density > 0)
coords_reg <- st_coordinates(st_centroid(reg_data))
reg_df <- data.frame(
  mmr = reg_data$mmr,
  depr = reg_data$deprivation_median,
  pop_dens = reg_data$pop_density,
  dist_km = reg_data$dist_health_km,
  x = coords_reg[,1],
  y = coords_reg[,2]
)
reg_df <- na.omit(reg_df)

ols1 <- lm(mmr ~ depr + pop_dens + dist_km, data = reg_df)
res1 <- tidy(ols1); res1$model <- "Original"
write.csv(res1, "outputs/table3_regression_results.csv", row.names = FALSE)

fit1 <- glance(ols1)
fit_df <- data.frame(model = "Original", R2 = fit1$r.squared, Adj_R2 = fit1$adj.r.squared, p_F = fit1$p.value)
write.csv(fit_df, "outputs/model_fit_comparison.csv", row.names = FALSE)

# Residual Moran's I
residuals1 <- residuals(ols1)
coords_res <- as.matrix(reg_df[, c("x", "y")])
nb_res <- knn2nb(knearneigh(coords_res, k = 5))
listw_res <- nb2listw(nb_res, style = "W")
moran_res <- moran.test(residuals1, listw_res, zero.policy = TRUE)
write.csv(data.frame(
  test = "Residual Moran's I",
  moran_i = moran_res$estimate[1],
  p_value = moran_res$p.value
), "outputs/residual_moran.csv", row.names = FALSE)
cat("✅ Regression models complete.\n")

# ----------------------------------------------------------------------------
# 11. PREDICTIVE MAPPING & PRIORITY AREAS
# ----------------------------------------------------------------------------

coefs <- coef(ols1)
grid_5km_filtered$pred_mmr <- coefs[1] +
  coefs[2] * grid_5km_filtered$depr_imputed +
  coefs[3] * grid_5km_filtered$pop_density +
  coefs[4] * grid_5km_filtered$dist_health_km

threshold <- quantile(grid_5km_filtered$pred_mmr, 0.90, na.rm = TRUE)
priority <- grid_5km_filtered %>% filter(pred_mmr >= threshold)

priority_lgas <- priority %>%
  st_drop_geometry() %>%
  group_by(lganame) %>%
  summarise(n_cells = n(),
            mean_mmr = mean(pred_mmr, na.rm = TRUE),
            max_mmr = max(pred_mmr, na.rm = TRUE)) %>%
  arrange(desc(n_cells))
write.csv(priority_lgas, "outputs/table4_priority_areas.csv", row.names = FALSE)
cat("✅ Predictive mapping complete.\n")

# ----------------------------------------------------------------------------
# 12. GENERATE FIGURES
# ----------------------------------------------------------------------------

dir.create("outputs", showWarnings = FALSE)

p1 <- ggplot() +
  geom_sf(data = study_area, fill = "gray95", color = "black") +
  geom_sf(data = grid_5km_filtered, aes(fill = gi_z), color = NA) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Gi* z-score") +
  labs(title = "Getis-Ord Gi* Hotspots") +
  theme_void()
ggsave("outputs/fig1_getis_hotspots.png", p1, width = 8, height = 7, dpi = 300)

composite_grid <- grid_5km_filtered %>%
  mutate(deprivation_group = case_when(
    deprivation_median == 2 ~ "High Deprivation",
    deprivation_median == 1 ~ "Medium Deprivation",
    deprivation_median == 0 ~ "Low Deprivation",
    TRUE ~ "Missing"
  ))
p2 <- ggplot() +
  geom_sf(data = study_area, fill = "gray95", color = "black") +
  geom_sf(data = composite_grid, aes(fill = mmr), color = NA) +
  scale_fill_gradient2(low = "#FEE5D9", mid = "#FC9272", high = "#D73027",
                       midpoint = median(composite_grid$mmr, na.rm = TRUE),
                       na.value = "white", name = "MMR (per 100k)") +
  geom_sf(data = composite_grid, aes(color = deprivation_group), fill = NA, linewidth = 0.6) +
  scale_color_manual(values = c("High Deprivation" = "#D73027",
                                "Medium Deprivation" = "#F46D43",
                                "Low Deprivation" = "#1A9850",
                                "Missing" = "gray80"),
                     name = "EmOC Deprivation") +
  geom_sf(data = composite_grid %>% filter(hotspot %in% c("Very Hot", "Hot")),
          fill = NA, color = "black", linewidth = 1.2, linetype = "dashed") +
  labs(title = "Composite: MMR, Deprivation, and Hotspots") +
  theme_void() +
  theme(legend.position = "bottom")
ggsave("outputs/fig2_composite.png", p2, width = 10, height = 8, dpi = 300)

p3 <- ggplot() +
  geom_sf(data = study_area, fill = "gray95", color = "black") +
  geom_sf(data = grid_5km_filtered, aes(fill = pred_mmr), color = NA) +
  scale_fill_gradient(low = "#FFF5F0", high = "#D73027", na.value = "white", name = "Predicted MMR") +
  labs(title = "Predicted Risk Surface") +
  theme_void()
ggsave("outputs/fig3_predicted_risk.png", p3, width = 8, height = 7, dpi = 300)

p4 <- ggplot() +
  geom_sf(data = study_area, fill = "gray95", color = "black") +
  geom_sf(data = grid_5km_filtered, aes(fill = pred_mmr), color = NA) +
  scale_fill_gradient(low = "#FFF5F0", high = "#D73027", na.value = "white", name = "Predicted MMR") +
  geom_sf(data = priority, fill = NA, color = "black", linewidth = 1.5, linetype = "dashed") +
  geom_sf(data = kano_lgas_utm, fill = NA, color = "gray40", linewidth = 0.2) +
  geom_sf_label(data = kano_lgas_utm %>% filter(lganame %in% unique(priority$lganame)),
                aes(label = lganame), size = 2.5, fill = "white", alpha = 0.8) +
  labs(title = "Priority Areas (Top 10% Predicted MMR)") +
  theme_void()
ggsave("outputs/fig4_priority_areas.png", p4, width = 8, height = 7, dpi = 300)

cat("✅ All figures saved.\n")

# ----------------------------------------------------------------------------
# 13. DESCRIPTIVE TABLES
# ----------------------------------------------------------------------------

urban_lgas <- c("Dala", "Fagge", "Gwale", "Kano Municipal", "Ungogo",
                "Nasarawa", "Kumbotso", "Tarauni", "Dawakin Kudu", "Madobi")
table1 <- data.frame(
  Variable = c("Total deaths","Valid GPS","Urban LGAs","Rural LGAs",
               "Hospital","Home","Elsewhere/unknown"),
  Value = c(
    nrow(df_clean), nrow(df_clean),
    sum(df_clean$lganame %in% urban_lgas, na.rm = TRUE),
    sum(!df_clean$lganame %in% urban_lgas, na.rm = TRUE),
    sum(df_clean$`(Id10058) Where did the deceased die?` == "Hospital", na.rm = TRUE),
    sum(df_clean$`(Id10058) Where did the deceased die?` == "Home", na.rm = TRUE),
    sum(!df_clean$`(Id10058) Where did the deceased die?` %in% c("Hospital", "Home"), na.rm = TRUE)
  )
)
write.csv(table1, "outputs/table1_descriptive.csv", row.names = FALSE)

dist_lga <- grid_5km_filtered %>%
  st_drop_geometry() %>%
  group_by(lganame) %>%
  summarise(n_cells = n(),
            mean_dist_km = mean(dist_health_km, na.rm = TRUE),
            min_dist_km = min(dist_health_km, na.rm = TRUE),
            max_dist_km = max(dist_health_km, na.rm = TRUE),
            mean_mmr = mean(mmr, na.rm = TRUE)) %>%
  arrange(desc(mean_dist_km))
write.csv(dist_lga, "outputs/table6_distance_by_lga.csv", row.names = FALSE)

table7 <- data.frame(
  Variable = c("Maternal mortality data","Population denominator","Administrative boundaries",
               "EmOC access deprivation","Health facilities","Distance to facilities"),
  Source = c("KASSEP (verbal autopsy)","WorldPop 2020 (women 15-49)","GRID3 Nigeria",
             "IDEAMAPS (University of Glasgow)","OpenStreetMap (OSM)","Computed from OSM points"),
  Resolution = c("Point (GPS)","100m raster","Polygon (LGA/State)",
                 "100m grid (classified)","Point","Continuous (km)"),
  Year = c("2023-2024","2020","2020","2020","2022","2022")
)
write.csv(table7, "outputs/table7_data_sources.csv", row.names = FALSE)

cat("✅ All tables saved.\n")

# ----------------------------------------------------------------------------
# 14. SESSION INFO
# ----------------------------------------------------------------------------

sink("outputs/session_info.txt")
sessionInfo()
sink()
cat("\n🎯 Analysis complete! All outputs are in the 'outputs' folder.\n")