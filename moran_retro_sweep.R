# Global Moran's I on death counts per death-containing cell, retrospective deaths only
# (same method as run_analysis_n116.R Table 2: Queen contiguity and KNN5, upper-tailed test).
setwd("/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR")
suppressMessages({library(sf); library(spdep); library(dplyr)})
st <- st_read("../NGA_State_Boundaries_V2_-6583428934468491568 (1)/grid3_nga_boundary_vaccstates.shp", quiet = TRUE)
sc <- grep("state|admin", names(st), value = TRUE, ignore.case = TRUE)[1]
kano <- st %>% filter(!!sym(sc) == "Kano") %>% st_transform(32632)
lgas <- st_read("../NGA_LGA_Boundaries_2_8280999904383118650/grid3_nga_boundary_vacclgas.shp", quiet = TRUE) %>% filter(!!sym(sc) == "Kano") %>% st_transform(32632)
raw <- read.csv("combined_points_n116.csv")
# One retrospective record has no GPS; as in run_analysis_n116.R it is placed at the Gwale LGA centroid.
gw <- lgas %>% filter(lganame == "Gwale") %>% st_centroid() %>% st_transform(4326) %>% st_coordinates()
miss <- is.na(raw$Lat); raw$Lng[miss] <- gw[1, "X"]; raw$Lat[miss] <- gw[1, "Y"]
ALL <- Sys.getenv("ALL") == "1"
pts <- raw %>% filter(ALL | wave == "retrospective") %>% st_as_sf(coords = c("Lng", "Lat"), crs = 4326) %>% st_transform(32632)
res <- data.frame()
for (km in c(2, 3, 5, 10)) {
  g <- st_sf(geometry = st_make_grid(kano, cellsize = km * 1000)) %>% st_intersection(kano)
  g$d <- lengths(st_intersects(g, pts)); s <- g %>% filter(d > 0)
  q <- moran.test(s$d, nb2listw(poly2nb(s, queen = TRUE), style = "W", zero.policy = TRUE), zero.policy = TRUE)
  nb <- knn2nb(knearneigh(st_coordinates(suppressWarnings(st_centroid(s))), k = 5))
  k <- moran.test(s$d, nb2listw(nb, style = "W"))
  res <- rbind(res, data.frame(grid = paste(km, "km"), weights = c("Queen", "KNN5"), cells = nrow(s),
    moran_i = round(c(q$estimate[1], k$estimate[1]), 3), p = signif(c(q$p.value, k$p.value), 3)))
}
print(res); dir.create("outputs_n116/aisha_hotspots", showWarnings = FALSE)
write.csv(res, if (ALL) "outputs_n116/aisha_hotspots/moran_all_deaths_check.csv" else "outputs_n116/aisha_hotspots/moran_retrospective.csv", row.names = FALSE)
