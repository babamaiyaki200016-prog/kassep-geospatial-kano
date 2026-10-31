# Aisha comment 1: Gi* hotspots shown separately for the retrospective deaths
# (n = 93) and for all recorded deaths (retrospective + prospective, n = 116).
# Same method as run_analysis_n116.R: 5 km grid, death-containing cells, KNN5 neighbours incl. self.
setwd("/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR")
suppressMessages({library(sf); library(spdep); library(dplyr)})
STATE_SHP <- "../NGA_State_Boundaries_V2_-6583428934468491568 (1)/grid3_nga_boundary_vaccstates.shp"
LGA_SHP   <- "../NGA_LGA_Boundaries_2_8280999904383118650/grid3_nga_boundary_vacclgas.shp"
st <- st_read(STATE_SHP, quiet = TRUE); sc <- grep("state|admin", names(st), value = TRUE, ignore.case = TRUE)[1]
kano <- st %>% filter(!!sym(sc) == "Kano") %>% st_transform(32632)
lgas <- st_read(LGA_SHP, quiet = TRUE) %>% filter(!!sym(sc) == "Kano") %>% st_transform(32632)
raw <- read.csv("combined_points_n116.csv")
# One retrospective record has no GPS; as in run_analysis_n116.R it is placed at the Gwale LGA centroid.
gw <- lgas %>% filter(lganame == "Gwale") %>% st_centroid() %>% st_transform(4326) %>% st_coordinates()
miss <- is.na(raw$Lat); raw$Lng[miss] <- gw[1, "X"]; raw$Lat[miss] <- gw[1, "Y"]
pts <- raw %>% st_as_sf(coords = c("Lng", "Lat"), crs = 4326) %>% st_transform(32632)
grid <- st_sf(geometry = st_make_grid(kano, cellsize = 5000)) %>% st_intersection(kano) %>%
  mutate(cell_id = row_number()) %>% select(cell_id)
THR <- as.numeric(Sys.getenv("THR", "1.96"))   # manuscript definition: Gi* z >= 1.96
run <- function(p, label) {
  g <- grid; g$deaths <- lengths(st_intersects(g, p)); g <- g %>% filter(deaths > 0)
  nb <- knn2nb(knearneigh(st_coordinates(st_centroid(g)), k = 5))
  g$z <- as.numeric(localG(g$deaths, nb2listw(include.self(nb), style = "W"), zero.policy = TRUE))
  g$hot <- g$z >= THR
  g <- suppressWarnings(st_join(g, lgas[, "lganame"], largest = TRUE)); g$lga <- g$lganame   # largest overlap, as in the manuscript
  cat("\n==", label, ": deaths", sum(g$deaths), " death cells", nrow(g), " hotspot cells (Gi* z >=", THR, "):", sum(g$hot), "\n")
  print(g %>% st_drop_geometry() %>% filter(hot) %>% count(lga, name = "cells") %>% arrange(desc(cells)))
  g %>% st_drop_geometry() %>% filter(hot) %>% pull(cell_id)
}
a <- run(pts %>% filter(wave == "retrospective"), "Retrospective only")
b <- run(pts, "All recorded deaths")
cat("\nHotspot cells in both:", length(intersect(a, b)), " retro-only:", length(setdiff(a, b)), " added by prospective:", length(setdiff(b, a)), "\n")
