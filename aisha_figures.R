# New/updated figures for the Aisha V3 revision (household-denominator analysis)
suppressMessages({library(sf);library(dplyr);library(ggplot2);library(readr)})
OUT <- "outputs_n116/aisha_figures"; dir.create(OUT, showWarnings = FALSE)
x <- readRDS("outputs_n116/aisha_option1/cells_for_figures.rds")
g5 <- x$g5; lgas <- x$lgas; kano <- x$kano; strat <- x$strat
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GREY <- "#d9d8d3"
thm <- theme_minimal(base_size = 11) + theme(
  panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
  plot.title = element_text(face = "bold", colour = INK, size = 12),
  plot.subtitle = element_text(colour = INK2, size = 9.5), legend.position = "bottom",
  legend.title = element_text(colour = INK2, size = 9), legend.text = element_text(colour = INK2, size = 9),
  plot.caption = element_text(colour = INK2, size = 8, hjust = 0))

# ---- Figure 4: household-denominator map (5 km) ----
lg <- lgas %>% left_join(strat %>% select(LGA, strategy), by = c("lganame" = "LGA"))
g5$dcat <- cut(g5$deaths, c(-1, 0, 1, 2, Inf), labels = c("0", "1", "2", "3 or more"))
f4 <- ggplot() +
  geom_sf(data = kano, fill = "#fcfcfb", colour = NA) +
  geom_sf(data = g5, aes(fill = dcat), colour = "white", linewidth = 0.15) +
  geom_sf(data = lg %>% filter(strategy == "Census"), fill = NA, colour = INK, linewidth = 0.35, linetype = "22") +
  geom_sf(data = lg %>% filter(is.na(strategy) | strategy != "Census"), fill = NA, colour = "#9b9a94", linewidth = 0.2) +
  scale_fill_manual(values = c("0" = GREY, "1" = "#9ec5f2", "2" = BLUE, "3 or more" = "#123f78"),
                    name = "Maternal deaths per 5 km cell") +
  labs(title = "Maternal deaths among 5 km cells with surveyed births, Kano State",
       subtitle = sprintf("%d cells with surveyed births; %d contain no death. Dashed outlines: LGAs surveyed by census.",
                          nrow(g5), sum(g5$deaths == 0)),
       caption = "Deaths and births from the KASSEP household listing (retrospective wave, n = 91 deaths).") + thm
ggsave(file.path(OUT, "fig4_household_denominator_map.png"), f4, width = 7.5, height = 7.5, dpi = 300, bg = "white")

# ---- Figure 5: forest plot ----
rd <- function(f, lab) {
  d <- read.csv(f); d <- d[d$term %in% c("dist_km", "census_share") &
        d$model %in% c("Poisson (LGA-cluster robust SE)", "Negative binomial (LGA-cluster robust SE)"), ]
  d$spec <- paste(lab, ifelse(grepl("Negative", d$model), "negative binomial", "Poisson"))
  d %>% transmute(term, spec, est = RR, lo, hi)
}
cm <- bind_rows(rd("outputs_n116/aisha_option1/cell_model_5km.csv", "5 km cells,"),
                rd("outputs_n116/aisha_option1/cell_model_10km.csv", "10 km cells,"))
hl <- read.csv("outputs_n116/aisha_option1/household_logistic.csv")
hh <- hl %>% filter(term %in% c("dist_km", "census")) %>%
  transmute(term = ifelse(term == "census", "census_share", term), spec = "Households, logistic", est = OR, lo, hi)
fd <- bind_rows(cm, hh) %>%
  mutate(term = factor(term, c("dist_km", "census_share"),
                       c("Distance to nearest facility (per km)", "Census strategy (vs sampling)")),
         spec = factor(spec, rev(c("5 km cells, Poisson", "5 km cells, negative binomial",
                                   "10 km cells, Poisson", "10 km cells, negative binomial", "Households, logistic"))),
         lab = sprintf("%.2f (%.2f to %.2f)", est, lo, hi))
f5 <- ggplot(fd, aes(est, spec)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = BLUE, linewidth = 0.7) +
  geom_point(colour = BLUE, size = 2.6) +
  geom_text(aes(x = hi, label = lab), hjust = -0.1, size = 2.9, colour = INK2) +
  facet_wrap(~term, ncol = 1, scales = "free_x") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.55))) +
  labs(x = "Rate ratio (cell models) or odds ratio (household model), log scale, 95% CI",
       y = NULL, title = "Distance was not associated with maternal deaths;\ncensus LGAs had more deaths",
       caption = "Cell models: log surveyed-births offset, LGA-clustered SE.\nHousehold model: all 41,007 completed households, EA-clustered SE.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = GRID), strip.text = element_text(face = "bold", colour = INK, hjust = 0),
        plot.title = element_text(face = "bold", colour = INK, size = 11.5), plot.title.position = "plot", plot.caption.position = "plot", axis.text = element_text(colour = INK2),
        axis.title.x = element_text(colour = INK2, size = 9), plot.caption = element_text(colour = INK2, size = 8, hjust = 0))
ggsave(file.path(OUT, "fig5_forest_plot.png"), f5, width = 8, height = 5.5, dpi = 300, bg = "white")

# ---- Supplementary Figure S1: IDEAMAPS extent ----
e <- read.csv("../dataset-outputs.csv") %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(32632)
hull <- st_convex_hull(st_union(e))
va <- read.csv("combined_points_n116.csv") %>% filter(!is.na(Lat)) %>%
  st_as_sf(coords = c("Lng", "Lat"), crs = 4326) %>% st_transform(32632)
va$inside <- lengths(st_intersects(va, hull)) > 0
fs1 <- ggplot() +
  geom_sf(data = kano, fill = "#fcfcfb", colour = "#9b9a94", linewidth = 0.3) +
  geom_sf(data = lgas, fill = NA, colour = "#cfcec8", linewidth = 0.15) +
  geom_sf(data = hull, fill = "#9ec5f2", colour = BLUE, alpha = 0.55, linewidth = 0.4) +
  geom_sf(data = va, aes(shape = inside), colour = INK, size = 1.4, fill = "white") +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1), labels = c("TRUE" = "Inside IDEAMAPS extent", "FALSE" = "Outside"),
                     name = "Recorded maternal death") +
  labs(title = "IDEAMAPS EmOC deprivation covers metropolitan Kano only",
       subtitle = sprintf("Blue area: extent of IDEAMAPS points (about 9%% of the state). %d of %d death locations fall inside it.",
                          sum(va$inside), nrow(va)),
       caption = "Death locations: verbal-autopsy GPS points (retrospective and prospective, n = 115 with coordinates).") + thm
ggsave(file.path(OUT, "figS1_ideamaps_extent.png"), fs1, width = 7.5, height = 7.5, dpi = 300, bg = "white")
cat("done\n")
