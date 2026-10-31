# ============================================================================
# Aisha (cell-level suggestion): Poisson model with a spatial random effect (BYM2),
# log surveyed-births offset, distance to nearest facility and census share.
# BYM2 handles overdispersion and spatial correlation together and smooths sparse cells.
# Retrospective wave, same cells as aisha_option1_household_denominator.R (cells with surveyed births).
# Requires INLA (install from https://inla.r-inla-download.org/R/stable).
# ============================================================================
setwd("/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR")
suppressMessages({library(INLA); library(sf); library(spdep); library(dplyr)})
OUT <- "outputs_n116/aisha_bym2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cells <- readRDS("outputs_n116/aisha_option1/cells_for_figures.rds")

fit_bym2 <- function(g, label) {
  g <- g %>% filter(births > 0) %>% mutate(id = row_number())
  nb <- poly2nb(g, queen = TRUE)
  nb2INLA(file.path(OUT, paste0("graph_", label, ".adj")), nb)
  gf <- file.path(OUT, paste0("graph_", label, ".adj"))
  d <- st_drop_geometry(g) %>% mutate(dist_c = dist_km, cen = census_share)
  cat("\n=====", label, ": cells", nrow(d), " deaths", sum(d$deaths), " cells with no neighbour:", sum(card(nb) == 0), "\n")
  ctl <- list(config = FALSE)
  base <- deaths ~ dist_c + cen
  fit <- function(extra) inla(update(base, extra), family = "poisson", data = d, offset = log(births),
          control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE), verbose = FALSE)
  m0 <- inla(base, family = "poisson", data = d, offset = log(births), control.compute = list(dic = TRUE, waic = TRUE))
  m1 <- inla(deaths ~ dist_c + cen + f(id, model = "bym2", graph = gf, scale.model = TRUE, constr = TRUE,
              hyper = list(phi = list(prior = "pc", param = c(0.5, 2/3)), prec = list(prior = "pc.prec", param = c(1, 0.01)))),
             family = "poisson", data = d, offset = log(births), control.compute = list(dic = TRUE, waic = TRUE))
  # BYM2 + LGA random intercept: census strategy is fixed within LGA, so LGA-level clustering is part of the design
  d$lga_id <- as.integer(factor(d$lga))
  m2 <- inla(deaths ~ dist_c + cen + f(lga_id, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
              f(id, model = "bym2", graph = gf, scale.model = TRUE, constr = TRUE,
                hyper = list(phi = list(prior = "pc", param = c(0.5, 2/3)), prec = list(prior = "pc.prec", param = c(1, 0.01)))),
             family = "poisson", data = d, offset = log(births), control.compute = list(dic = TRUE, waic = TRUE))
  f2 <- m2$summary.fixed
  res2 <- data.frame(grid = label, term = rownames(f2)[-1], RR = exp(f2$mean[-1]), lo = exp(f2[["0.025quant"]][-1]), hi = exp(f2[["0.975quant"]][-1]))
  cat("BYM2 + LGA random intercept (rate ratio, 95% CrI):\n"); print(res2 %>% mutate(across(where(is.numeric), ~round(.x, 3))))
  cat(sprintf("DIC BYM2+LGA %.1f ; WAIC %.1f\n", m2$dic$dic, m2$waic$waic))
  fx <- m1$summary.fixed
  res <- data.frame(grid = label, term = rownames(fx)[-1], RR = exp(fx$mean[-1]),
                    lo = exp(fx[["0.025quant"]][-1]), hi = exp(fx[["0.975quant"]][-1]))
  hy <- m1$summary.hyperpar
  cat("Fixed effects (rate ratio, 95% credible interval):\n"); print(res %>% mutate(across(where(is.numeric), ~round(.x, 3))))
  cat("Hyperparameters (precision, mixing phi):\n"); print(round(hy[, c("mean", "0.025quant", "0.975quant")], 3))
  cat(sprintf("DIC: Poisson %.1f | BYM2 %.1f ;  WAIC: Poisson %.1f | BYM2 %.1f\n", m0$dic$dic, m1$dic$dic, m0$waic$waic, m1$waic$waic))
  cat("Effective number of parameters (pD) BYM2:", round(m1$dic$p.eff, 1), "\n")
  res$model <- "BYM2"; res2$model <- "BYM2 + LGA intercept"
  hy2 <- m2$summary.hyperpar
  res2$phi_mean <- hy2["Phi for id", "mean"]; res2$phi_lo <- hy2["Phi for id", "0.025quant"]; res2$phi_hi <- hy2["Phi for id", "0.975quant"]
  res2$dic_poisson <- m0$dic$dic; res2$dic_bym2 <- m2$dic$dic; res2$waic_poisson <- m0$waic$waic; res2$waic_bym2 <- m2$waic$waic; res2$prec_mean <- hy2["Precision for id", "mean"]
  rbind(data.frame(res, dic_poisson = m0$dic$dic, dic_bym2 = m1$dic$dic, waic_poisson = m0$waic$waic, waic_bym2 = m1$waic$waic,
             phi_mean = hy["Phi for id", "mean"], phi_lo = hy["Phi for id", "0.025quant"], phi_hi = hy["Phi for id", "0.975quant"],
             prec_mean = hy["Precision for id", "mean"]), res2[, c("grid","term","RR","lo","hi","model","dic_poisson","dic_bym2","waic_poisson","waic_bym2","phi_mean","phi_lo","phi_hi","prec_mean")])
}
out <- rbind(fit_bym2(cells$g5, "5km"), fit_bym2(cells$g10, "10km"))
write.csv(out, file.path(OUT, "bym2_results.csv"), row.names = FALSE)
