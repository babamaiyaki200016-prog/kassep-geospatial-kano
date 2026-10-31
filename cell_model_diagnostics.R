# Aisha: Poisson first; test overdispersion (Pearson chi2/df; LRT Poisson vs NB, p halved);
# check zeros against the fitted model before considering zero-inflation.
setwd("/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR")
suppressMessages({library(MASS)})
for (km in c("5km","10km")) {
  d <- read.csv(sprintf("outputs_n116/aisha_option1/cell_data_%s.csv", km))
  f <- deaths ~ dist_km + census_share + offset(log(births))
  p <- glm(f, family = poisson, data = d)
  nb <- glm.nb(f, data = d)
  chi <- sum(residuals(p, "pearson")^2); df <- p$df.residual
  lr <- 2 * (as.numeric(logLik(nb)) - as.numeric(logLik(p)))
  cat("\n==", km, ": cells", nrow(d), " deaths", sum(d$deaths), "\n")
  cat("Pearson chi2/df =", round(chi/df, 2), " (chi2", round(chi,1), "df", df, " p =", signif(pchisq(chi, df, lower.tail = FALSE), 3), ")\n")
  cat("NB theta =", round(nb$theta, 2), " LRT =", round(lr, 2), " halved p =", signif(0.5 * pchisq(max(lr,0), 1, lower.tail = FALSE), 3), "\n")
  cat("AIC Poisson", round(AIC(p), 1), " NB", round(AIC(nb), 1), "\n")
  cat("Zeros: observed", sum(d$deaths == 0), " expected under Poisson", round(sum(dpois(0, fitted(p))), 1),
      " under NB", round(sum(dnbinom(0, mu = fitted(nb), size = nb$theta)), 1), "\n")
  cat("Max cell deaths:", max(d$deaths), " cells with >=2 deaths:", sum(d$deaths >= 2), "\n")
}
cat("\nINLA installed:", requireNamespace("INLA", quietly = TRUE), " | brms:", requireNamespace("brms", quietly = TRUE), " | CARBayes:", requireNamespace("CARBayes", quietly = TRUE), "\n")
