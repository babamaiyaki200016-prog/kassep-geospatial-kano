# ============================================================================
# Aisha (part B): individual-level model, death or survival for each woman.
# Women = women of reproductive age listed in the household rosters of completed
# households (survivors) plus the deceased women (one per recorded death, age from
# Age_maternal_death). Design: EA = cluster, LGA = stratum. KASSEP weights are NOT
# available (EA sizes missing), so the design is unweighted (equal weights).
# Exposures: distance home -> nearest facility; IDEAMAPS category of the 100 m cell
# (covered households only). Covariates available for BOTH survivors and deceased:
# age, settlement type, survey strategy. Education and parity are recorded for survivors
# only in the roster, so they cannot be adjusted without linking the verbal-autopsy file.
# ============================================================================
setwd("/Users/babamusa/Documents/KASSEP GEOSPATIAL/kassep_MMR")
HH_XLSX <- "/Users/babamusa/Downloads/Sharable_Final KASSEP BASELINE DATA.xlsx"
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(sf); library(survey)})
OUT <- "outputs_n116/aisha_individual"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

raw <- suppressMessages(read_excel(HH_XLSX, sheet = 1, .name_repair = "unique_quiet"))
raw <- raw %>% mutate(hh_id = row_number()) %>% filter(`survey_type...654` == 1)
hh <- raw %>% transmute(hh_id, LGA = trimws(LGA), EA, Settlement_type, strategy = `survey_type...10`,
        lat = as.numeric(latitude), lon = as.numeric(longitude),
        death = suppressWarnings(as.numeric(No_maternal_death)), age_death = suppressWarnings(as.numeric(Age_maternal_death)))
hh$LGA[hh$LGA == "Ungoggo"] <- "Ungogo"; hh$LGA[hh$LGA == "Garum Mallam"] <- "Garun Mallam"
hh <- hh %>% filter(!is.na(lat), !is.na(lon), lat > 10, lat < 13, lon > 7, lon < 10, !is.na(EA))
hh$death[is.na(hh$death)] <- 0

# ---- survivors: women of reproductive age in the roster ----------------------
slots <- 0:38
ro <- bind_rows(lapply(slots, function(k) {
  a <- paste0("Member_", k, "_age"); g <- paste0("Member_", k, "_gender"); r <- paste0("Member_", k, "_repro_age")
  data.frame(hh_id = raw$hh_id, age = suppressWarnings(as.numeric(raw[[a]])), sex = raw[[g]], repro = raw[[r]])
}))
surv <- ro %>% filter(sex == "Female", repro == "Yes", !is.na(age), age >= 15, age <= 49) %>%
  inner_join(hh %>% select(hh_id), by = "hh_id") %>% transmute(hh_id, age, died = 0L)
cases <- hh %>% filter(death > 0) %>% select(hh_id, death, age_death) %>%
  tidyr::uncount(death) %>% transmute(hh_id, age = age_death, died = 1L)
cat("Surviving women 15-49 in roster:", nrow(surv), " | deaths:", nrow(cases),
    " | deaths with age:", sum(!is.na(cases$age)), "\n")
cat("Death ages:", paste(range(cases$age, na.rm = TRUE), collapse = "-"), "\n")
w <- bind_rows(surv, cases) %>% filter(!is.na(age), age >= 15, age <= 49) %>% inner_join(hh, by = "hh_id")

# ---- exposures --------------------------------------------------------------
STATE_SHP <- "../NGA_State_Boundaries_V2_-6583428934468491568 (1)/grid3_nga_boundary_vaccstates.shp"
OSM_SHP <- "nga_health_facilities_osm/health_facilities_points.shp"
fac <- st_read(OSM_SHP, quiet = TRUE) %>% st_transform(32632)
hs <- hh %>% distinct(hh_id, lon, lat) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(32632)
hs$dist_km <- as.numeric(st_distance(hs, fac[st_nearest_feature(hs, fac), ], by_element = TRUE)) / 1000
em <- read.csv("../dataset-outputs.csv") %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(32632)
nn <- st_nearest_feature(hs, em); dd <- as.numeric(st_distance(hs, em[nn, ], by_element = TRUE))
hs$depr_cat <- ifelse(dd <= 80, em$result[nn], NA)       # own 100 m cell (half-diagonal ~ 80 m)
hs$depr_cat <- factor(hs$depr_cat, 0:2, c("low", "moderate", "high"))
w <- w %>% left_join(st_drop_geometry(hs), by = "hh_id")
w$agegrp <- cut(w$age, c(14, 19, 24, 29, 34, 39, 49), labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-49"))
w$census <- as.integer(w$strategy == "Census"); w$urban <- as.integer(w$Settlement_type == "Urban")
w$semiurban <- as.integer(w$Settlement_type == "Semi-Urban")
w <- w %>% filter(!is.na(urban))
cat("Women in analysis:", nrow(w), " deaths:", sum(w$died), " (", round(100 * mean(w$died), 2), "% )\n")
cat("With IDEAMAPS 100 m cell:", sum(!is.na(w$depr_cat)), " deaths among them:", sum(w$died[!is.na(w$depr_cat)]), "\n")
print(table(w$depr_cat, w$died, useNA = "ifany"))

# ---- survey-design models (EA cluster, LGA strata, equal weights) -----------
des <- svydesign(ids = ~EA, strata = ~LGA, data = w, nest = TRUE)
options(survey.lonely.psu = "adjust")
tidyo <- function(m, nm, expo = TRUE) { s <- summary(m)$coefficients; ci <- confint(m)
  data.frame(model = nm, term = rownames(s), est = exp(s[, 1]), lo = exp(ci[, 1]), hi = exp(ci[, 2]), p = s[, 4], row.names = NULL) }
m1 <- svyglm(died ~ dist_km + agegrp + urban + semiurban + census, design = des, family = quasibinomial)
m1r <- svyglm(died ~ dist_km + agegrp + urban + semiurban + census, design = des, family = quasipoisson)
res <- rbind(tidyo(m1, "Logistic (OR)"), tidyo(m1r, "Poisson robust (RR)"))
cov_des <- subset(des, !is.na(depr_cat))
if (sum(w$died[!is.na(w$depr_cat)]) >= 10) {
  m2 <- svyglm(died ~ dist_km + depr_cat + agegrp + urban + semiurban + census, design = cov_des, family = quasibinomial)
  res <- rbind(res, tidyo(m2, "Logistic + deprivation category (IDEAMAPS-covered women)"))
}
print(res %>% mutate(across(where(is.numeric), ~round(.x, 3))))
write.csv(res, file.path(OUT, "individual_models.csv"), row.names = FALSE)
write.csv(data.frame(women = nrow(w), deaths = sum(w$died), survivors = sum(w$died == 0),
  covered = sum(!is.na(w$depr_cat)), covered_deaths = sum(w$died[!is.na(w$depr_cat)])), file.path(OUT, "individual_counts.csv"), row.names = FALSE)
