# kassep-geospatial-kano

Complete R script and outputs for the geospatial analysis of maternal mortality in Kano State, Nigeria.

## Contents

- `run_analysis.R` – full reproducible workflow, retrospective cohort only (n=93)
- `outputs/` – all tables and figures generated from the n=93 analysis
- `run_analysis_n116.R` – combined retrospective + prospective workflow (n=116: 93 retrospective + 23 prospective deaths); adds a 10 km hotspot comparison (§8b) not present in the n=93 script
- `outputs_n116/` – all tables and figures generated from the n=116 analysis. Hotspots use Getis-Ord Gi* (own count included, via `include.self()`); earlier versions of the script computed Gi.
- `aisha_option1_household_denominator.R` – area level of the two-level analysis (primary analysis in the revised manuscript). Deaths and surveyed births come from the same households in the KASSEP household listing, so cells without deaths are kept. Poisson model with a log surveyed-births offset and LGA-clustered SEs, at 5 km (primary) and 10 km. Set `BIRTHS=live` to approximate live births (babies delivered minus stillbirths). Retrospective wave only (91 of 93 deaths; the 23 prospective deaths have no household denominator). Outputs: `outputs_n116/aisha_option1/`, `outputs_n116/aisha_option1_livebirths/`.
- `cell_model_diagnostics.R` – overdispersion (Pearson chi-square/df, likelihood-ratio test against negative binomial with halved p) and zero-count checks for the area-level model.
- `aisha_individual_level.R` – individual level: death or survival of each woman aged 15-49 in the household rosters, logistic regression with the survey design (EA cluster, LGA stratum, unweighted because EA sizes are unavailable). Outputs: `outputs_n116/aisha_individual/`.
- `hotspot_retro_vs_combined.R`, `moran_retro_sweep.R` – Gi* hotspots and global Moran's I for the retrospective deaths alone and for all recorded deaths (outputs in `outputs_n116/aisha_hotspots/`). Setting `ALL=1` for `moran_retro_sweep.R` reproduces the all-deaths table in the manuscript.
- `aisha_figures.R` – figures for the two-level analysis (Figures 4 and 5, Supplementary Figure S1).
- `Rplot.pdf` – supplementary plots

Individual-level GPS coordinates are never committed to this repository (see Data Sources below) — only grid-aggregated tables and maps.

## Data Sources

Data are not included due to size, and individual-level verbal autopsy records (including GPS coordinates) are not publicly available due to ethical restrictions. See the script for download instructions:
- KASSEP household listing (GPS, births, deaths, survey strategy; restricted, same ethics terms as the verbal autopsy data)
- KASSEP (verbal autopsy data) — retrospective (n=93) and, for `run_analysis_n116.R`, prospective (n=23) ascertainment components
- WorldPop 2020 (women 15-49) — DOI [10.5258/SOTON/WP00699](https://hub.worldpop.org/geodata/summary?id=50493)
- IDEAMAPS EmOC deprivation
- OpenStreetMap health facilities

## How to Run

Set your working directory and run `source("run_analysis.R")`.

## Citation

Please cite original data sources: KASSEP, WorldPop, IDEAMAPS, and GRID3.

## License

[Choose a license – MIT recommended for code]
