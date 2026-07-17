# kassep-geospatial-kano

Complete R script and outputs for the geospatial analysis of maternal mortality in Kano State, Nigeria.

## Contents

- `run_analysis.R` – full reproducible workflow, retrospective cohort only (n=93)
- `outputs/` – all tables and figures generated from the n=93 analysis
- `run_analysis_n116.R` – combined retrospective + prospective workflow (n=116: 93 retrospective + 23 prospective deaths); adds a 10 km hotspot comparison (§8b) not present in the n=93 script
- `outputs_n116/` – all tables and figures generated from the n=116 analysis
- `Rplot.pdf` – supplementary plots

Individual-level GPS coordinates are never committed to this repository (see Data Sources below) — only grid-aggregated tables and maps.

## Data Sources

Data are not included due to size, and individual-level verbal autopsy records (including GPS coordinates) are not publicly available due to ethical restrictions. See the script for download instructions:
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
