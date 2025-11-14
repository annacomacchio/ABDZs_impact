# ABDZs_impact

This repository contains code and input data used in the paper under review:

**Comacchio et al. (2025)** – *Area-based community management of agrobiodiversity zones reduces agricultural expansion and natural cover loss.*

The repository is structured to allow full replication of the **matching and statistical analyses** without requiring users to redo all raw GIS preprocessing.

---

## Repository structure

### 1. Analysis-ready datasets

These files are sufficient to reproduce all analyses in this repo.

#### Pre-matching datasets

- `COVARIATES_ALL_NO_COTAHUASI.zip`
  - Contains the pre-matching grid (e.g. `COVARIATES_ALL_NO_COTAHUASI.*`).
  - Includes:
    - Unique cell identifier (`ID`)
    - Treatment indicator (`TREATMENT`)
    - Land-cover summaries
    - Environmental covariates used for matching:
      - precipitation, temperature, elevation, slope  
      - distance to road, accessibility
      - population density, income  
      - baseline land cover / agriculture shares

This is the **input** for all land-use matching scripts and permutation/specification-curve scripts below.

- 'GPP_ALL_WITH_COVARIATES_part1.zip'+'GPP_ALL_WITH_COVARIATES_part2.zip'
  - same as above +
    - GPP-based metrics for each year 2001-2021

#### Matched datasets

Each zip contains the shapefile(s)/CSV(s) corresponding to the final matched samples used in the paper.

- `Agri_exp_nat_cover_loss_matched_dataset_.zip`
  - Matched samples for agricultural expansion and natural cover loss analyses.

- `GPP_CV_matched_dataset.zip`
  - Matched samples for GPP interannual variability (CV) analyses.

- `GPP_mean_matched_dataset.zip`
  - Matched samples for mean GPP (2001–2003 vs 2020–2022 / 2019–2021) analyses.

- `Mean_field_size_matched_dataset.zip`
  - Matched samples for mean field size (MFS) analyses (4-km lattice etc.).

These matched datasets are the **exact inputs** consumed by the analysis. 

---

## 2. Code

All core scripts are numbered in the order they are typically run.

### Matching & main analyses

- `01_agri_expansion_natural_cover_loss_matching_spec.R`  
  Matching specifications for agricultural expansion & natural cover loss.

- `02_analysis_agri_expansion_natural_cover_loss.R`  
  Main analysis for agricultural expansion & natural cover loss on matched data.

- `03_GPP_CV_matching_spec.R`  
  Matching specifications for GPP CV (interannual variability).

- `04_analysis_GPP_CV.R`  
  Main GPP CV analysis on matched data.

- `05_gpp_mean_matching_spec.R`  
  Matching specifications for Δ GPP mean between time windows.

- `06_analysis_GPP_mean.R`  
  Main Δ GPP mean analysis on matched data.

- `07_mean_field_size_matching_spec.R`  
  Matching specifications for mean field size (MFS) with 4-km lattice.

- `08_analysis_mfs_mean_field_size.R`  
  Main MFS analysis on matched data.

### Specification-curve / permutations

- `09_Agri_expansion_nat_cover_loss_permutations.R`  
  Specification curve / permutation analysis for agricultural expansion & natural cover loss.

- `10_gpp_cv_permutations.R`  
  Specification curve / permutation analysis for GPP CV.

- `11_GPP_mean_permutations.R`  
  Specification curve / permutation analysis for Δ GPP mean.

- `12_Mean_field_size_permutations.R`  
  Specification curve / permutation analysis for mean field size.

### Extraction & preprocessing helpers

- `13_Mean_field_size_extraction_loop.R`  
  Loop/scripts used to compute mean field size from raster classifications.

- `Averaging_land_cover_and_reclassification.R`  
  Scripts to:
  - compute modal land cover for 2000–2002 and 2020–2022 windows,
  - reclassify MapBiomas classes into aggregated groups used in analysis.

- `GPP_data_export_GEE/`  
  Google Earth Engine export (e.g. polygon-level GPP extraction).  
  Includes logic used to create `GPP_*` features in the pre-matching dataset.

- `PCA.R`  
  Code and outputs for PCA of environmental covariates (used for descriptive analysis; underlying raw rasters not included here).

---

## 3. What you can reproduce with this repo

Using only the files in this repository, you can:

1. Inspect the **analysis-ready pre-matching grid**.
2. Inspect and verify all **matched datasets**.
3. Re-run:
   - Matching outcomes (using the pre-matching dataset),
   - Main outcome models (OLS, sensitivity checks,figures),
   - Specification-curve / permutation exercises,
   - Balance diagnostics and plots.

All statistical results shown in the paper that depend on matching and regression can be regenerated from these materials.

---

## 4. Not included (by design)

The repository does **not** include:

- Raw remote sensing and ancillary layers, such as:
  - MapBiomas Peru land-cover rasters,
  - MODIS GPP time series,
  - Elevation, slope, accessibility, roads,
  - Population density, income surfaces,
  - Other large rasters or country-wide covariate layers.
- Intermediate QGIS project files and scratch layers.
- Full-resolution inputs used to build the pre-matching dataset from raw sources.

---

## 5. How to rebuild from raw data (optional)

If users wish to reconstruct the pre-matching dataset from scratch, they can follow:

- **Google Earth Engine scripts** in `GPP_data_export_GEE/`  
  for cropland masking and GPP aggregation.
- **Land-cover processing scripts** in `Averaging_land_cover_and_reclassification/`  
  for temporal averaging and reclassification.
- The methodological description in the paper (and supplements).

These steps will reproduce the covariates and outcomes that are already bundled in:

- `COVARIATES_ALL_NO_COTAHUASI.*` (pre-matching grid for land use outcomes), 

---


- `PCA.R`  
  Code and outputs for PCA of environmental covariates (used for descriptive analysis; underlying raw rasters not included here).  
  The script expects the following **input layers**, which are **not** bundled in this repository:

  - Environmental rasters (1 km):
    - `PCA/processed/bio2_1km.tif` (temperature)
    - `PCA/processed/bio12_1km.tif` (precipitation)
    - `PCA/processed/radiation_1km.tif` (solar radiation)
    - `PCA/processed/elevation_1km.tif` (elevation)
    - `PCA/processed/slope_1km.tif` (slope)
    - `PCA/processed/swb_1km.tif` (soil water balance)
    - `PCA/processed/humidity_1km.tif` (humidity)
    - `PCA/processed/gsl_1km.tif` (growing season length)

  These layers are omitted due to size/licensing constraints but can be reconstructed from the CHELSA products (CHELSA climatologies & bioclim variables: https://chelsa-climate.org/datasets/chelsa_climatologies and https://chelsa-climate.org/datasets/chelsa_bioclim), except for slope and elevation 


## 6. Contact

For questions about replication or data access, please contact:

- **Repository owner:** `annacomacchio`
- Or use the issue tracker in this repository.
