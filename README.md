# ABDZs_impact
This repository contains code and input data used in the submitted paper: Area-based community management of agrobiodiversity zones reduces agricultural expansion and natural cover loss. (Comacchio et al., 2025)


## Data availability & repository contents

This repository includes:

- **Pre-matching analysis dataset**
  - `data/prematch_ABDZ_grid.*` (e.g. GeoPackage or Shapefile)
  - Contains:
    - Unique cell identifier (`ID`)
    - Treatment indicator (`TREATMENT`)
    - Summary outcomes (e.g. mean field size, GPP stats)
    - All covariates used for matching (already extracted & joined)

- **Matched datasets**
  - `data/matched_MFS_4km_*.shp` / `.csv`
  - `data/matched_GPP_*.shp` / `.csv`
  - One file per matching specification (e.g. caliper, Mahalanobis, etc.)
  - These are the exact inputs used by the analysis / spec-curve scripts.

These files are sufficient to:
- reproduce all **statistical analyses** in this repo, and
- inspect balance, matched samples, and estimates.

### Not included (by design)

The repository does **not** include:

- Raw layers (e.g. MapBiomas rasters, MODIS GPP, all environmental covariates processed and used for PCA analysis  
- Intermediate QGIS projects and scratch layers
- The full stack of rasters and vectors used to build the pre-matching dataset (temperature, precipitation elevation, slop, distance to road, accessibility, population density, income, etc.)

These are large, sometimes license-restricted, and available from public sources:

- **GPP:** MODIS/061/MOD17A2HGF (via Google Earth Engine)
- **Land cover:** MapBiomas Peru collections
- **Population, income, roads, etc.:** from standard global / national datasets - see covariates table S2 in the publication

We provide:

1. **Google Earth Engine scripts** (`gee/`) that show how polygon-level GPP was computed.
2. **R scripts** (`code/`) documenting all matching and analysis steps.
3. The **pre-matching dataset** and **post-matching dataset** , so users can rerun the matching & analysis *without* redoing all GIS work.

Researchers who wish to rebuild the pre-matching dataset from raw sources can follow the documented steps in:
- `gee/gpp_cropland_extraction_polygons.js`
- `code/08_landcover_mode_mapbiomas_2000_2022.R`
- `code/09_landcover_grouping_to_5classes.R`
and the method description in the publication.


