# =========================================================
# Script name: 01_matching_spec_agri_expansion_natural_cover_loss.R
# Project:     ABDZs_impact
#
# Purpose:
#   Construct a stable lattice of grid cells, select 1 cell per lattice
#   cell (prioritizing ABDZ core and treated cells), perform matching
#   on environmental and socio-economic covariates, and export matched
#   samples. Additionally, compute and visualize the semivariogram of
#   agricultural expansion for matched units.
#
# Inputs (expected in working directory): 
#   - COVARIATES_ALL_NO_COTAHUASI.shp
#   - OUTER_BUFFER_10km.shp
#
# Outputs:
#   - matched_caliper_4km_NO_COTA_NO_10KM.shp
#   - matched_nearest_4km_NO_COTA_NO_10KM.shp
#   - matched_caliperR_4km_NO_COTA_NO_10KM.shp
#   - matched_mahalanobis_4km_NO_COTA_NO_10KM.shp
#
# Notes:
#   - Coordinates are projected to EPSG:32718 (meters).
#   - spacing_m controls lattice spacing (currently 2000 m).
#   - Update file paths or CRS here if applying to other regions.
# =========================================================

# ---- 0. Packages & global options ----

library(sf)
library(dplyr)
library(MatchIt)
library(cobalt)
library(ggplot2)
library(gstat)

set.seed(42)

# ---- 1. Parameters ----

epsg_utm           <- 32718   # projected CRS in meters
pixel_size_m       <- 1000    # base pixel: 1 km
pixel_area_full    <- pixel_size_m^2
min_pixel_fullness <- 1.0     # require full 1x1 km pixels
spacing_m          <- 2000    # lattice spacing in meters (adjust as needed)
cal_sd             <- 0.2     # matching caliper in SD units (std.caliper = TRUE)

# ---- 2. Helper functions ----

# Convert various logical/label encodings to 0/1 integer
to_int01 <- function(x) {
  if (is.factor(x))   x <- as.character(x)
  if (is.character(x)) {
    xl <- tolower(trimws(x))
    x <- ifelse(
      xl %in% c("1", "true", "t", "yes", "y", "treated"), 1L,
      ifelse(
        xl %in% c("0", "false", "f", "no", "n", "control"), 0L,
        NA_integer_
      )
    )
  }
  if (is.logical(x)) x <- as.integer(x)
  if (is.numeric(x)) x <- as.integer(x != 0 & !is.na(x))
  as.integer(x)
}

as_df_no_geom <- function(x) {
  if (inherits(x, "sf")) sf::st_drop_geometry(x) else x
}

practical_range <- function(model) {
  m <- model[model$model != "Nug", ][1, ]
  switch(
    m$model,
    "Sph" = m$range,
    "Exp" = 3 * m$range,
    "Gau" = sqrt(3) * m$range,
    "Mat" = 3 * m$range,
    m$range
  )
}

# ---- 3. Load covariates & remove 10-km buffered area ----

all_gridcells <- st_read("COVARIATES_ALL_NO_COTAHUASI.shp") %>%
  st_transform(epsg_utm)

if (!"ID" %in% names(all_gridcells)) {
  all_gridcells$ID <- seq_len(nrow(all_gridcells))
}

outer10   <- st_read("OUTER_BUFFER_10km.shp")
outer10_u <- st_union(outer10)

n0 <- nrow(all_gridcells)

# Keep only cells that do NOT intersect the outer 10-km buffer
keep_idx      <- lengths(st_intersects(all_gridcells, outer10_u)) == 0L
all_gridcells <- all_gridcells[keep_idx, , drop = FALSE]

message(
  "Removed ", n0 - nrow(all_gridcells),
  " cells intersecting OUTER_BUFFER_10km; kept ", nrow(all_gridcells)
)

# Sanity check: confirm no intersections remain
stopifnot(all(lengths(st_intersects(all_gridcells, outer10_u)) == 0L))

# ---- 4. Standardize variable names & treatment indicators ----

all_gridcells <- all_gridcells %>%
  rename(
    precipitation      = PRECIPITAT,
    temperature        = TEMP_MEAN,
    population_density = POP_TOTAL,
    accessibility      = ACCESSIBIL,
    elevation          = ELEVATIONM,
    distance_to_road   = Hub.distan,
    slope              = SLOPEMEAN,
    income             = ing_pc_,
    Forest_2000        = FOREST_200,
    Grassland_2000     = GRASSLAND_,
    Agriculture_2000   = AGRICULTUR,
    agri_expansion     = AGRICULTU3
  ) %>%
  mutate(
    TREATMENT  = to_int01(TREATMENT),
    ABDZS_CORE = to_int01(ABDZS_CORE),
    # Ensure all ABDZS core cells are treated
    TREATMENT  = ifelse(ABDZS_CORE == 1L, 1L, TREATMENT) %>% as.integer()
  )

# ---- 5. Keep only full 1x1 km pixels ----

all_gridcells <- all_gridcells %>%
  mutate(
    pixel_area_m2  = as.numeric(st_area(geometry)),
    pixel_fullness = pmin(pixel_area_m2 / pixel_area_full, 1)
  )

message("Pixels before fullness filter: ", nrow(all_gridcells))

all_gridcells <- all_gridcells %>%
  filter(pixel_fullness >= min_pixel_fullness)

message("Pixels kept (fullness >= ", min_pixel_fullness, "): ", nrow(all_gridcells))

# ---- 6. Build lattice & assign lattice cell IDs ----

centroids <- st_centroid(all_gridcells)
coords    <- st_coordinates(centroids)

all_gridcells <- all_gridcells %>%
  mutate(
    x_coord = coords[, 1],
    y_coord = coords[, 2]
  )

x0   <- min(all_gridcells$x_coord, na.rm = TRUE)
y0   <- min(all_gridcells$y_coord, na.rm = TRUE)
x_off <- x0 %% spacing_m
y_off <- y0 %% spacing_m

all_gridcells <- all_gridcells %>%
  mutate(
    cell_x = floor((x_coord - x_off) / spacing_m),
    cell_y = floor((y_coord - y_off) / spacing_m)
  )

# ---- 7. Global thinning: 1 unit per lattice cell ----
# Priority: ABDZS core (ABDZS_CORE==1) > treated (TREATMENT==1) > nearest to cell center

sample_global_pref_core_treated <- function(df) {
  df %>%
    group_by(cell_x, cell_y) %>%
    mutate(
      cx    = cell_x * spacing_m + spacing_m / 2 + x_off,
      cy    = cell_y * spacing_m + spacing_m / 2 + y_off,
      d2    = (x_coord - cx)^2 + (y_coord - cy)^2,
      p_core = if_else(ABDZS_CORE == 1L, 0L, 1L),
      p_trt  = if_else(TREATMENT  == 1L, 0L, 1L)
    ) %>%
    arrange(p_core, p_trt, d2, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-cx, -cy, -d2, -p_core, -p_trt)
}

sf_sampled <- sample_global_pref_core_treated(all_gridcells)

stopifnot(
  !any(duplicated(
    sf::st_drop_geometry(sf_sampled)[, c("cell_x", "cell_y")]
  ))
)

# ---- 8. Matching setup ----

sampled_gridcells_clean <- sf_sampled
geo2 <- sampled_gridcells_clean %>%
  dplyr::select(ID, geometry)

covariates <- c(
  "precipitation", "temperature", "population_density",
  "elevation", "slope", "distance_to_road",
  "accessibility", "income",
  "Forest_2000", "Grassland_2000", "Agriculture_2000"
)

# Ensure covariates are numeric
sampled_gridcells_clean <- sampled_gridcells_clean %>%
  mutate(across(all_of(covariates), ~ suppressWarnings(as.numeric(.))))

data_match <- sampled_gridcells_clean %>%
  filter(if_all(all_of(covariates), ~ is.finite(.) & !is.na(.))) %>%
  mutate(TREATMENT = as.integer(TREATMENT))

covariate_formula <- as.formula(
  paste("TREATMENT ~", paste(covariates, collapse = " + "))
)

# ---- 9. Matching specifications ----

m_nn <- MatchIt::matchit(
  covariate_formula, data = data_match,
  method = "nearest"
)

m_cal <- MatchIt::matchit(
  covariate_formula, data = data_match,
  method = "nearest",
  caliper = cal_sd, std.caliper = TRUE
)

m_calr <- MatchIt::matchit(
  covariate_formula, data = data_match,
  method = "nearest",
  caliper = cal_sd, std.caliper = TRUE,
  replace = TRUE
)

m_mah <- MatchIt::matchit(
  covariate_formula, data = data_match,
  method = "nearest",
  distance = "mahalanobis"
)

d_nn   <- MatchIt::match.data(m_nn)
d_cal  <- MatchIt::match.data(m_cal)
d_calr <- MatchIt::match.data(m_calr)
d_mah  <- MatchIt::match.data(m_mah)

# ---- 10. Balance diagnostics (Love plots & summaries) ----

cobalt::love.plot(
  m_nn,
  title = "Nearest Neighbor",
  thresholds = c(m = 0.1)
)

cobalt::love.plot(
  m_cal,
  title = "Caliper 0.2 SD (No Replacement)",
  thresholds = c(m = 0.1)
)

cobalt::love.plot(
  m_calr,
  title = "Caliper 0.2 SD (With Replacement)",
  thresholds = c(m = 0.1)
)

cobalt::love.plot(
  m_mah,
  title = "Mahalanobis",
  thresholds = c(m = 0.1)
)

print(summary(m_nn))
print(summary(m_cal))
print(summary(m_calr))
print(summary(m_mah))

# ---- 11. Join matched units to geometry ----

stopifnot(inherits(geo2, "sf"))
stopifnot(anyDuplicated(geo2$ID) == 0)

sf_nn   <- geo2 %>% dplyr::inner_join(as_df_no_geom(d_nn),   by = "ID")
sf_cal  <- geo2 %>% dplyr::inner_join(as_df_no_geom(d_cal),  by = "ID")
sf_calr <- geo2 %>% dplyr::inner_join(as_df_no_geom(d_calr), by = "ID")
sf_mah  <- geo2 %>% dplyr::inner_join(as_df_no_geom(d_mah),  by = "ID")

stopifnot(setequal(unique(d_cal$ID), unique(sf_cal$ID)))

# ---- 12. Save matched samples ----

st_write(sf_cal,  "matched_caliper_4km_NO_COTA_NO_10KM.shp",    delete_dsn = TRUE)
st_write(sf_nn,   "matched_nearest_4km_NO_COTA_NO_10KM.shp",    delete_dsn = TRUE)
st_write(sf_calr, "matched_caliperR_4km_NO_COTA_NO_10KM.shp",   delete_dsn = TRUE)
st_write(sf_mah,  "matched_mahalanobis_4km_NO_COTA_NO_10KM.shp", delete_dsn = TRUE)
