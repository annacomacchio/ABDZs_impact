# =========================================================
# Script name: 05_gpp_mean_matching_spec.R
# Project:     ABDZs_impact
#
# Purpose:
#   - Load GPP + covariates shapefile
#   - Filter out outer 10-km buffer
#   - Keep only full 1 km² pixels
#   - Apply 2-km lattice thinning
#   - For GPP means:
#       * W1 = 2001–2003
#       * W2 = 2020–2022
#     require FULL coverage in BOTH windows (all years present)
#   - Construct:
#       * mean_w1, mean_w2
#       * delta_mean = mean_w2 - mean_w1
#   - Run matching on standard baseline covariates
#     (NO mean_w1 / mean_w2 / delta_mean in covariates)
#   - Export matched shapefiles for downstream analysis
#
# Inputs:
#   - GPP_exports/GPP_ALL_WITH_COVARIATES.shp
#   - OUTER_BUFFER_10km.shp
#
# Outputs (to GPP_exports/):
#   - GPPmean_matched_NN.shp
#   - GPPmean_matched_CAL.shp   (caliper 0.2 SD, no replacement)
#   - GPPmean_matched_MAH.shp
#
# Notes:
#   - This script ONLY does preparation + matching for mean GPP.
#   - Analysis script (e.g. OLS on delta_mean) should be separate.
# =========================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(MatchIt)
  library(cobalt)
  library(tibble)
})

set.seed(42)

# -----------------------------
# Parameters
# -----------------------------
in_shp            <- "GPP_exports/GPP_ALL_WITH_COVARIATES.shp"
outer_buffer_path <- "OUTER_BUFFER_10km.shp"

out_dir           <- "GPP_exports"
out_mat_nn_shp    <- file.path(out_dir, "GPPmean_matched_NN.shp")
out_mat_cal_shp   <- file.path(out_dir, "GPPmean_matched_CAL.shp")
out_mat_mah_shp   <- file.path(out_dir, "GPPmean_matched_MAH.shp")

epsg_utm          <- 32718
pixel_size_m      <- 1000
pixel_area_full   <- pixel_size_m^2
min_pixel_fullness<- 1.0      # keep only full 1 km² pixels
spacing_m         <- 2000     # 2-km lattice spacing

years_w1          <- 2001:2003
years_w2          <- 2020:2022

treatment_col     <- "TREATMENT"
cal_sd            <- 0.2      # caliper (std. dev.) for nearest neighbor matching

# -----------------------------
# Helper functions
# -----------------------------
to_int01 <- function(x){
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    xl <- tolower(trimws(x))
    x <- ifelse(
      xl %in% c("1","true","t","yes","y","treated","treatment","abdz"), 1L,
      ifelse(xl %in% c("0","false","f","no","n","control"), 0L, NA_integer_)
    )
  }
  if (is.logical(x)) x <- as.integer(x)
  if (is.numeric(x)) x <- as.integer(x != 0 & !is.na(x))
  as.integer(x)
}

summarize_mean_full <- function(df, cols){
  stopifnot(length(cols) > 0)
  # Require ALL years present (no NAs) for that window
  n_vals <- df %>%
    transmute(n = rowSums(!is.na(across(all_of(cols)))))
  means  <- df %>%
    transmute(mean = rowMeans(across(all_of(cols)), na.rm = FALSE))
  keep_full <- n_vals$n == length(cols)
  means$mean[!keep_full] <- NA_real_
  tibble(n = n_vals$n, mean = means$mean)
}

as_df_no_geom <- function(z) {
  if (inherits(z, "sf")) sf::st_drop_geometry(z) else z
}

# -----------------------------
# 1) Load & project
# -----------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_shp))

x <- sf::st_read(in_shp, quiet = TRUE) |>
  sf::st_transform(epsg_utm)

if (!"ID" %in% names(x)) x$ID <- seq_len(nrow(x))
x$ID <- as.character(x$ID)

x[[treatment_col]] <- to_int01(x[[treatment_col]])

# -----------------------------
# 2) Remove cells intersecting OUTER_BUFFER_10km
# -----------------------------
stopifnot(file.exists(outer_buffer_path))

outer10 <- sf::st_read(outer_buffer_path, quiet = TRUE) |>
  sf::st_transform(epsg_utm) |>
  sf::st_make_valid()

outer10_u <- sf::st_union(outer10)

n0 <- nrow(x)
x  <- x[lengths(sf::st_intersects(x, outer10_u)) == 0L, ]
message(
  "Removed ", n0 - nrow(x),
  " cells intersecting OUTER_BUFFER_10km; kept ", nrow(x), "."
)

# -----------------------------
# 3) Keep ONLY full pixels & thin to 2-km lattice
# -----------------------------
x <- x |>
  mutate(
    pixel_area_m2  = as.numeric(sf::st_area(geometry)),
    pixel_fullness = pmin(pixel_area_m2 / pixel_area_full, 1)
  ) |>
  filter(pixel_fullness >= min_pixel_fullness)

cent <- sf::st_centroid(x)
xy   <- sf::st_coordinates(cent)

x <- x %>%
  mutate(
    x_coord = xy[,1],
    y_coord = xy[,2]
  )

x0    <- min(x$x_coord, na.rm = TRUE)
y0    <- min(x$y_coord, na.rm = TRUE)
x_off <- x0 %% spacing_m
y_off <- y0 %% spacing_m

x <- x %>%
  mutate(
    cell_x = floor((x_coord - x_off) / spacing_m),
    cell_y = floor((y_coord - y_off) / spacing_m)
  )

sample_2km_nearest_center <- function(df) {
  df %>%
    group_by(cell_x, cell_y) %>%
    mutate(
      cx = cell_x * spacing_m + spacing_m/2 + x_off,
      cy = cell_y * spacing_m + spacing_m/2 + y_off,
      d2 = (x_coord - cx)^2 + (y_coord - cy)^2
    ) %>%
    arrange(d2, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-cx, -cy, -d2)
}

sf_sampled <- sample_2km_nearest_center(x)

message("After 2-km thinning: n = ", nrow(sf_sampled))
stopifnot(
  nrow(
    sf_sampled %>%
      count(cell_x, cell_y) %>%
      filter(n > 1)
  ) == 0L
)

# -----------------------------
# 4) Compute means & require FULL coverage in BOTH windows
# -----------------------------
gpp_cols <- names(sf_sampled)[stringr::str_detect(names(sf_sampled), "^GPP_\\d{4}$")]
if (length(gpp_cols) == 0L) {
  stop("No GPP_YYYY columns found in the data.")
}

gpp_years_present <- as.integer(stringr::str_remove(gpp_cols, "^GPP_"))

w1_cols <- paste0("GPP_", intersect(years_w1, gpp_years_present))
w2_cols <- paste0("GPP_", intersect(years_w2, gpp_years_present))

if (length(w1_cols) < length(years_w1) || length(w2_cols) < length(years_w2)) {
  stop(
    "Full-coverage requirement failed: not all requested years exist in the file.\n",
    "Available GPP years: ", paste(sort(gpp_years_present), collapse = ", ")
  )
}

stats_w1 <- summarize_mean_full(sf_sampled, w1_cols)
stats_w2 <- summarize_mean_full(sf_sampled, w2_cols)

ok_idx <- which(is.finite(stats_w1$mean) & is.finite(stats_w2$mean))
if (length(ok_idx) == 0L) {
  stop("No rows have FULL GPP coverage in BOTH windows (2001–2003 AND 2020–2022).")
}

xv <- sf_sampled[ok_idx, ] %>%
  bind_cols(
    tibble(
      n_w1    = stats_w1$n[ok_idx],
      mean_w1 = stats_w1$mean[ok_idx],
      n_w2    = stats_w2$n[ok_idx],
      mean_w2 = stats_w2$mean[ok_idx]
    )
  ) %>%
  mutate(
    delta_mean = mean_w2 - mean_w1
  )

message("Kept after FULL coverage filter (means): n = ", nrow(xv))

# -----------------------------
# 5) Matching setup (rename covariates; NO mean_w1/mean_w2/delta_mean)
# -----------------------------
geo_only <- xv %>% select(ID, geometry)

# Standardize covariate names (only if present)
rename_map <- c(
  "PRECIPITAT" = "precipitation",
  "TEMP_MEAN"  = "temperature",
  "POPDENS_ME" = "population_density",
  "ELEVATIONM" = "elevation",
  "SLOPEMEAN"  = "slope",
  "Hub_distan" = "distance_to_road",
  "ACCESSIBIL" = "accessibility",
  "AGRICULTUR" = "crop_2000s",
  "ing_pc_"    = "income"
)

present_old <- intersect(names(xv), names(rename_map))
xv <- dplyr::rename(xv, !!!setNames(rename_map[present_old], present_old))

covariates <- c(
  "precipitation","temperature","population_density",
  "elevation","income",
  "slope","distance_to_road","accessibility","crop_2000s"
)
covariates <- intersect(covariates, names(xv))

if (length(covariates) == 0L) {
  stop("No matching covariate columns found after renaming.")
}

# Build dataset for MatchIt
data_match <- xv %>%
  mutate(
    !!treatment_col := to_int01(.data[[treatment_col]]),
    across(all_of(covariates), ~ suppressWarnings(as.numeric(.)))
  ) %>%
  filter(.data[[treatment_col]] %in% c(0L, 1L)) %>%
  filter(if_all(all_of(covariates), ~ is.finite(.) & !is.na(.)))

if (nrow(data_match) == 0L) {
  stop("No rows left for matching after covariate cleaning.")
}

covariate_formula <- stats::as.formula(
  paste(treatment_col, "~", paste(covariates, collapse = " + "))
)

# -----------------------------
# 6) Run MatchIt specifications
# -----------------------------
m_nn <- MatchIt::matchit(
  covariate_formula,
  data   = data_match,
  method = "nearest"
)

m_cal <- MatchIt::matchit(
  covariate_formula,
  data        = data_match,
  method      = "nearest",
  caliper     = cal_sd,
  std.caliper = TRUE,
  replace     = FALSE
)

m_mah <- MatchIt::matchit(
  covariate_formula,
  data     = data_match,
  method   = "nearest",
  distance = "mahalanobis"
)

d_nn  <- MatchIt::match.data(m_nn)
d_cal <- MatchIt::match.data(m_cal)
d_mah <- MatchIt::match.data(m_mah)

# -----------------------------
# 7) Balance diagnostics (Love plots & summaries)
# -----------------------------
cobalt::love.plot(
  m_nn,
  title      = "GPP mean: Nearest Neighbor",
  thresholds = c(m = 0.1)
)

cobalt::love.plot(
  m_cal,
  title      = "GPP mean: Caliper 0.2 SD (no replacement)",
  thresholds = c(m = 0.1)
)

cobalt::love.plot(
  m_mah,
  title      = "GPP mean: Mahalanobis",
  thresholds = c(m = 0.1)
)

cat("\n--- Matching summaries (GPP mean) ---\n")
print(summary(m_nn))
print(summary(m_cal))
print(summary(m_mah))

# -----------------------------
# 8) Join matched sets back to geometry & write shapefiles
# -----------------------------
sf_nn  <- geo_only %>% inner_join(as_df_no_geom(d_nn),  by = "ID")
sf_cal <- geo_only %>% inner_join(as_df_no_geom(d_cal), by = "ID")
sf_mah <- geo_only %>% inner_join(as_df_no_geom(d_mah), by = "ID")

sf::st_write(sf_nn,  out_mat_nn_shp,  delete_layer = TRUE, quiet = TRUE)
sf::st_write(sf_cal, out_mat_cal_shp, delete_layer = TRUE, quiet = TRUE)
sf::st_write(sf_mah, out_mat_mah_shp, delete_layer = TRUE, quiet = TRUE)

message("✅ Wrote matched shapefiles (GPP mean, full coverage, 2-km grid):")
message("  - ", out_mat_nn_shp,  "  [Nearest neighbor]")
message("  - ", out_mat_cal_shp, "  [Caliper 0.2 SD]")
message("  - ", out_mat_mah_shp, "  [Mahalanobis]")
