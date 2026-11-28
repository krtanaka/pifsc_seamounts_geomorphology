# ==============================================================================
# process_knolls_seamounts.R
#
# Builds global distance-to-seamount rasters (NetCDF) for use as predictors in
# Pacific longline CPUE / catch models, plus diagnostic maps.
#
# Workflow
#   0. USER OPTIONS      - data paths, filter scenarios, compute/plot switches
#   1. Helpers           - depth standardisation, filtering, file naming
#   2. Load features     - Yesson et al. (2011) seamounts + knolls
#   3. Feature maps      - optional exploratory maps / hex plots (unfiltered)
#   4. Distance rasters  - one NetCDF per scenario x height class
#   5. Distance maps     - one PNG per NetCDF
#
# Re-running: edit the OPTIONS block only. Existing NetCDFs are skipped unless
# `recalculate = "all"`. Every output file name encodes the filters used, so
# different scenarios never overwrite each other. See README.md for the
# literature rationale behind the default filters.
# ==============================================================================

library(sf)
library(maps)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(terra)
library(parallel)
library(pbapply)

# ==============================================================================
# 0. USER OPTIONS
# ==============================================================================

base_path  <- "A:/My Drive/Data/Seamounts & Knolls/Modelled_Seamounts_2011_v1/01_Data/"
base_path  <- "D:/Environmental_Data_Summary/deepsea/data_download/Seamounts_Knolls/Modelled_Seamounts_2011_v1/01_Data/"

nc_dir     <- "."        # where NetCDF rasters are written (repo root, as before)
output_dir <- "output"   # where PNGs and the manifest CSV are written

# Filter scenarios. Each produces its own NetCDF(s). Fields not given fall back
# to `scenario_defaults`. Comment out scenarios you don't need.
#   include_seamounts / include_knolls : which Yesson et al. layers to use
#   max_summit_depth_m : keep features with summit depth <= this (m, positive
#                        down). NULL = no depth filter.
#   min_height_m       : keep features with HEIGHT >= this (m). NULL = none.
#   height_tertiles    : TRUE also writes Small/Medium/Big rasters, split by
#                        HEIGHT tertiles *within* the filtered set.
scenario_defaults <- list(
  include_seamounts  = TRUE,
  include_knolls     = TRUE,
  max_summit_depth_m = NULL,
  min_height_m       = NULL,
  height_tertiles    = FALSE
)

scenarios <- list(
  # Primary predictor: shallow summits, within reach of epipelagic/mesopelagic
  # habitat and longline gear (colleague's 400 m threshold; see README).
  primary_shallow = list(max_summit_depth_m = 400),

  # Sensitivity: deeper summits (proposed, not a published cutoff).
  sensitivity_1000m = list(max_summit_depth_m = 1000),

  # Baseline: every seamount and knoll, no depth filter.
  baseline_all = list()

  # Legacy run (previous version of this script: seamounts only, HEIGHT >= 2500,
  # with height tertiles). Uncomment to reproduce.
  # legacy_tall = list(include_knolls = FALSE, min_height_m = 2500, height_tertiles = TRUE)
)

# How the shapefile DEPTH field is signed. "auto" infers from the median
# (negative median => negative-down elevation). Override with
# "negative_down" (e.g. -350 = 350 m deep) or "positive_down" (350 = 350 m deep).
depth_sign_convention <- "auto"

# Drop duplicate summit locations (identical coordinates to ~1 m) before
# computing distances.
drop_duplicate_locations <- TRUE

# Compute switches
#   "missing" : compute only NetCDFs that do not exist yet (default)
#   "all"     : recompute and overwrite every scenario NetCDF
#   "none"    : never compute; only plot files that already exist
recalculate <- "missing"
res_km      <- 5                                         # grid cell size (km at equator)
n_chunks    <- 60                                        # latitude strips for parallel work
# Respect cloud/container CPU quotas; reserve two cores, with at least one worker.
n_cores     <- min(n_chunks, parallelly::availableCores(omit = 2))

# Plot switches
make_feature_maps  <- FALSE  # exploratory maps/hexplots of the raw feature layers
make_distance_maps <- TRUE   # one log10 distance map per NetCDF

# Figure style
dark_background <- "black"
dark_panel      <- "black"
dark_text       <- "#e6edf3"
muted_text      <- "#a9b8c9"
slide_width     <- 16
slide_height    <- 9

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(nc_dir)) dir.create(nc_dir, recursive = TRUE)
stopifnot(recalculate %in% c("missing", "all", "none"))

metadata_caption <- paste0(
  "Data Provider: Zoological Society of London (2011) | Version: 1.0 | License: CC BY 3.0\n",
  "Citation: Yesson C, Clark MR, Taylor M, Rogers AD (2011). The global distribution of seamounts based on 30-second bathymetry data.\n",
  "Deep Sea Research Part I 58: 442-453. doi: 10.1016/j.dsr.2011.02.004 | Data URL: http://data.unep-wcmc.org/datasets/41"
)

# ==============================================================================
# 1. HELPERS
# ==============================================================================

# Converts the raw DEPTH field to positive-down summit depth in metres.
standardise_depth <- function(depth, convention = "auto") {
  if (convention == "auto") {
    convention <- if (median(depth, na.rm = TRUE) < 0) "negative_down" else "positive_down"
    message(sprintf("DEPTH sign convention inferred as '%s' (median raw DEPTH = %.0f).",
                    convention, median(depth, na.rm = TRUE)))
  }
  switch(convention,
         negative_down = -depth,
         positive_down = depth,
         stop("depth_sign_convention must be 'auto', 'negative_down' or 'positive_down'"))
}

resolve_scenario <- function(s) modifyList(scenario_defaults, s, keep.null = TRUE)

# Builds the file-name tag that encodes every active filter, e.g.
#   SeamountsKnolls_DepthLE400m_HgtAll_All_5km
scenario_tag <- function(s, size_class = "All") {
  feat <- paste0(if (s$include_seamounts) "Seamounts" else "",
                 if (s$include_knolls) "Knolls" else "")
  depth <- if (is.null(s$max_summit_depth_m)) "DepthAll" else sprintf("DepthLE%gm", s$max_summit_depth_m)
  hgt   <- if (is.null(s$min_height_m)) "HgtAll" else sprintf("HgtGE%gm", s$min_height_m)
  paste(feat, depth, hgt, size_class, sprintf("%gkm", res_km), sep = "_")
}

scenario_description <- function(s, size_class = "All") {
  feat <- paste(c(if (s$include_seamounts) "seamounts", if (s$include_knolls) "knolls"), collapse = " + ")
  depth <- if (is.null(s$max_summit_depth_m)) "any summit depth" else sprintf("summit depth <= %g m", s$max_summit_depth_m)
  hgt   <- if (is.null(s$min_height_m)) "any height" else sprintf("height >= %g m", s$min_height_m)
  cls   <- if (size_class == "All") "" else sprintf("; %s height tertile", tolower(size_class))
  sprintf("%s; %s; %s%s", feat, depth, hgt, cls)
}

filter_features <- function(features, s) {
  types <- c(if (s$include_seamounts) "Seamount", if (s$include_knolls) "Knoll")
  if (length(types) == 0) stop("Scenario must include seamounts and/or knolls.")
  out <- features %>% filter(feature_type %in% types)
  if (!is.null(s$max_summit_depth_m)) {
    out <- out %>% filter(!is.na(summit_depth_m), summit_depth_m <= s$max_summit_depth_m)
  }
  if (!is.null(s$min_height_m)) {
    out <- out %>% filter(!is.na(HEIGHT), HEIGHT >= s$min_height_m)
  }
  if (nrow(out) == 0) stop("Scenario filters removed every feature: ", scenario_description(s))
  out
}

# Returns a named list of point sets: "All", plus Small/Medium/Big if requested.
split_size_classes <- function(pts, s) {
  classes <- list(All = pts)
  if (isTRUE(s$height_tertiles)) {
    pts <- pts %>% mutate(size_group = ntile(HEIGHT, 3))
    classes$Small  <- pts %>% filter(size_group == 1)
    classes$Medium <- pts %>% filter(size_group == 2)
    classes$Big    <- pts %>% filter(size_group == 3)
  }
  classes
}

# Runs on each worker: geodesic distance (km) from every cell in one latitude
# strip to the nearest point. All points are passed so the nearest feature is
# found even when it lies in another strip.
process_strip <- function(i, xy, lat_breaks, res_deg) {
  strip_rast <- terra::rast(terra::ext(0, 360, lat_breaks[i], lat_breaks[i + 1]),
                            resolution = res_deg, crs = "EPSG:4326")
  v_pts <- terra::vect(xy, crs = "EPSG:4326")
  d_km <- terra::distance(strip_rast, v_pts) / 1000
  names(d_km) <- "Distance_km"
  tmp_file <- tempfile(fileext = ".tif")
  terra::writeRaster(d_km, tmp_file, overwrite = TRUE)
  tmp_file
}

compute_distance_nc <- function(pts, out_file, longname, atts, cl) {
  res_deg    <- res_km / 111.12
  lat_breaks <- seq(-90, 90, length.out = n_chunks + 1)
  xy <- st_coordinates(pts)[, c("X", "Y"), drop = FALSE]

  tif_files <- pblapply(seq_len(n_chunks), process_strip,
                        xy = xy, lat_breaks = lat_breaks, res_deg = res_deg, cl = cl)
  on.exit(unlink(unlist(tif_files)), add = TRUE)

  final_raster <- mosaic(sprc(unlist(tif_files)))
  writeCDF(final_raster, out_file, overwrite = TRUE,
           varname = "distance", longname = longname, unit = "km", atts = atts)
}

# ==============================================================================
# 2. LOAD FEATURES (Yesson et al. 2011)
# ==============================================================================

read_layer <- function(path, type) {
  st_read(path, quiet = TRUE) %>%
    select(HEIGHT, DEPTH) %>%
    mutate(feature_type = type)
}

features <- rbind(
  read_layer(file.path(base_path, "Seamounts/Seamounts.shp"), "Seamount"),
  read_layer(file.path(base_path, "Knolls/Knolls.shp"), "Knoll")
) %>%
  st_transform(4326) %>%
  st_shift_longitude() %>%   # 0-360 longitudes, Pacific-centred
  mutate(summit_depth_m = standardise_depth(DEPTH, depth_sign_convention))

# Features whose summit is above sea level are islands, not seamounts
# (Morato et al. 2010 definition: does not break the surface).
n_above <- sum(features$summit_depth_m < 0, na.rm = TRUE)
if (n_above > 0) message(sprintf("Dropping %d features with summit above sea level.", n_above))
features <- features %>% filter(is.na(summit_depth_m) | summit_depth_m >= 0)

if (drop_duplicate_locations) {
  key <- apply(round(st_coordinates(features), 5), 1, paste, collapse = "_")
  n_dup <- sum(duplicated(key))
  if (n_dup > 0) message(sprintf("Dropping %d duplicate summit locations.", n_dup))
  features <- features[!duplicated(key), ]
}

print(features %>% st_drop_geometry() %>% group_by(feature_type) %>%
        summarise(n = n(), n_missing_depth = sum(is.na(summit_depth_m)),
                  n_le_400m = sum(summit_depth_m <= 400, na.rm = TRUE),
                  n_le_1000m = sum(summit_depth_m <= 1000, na.rm = TRUE)))

# ==============================================================================
# 3. EXPLORATORY FEATURE MAPS (unfiltered layers)
# ==============================================================================

map_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = dark_background, color = NA),
    panel.background = element_rect(fill = dark_background, color = NA),
    plot.title = element_text(face = "bold", size = 18, color = dark_text),
    plot.caption = element_text(size = 9, color = muted_text, hjust = 0, lineheight = 1.1),
    axis.title = element_text(size = 12, color = dark_text),
    axis.text = element_text(size = 10, color = muted_text),
    legend.background = element_rect(fill = dark_panel, color = "#415166", linewidth = 0.3),
    legend.key = element_blank(),
    legend.title = element_text(size = 12, face = "bold", color = dark_text),
    legend.text = element_text(size = 10, color = dark_text),
    legend.key.size = unit(0.55, "cm")
  )

if (make_feature_maps) {
  sf_use_s2(FALSE)
  world_base <- st_as_sf(map("world2", plot = FALSE, fill = TRUE)) %>% st_make_valid()
  st_crs(world_base) <- 4326
  world_cropped <- st_crop(world_base, c(xmin = 110, ymin = -50, xmax = 300, ymax = 70))
  sf_use_s2(TRUE)

  pacific_view <- coord_sf(crs = "+proj=robin +lon_0=180",
                           xlim = c(120, 290), ylim = c(-45, 60),
                           default_crs = 4326, expand = FALSE)
  spatial_theme <- map_theme + theme(panel.grid = element_blank(), legend.position = c(0.86, 0.85))
  hex_theme     <- map_theme + theme(panel.grid = element_line(color = "#263445"), legend.position = c(0.12, 0.78))

  depth_thresholds <- unique(unlist(lapply(scenarios, function(s) resolve_scenario(s)$max_summit_depth_m)))

  for (type in c("Seamount", "Knoll")) {
    pts <- features %>% filter(feature_type == type)
    lower <- tolower(type)

    for (var in c("HEIGHT", "summit_depth_m")) {
      is_height <- var == "HEIGHT"
      p <- ggplot() +
        geom_sf(data = world_cropped, fill = "#202b38", color = "#526170", linewidth = 0.1) +
        geom_sf(data = pts, aes(color = .data[[var]]), size = if (type == "Seamount") 1.2 else 0.5, alpha = 0.8) +
        scale_color_viridis_c(option = if (is_height) "plasma" else "mako",
                              name = if (is_height) "Height (m)" else "Summit depth (m)",
                              trans = if (is_height) "log10" else "identity",
                              direction = if (is_height) 1 else -1) +
        pacific_view + spatial_theme +
        labs(title = sprintf("Pacific %ss by %s", type, if (is_height) "Height" else "Summit Depth"),
             x = "Longitude", y = "Latitude", caption = metadata_caption)
      ggsave(file.path(output_dir, sprintf("%ss_by_%s.png", lower, if (is_height) "height" else "depth")),
             p, width = slide_width, height = slide_height, dpi = 300, bg = dark_background)
    }

    # Height vs summit depth; dashed lines mark the depth thresholds in use.
    p <- ggplot(pts, aes(x = summit_depth_m, y = HEIGHT)) +
      geom_hex(bins = 50, alpha = 0.9) +
      geom_vline(xintercept = depth_thresholds, linetype = "dashed", color = "#ffcc00") +
      scale_fill_viridis_c(option = "plasma", name = "Count") +
      hex_theme +
      labs(title = sprintf("%ss: Morphological Density Structure", type),
           x = "Summit depth (m)", y = "Height (m)", caption = metadata_caption)
    ggsave(file.path(output_dir, sprintf("%ss_morphology_hex.png", lower)),
           p, width = slide_width, height = slide_height, dpi = 300, bg = dark_background)
  }
}

# ==============================================================================
# 4. DISTANCE RASTERS
# ==============================================================================

# Enumerate every (scenario x height class) job and its output file.
jobs <- list()
for (sc_name in names(scenarios)) {
  s <- resolve_scenario(scenarios[[sc_name]])
  classes <- split_size_classes(filter_features(features, s), s)
  for (cls in names(classes)) {
    tag <- scenario_tag(s, cls)
    jobs[[length(jobs) + 1]] <- list(
      scenario = sc_name, settings = s, size_class = cls, pts = classes[[cls]], tag = tag,
      nc_file  = file.path(nc_dir, sprintf("Distance_%s.nc", tag)),
      png_file = file.path(output_dir, sprintf("Map_%s_Log10.png", tag))
    )
  }
}

manifest <- bind_rows(lapply(jobs, function(j) data.frame(
  scenario           = j$scenario,
  size_class         = j$size_class,
  description        = scenario_description(j$settings, j$size_class),
  n_features         = nrow(j$pts),
  n_seamounts        = sum(j$pts$feature_type == "Seamount"),
  n_knolls           = sum(j$pts$feature_type == "Knoll"),
  max_summit_depth_m = if (is.null(j$settings$max_summit_depth_m)) NA else j$settings$max_summit_depth_m,
  min_height_m       = if (is.null(j$settings$min_height_m)) NA else j$settings$min_height_m,
  res_km             = res_km,
  nc_file            = j$nc_file
)))
print(manifest[, c("scenario", "size_class", "n_features", "nc_file")])
write.csv(manifest, file.path(output_dir, "distance_raster_manifest.csv"), row.names = FALSE)

to_compute <- Filter(function(j) recalculate == "all" || (recalculate == "missing" && !file.exists(j$nc_file)), jobs)

if (length(to_compute) > 0) {
  cat(sprintf("Computing %d distance raster(s) on %d cores...\n", length(to_compute), n_cores))
  cl <- makeCluster(n_cores)
  invisible(clusterEvalQ(cl, library(terra)))

  tryCatch({
    for (j in to_compute) {
      cat(sprintf("\n=== %s [%s] : %d features -> %s\n",
                  j$scenario, j$size_class, nrow(j$pts), j$nc_file))
      s <- j$settings
      atts <- c(
        sprintf("scenario=%s", j$scenario),
        sprintf("size_class=%s", j$size_class),
        # writeCDF splits atts on "=", so spell out comparison operators.
        sprintf("filter_description=%s",
                gsub(">=", "GE", gsub("<=", "LE", scenario_description(s, j$size_class)))),
        sprintf("include_seamounts=%s", s$include_seamounts),
        sprintf("include_knolls=%s", s$include_knolls),
        sprintf("max_summit_depth_m=%s", if (is.null(s$max_summit_depth_m)) "none" else s$max_summit_depth_m),
        sprintf("min_height_m=%s", if (is.null(s$min_height_m)) "none" else s$min_height_m),
        sprintf("n_features=%d", nrow(j$pts)),
        "source=Yesson et al. 2011, doi:10.1016/j.dsr.2011.02.004",
        "longitude_convention=0-360"
      )
      compute_distance_nc(
        j$pts, j$nc_file,
        longname = sprintf("Geodesic distance to nearest feature (%s)", scenario_description(s, j$size_class)),
        atts = atts, cl = cl
      )
    }
  }, finally = stopCluster(cl))
} else {
  cat("All requested distance rasters already exist; nothing to compute.\n")
}

# ==============================================================================
# 5. DISTANCE MAPS (one PNG per NetCDF)
# ==============================================================================

if (make_distance_maps) {
  robinson_crs <- "+proj=robin +lon_0=180 +datum=WGS84 +units=m +no_defs"

  world_map <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf") %>%
    st_geometry() %>%
    # A Pacific-centred projection has its seam at Greenwich, not at 180 degrees.
    st_break_antimeridian(lon_0 = 180) %>%
    st_collection_extract("POLYGON") %>%
    st_transform(crs = robinson_crs)

  for (j in jobs) {
    if (!file.exists(j$nc_file)) {
      message("Skipping map, NetCDF not found: ", j$nc_file)
      next
    }
    cat(sprintf("Mapping %s...\n", j$nc_file))

    r <- aggregate(rast(j$nc_file), fact = 2, fun = mean)
    # GDAL expects geographic longitudes in [-180, 180] when reprojecting.
    r_rob <- project(crop(rotate(r), ext(-180, 180, -75, 75)), robinson_crs, method = "bilinear")

    r_df <- as.data.frame(r_rob, xy = TRUE, na.rm = TRUE)
    names(r_df)[3] <- "Distance"
    # Cap at 1 km so log10(0) does not produce -Inf / blank pixels.
    r_df$Distance <- pmax(r_df$Distance, 1)

    pts_xy <- as.data.frame(st_coordinates(st_transform(j$pts, robinson_crs)))

    p <- ggplot() +
      geom_tile(data = r_df, aes(x = x, y = y, fill = Distance),
                width = res(r_rob)[1], height = res(r_rob)[2]) +
      geom_sf(data = world_map, fill = "grey15", color = "grey35", linewidth = 0.1) +
      geom_point(data = pts_xy, aes(x = X, y = Y), color = "#ffcc00", size = 0.3, alpha = 0.7) +
      scale_fill_viridis_c(
        option = "mako", direction = -1, trans = "log10",
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)),
        name = "Distance (km) [Log10 Scale]", na.value = "#0b141f"
      ) +
      coord_sf(crs = robinson_crs, default_crs = robinson_crs, datum = 4326,
               xlim = c(xmin(r_rob), xmax(r_rob)), ylim = c(ymin(r_rob), ymax(r_rob)),
               expand = FALSE) +
      labs(
        title = sprintf("Distance to Nearest Feature: %s", j$scenario),
        subtitle = sprintf("Filters: %s (n = %s)",
                           scenario_description(j$settings, j$size_class),
                           format(nrow(j$pts), big.mark = ",")),
        caption = sprintf("File: %s\nSource: Yesson et al. (2011) doi:10.1016/j.dsr.2011.02.004",
                          basename(j$nc_file))
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = dark_background, color = NA),
        plot.background = element_rect(fill = dark_background, color = NA),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 20, color = dark_text, margin = margin(b = 8)),
        plot.subtitle = element_text(size = 13, color = muted_text, margin = margin(b = 10)),
        plot.caption = element_text(size = 11, face = "italic", color = muted_text, hjust = 1, margin = margin(t = 15)),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(4, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title.position = "top",
        legend.background = element_rect(fill = dark_panel, color = NA),
        legend.title = element_text(hjust = 0.5, face = "bold", size = 13, color = dark_text),
        legend.text = element_text(size = 11, color = dark_text)
      )

    ggsave(j$png_file, plot = p, width = slide_width, height = slide_height, dpi = 300, bg = dark_background)
  }
  cat("Distance maps saved to", output_dir, "\n")
}
