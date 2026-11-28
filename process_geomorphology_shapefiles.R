# ==============================================================================
# MASTER SEAFLOOR GEOMORPHOLOGY WORKFLOW
# Model-ready GIS exports + optional Pacific visualization
#
# Source:
# Harris, P.T., Macmillan-Lawler, M., Rupp, J., & Baker, E.K. (2014).
# Geomorphology of the oceans. Marine Geology, 352, 4-24.
# https://doi.org/10.1016/j.margeo.2014.01.011
# ==============================================================================

library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(scales)
library(foreach)
library(doParallel)
library(dplyr)
library(readr)

rm(list = ls())

# ==============================================================================
# 1. PATHS & PARAMETERS
# ==============================================================================

data_dir <- "A:/My Drive/Data/Seamounts & Knolls/global seafloor geomorphic features map/"
out_dir  <- "A:/My Drive/Data/Seamounts & Knolls/rasterized_geomorph/"

# REVISED: keep model-ready rasters separate from projected map products.
model_out_dir <- file.path(out_dir, "model_ready_0p05deg")
map_out_dir   <- file.path(out_dir, "maps")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(map_out_dir, showWarnings = FALSE, recursive = TRUE)

pacific_crs <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# REVISED:
# Model-ready grid matches the current DSC environmental domain/resolution.
# These are the files to import into 3. prep_prediction_layers_baseline.R.
model_template <- rast(
  ext(-180, 180, -60, 60),
  resolution = 0.05,
  crs = "EPSG:4326"
)

# Fractional coverage is recommended for modeling.
# Binary presence rasters are also exported for QA/QC and sensitivity tests.
export_binary <- TRUE
make_maps <- TRUE

ordered_layers <- c(
  "Abyss", "Abyssal_Classification", "Basins", "Shelf", "Shelf_Classification",
  "Slope", "Plateaus", "Terraces", "Rises", "Fans", "Bridges", "Sills",
  "Troughs", "Glacial_troughs", "Ridges", "Spreading_ridges", "Escarpments",
  "Canyons", "Shelf_valleys", "Rift_valleys", "Guyots", "Seamounts",
  "Trenches", "Hadal"
)

# ==============================================================================
# 2. FILE DISCOVERY
# ==============================================================================

shp_paths <- list.files(
  data_dir,
  pattern = "\\.shp$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(shp_paths) == 0) {
  stop("No shapefiles found in: ", data_dir)
}

source_layer_names <- tools::file_path_sans_ext(basename(shp_paths))

source_lookup <- setNames(
  shp_paths,
  tolower(source_layer_names)
)

missing_layers <- ordered_layers[
  !tolower(ordered_layers) %in% names(source_lookup)
]

if (length(missing_layers) > 0) {
  message(
    "Expected layers not found and skipped: ",
    paste(missing_layers, collapse = ", ")
  )
}

layers_to_process <- ordered_layers[
  tolower(ordered_layers) %in% names(source_lookup)
]

if (length(layers_to_process) == 0) {
  stop("None of the expected geomorphology layers were found.")
}

# ==============================================================================
# PHASE A: MODEL-READY RASTERIZATION
# ==============================================================================

message(
  "Rasterizing ",
  length(layers_to_process),
  " geomorphology layers to 0.05-degree WGS84..."
)

# Respect cloud/container CPU quotas; reserve two cores, with at least one worker.
num_cores <- min(length(layers_to_process), parallelly::availableCores(omit = 2))

cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

rasterization_results <- foreach(
  layer = layers_to_process,
  .packages = "terra",
  .export = c(
    "source_lookup",
    "model_out_dir",
    "export_binary"
  )
) %dopar% {
  
  # REVISED:
  # Create the terra raster independently inside each worker.
  # SpatRaster objects should not be exported across PSOCK workers.
  model_template <- terra::rast(
    terra::ext(-180, 180, -60, 60),
    resolution = 0.05,
    crs = "EPSG:4326"
  )

  source_file <- source_lookup[[tolower(layer)]]

  fraction_file <- file.path(
    model_out_dir,
    paste0("geomorph_", layer, "_fraction.tif")
  )

  binary_file <- file.path(
    model_out_dir,
    paste0("geomorph_", layer, "_presence.tif")
  )

  fraction_ok <- file.exists(fraction_file) &&
    file.size(fraction_file) > 0

  binary_ok <- !export_binary ||
    (file.exists(binary_file) && file.size(binary_file) > 0)

  if (fraction_ok && binary_ok) {
    return(
      data.frame(
        layer = layer,
        source_file = source_file,
        fraction_file = fraction_file,
        binary_file = if (export_binary) binary_file else NA_character_,
        status = "existing"
      )
    )
  }

  v <- terra::vect(source_file)

  # REVISED: never silently assign CRS to a source vector.
  if (is.na(terra::crs(v)) || terra::crs(v) == "") {
    stop("Missing CRS in source shapefile: ", source_file)
  }

  if (!terra::same.crs(v, model_template)) {
    v <- terra::project(v, terra::crs(model_template))
  }

  v <- terra::crop(v, terra::ext(model_template))

  if (nrow(v) == 0) {

    fraction_r <- model_template
    terra::values(fraction_r) <- 0

  } else {

    geometry_types <- unique(as.character(terra::geomtype(v)))

    if (!all(grepl("polygon", geometry_types, ignore.case = TRUE))) {
      stop(
        "Expected polygon geometry for ",
        layer,
        "; found: ",
        paste(geometry_types, collapse = ", ")
      )
    }

    # REVISED:
    # cover=TRUE returns the estimated fraction of each 0.05-degree cell
    # covered by the geomorphic polygon. Outside-feature cells are explicitly 0.
    fraction_r <- terra::rasterize(
      v,
      model_template,
      cover = TRUE,
      background = 0
    )

    fraction_r <- terra::clamp(
      fraction_r,
      lower = 0,
      upper = 1,
      values = TRUE
    )
  }

  names(fraction_r) <- paste0(
    "geomorph_",
    tolower(layer),
    "_fraction"
  )

  terra::writeRaster(
    fraction_r,
    fraction_file,
    overwrite = TRUE,
    datatype = "FLT4S",
    gdal = c("COMPRESS=LZW", "PREDICTOR=2")
  )

  if (export_binary) {

    binary_r <- terra::ifel(fraction_r > 0, 1, 0)

    names(binary_r) <- paste0(
      "geomorph_",
      tolower(layer),
      "_presence"
    )

    terra::writeRaster(
      binary_r,
      binary_file,
      overwrite = TRUE,
      datatype = "INT1U",
      gdal = "COMPRESS=LZW"
    )
  }

  data.frame(
    layer = layer,
    source_file = source_file,
    fraction_file = fraction_file,
    binary_file = if (export_binary) binary_file else NA_character_,
    status = "created"
  )
}

parallel::stopCluster(cl)

rasterization_results <- bind_rows(rasterization_results)

# ==============================================================================
# 3. VALIDATE AND EXPORT STACK/MANIFEST
# ==============================================================================

fraction_stack <- rast(rasterization_results$fraction_file)

if (!compareGeom(
  fraction_stack,
  model_template,
  stopOnError = FALSE
)) {
  stop("Exported geomorphology rasters do not match the model template.")
}

fraction_ranges <- global(
  fraction_stack,
  range,
  na.rm = TRUE
)

fraction_ranges$layer <- rasterization_results$layer

if (
  any(fraction_ranges$min < 0, na.rm = TRUE) ||
  any(fraction_ranges$max > 1, na.rm = TRUE)
) {
  stop("One or more fractional rasters have values outside 0-1.")
}

# REVISED:
# Export all source categories, but flag a smaller set for initial model testing.
# We will still evaluate these against the existing bathymetry/slope/TRI/TPI
# predictors before final inclusion.
initial_candidate_layers <- c(
  "Plateaus",
  "Ridges",
  "Spreading_ridges",
  "Escarpments",
  "Canyons",
  "Guyots",
  "Seamounts"
)

manifest <- rasterization_results %>%
  mutate(
    model_variable = paste0(
      "geomorph_",
      tolower(layer),
      "_fraction"
    ),
    raster_type = "fractional polygon coverage",
    units = "proportion of 0.05-degree cell",
    minimum = fraction_ranges$min,
    maximum = fraction_ranges$max,
    candidate_for_initial_test =
      layer %in% initial_candidate_layers,
    source_citation = paste0(
      "Harris PT, Macmillan-Lawler M, Rupp J, Baker EK (2014) ",
      "Geomorphology of the oceans. Marine Geology 352:4-24. ",
      "doi:10.1016/j.margeo.2014.01.011"
    )
  )

write_csv(
  manifest,
  file.path(
    model_out_dir,
    "geomorphology_model_layer_manifest.csv"
  )
)

writeRaster(
  fraction_stack,
  file.path(
    model_out_dir,
    "geomorphology_fraction_stack.tif"
  ),
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW", "PREDICTOR=2")
)

# ==============================================================================
# PHASE B: OPTIONAL MAPPING
# ==============================================================================

if (make_maps) {

  message("Preparing optional Pacific geomorphology maps...")

  # Categorical display raster only. NOT for MaxEnt.
  map_class <- model_template
  values(map_class) <- NA

  for (i in seq_along(layers_to_process)) {

    layer <- layers_to_process[i]

    fraction_file <- rasterization_results$fraction_file[
      rasterization_results$layer == layer
    ][1]

    feature_r <- rast(fraction_file)

    # Preserve original stacking order: later layers overwrite earlier layers.
    map_class <- ifel(
      feature_r > 0,
      i,
      map_class
    )
  }

  names(map_class) <- "geomorphology_display_class"

  writeRaster(
    map_class,
    file.path(
      map_out_dir,
      "geomorphology_display_class_WGS84.tif"
    ),
    overwrite = TRUE,
    datatype = "INT2U",
    gdal = "COMPRESS=LZW"
  )

  map_class_pac <- project(
    map_class,
    pacific_crs,
    method = "near"
  )

  df_all_pac <- as.data.frame(
    map_class_pac,
    xy = TRUE,
    na.rm = TRUE
  )

  df_all_pac$Feature <- factor(
    df_all_pac$geomorphology_display_class,
    levels = seq_along(layers_to_process),
    labels = gsub("_", " ", layers_to_process)
  )

  world_land <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf"
  )

  world_land_cropped <- suppressWarnings(
    sf::st_crop(
      world_land,
      xmin = -180,
      xmax = 180,
      ymin = -85,
      ymax = 85
    )
  )

  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)

  land_pac <- world_land_cropped %>%
    sf::st_make_valid() %>%
    sf::st_transform(crs = pacific_crs)

  sf::sf_use_s2(old_s2)

  vibrant_colors <- c(
    "Abyss" = "#122643",
    "Abyssal Classification" = "#1a365d",
    "Basins" = "#22467b",
    "Rises" = "#2e5b9a",
    "Shelf" = "#2d6a4f",
    "Shelf Classification" = "#40916c",
    "Shelf valleys" = "#52b788",
    "Slope" = "#1b4332",
    "Plateaus" = "#74c69d",
    "Terraces" = "#95d5b2",
    "Fans" = "#b7e4c7",
    "Ridges" = "#d68c45",
    "Spreading ridges" = "#ffb703",
    "Rift valleys" = "#fb8500",
    "Escarpments" = "#ee9b00",
    "Guyots" = "#ca6702",
    "Seamounts" = "#ae2012",
    "Bridges" = "#9b5de5",
    "Sills" = "#f15bb5",
    "Troughs" = "#00bbf9",
    "Glacial troughs" = "#00f5d4",
    "Canyons" = "#00f5d4",
    "Trenches" = "#00ffff",
    "Hadal" = "#ffffff"
  )

  source_text <- paste0(
    "Source: Harris, P.T., Macmillan-Lawler, M., Rupp, J. & Baker, E.K. ",
    "(2014). Geomorphology of the oceans. Marine Geology 352:4-24. ",
    "doi:10.1016/j.margeo.2014.01.011."
  )

  message("Rendering Global Pacific map...")

  p_global <- ggplot() +
    geom_raster(
      data = df_all_pac,
      aes(x = x, y = y, fill = Feature)
    ) +
    geom_sf(
      data = land_pac,
      fill = "#0a0d14",
      color = "#1c2331",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    coord_sf(
      crs = pacific_crs,
      datum = NA
    ) +
    scale_fill_manual(
      values = vibrant_colors,
      guide = guide_legend(nrow = 2),
      drop = FALSE
    ) +
    labs(
      title = "Pacific Ocean Seafloor Geomorphology",
      caption = source_text
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "#020408", color = NA),
      plot.background = element_rect(fill = "#020408", color = NA),
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0.5,
        color = "#ffffff"
      ),
      plot.caption = element_text(
        size = 10,
        color = "#9aabc2",
        hjust = 0
      ),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.background = element_rect(
        fill = "#0a0d14",
        color = "#1c2331",
        linewidth = 0.5
      ),
      legend.key = element_blank(),
      legend.text = element_text(
        color = "#ffffff",
        size = 11,
        face = "bold"
      ),
      legend.title = element_blank(),
      legend.key.size = grid::unit(0.5, "cm")
    )

  ggsave(
    filename = file.path(
      map_out_dir,
      "Pacific_Geomorphology.png"
    ),
    plot = p_global,
    width = 16,
    height = 9,
    dpi = 300,
    bg = "#020408"
  )

  message("Rendering Focused Basin map...")

  p_focused <- ggplot() +
    geom_raster(
      data = df_all_pac,
      aes(x = x, y = y, fill = Feature)
    ) +
    geom_sf(
      data = land_pac,
      fill = "#0a0d14",
      color = "#1c2331",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    coord_sf(
      crs = pacific_crs,
      xlim = c(-6500000, 11500000),
      ylim = c(-6000000, 6500000),
      datum = NA,
      expand = FALSE
    ) +
    scale_fill_manual(
      values = vibrant_colors,
      guide = guide_legend(ncol = 1),
      drop = FALSE
    ) +
    labs(
      title = "Pacific Ocean Basin Seafloor Geomorphology",
      subtitle = "Global Seafloor Geomorphic Features Map",
      caption = source_text
    ) +
    theme_minimal(base_size = 30) +
    theme(
      panel.background = element_rect(fill = "#020408", color = NA),
      plot.background = element_rect(fill = "#020408", color = NA),
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0.5,
        color = "#ffffff"
      ),
      plot.subtitle = element_text(
        size = 16,
        hjust = 0.5,
        color = "#8b9bb4"
      ),
      plot.caption = element_text(
        size = 10,
        color = "#9aabc2",
        hjust = 0
      ),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.background = element_rect(
        fill = "#0a0d14",
        color = "#1c2331",
        linewidth = 0.5
      ),
      legend.key = element_blank(),
      legend.text = element_text(
        color = "#ffffff",
        size = 11,
        face = "bold"
      ),
      legend.title = element_blank(),
      legend.key.size = grid::unit(0.5, "cm")
    )

  ggsave(
    filename = file.path(
      map_out_dir,
      "Pacific_Basin_Focused_Map.png"
    ),
    plot = p_focused,
    width = 16,
    height = 9,
    dpi = 300,
    bg = "#020408"
  )
}

message(
  "Processing complete. Model-ready layers are in: ",
  model_out_dir
)
