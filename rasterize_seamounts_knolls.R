# Yesson et al. (2011): distance to summit points, not polygon footprints.
library(terra)

data_root <- Sys.getenv("PIFSC_GIS_DATA_DIR", unset = path.expand("~/data/deepsea/data_download"))
input_dir <- file.path(data_root, "Seamounts_Knolls", "Modelled_Seamounts_2011_v1", "01_Data")
output_dir <- "output/rasters/seamounts_knolls"
max_depth_m <- 400 # Use 1000 for sensitivity, or NULL for all depths.
depth_sign <- "auto" # Or "positive_down" / "negative_down"; DEPTH must be in metres.
overwrite <- FALSE

template <- rast(ext(-180, 180, -60, 60), resolution = 0.05, crs = "EPSG:4326")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

layers <- lapply(c("Seamounts", "Knolls"), function(layer) {
  points <- vect(file.path(input_dir, layer, paste0(layer, ".shp")))
  if (crs(points) == "") stop("Missing CRS: ", layer)
  if (geomtype(points) != "points") stop("Expected summit points: ", layer)
  if (!"DEPTH" %in% names(points)) stop("Missing DEPTH field: ", layer)
  project(points[, "DEPTH"], "EPSG:4326")
})
points <- do.call(rbind, layers)
depth <- points$DEPTH
if (!is.numeric(depth)) stop("DEPTH must be numeric, in metres.")
if (depth_sign == "auto") {
  if (!any(is.finite(depth))) stop("Cannot infer DEPTH sign: no finite depths.")
  depth_sign <- if (median(depth[is.finite(depth)]) < 0) "negative_down" else "positive_down"
}
if (!depth_sign %in% c("positive_down", "negative_down")) stop("Invalid depth_sign.")
message("DEPTH convention: ", depth_sign, ". Check against source metadata.")
if (depth_sign == "negative_down") depth <- -depth
keep <- !is.finite(depth) | depth > 0
if (!is.null(max_depth_m)) {
  if (length(max_depth_m) != 1 || !is.finite(max_depth_m) || max_depth_m <= 0) {
    stop("max_depth_m must be NULL or a positive number.")
  }
  keep <- is.finite(depth) & depth > 0 & depth <= max_depth_m
}
message("Depth screening removed ", sum(!keep), " features.")
points <- points[keep, ]
if (nrow(points) == 0) stop("No features remain after depth screening.")
xy <- crds(points)
xy[, 1] <- (xy[, 1] + 180) %% 360 - 180
if (any(!is.finite(xy))) stop("Non-finite summit coordinates.")
duplicates <- duplicated(data.frame(round(xy, 5)))
message("Removing ", sum(duplicates), " duplicate summit locations.")
points <- vect(xy[!duplicates, , drop = FALSE], crs = "EPSG:4326")

# Keep all global points, including those outside the output grid, to avoid edge effects.
tag <- if (is.null(max_depth_m)) "DepthAll" else paste0("DepthLE", max_depth_m, "m")
output_file <- file.path(output_dir, paste0("Distance_SeamountsKnolls_", tag, "_0p05deg_km.tif"))
if (!file.exists(output_file) || overwrite) {
  message("Calculating distance to ", nrow(points), " summit points...")
  distance(template, points, unit = "km", method = "geo", rasterize = FALSE,
           filename = output_file, overwrite = overwrite)
} else {
  message("Using existing raster; set overwrite = TRUE to recalculate.")
}
message("Output: ", output_file)
