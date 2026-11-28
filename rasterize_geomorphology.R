# Harris et al. (2014): one fractional-coverage raster per polygon layer.
library(terra)

data_root <- Sys.getenv("PIFSC_GIS_DATA_DIR", unset = path.expand("~/data/deepsea/data_download"))
input_dir <- file.path(data_root, "Geomorphology", "Global_Seafloor_Geomorphic_Features_Map")
output_dir <- "output/rasters/geomorphology"
overwrite <- FALSE

template <- rast(ext(-180, 180, -60, 60), resolution = 0.05, crs = "EPSG:4326")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(input_dir, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("No shapefiles found in: ", input_dir)

for (file in files) {
  layer <- tools::file_path_sans_ext(basename(file))
  output_file <- file.path(output_dir, paste0("geomorph_", layer, "_fraction_0p05deg.tif"))
  if (file.exists(output_file) && !overwrite) {
    message("Using existing raster: ", output_file)
    next
  }
  message("Rasterizing: ", layer)
  polygons <- vect(file)
  if (crs(polygons) == "") stop("Missing CRS: ", file)
  if (geomtype(polygons) != "polygons") stop("Expected polygons: ", file)
  polygons <- project(polygons, crs(template))
  polygons <- crop(polygons, ext(template))
  if (nrow(polygons) == 0) {
    message("No features inside the grid for ", layer, "; writing zero coverage.")
    coverage <- init(template, 0)
  } else {
    # cover = TRUE estimates the fraction of each cell occupied by polygons.
    coverage <- rasterize(polygons, template, field = 1, background = 0, cover = TRUE)
  }
  writeRaster(coverage, output_file, overwrite = overwrite)
}
message("Outputs: ", output_dir)
