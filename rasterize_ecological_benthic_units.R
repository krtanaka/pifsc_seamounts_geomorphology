# EBU already contains rasters: align them, rather than duplicate vector products.
library(terra)

data_root <- Sys.getenv("PIFSC_GIS_DATA_DIR", unset = path.expand("~/data/deepsea/data_download"))
input_dir <- file.path(data_root, "Ecological_Benthic_Units", "EBU_GIS")
output_dir <- "output/rasters/ecological_benthic_units"
overwrite <- FALSE

template <- rast(ext(-180, 180, -60, 60), resolution = 0.05, crs = "EPSG:4326")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
folders <- c("01_classification", "02_environmental_predictors", "04_supporting")
files <- unlist(lapply(folders, function(folder) {
  list.files(file.path(input_dir, folder), pattern = "\\.tif$", full.names = TRUE,
             ignore.case = TRUE)
}), use.names = FALSE)
files <- files[!grepl("\\.(part|partial|tmp)(\\.|$)", basename(files), ignore.case = TRUE)]
if (length(files) == 0) stop("No EBU GeoTIFFs found in: ", input_dir)
if (anyDuplicated(basename(files))) stop("Duplicate EBU filenames; output names would collide.")

for (file in files) {
  layer <- tools::file_path_sans_ext(basename(file))
  output_file <- file.path(output_dir, paste0(layer, "_0p05deg.tif"))
  if (file.exists(output_file) && !overwrite) {
    message("Using existing raster: ", output_file)
    next
  }
  source <- rast(file)
  if (crs(source) == "") stop("Missing CRS: ", file)
  categorical <- basename(dirname(file)) == "01_classification" ||
    grepl("lithology|slope_hemisphere_class", layer, ignore.case = TRUE)
  method <- if (categorical) "near" else "bilinear"
  if (grepl("bathymetry_0bd2bd53$", layer)) {
    warning("Source metadata reports a 0-360 range for ", layer,
            "; verify its meaning before using it as bathymetry.")
  }
  message("Aligning ", layer, " (", method, ")")
  if (is.lonlat(source) && xmax(source) > 180) source <- rotate(source)
  project(source, template, method = method, filename = output_file, overwrite = overwrite)
}
message("Outputs: ", output_dir)
