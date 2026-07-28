# ---------------------------------------------------------
# The Bulletproof Raster Geomorphology Workflow
# ---------------------------------------------------------
library(terra)

# 1. Define Directories
data_dir <- "C://Users/Kisei.Tanaka/Desktop/global seafloor geomorphic features map/"
out_dir <- "C://Users/Kisei.Tanaka/Desktop/rasterized_geomorph_10km/"

# Create the output folder if it doesn't exist
dir.create(out_dir, showWarnings = FALSE)

# Find all shapefiles
shp_paths <- list.files(data_dir, pattern = "\\.shp$", full.names = TRUE)
layer_names <- tools::file_path_sans_ext(basename(shp_paths))

# 2. Create the 10km Global Raster Template
# Extent is global (-180 to 180). Resolution is 0.09 degrees (~10km)
template_rast <- rast(ext(-180, 180, -90, 90), resolution = 0.05, crs = "EPSG:4326")

# 3. Batch Rasterize and Export
message("Rasterizing shapefiles. This completely bypasses vector topology bugs...")

for (i in seq_along(shp_paths)) {
  name <- layer_names[i]
  cat("Converting & Saving:", name, "\n")
  
  # Load vector with terra
  v <- vect(shp_paths[i])
  
  # Create a dummy column so every pixel gets a value of 1
  v$present <- 1
  
  # Rasterize: If the feature exists in the pixel, assign 1. Otherwise, NA.
  r <- rasterize(v, template_rast, field = "present", background = NA)
  
  # Export as GeoTIFF
  out_file <- paste0(out_dir, name, "_10km.tif")
  writeRaster(r, out_file, overwrite = TRUE)
}

library(terra)
library(ggplot2)

# 1. Define Directories and Projections
out_dir <- "C://Users/Kisei.Tanaka/Desktop/rasterized_geomorph_10km/"
pacific_crs <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

message("Loading rasters and projecting to the Pacific Robinson...")

# 2. Load and Project Global Features
# (We use method = "near" for categorical/binary data like geomorphology)
abyss_pac     <- project(rast(paste0(out_dir, "Abyss_10km.tif")), pacific_crs, method = "near")
ridges_pac    <- project(rast(paste0(out_dir, "Ridges_10km.tif")), pacific_crs, method = "near")
trenches_pac  <- project(rast(paste0(out_dir, "Trenches_10km.tif")), pacific_crs, method = "near")
seamounts_pac <- project(rast(paste0(out_dir, "Seamounts_10km.tif")), pacific_crs, method = "near")

message("Converting rasters to data frames for ggplot2...")

# 3. Convert Global Rasters to Data Frames
# na.rm = TRUE is critical so we don't plot empty pixels
df_abyss_pac     <- as.data.frame(abyss_pac, xy = TRUE, na.rm = TRUE)
df_ridges_pac    <- as.data.frame(ridges_pac, xy = TRUE, na.rm = TRUE)
df_trenches_pac  <- as.data.frame(trenches_pac, xy = TRUE, na.rm = TRUE)
df_seamounts_pac <- as.data.frame(seamounts_pac, xy = TRUE, na.rm = TRUE)

# 4. Build the Global Pacific ggplot
message("Building global Pacific map...")
p_global <- ggplot() +
  # Note: The fill values here map to the 'present' column (1) we created during rasterization.
  # We force the fill color directly since these are effectively binary presence layers.
  geom_tile(data = df_abyss_pac, aes(x = x, y = y), fill = "#1a2a4a", alpha = 0.8) +
  geom_tile(data = df_ridges_pac, aes(x = x, y = y), fill = "#d68c45") +
  geom_tile(data = df_trenches_pac, aes(x = x, y = y), fill = "#00ffff") +
  geom_tile(data = df_seamounts_pac, aes(x = x, y = y), fill = "#e63946") +
  
  # Ensure the Robinson projection scales correctly
  coord_fixed() + 
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#0b1325", color = NA),
    plot.background = element_rect(fill = "#ffffff", color = NA),
    panel.grid.major = element_line(color = "#ffffff1a", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.text = element_blank(), # Hide arbitrary meter coordinates
    axis.title = element_blank()
  ) +
  labs(
    title = "Pacific Ocean Seafloor Geomorphology",
    subtitle = "Abyss (Dark Blue), Ridges (Orange), Trenches (Cyan), Seamounts (Red)"
  )

# 5. Save the Global Map
ggsave(filename = "C://Users/Kisei.Tanaka/Desktop/Pacific_Geomorphology_Global.png", 
       plot = p_global, width = 12, height = 7, dpi = 300)
message("Global Pacific map saved.")

# ---------------------------------------------------------
# REGIONAL: Main Hawaiian Islands 
# ---------------------------------------------------------

http://googleusercontent.com/map_location_reference/1
message("Processing regional [Island of Hawai'i](http://googleusercontent.com/map_location_reference/0) map...")

# 6. Define MHI Extent (Unprojected WGS84)
mhi_ext <- ext(-161.5, -154.0, 18.5, 22.5)

# 7. Load and Crop Unprojected Rasters
mhi_abyss     <- crop(rast(paste0(out_dir, "Abyss_10km.tif")), mhi_ext)
mhi_slope     <- crop(rast(paste0(out_dir, "Slope_10km.tif")), mhi_ext)
mhi_shelf     <- crop(rast(paste0(out_dir, "Shelf_10km.tif")), mhi_ext)
mhi_seamounts <- crop(rast(paste0(out_dir, "Seamounts_10km.tif")), mhi_ext)

# 8. Convert MHI Rasters to Data Frames
df_mhi_abyss     <- as.data.frame(mhi_abyss, xy = TRUE, na.rm = TRUE)
df_mhi_slope     <- as.data.frame(mhi_slope, xy = TRUE, na.rm = TRUE)
df_mhi_shelf     <- as.data.frame(mhi_shelf, xy = TRUE, na.rm = TRUE)
df_mhi_seamounts <- as.data.frame(mhi_seamounts, xy = TRUE, na.rm = TRUE)

# 9. Build the MHI ggplot
p_mhi <- ggplot() +
  geom_tile(data = df_mhi_abyss, aes(x = x, y = y), fill = "#1a2a4a") +
  geom_tile(data = df_mhi_slope, aes(x = x, y = y), fill = "#2c4c7c") +
  geom_tile(data = df_mhi_shelf, aes(x = x, y = y), fill = "#8b9bb4") +
  geom_tile(data = df_mhi_seamounts, aes(x = x, y = y), fill = "#e63946") +
  
  # Keeps lat/long aspect ratio accurate
  coord_quickmap() +
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#0b1325", color = NA), 
    plot.background = element_rect(fill = "#ffffff", color = NA),
    panel.grid.major = element_line(color = "#ffffff1a", linewidth = 0.2), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Seafloor Geomorphology: Main Hawaiian Islands",
    subtitle = "Abyss (Dark Blue), Slope (Mid Blue), Shelf (Light Blue), Seamounts (Red)",
    x = "Longitude",
    y = "Latitude"
  )

# 10. Save the MHI Map
ggsave(filename = "C://Users/Kisei.Tanaka/Desktop/MHI_Geomorphology_Regional.png", 
       plot = p_mhi, width = 10, height = 6, dpi = 300)
message("Regional MHI map saved.")