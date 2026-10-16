library(sf)
library(maps)
library(dplyr)
library(ggplot2)

# Ensure output directory exists
if(!dir.exists("output")) dir.create("output")

# ==============================================================================
# 1. SETUP GLOBAL BASE MAP
# ==============================================================================
sf_use_s2(FALSE)
world_base <- st_as_sf(map("world2", plot = FALSE, fill = TRUE)) %>% st_make_valid()
st_crs(world_base) <- 4326
world_cropped <- st_crop(world_base, c(xmin = 110, ymin = -50, xmax = 300, ymax = 70))
sf_use_s2(TRUE)

# ==============================================================================
# 2. LOAD AND PREPARE DATASET LAYERS
# ==============================================================================
base_path <- "A:/My Drive/Data/Seamounts & Knolls/ModelledSeamounts2011_v1/01_Data/"

seamounts_sf  <- st_read(file.path(base_path, "Seamounts/Seamounts.shp"))
seamounts_sf = seamounts_sf %>% filter(HEIGHT >= 2500)
seamounts_360 <- st_shift_longitude(seamounts_sf)

knolls_sf  <- st_read(file.path(base_path, "Knolls/Knolls.shp"))
knolls_360 <- st_shift_longitude(knolls_sf)

# Common view constraints for spatial maps
pacific_view <- coord_sf(crs = "+proj=robin +lon_0=180", 
                         xlim = c(120, 290), ylim = c(-45, 60), 
                         default_crs = 4326, expand = FALSE)

# Text caption block for all files
metadata_caption <- paste0(
  "Data Provider: Zoological Society of London (2011) | Version: 1.0 | License: CC BY 3.0\n",
  "Citation: Yesson C, Clark MR, Taylor M, Rogers AD (2011). The global distribution of seamounts based on 30-second bathymetry data.\n",
  "Deep Sea Research Part I 58: 442-453. doi: 10.1016/j.dsr.2011.02.004 | Data URL: http://data.unep-wcmc.org/datasets/41"
)

dark_background <- "black"
dark_panel <- "black"
dark_text <- "#e6edf3"
muted_text <- "#a9b8c9"
slide_width <- 16
slide_height <- 9
recalculate_distance_rasters <- FALSE # Set TRUE to recompute the 5 km distance rasters.

# Shared theme adjustments for spatial maps (Legend Top-Right / Over USA)
map_theme_adjust <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = dark_background, color = NA),
    panel.background = element_rect(fill = dark_background, color = NA),
    plot.title = element_text(face = "bold", size = 18, color = dark_text),
    plot.caption = element_text(size = 9, color = muted_text, hjust = 0, lineheight = 1.1),
    axis.title = element_text(size = 12, color = dark_text),
    axis.text = element_text(size = 10, color = muted_text),
    panel.grid = element_blank(),
    legend.position = c(0.86, 0.85),             
    legend.background = element_rect(fill = dark_panel, color = "#415166", linewidth = 0.3),
    legend.key = element_blank(),
    legend.title = element_text(size = 12, face = "bold", color = dark_text),
    legend.text = element_text(size = 10, color = dark_text),
    legend.key.size = unit(0.55, "cm")
  )

# Shared theme adjustments for hex plots (Legend Top-Left)
hex_theme_adjust <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = dark_background, color = NA),
    panel.background = element_rect(fill = dark_background, color = NA),
    plot.title = element_text(face = "bold", size = 18, color = dark_text),
    plot.caption = element_text(size = 9, color = muted_text, hjust = 0, lineheight = 1.1),
    axis.title = element_text(size = 12, color = dark_text),
    axis.text = element_text(size = 10, color = muted_text),
    panel.grid = element_line(color = "#263445"),
    legend.position = c(0.12, 0.78),             
    legend.background = element_rect(fill = dark_panel, color = "#415166", linewidth = 0.3),
    legend.key = element_blank(),
    legend.title = element_text(size = 12, face = "bold", color = dark_text),
    legend.text = element_text(size = 10, color = dark_text),
    legend.key.size = unit(0.55, "cm")
  )


# ==============================================================================
# 3. GENERATE & SAVE MAPS INDIVIDUALLY
# ==============================================================================

# Map 1: Seamounts by Height
ggplot() +
  geom_sf(data = world_cropped, fill = "#202b38", color = "#526170", linewidth = 0.1) +
  geom_sf(data = seamounts_360, aes(fill = HEIGHT), color = "gray", stroke = 0.5, size = 2, alpha = 0.8, shape = 21) +
  scale_fill_viridis_c(option = "plasma", name = "Height (m)", trans = "log10") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Seamounts by Height", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/seamounts_by_height.png", width = slide_width, height = slide_height, dpi = 300, bg = dark_background)

# Map 2: Knolls by Height
ggplot() +
  geom_sf(data = world_cropped, fill = "#202b38", color = "#526170", linewidth = 0.1) +
  geom_sf(data = knolls_360, aes(color = HEIGHT), size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Height (m)", trans = "log10") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Knolls by Height", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/knolls_by_height.png", width = slide_width, height = slide_height, dpi = 300, bg = dark_background)

# Map 3: Seamounts by Summit Depth
ggplot() +
  geom_sf(data = world_cropped, fill = "#202b38", color = "#526170", linewidth = 0.1) +
  geom_sf(data = seamounts_360, aes(fill = DEPTH), color = "gray", stroke = 0.5, size = 2, alpha = 0.8, shape = 21) +
  scale_fill_viridis_c(option = "mako", name = "Depth (m)") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Seamounts by Summit Depth", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/seamounts_by_depth.png", width = slide_width, height = slide_height, dpi = 300, bg = dark_background)

# Map 4: Knolls by Summit Depth
ggplot() +
  geom_sf(data = world_cropped, fill = "#202b38", color = "#526170", linewidth = 0.1) +
  geom_sf(data = knolls_360, aes(color = DEPTH), size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "mako", name = "Depth (m)") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Knolls by Summit Depth", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/knolls_by_depth.png", width = slide_width, height = slide_height, dpi = 300, bg = dark_background)


# ==============================================================================
# 4. GENERATE & SAVE MORPHOLOGY HEXPLOTS (LEGEND TOP-LEFT)
# ==============================================================================

# Plot 5: Seamounts Morphological Density
ggplot(seamounts_360, aes(x = DEPTH, y = HEIGHT)) +
  geom_hex(bins = 50, alpha = 0.9) +
  scale_fill_viridis_c(option = "plasma", name = "Count") +
  hex_theme_adjust + 
  labs(title = "Seamounts: Morphological Density Structure", 
       x = "Depth", y = "Height", caption = metadata_caption)
ggsave("output/seamounts_morphology_hex.png", width = slide_width, height = slide_height, dpi = 300, bg = dark_background)

# Plot 6: Knolls Morphological Density
ggplot(knolls_360, aes(x = DEPTH, y = HEIGHT)) +
  geom_hex(bins = 50, alpha = 0.9) +
  scale_fill_viridis_c(option = "plasma", name = "Count") +
  hex_theme_adjust + 
  labs(title = "Knolls: Morphological Density Structure", 
       x = "Depth", y = "Height", caption = metadata_caption)
ggsave("output/knolls_morphology_hex.png", width = slide_width, height = slide_height, dpi = 300, bg = dark_background)

# 1. Load the necessary libraries
library(sf)
library(dplyr)
library(terra)
library(parallel)
library(pbapply) # For the progress bar!

# 2. Categorize the seamounts
# (Assuming your dataset is currently an sf object named 'seamounts_360')

# doi.org/10.5670/oceanog.2010.60
# doi.org/10.5194/bg-10-3039-2013
# doi.org/10.3389/fmars.2020.00238

# Create a new column dividing them into 3 equal-sized groups based on HEIGHT
# 1 = Smallest 33%, 2 = Middle 33%, 3 = Tallest 33%
seamounts_360 <- seamounts_360 %>% mutate(size_group = ntile(HEIGHT, 3))

if (recalculate_distance_rasters) {
  # 2. Categorize the seamounts based on Statistical Tertiles (33.3% splits)
  cat("Filtering seamount categories by tertiles...\n")

categories <- list(
  "All"    = seamounts_360,
  "Small"  = seamounts_360 %>% filter(size_group == 1),
  "Medium" = seamounts_360 %>% filter(size_group == 2),
  "Big"    = seamounts_360 %>% filter(size_group == 3)
)

# 3. Define Grid Parameters
res_5km <- 5 / 111.12  # ~5km resolution

# Create 60 latitudinal boundaries (chunks) from -90 to 90
n_chunks <- 60
lat_breaks <- seq(-90, 90, length.out = n_chunks + 1)

# 4. Setup the Parallel Cluster
n_cores <- 20 # Using 20 cores, leaving 2 for the OS
cat(sprintf("Firing up cluster with %d cores...\n", n_cores))
cl <- makeCluster(n_cores)

# Load required packages on all 20 worker cores
invisible(clusterEvalQ(cl, {
  library(terra)
  library(sf)
}))

# 5. Define the Worker Function
# This runs independently on every core. It calculates the distance for just 
# one latitudinal strip, saves it as a temporary file, and returns the path.
process_strip <- function(i, sf_pts, lat_breaks, res_5km) {
  
  ymin <- lat_breaks[i]
  ymax <- lat_breaks[i + 1]
  
  # Create a blank raster for just this strip
  strip_rast <- rast(ext(0, 360, ymin, ymax), resolution = res_5km, crs = "EPSG:4326")
  
  # Note: Even though we calculate distance for a small strip, we pass the 
  # ENTIRE seamount dataset (v_pts) so it accurately finds the closest point, 
  # even if that point is physically located in a different strip!
  v_pts <- vect(sf_pts)
  
  # Calculate geodesic distance (meters), then convert to km
  d_meters <- distance(strip_rast, v_pts)
  d_km <- d_meters / 1000
  names(d_km) <- "Distance_km"
  
  # Write to temporary file (safest way to pass rasters back from parallel workers)
  tmp_file <- tempfile(fileext = ".tif")
  writeRaster(d_km, tmp_file, overwrite = TRUE)
  
  return(tmp_file)
}

# 6. Execute the processing loop for each category
for (cat_name in names(categories)) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Calculating distances for %s Seamounts...\n", toupper(cat_name)))
  
  sf_data <- categories[[cat_name]]
  
  # pblapply acts like a standard loop, but runs across your 20 cores AND 
  # displays a live progress bar in the console!
  tif_files <- pblapply(1:n_chunks, process_strip, 
                        sf_pts = sf_data, 
                        lat_breaks = lat_breaks, 
                        res_5km = res_5km, 
                        cl = cl)
  
  cat("Stitching chunks together...\n")
  # Read the temporary files as a SpatRasterCollection, then mosaic them
  raster_collection <- sprc(unlist(tif_files))
  final_raster <- mosaic(raster_collection)
  
  # Export to NetCDF
  out_name <- sprintf("Distance_%s_Seamounts_5km.nc", cat_name)
  cat(sprintf("Exporting to %s...\n", out_name))
  
  writeCDF(final_raster, out_name, overwrite = TRUE,
           varname = "distance", 
           longname = sprintf("Geodesic distance to nearest %s seamount", tolower(cat_name)), 
           unit = "km")
  
  # Clean up the temporary tif files from your hard drive
  unlink(unlist(tif_files))
}

# 7. Shut down the cluster to free up memory
stopCluster(cl)
cat("\nAll 4 NetCDF files generated successfully!\n")
}

# Load the required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(maps)
library(scales)

# =============================================================================
# 1. PREPARE THE DATA (Same as before)
# =============================================================================
nc_files <- c(
  "Distance_All_Seamounts_5km.nc", 
  "Distance_Big_Seamounts_5km.nc",
  "Distance_Medium_Seamounts_5km.nc", 
  "Distance_Small_Seamounts_5km.nc"
)
robinson_crs <- "+proj=robin +lon_0=180 +datum=WGS84 +units=m +no_defs"

r_stack <- rast(nc_files)
names(r_stack) <- c("All Seamounts", "Big Seamounts", "Medium Seamounts", "Small Seamounts")

cat("Aggregating raster for faster plotting...\n")
r_agg <- aggregate(r_stack, fact = 2, fun = mean)
r_robinson <- project(
  crop(r_agg, ext(0, 360, -75, 75)),
  robinson_crs,
  method = "bilinear"
)

r_df <- as.data.frame(r_robinson, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = -c(x, y), names_to = "Category", values_to = "Distance")

# To prevent log10(0) returning -Inf (which turns into empty NA pixels), 
# we cap the minimum distance at 1 km. 
r_df$Distance <- ifelse(r_df$Distance < 1, 1, r_df$Distance)

# Create overlaid points (Assuming seamounts_360 already has 'size_group' from ntile)
pts_all <- seamounts_360 %>% mutate(Category = "All Seamounts")
pts_big <- seamounts_360 %>% filter(size_group == 3) %>% mutate(Category = "Big Seamounts")
pts_med <- seamounts_360 %>% filter(size_group == 2) %>% mutate(Category = "Medium Seamounts")
pts_small <- seamounts_360 %>% filter(size_group == 1) %>% mutate(Category = "Small Seamounts")

pts_combined <- bind_rows(pts_all, pts_big, pts_med, pts_small)
pts_combined <- st_transform(pts_combined, crs = robinson_crs)
coords <- st_coordinates(pts_combined)
pts_df <- pts_combined %>%
  st_drop_geometry() %>%
  bind_cols(x = coords[, "X"], y = coords[, "Y"])

# =============================================================================
# 2. DEFINE THE CAPTIONS
# =============================================================================
# These will be dynamically injected into the plot captions
captions <- list(
  "All Seamounts" = "Category Definition: All documented seamounts in the dataset.",
  "Big Seamounts" = "Category Definition: Top 33.3% of seamounts by height (Statistical Tertile).",
  "Medium Seamounts" = "Category Definition: Middle 33.3% of seamounts by height (Statistical Tertile).",
  "Small Seamounts" = "Category Definition: Bottom 33.3% of seamounts by height (Statistical Tertile)."
)

# =============================================================================
# 3. LOOP THROUGH CATEGORIES AND SAVE INDIVIDUAL MAPS
# =============================================================================
world_map <- st_as_sf(map("world2", plot = FALSE, fill = TRUE)) %>%
  st_set_crs(4326) %>%
  st_transform(crs = robinson_crs)
categories <- names(captions)

for (cat_name in categories) {
  
  cat(sprintf("Generating and saving map for %s...\n", cat_name))
  
  # Subset data for the current category
  r_sub <- r_df %>% filter(Category == cat_name)
  pts_sub <- pts_df %>% filter(Category == cat_name)
  
  # Build the individual plot
  p <- ggplot() +
    
    # Raster Heatmap
    geom_tile(data = r_sub, aes(x = x, y = y, fill = Distance)) +
    
    # Landmass
    geom_sf(data = world_map, fill = "grey15", color = "grey35", linewidth = 0.1) +
    
    # Points
    geom_point(data = pts_sub, aes(x = x, y = y), 
               color = "#ffcc00", size = 0.3, alpha = 0.7) +
    
    # Log10 Color Scale
    scale_fill_viridis_c(
      option = "mako", 
      direction = -1, 
      trans = "log10",                   # Applies the Log10 transformation
      breaks = trans_breaks("log10", function(x) 10^x), 
      labels = trans_format("log10", math_format(10^.x)), 
      name = "Distance (km) [Log10 Scale]",
      na.value = "#0b141f"
    ) +
    
    # Limits
    coord_sf(
      crs = robinson_crs,
      default_crs = robinson_crs,
      datum = 4326,
      expand = FALSE
    ) +
    
    # Labels (Title + Subtitle + Caption)
    labs(
      title = sprintf("Global Oceanic Proximity: %s", cat_name),
      subtitle = "Geodesic distance to the nearest seamount with specific locations overlaid",
      caption = captions[[cat_name]],
      x = "Longitude", 
      y = "Latitude"
    ) +
    
    # Professional Theme
    theme_minimal(base_size = 14) +
    theme(
    panel.background = element_rect(fill = dark_background, color = NA),
    plot.background = element_rect(fill = dark_background, color = NA),
      panel.grid = element_blank(),
      
    plot.title = element_text(face = "bold", size = 20, color = dark_text, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 13, color = muted_text, margin = margin(b = 10)),
    plot.caption = element_text(size = 11, face = "italic", color = muted_text, hjust = 1, margin = margin(t = 15)),
      
    axis.text = element_text(color = muted_text, size = 10),
    axis.title = element_text(color = dark_text, face = "bold"),
      
      legend.position = "bottom",
      legend.key.width = unit(4, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.title.position = "top",
    legend.background = element_rect(fill = dark_panel, color = NA),
    legend.title = element_text(hjust = 0.5, face = "bold", size = 13, color = dark_text),
    legend.text = element_text(size = 11, color = dark_text)
    )
  
  # Format file name (e.g., "Map_Big_Seamounts.png")
  safe_name <- gsub(" ", "_", cat_name)
  file_out <- sprintf("output/Map_%s_Log10.png", safe_name)
  
  # Save the individual map
  ggsave(file_out, plot = p, width = slide_width, height = slide_height, dpi = 300, bg = dark_background)
}

cat("All 4 individual Log10 maps have been saved successfully to the output folder!\n")
