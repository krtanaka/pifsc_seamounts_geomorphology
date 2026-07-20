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
base_path <- "C:/Users/Kisei.Tanaka/Desktop/ModelledSeamounts2011_v1/01_Data"

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

# Shared theme adjustments for spatial maps (Legend Top-Right / Over USA)
map_theme_adjust <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.caption = element_text(size = 7.5, color = "gray30", hjust = 0, lineheight = 1.1),
    panel.background = element_rect(fill = alpha("black", 0.8)), 
    panel.grid = element_blank(),
    legend.position = c(0.86, 0.85),             
    legend.background = element_rect(fill = alpha("white", 0.8), color = "white", linewidth = 0.3), 
    legend.key = element_blank(), # Ensures individual legend keys are transparent
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

# Shared theme adjustments for hex plots (Legend Top-Left)
hex_theme_adjust <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.caption = element_text(size = 7.5, color = "gray30", hjust = 0, lineheight = 1.1),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = c(0.12, 0.78),             
    legend.background = element_rect(fill = alpha("white", 0.6), color = "white", linewidth = 0.3),
    legend.key = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# ==============================================================================
# 3. GENERATE & SAVE MAPS INDIVIDUALLY
# ==============================================================================

# Map 1: Seamounts by Height
ggplot() +
  geom_sf(data = world_cropped,fill = "gray90", color = "gray20", linewidth = 0.1) +
  geom_sf(data = seamounts_360, aes(fill = HEIGHT), color = "gray", stroke = 0.5, size = 2, alpha = 0.8, shape = 21) +
  scale_fill_viridis_c(option = "plasma", name = "Height (m)", trans = "log10") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Seamounts by Height", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/seamounts_by_height.png", width = 8, height = 9, dpi = 300, bg = "white")

# Map 2: Knolls by Height
ggplot() +
  geom_sf(data = world_cropped, fill = "gray90", color = "gray20", linewidth = 0.1) +
  geom_sf(data = knolls_360, aes(color = HEIGHT), size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Height (m)", trans = "log10") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Knolls by Height", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/knolls_by_height.png", width = 8, height = 9, dpi = 300, bg = "white")

# Map 3: Seamounts by Summit Depth
ggplot() +
  geom_sf(data = world_cropped, fill = "gray90", color = "gray20", linewidth = 0.1) +
  geom_sf(data = seamounts_360, aes(fill = DEPTH), color = "gray", stroke = 0.5, size = 2, alpha = 0.8, shape = 21) +
  scale_fill_viridis_c(option = "mako", name = "Depth (m)") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Seamounts by Summit Depth", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/seamounts_by_depth.png", width = 8, height = 9, dpi = 300, bg = "white")

# Map 4: Knolls by Summit Depth
ggplot() +
  geom_sf(data = world_cropped,fill = "gray90", color = "gray20", linewidth = 0.1) +
  geom_sf(data = knolls_360, aes(color = DEPTH), size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "mako", name = "Depth (m)") +
  pacific_view + map_theme_adjust +
  labs(title = "Pacific Knolls by Summit Depth", x = "Longitude", y = "Latitude", caption = metadata_caption)
ggsave("output/knolls_by_depth.png", width = 11, height = 7, dpi = 300, bg = "white")


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
ggsave("output/seamounts_morphology_hex.png", width = 7, height = 7, dpi = 300, bg = "white")

# Plot 6: Knolls Morphological Density
ggplot(knolls_360, aes(x = DEPTH, y = HEIGHT)) +
  geom_hex(bins = 50, alpha = 0.9) +
  scale_fill_viridis_c(option = "plasma", name = "Count") +
  hex_theme_adjust + 
  labs(title = "Knolls: Morphological Density Structure", 
       x = "Depth", y = "Height", caption = metadata_caption)
ggsave("output/knolls_morphology_hex.png", width = 7, height = 7, dpi = 300, bg = "white")

# 1. Load the necessary libraries
library(sf)
library(terra)

# (Assuming your dataset is currently an sf object named 'seamounts_360')

# 2. Convert the sf point object to a terra SpatVector
seamounts_vect <- vect(seamounts_360)

# 3. Create a template global raster
# Since your data uses a 0 to 360 longitude range (based on the bounding box),
# we define the extent as xmin=0, xmax=360, ymin=-90, ymax=90.
# The resolution is set to 0.5 degrees (~55km at the equator). 
# Note: Decrease this number (e.g., 0.1) for a smoother, higher-resolution map, 
# but it will take longer to compute.
template_raster <- rast(ext(0, 360, -90, 90), resolution = 0.5, crs = "EPSG:4326")

# 4. Calculate the distance from all raster cells to the nearest point
# Because the CRS is unprojected (WGS 84), terra::distance automatically 
# calculates geodesic distances in meters.
dist_raster <- distance(template_raster, seamounts_vect)

# 5. Convert the distance from meters to kilometers for easier interpretation
dist_raster_km <- dist_raster / 1000
names(dist_raster_km) <- "Distance_to_Seamount_km"

# 6. Plot the results to visualize
plot(dist_raster_km, 
     main = "Distance to Nearest Seamount (km)",
     col = mapview::mapviewGetOption("raster.palette")(100)) # Optional: smooth color palette
plot(seamounts_vect, add = TRUE, col = "black", cex = 0.1, pch = 16)

library(ggplot2)
library(maps)   # For the world2 base map
library(scales) # For clean number formatting in the legend

# 6. Convert the SpatRaster to a dataframe so ggplot2 can read it
# xy = TRUE keeps the longitude and latitude coordinates
dist_df <- as.data.frame(dist_raster_km, xy = TRUE)
colnames(dist_df) <- c("lon", "lat", "distance")

# Extract the 0-360 world map for the landmass overlay
world_map <- map_data("world2")

# 7. Build the Professional Plot
p_dist <- ggplot() +
  
  # A. The Raster Heatmap
  geom_raster(data = dist_df, aes(x = lon, y = lat, fill = distance)) +
  
  # B. The Landmass Overlay (Drawn on top of the raster)
  geom_polygon(
    data = world_map, 
    aes(x = long, y = lat, group = group),
    fill = "grey15", color = "grey30", linewidth = 0.2
  ) +
  
  # C. The Color Scale
  # 'mako' is a stunning dark-blue to bright-cyan palette perfect for ocean data.
  # direction = -1 makes distances close to 0 bright, and far distances dark.
  scale_fill_viridis_c(
    option = "mako", 
    direction = -1, 
    name = "Distance (km)", 
    labels = label_comma()
  ) +
  
  # D. Projections and Limits
  # Crops out the extreme poles (Antarctica) to focus on the global oceans
  coord_cartesian(xlim = c(0, 360), ylim = c(-75, 75), expand = FALSE) +
  
  # E. Labels
  labs(
    title = "Global Proximity to Seamounts",
    subtitle = "Geodesic distance to the nearest seamount (0.5° resolution grid)",
    x = "Longitude", 
    y = "Latitude"
  ) +
  
  # F. Dark Mode Theme
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#0b141f", color = NA),
    panel.grid = element_blank(), # Removes gridlines for a cleaner look
    
    plot.title = element_text(face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "grey40", margin = margin(b = 15)),
    axis.text = element_text(color = "grey50"),
    
    # Move the legend to the bottom and make it a wide color bar
    legend.position = "bottom",
    legend.key.width = unit(3, "cm"),
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5, face = "bold")
  )

# View the plot
print(p_dist)

# Optional: Save the high-resolution output
ggsave("output/Global_Seamount_Proximity.png", plot = p_dist, width = 16, height = 9, dpi = 300)
