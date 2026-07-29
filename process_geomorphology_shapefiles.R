# ==============================================================================
# MASTER SEAFLOOR GEOMORPHOLOGY WORKFLOW (WITH LAYER CHECKING & SOURCE CITATION)
# Bypasses Vector Bugs, Parallelizes Spatial Data, and Exports Multiple Maps
# ==============================================================================

# 1. LOAD ALL LIBRARIES ONCE
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth) 
library(scales)
library(foreach)
library(doParallel)

# 2. GLOBAL PATHS & PARAMETERS
data_dir    <- "C://Users/Kisei.Tanaka/Desktop/global seafloor geomorphic features map/"
out_dir     <- "C://Users/Kisei.Tanaka/Desktop/rasterized_geomorph/"
pacific_crs <- "+proj=robin +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

dir.create(out_dir, showWarnings = FALSE)

# Master Stacking Order (From Background Abyss to Foreground Features)
ordered_layers <- c(
  "Abyss", "Abyssal_Classification", "Basins", "Shelf", "Shelf_Classification", 
  "Slope", "Plateaus", "Terraces", "Rises", "Fans", "Bridges", "Sills", 
  "Troughs", "Glacial_troughs", "Ridges", "Spreading_ridges", "Escarpments", 
  "Canyons", "Shelf_valleys", "Rift_valleys", "Guyots", "Seamounts", 
  "Trenches", "Hadal"
)

# ==============================================================================
# PHASE A: BATCH RASTERIZATION (Vector to GeoTIFF - With Smart Skip Check)
# ==============================================================================
shp_paths <- list.files(data_dir, pattern = "\\.shp$", full.names = TRUE)
layer_names <- tools::file_path_sans_ext(basename(shp_paths))

if (length(shp_paths) > 0) {
  message("Checking and rasterizing shapefiles to 0.05° (~5.5km) template...")
  template_rast <- rast(ext(-180, 180, -90, 90), resolution = 0.05, crs = "EPSG:4326")
  
  for (i in seq_along(shp_paths)) {
    name <- layer_names[i]
    out_file <- paste0(out_dir, name, ".tif")
    
    # SMART SKIP: Bypass processing if the output GeoTIFF already exists
    if (file.exists(out_file)) {
      cat("Skipping (already rasterized):", name, "\n")
      next
    }
    
    cat("Processing raw vector:", name, "\n")
    v <- vect(shp_paths[i])
    v$present <- 1
    r <- rasterize(v, template_rast, field = "present", background = NA)
    writeRaster(r, out_file, overwrite = TRUE)
  }
}

# ==============================================================================
# PHASE B: SPATIAL DATA PREPARATION (Run Once)
# ==============================================================================

# 1. Clean Landmass Boundaries with Full Polar Corrections
message("Preparing landmass boundaries...")
world_land <- ne_countries(scale = "medium", returnclass = "sf")
world_land_cropped <- st_crop(world_land, xmin = -180, xmax = 180, ymin = -85, ymax = 85)

sf_use_s2(FALSE) 
land_pac <- world_land_cropped %>% 
  st_make_valid() %>% 
  st_transform(crs = pacific_crs) %>% 
  st_collection_extract("POLYGON") 
sf_use_s2(TRUE)

# 2. Turbocharged Parallel Raster Processing Engine
num_cores <- parallel::detectCores() - 1
message("Spawning parallel cluster using ", num_cores, " CPU cores...")

cl <- makeCluster(num_cores)
registerDoParallel(cl)

message("Processing 24 layers in parallel (Reprojecting & Downsampling)...")
df_list <- foreach(
  layer = ordered_layers, 
  .packages = c("terra"), 
  .export = c("out_dir", "pacific_crs")
) %dopar% {
  file_path <- paste0(out_dir, layer, ".tif")
  
  if (file.exists(file_path)) {
    r       <- rast(file_path)
    r_proj  <- project(r, pacific_crs, method = "near")
    r_agg   <- aggregate(r_proj, fact = 5, fun = "max") 
    df      <- as.data.frame(r_agg, xy = TRUE, na.rm = TRUE)
    
    if (nrow(df) > 0) {
      df$Feature <- gsub("_", " ", layer)
      return(df[, c("x", "y", "Feature")])
    }
  }
  return(NULL)
}

stopCluster(cl) 
message("Parallel extraction complete. Assembling master dataset...")

df_list    <- Filter(Negate(is.null), df_list)
df_all_pac <- do.call(rbind, df_list)

# 3. Premium High-Contrast Thematic Palette
vibrant_colors <- c(
  # Deep Ocean Floor Family (Quiet Blues)
  "Abyss"                  = "#122643", "Abyssal Classification" = "#1a365d",
  "Basins"                 = "#22467b", "Rises"                  = "#2e5b9a",
  
  # Shallow & Coastlines Family (Muted Greens/Teals)
  "Shelf"                  = "#2d6a4f", "Shelf Classification"   = "#40916c",
  "Shelf valleys"          = "#52b788", "Slope"                  = "#1b4332",
  "Plateaus"               = "#74c69d", "Terraces"               = "#95d5b2", "Fans" = "#b7e4c7",
  
  # Volcanic High-Relief Family (Warm Accents)
  "Ridges"                 = "#d68c45", "Spreading ridges"       = "#ffb703",
  "Rift valleys"           = "#fb8500", "Escarpments"            = "#ee9b00",
  "Guyots"                 = "#ca6702", "Seamounts"              = "#ae2012",
  
  # Dynamic Structural Cuts (Neon Cyans, Purples & Whites)
  "Bridges"                = "#9b5de5", "Sills"                  = "#f15bb5",
  "Troughs"                = "#00bbf9", "Glacial troughs"        = "#00f5d4",
  "Canyons"                = "#00f5d4", "Trenches"               = "#00ffff", "Hadal" = "#ffffff"
)

# 4. Map Citation text from Document Source
source_text <- "Source: Harris, P.T. & McMillan-Lawler, M. (GRID-Arendal). Global Seafloor Geomorphic Features Map (GSFM), interpreted from the SRTM 30 arc-second (~1 km) global bathymetry grid (Harris et al., 2014)."

# ==============================================================================
# PHASE C: MAP RENDERING & EXPORTS
# ==============================================================================

# --- MAP 1: FULL GLOBAL PACIFIC VIEW ---
message("Rendering Global Pacific map...")
p_global <- ggplot() +
  geom_raster(data = df_all_pac, aes(x = x, y = y, fill = Feature)) +
  geom_sf(data = land_pac, fill = "#0a0d14", color = "#1c2331", linewidth = 0.2, inherit.aes = FALSE) +
  coord_sf(crs = pacific_crs, datum = NA) + 
  scale_fill_manual(values = vibrant_colors, guide = guide_legend(nrow = 2)) +
  base_map_theme +
  labs(
    title    = "Pacific Ocean Seafloor Geomorphology",
    caption  = source_text
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#020408", color = NA),
    plot.background  = element_rect(fill = "#020408", color = NA),
    plot.title       = element_text(face = "bold", size = 20, hjust = 0.5, color = "#ffffff"),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "#8b9bb4", margin = margin(b = 15)),
    plot.caption     = element_text(size = 8, color = "#54657e", hjust = 0, margin = margin(t = 15, l = 20)),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    
    legend.position   = "bottom",
    legend.background = element_rect(fill = "#0a0d14", color = "#1c2331", linewidth = 0.5),
    legend.key        = element_blank(),
    legend.text       = element_text(color = "#ffffff", size = 8, face = "bold"),
    legend.title      = element_blank(),
    legend.key.size   = unit(0.35, "cm") 
  )

ggsave(filename = "output/Pacific_Geomorphology.png", 
       plot = p_global, width = 16, height = 9, dpi = 300)

# --- MAP 2: TIGHTLY FOCUSED PACIFIC BASIN VIEW ---
message("Rendering Focused Basin map...")
p_focused <- ggplot() +
  geom_raster(data = df_all_pac, aes(x = x, y = y, fill = Feature)) +
  geom_sf(data = land_pac, fill = "#0a0d14", color = "#1c2331", linewidth = 0.2, inherit.aes = FALSE) +
  coord_sf(
    crs = pacific_crs, 
    xlim = c(-6500000, 11500000), 
    ylim = c(-6000000, 6500000), 
    datum = NA, 
    expand = FALSE
  ) + 
  scale_fill_manual(values = vibrant_colors, guide = guide_legend(ncol = 1)) +
  base_map_theme +
  labs(
    title    = "Pacific Ocean Basin Seafloor Geomorphology",
    subtitle = "Focused 24-feature bathymetric classification of the Pacific Rim",
    caption  = source_text
  ) + 
  # Shared Core Theme Engine (With Single-Column Legend Adaptations)
  theme_minimal(base_size = 30) +
  theme(
    panel.background = element_rect(fill = "#020408", color = NA),
    plot.background  = element_rect(fill = "#020408", color = NA),
    plot.title       = element_text(face = "bold", size = 20, hjust = 0.5, color = "#ffffff"),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "#8b9bb4", margin = margin(b = 15)),
    plot.caption     = element_text(size = 8, color = "#54657e", hjust = 0, margin = margin(t = 15, l = 20)),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    
    legend.position   = "right",
    legend.background = element_rect(fill = "#0a0d14", color = "#1c2331", linewidth = 0.5),
    legend.key        = element_blank(),
    legend.text       = element_text(color = "#ffffff", size = 8, face = "bold"),
    legend.title      = element_blank(),
    legend.key.size   = unit(0.35, "cm") 
  )

ggsave(filename = "output/Pacific_Basin_Focused_Map.png", 
       plot = p_focused, width = 16, height = 9, dpi = 300)

message("Process complete! Both maps updated and saved directly to your Desktop.")
