#install.packages("s3")
#-------------------------------------------------------------------------------------
library(terra)
library(s3)
library(sf)
library(dplyr)
library(stringr)
library(purrr)
library(sp)
library(DBI)
library(duckdb)
#-------------------------------------------------------------------------------------

#iso3 <- "BRA"

#-------------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) { 
  iso3 <- args[1]  #country to process
}
#-------------------------------------------------------------------------------

f.path <- "/projects/my-public-bucket/GEDI_global_change/"
s3.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_change/"
path2gedi <- "s3://maap-ops-workspace/shared/ameliah/gedi-test/brazil_tiles/data/"

s3_get_files(c(paste0(s3.path, "INPUT_countries/shp/", iso3, ".shp"),
               paste0(s3.path, "INPUT_countries/shp/", iso3, ".shx"),
               paste0(s3.path, "INPUT_countries/shp/", iso3, ".prj"),
               paste0(s3.path, "INPUT_countries/shp/", iso3, ".dbf")), confirm = FALSE)
#adm <- st_read(s3_get(paste(f.path,"WDPA_countries/shp/",iso3,".shp",sep=""), force=TRUE))
#adm_prj <- project(vect(adm), "epsg:6933")
adm <- vect(s3_get(paste0(s3.path, "INPUT_countries/shp/", iso3, ".shp")))
crs(adm) <- "EPSG:4326"
adm_prj <- project(adm, "EPSG:6933")

grid_rast <- rast(s3_get(paste0(s3.path, "GEDI04_B_MW019MW223_02_002_02_R01000M_MU.tif")))
grid_crop <- crop(grid_rast, adm_prj)
grid_mask <- mask(grid_crop, adm_prj)
#grid_mask

GRID.lats <- rast(s3_get(paste0(s3.path,"EASE2_M01km_lats.tif"), force=TRUE))
GRID.lons <- rast(s3_get(paste0(s3.path,"EASE2_M01km_lons.tif"), force=TRUE))
GRID.lats.adm   <- crop(GRID.lats, adm_prj)
GRID.lats.adm.m <- mask(GRID.lats.adm, adm_prj)
GRID.lons.adm   <- crop(GRID.lons, adm_prj)
GRID.lons.adm.m <- mask(GRID.lons.adm, adm_prj)
template_rast <- GRID.lons.adm.m

allPAs <- readRDS(s3_get(paste0(s3.path, "INPUT_shapefiles/", iso3, "_PA_poly.rds"), force=TRUE))

matching_tifs <- c("d2roads", "dcities", "dem", "slope",
                   "pop_cnt_2020", "pop_den_2020", "tt2cities_2015",
                   "wc_prec_2010-2018", "wc_tavg_2010-2018",
                   "wc_tmax_2010-2018", "wc_tmin_2010-2018",
                   "MapBiomas_brasil_coverage_2020", "gedi_l4b")

#-------------------------------------------------------------------------------------
# Get country adm bbox & extract bounds
bbox <- st_bbox(adm)

# Generate all possible intersecting tile IDs from coordinates
generate_tile_id <- function(lat, lon) {
  lat_prefix <- ifelse(lat >= 0, "N", "S")
  lon_prefix <- ifelse(lon >= 0, "E", "W")
  lat_val <- sprintf("%02d", abs(floor(lat)))
  lon_val <- sprintf("%03d", abs(floor(lon)))
  paste0(lat_prefix, lat_val, "_", lon_prefix, lon_val)
}

# Candidate tile ranges
lat_range <- seq(floor(bbox["ymin"]), ceiling(bbox["ymax"]))
lon_range <- seq(floor(bbox["xmin"]), ceiling(bbox["xmax"]))

# Generate all tile IDs in bbox
tile_ids <- expand.grid(lat = lat_range, lon = lon_range) %>%
  mutate(tile_id = generate_tile_id(lat, lon)) %>%
  pull(tile_id)

# Function to create a 1° tile polygon from lat/lon with tile_id attribute
tile_to_polygon <- function(tile_id) {
  # Extract latitude
  lat_str <- str_extract(tile_id, "[NS]\\d+")
  lat <- as.numeric(str_remove(lat_str, "[NS]"))
  if (str_starts(lat_str, "S")) lat <- -lat
  
  # Extract longitude
  lon_str <- str_extract(tile_id, "[EW]\\d+")
  lon <- as.numeric(str_remove(lon_str, "[EW]"))
  if (str_starts(lon_str, "W")) lon <- -lon
  
  # Create matrix of coordinates (must close polygon)
  coords <- rbind(
    c(lon, lat),         # top-left
    c(lon + 1, lat),     # top-right
    c(lon + 1, lat - 1), # bottom-right
    c(lon, lat - 1),     # bottom-left
    c(lon, lat)          # close polygon
  )
  
  # Return an sfg object
  st_polygon(list(coords))
}

# Convert each tile_id to an sfg polygon
polygons <- lapply(tile_ids, tile_to_polygon)

# Make sure st_sfc gets a list of sfg objects
tiles_sf <- st_sf(
  tile_id = tile_ids,
  geometry = st_sfc(polygons, crs = 4326)
)

#adm_sf <- st_as_sf(adm)
sf_use_s2(FALSE) # Temporarily switch off S2 for planar 
adm_sf <- st_as_sf(adm)
adm_sf <- st_make_valid(adm_sf)
adm_sf <- st_union(adm_sf)
adm_sf <- st_collection_extract(adm_sf, "POLYGON")
selected_tiles <- st_intersection(tiles_sf, adm_sf)
selected_tile_ids <- selected_tiles$tile_id
sf_use_s2(TRUE) # Restore normal behavior

cat(sprintf("Found %d overlapping tiles with adm boundary\n", length(selected_tile_ids)))

# --- PARAMETERS -----------------------------------------------------------
path2gedi <- "s3://maap-ops-workspace/shared/ameliah/gedi-test/brazil_tiles/data/"
cols_to_read <- c("lat_lowestmode", "lon_lowestmode")

yearsT1 <- c(2019, 2020)
yearsT2 <- c(2022, 2023)

min_count <- 1  # Minimum GEDI points per grid cell in both periods

coords_list <- vector("list", length(selected_tile_ids))
list_idx <- 1

# --- duckdb -----------------------------------------------------------------
con <- dbConnect(duckdb(), dbdir = ":memory:")

safe_read_duck <- function(p, cols) {
  tryCatch({
    # Download to a local temporary path
    local_path <- s3_get(p, force = TRUE)
    col_list <- paste(cols, collapse = ", ")
    q <- sprintf("SELECT %s FROM read_parquet('%s')", col_list, local_path)
    df <- dbGetQuery(con, q)
    if (nrow(df) == 0) {
      cat(sprintf("  → %s: empty parquet file\n", p))
      return(NULL)
    }
    return(df)
  }, error = function(e) {
    cat(sprintf("  → Could not read %s: %s\n", p, e$message))
    return(NULL)
  })
}

# --- LOOP -----------------------------------------------------------------
for (i in seq_along(selected_tile_ids)) {
  tile_id <- selected_tile_ids[i]
  cat(sprintf("Processing tile %d of %d : %s\n", i, length(selected_tile_ids), tile_id))
    
  # --- Read GEDI points for T1 ---
  paths_T1 <- paste0(path2gedi, "tile_id=", tile_id, "/year=", yearsT1, "/data_0.parquet")
  tile_T1 <- purrr::map_dfr(paths_T1, safe_read_duck, cols = cols_to_read)
  if (is.null(tile_T1) || nrow(tile_T1) == 0) next
  
  gedi_T1 <- vect(tile_T1, geom = c("lon_lowestmode", "lat_lowestmode"), crs = "EPSG:4326")
  gedi_T1_prj <- project(gedi_T1, "EPSG:6933")
  
  # --- Read GEDI points for T2 ---
  paths_T2 <- paste0(path2gedi, "tile_id=", tile_id, "/year=", yearsT2, "/data_0.parquet")
  tile_T2 <- purrr::map_dfr(paths_T2, safe_read_duck, cols = cols_to_read)
  if (is.null(tile_T2) || nrow(tile_T2) == 0) next
  
  gedi_T2 <- vect(tile_T2, geom = c("lon_lowestmode", "lat_lowestmode"), crs = "EPSG:4326")
  gedi_T2_prj <- project(gedi_T2, "EPSG:6933")
  
  # --- Rasterize points on template raster ---
  gcount_T1 <- rasterize(gedi_T1_prj, template_rast, fun = "count", background = NA)
  gcount_T2 <- rasterize(gedi_T2_prj, template_rast, fun = "count", background = NA)
  
  # --- Identify cells with >= min_count points in both periods ---
  r_both <- mask(gcount_T1 >= min_count & gcount_T2 >= min_count, template_rast)
  
  # --- Check for valid cells ---
  valid_cells <- global(!is.na(r_both), "sum", na.rm = TRUE)[1,1]
  if (is.na(valid_cells) || valid_cells == 0) {
    cat(sprintf("  → Tile %s: no valid overlap cells.\n", tile_id))
    next
  }
  
  # --- Convert raster to points (cell centroids) ---
  pts <- try(as.points(r_both, values = FALSE), silent = TRUE)
  if (inherits(pts, "try-error") || nrow(pts) == 0) {
    cat(sprintf("  → Tile %s: empty raster, skipping.\n", tile_id))
    next
  }
  
  # --- Convert back to lat/lon and store coordinates ---
  pts_ll <- project(pts, "EPSG:4326")
  coords <- crds(pts_ll)
  xy_overlap <- as.data.frame(coords[,1:2, drop = FALSE])
  colnames(xy_overlap) <- c("x","y")
  
  coords_list[[list_idx]] <- xy_overlap
  list_idx <- list_idx + 1
  
  cat(sprintf("  → Tile %s: added %d cells.\n", tile_id, nrow(xy_overlap)))
}

dbDisconnect(con, shutdown = TRUE)

# --- Combine all points and deduplicate ---
GRID.coords <- unique(do.call(rbind, coords_list))
cat(sprintf("\n✅ Finished: %d grid cells meet min_count = %d in both periods.\n", nrow(GRID.coords), min_count))

# --- Save output to file -------------------------------------------------
GRID.for.matching <- vect(GRID.coords, geom=c("x","y"), crs = "epsg:4326")

filename_out <- paste0("output/", iso3, "_grid.RDS")
#filename_out <- paste0(f.path, "INPUT_grids/", iso3, "_grid.RDS")
print(filename_out)

saveRDS(GRID.for.matching, file = filename_out)

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# STEP2. Clip sampling grid to nonPA areas & sample raster layers
#-------------------------------------------------------------------------------------
#GRID.for.matching <- readRDS(s3_get(paste0(s3.path, "INPUT_grids/", iso3, "_grid.RDS"), force=TRUE))

GRID.pts.nonPA <- project(GRID.for.matching, "epsg:4326")

# Project all PAs and buffer once
allPAs_vect <- vect(allPAs)
allPAs_prj  <- project(allPAs_vect, "EPSG:6933")
allPAs_buff <- buffer(allPAs_prj, width = 10000)  # 10 km buffer

# Combine all PAs into a single polygon
allPAs_union <- st_union(st_as_sf(allPAs_buff))
allPAs_union <- st_make_valid(allPAs_union)

# Project back to lat/lon
allPAs_union_ll <- st_transform(allPAs_union, 4326)

# Remove all nonPA points in one operation
GRID.pts.nonPA_sf <- st_as_sf(GRID.pts.nonPA)
GRID.pts.nonPA_sf <- st_difference(GRID.pts.nonPA_sf, allPAs_union_ll)
GRID.pts.nonPA <- vect(GRID.pts.nonPA_sf)

# Convert to XY dataframe
nonPA_xy <- geom(GRID.pts.nonPA)[,c("x","y")]
colnames(nonPA_xy) <- c("x","y")
nonPA_spdf <- tryCatch(vect(nonPA_xy, crs="EPSG:4326"),      
                          error=function(cond){
                            cat("Country too small - quit processing ", iso3, dim(nonPA_xy),"\n")
                            return(quit(save="no"))})

# Extract raster values for nonPA points
for (j in seq_along(matching_tifs)){
  ras <- rast(s3_get(paste0(s3.path, "INPUT_covariates_2020/", matching_tifs[j], ".tif")))
  ras_ex <- extract(ras, nonPA_spdf, method="simple", factors=FALSE)
  nonPA_spdf[[matching_tifs[j]]] <- ras_ex[, matching_tifs[j]]
}

# Add coordinates
nonPA_spdf$x <- geom(nonPA_spdf)[,"x"]
nonPA_spdf$y <- geom(nonPA_spdf)[,"y"]

# Convert to dataframe and rename columns
d_control <- data.frame(nonPA_spdf)
d_control$status <- FALSE
names(d_control) <- make.names(names(d_control), allow_ = TRUE)

# Rename variables
d_control <- data.frame(d_control) %>%
    dplyr::rename(
      ### land_cover = lc2000,
      gedi_l4b = gedi_l4b,
      land_cover = MapBiomas_brasil_coverage_2020,
      slope = slope,
      elevation = dem,
      popden = pop_den_2020,
      popcnt = pop_cnt_2020,
      min_temp = wc_tmin_2010.2018,
      max_temp = wc_tmax_2010.2018,
      mean_temp = wc_tavg_2010.2018,
      prec = wc_prec_2010.2018,
      tt2city = tt2cities_2015,
      ### wwfbiom = wwf.biomes,
      ### wwfecoreg = wwf.ecoreg,
      d2city = dcities,
      d2road = d2roads,
      lon = x,
      lat = y)
# Factor land cover
d_control$land_cover <- factor(d_control$land_cover, levels=sequence(10),
                                 labels = c("l1_forest",
                                            "l2_savanna",
                                            "l3_mangrove",
                                            "l4_floodedforest",
                                            "l5_plantation",
                                            "l6_wetland",
                                            "l7_grassland",
                                            "l8_agriculture",
                                            "l9_nonvegetated",
                                            "l10_water"))
  
d_control$UID <- seq.int(nrow(d_control))
   
# Save
filename_out <- paste("output/", iso3, "_prepped_control.RDS")
#filename_out <- paste0(f.path, "/MATCHING_points/", iso3, "_prepped_control.RDS")
print(filename_out)

saveRDS(d_control, file = filename_out)  

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#STEP3. Loop through all PAs in iso3 country:
# - clip sampling grid to each PA
# - sample raster layers on each PA grid
# - save each PA sample into prepped_pa_##.RDS file
#-------------------------------------------------------------------------------------
#GRID.for.matching <- readRDS(s3_get(paste0(s3.path, "INPUT_grids/", iso3, "_grid.RDS"), force=TRUE))

for(i in 1:length(allPAs)){
    
  testPA <- vect(allPAs[i,])
  testPA <- project(testPA, "EPSG:4326")
    
  # Select points inside PA
  GRID.pts.testPA <- GRID.for.matching[testPA]
  if(length(GRID.pts.testPA) <= 1) next  # skip empty PAs

  # Get XY
  testPA_xy <- geom(GRID.pts.testPA)[,c("x","y")]
  colnames(testPA_xy) <- c("x","y")
  testPA_spdf <- vect(testPA_xy, crs="EPSG:4326")
  
  # Extract rasters
  for(j in seq_along(matching_tifs)){
    ras <- rast(s3_get(paste0(s3.path, "INPUT_covariates_2020/", matching_tifs[j], ".tif")))
    ras_crop <- crop(ras, testPA)
    ras_ex <- extract(ras_crop, testPA_spdf, method="simple", factors=FALSE)
    testPA_spdf[[matching_tifs[j]]] <- ras_ex[, matching_tifs[j]]
  }
  
  # Add coordinates
  testPA_spdf$x <- geom(testPA_spdf)[,"x"]
  testPA_spdf$y <- geom(testPA_spdf)[,"y"]
  
  # Convert to dataframe
  d_pa <- data.frame(testPA_spdf)
  d_pa$status <- TRUE
  
  # Add PA attributes
  d_pa$DESIG_ENG <- testPA$DESIG_ENG
  d_pa$REP_AREA <- testPA$REP_AREA
  d_pa$PA_STATUS <- testPA$STATUS
  d_pa$PA_STATUSYR <- testPA$STATUS_YR
  d_pa$GOV_TYPE <- testPA$GOV_TYPE
  d_pa$OWN_TYPE <- testPA$OWN_TYPE
  d_pa$MANG_AUTH <- testPA$MANG_AUTH
  
  # Rename columns
  names(d_pa) <- make.names(names(d_pa), allow_ = TRUE)
  d_pa <- d_pa %>%
    dplyr::rename(
      ### land_cover = lc2000,
      gedi_l4b = gedi_l4b,
      land_cover = MapBiomas_brasil_coverage_2020,
      slope = slope,
      elevation = dem,
      popden = pop_den_2020,
      popcnt = pop_cnt_2020,
      min_temp = wc_tmin_2010.2018,
      max_temp = wc_tmax_2010.2018,
      mean_temp = wc_tavg_2010.2018,
      prec = wc_prec_2010.2018,
      tt2city = tt2cities_2015,
      ### wwfbiom = wwf.biomes,
      ### wwfecoreg = wwf.ecoreg,
      d2city = dcities,
      d2road = d2roads,
      lon = x,
      lat = y)
  
  d_pa$land_cover <- factor(d_pa$land_cover, levels=1:10,
                            labels=c("l1_forest","l2_savanna","l3_mangrove",
                                     "l4_floodedforest","l5_plantation","l6_wetland",
                                     "l7_grassland","l8_agriculture","l9_nonvegetated","l10_water"))
  
  d_pa$UID <- seq.int(nrow(d_pa))

  # Save
  filename_out <- paste0("output/", iso3, "_prepped_pa_", testPA$WDPAID, ".RDS")
#  filename_out <- paste0(f.path, "/MATCHING_points/", iso3, "/", iso3, "_prepped_pa_", testPA$WDPAID, ".RDS")
  saveRDS(d_pa, file=filename_out)  
}

#-------------------------------------------------------------------------------------
png(paste0("output/", iso3, "_matching_points_map.png"), width = 1000, height = 1000, res = 300)
#png(paste0(f.path, iso3, "_matching_points_map.png"), width = 1000, height = 1000, res = 300)
plot(allPAs)
plot(GRID.pts.nonPA, pch=".", col="blue", add=T)
plot(adm, border="red", add=T)
dev.off()

