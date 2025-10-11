#install.packages("arrow")
#install.packages("s3")
#-------------------------------------------------------------------------------------
library(terra)
library(s3)
library(sf)
library(dplyr)
library(arrow)
library(stringr)
library(purrr)
library(sp)
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

adm_sf <- st_as_sf(adm)
selected_tiles <- st_intersection(tiles_sf, adm_sf)
selected_tile_ids <- selected_tiles$tile_id

cat(sprintf("Found %d overlapping tiles with adm boundary\n", length(selected_tile_ids)))


# --- PARAMETERS -----------------------------------------------------------
path2gedi <- "s3://maap-ops-workspace/shared/ameliah/gedi-test/brazil_tiles/data/"
cols_to_read <- c("lat_lowestmode", "lon_lowestmode")

yearsT1 <- c(2019, 2020)
yearsT2 <- c(2022, 2023)

min_count <- 1  # Minimum GEDI points per grid cell in both periods
GRID.coords <- data.frame()

# --- LOOP -----------------------------------------------------------------
for (i in seq_along(selected_tile_ids)) {
  tile_id <- selected_tile_ids[i]
  cat(sprintf("Processing tile %d of %d : %s\n", i, length(selected_tile_ids), tile_id))

  # --- Helper to safely read parquet files from S3 ---
  safe_read <- function(p, cols) {
  tryCatch(
    read_parquet(p, col_select = all_of(cols)),
    error = function(e) {
      cat(sprintf("  → Could not read %s: %s\n", p, e$message))
      return(NULL)
      }
    )
  }

  # --- Read GEDI points for T1 ---
  paths_T1 <- paste0(path2gedi, "tile_id=", tile_id, "/year=", yearsT1, "/data_0.parquet")
  tile_T1 <- purrr::map_dfr(paths_T1, safe_read, cols = cols_to_read)

  if (nrow(tile_T1) == 0) {
    cat(sprintf("  → Tile %s: no points for T1, skipping.\n", tile_id))
    next
  }
    
  gedi_T1 <- vect(tile_T1, geom = c("lon_lowestmode", "lat_lowestmode"),
                  crs = "EPSG:4326", keepgeom = FALSE)
  gedi_T1_prj <- project(gedi_T1, "EPSG:6933")
  gcount_T1 <- rasterize(geom(gedi_T1_prj)[, c("x", "y")],
                         GRID.lons.adm.m, fun = "count", background = NA)

  # --- Read GEDI points for T2 ---
  paths_T2 <- paste0(path2gedi, "tile_id=", tile_id, "/year=", yearsT2, "/data_0.parquet")
  tile_T2 <- purrr::map_dfr(paths_T2, safe_read, cols = cols_to_read)

  if (nrow(tile_T2) == 0) {
    cat(sprintf("  → Tile %s: no points for T2, skipping.\n", tile_id))
    next
  }
  
  gedi_T2 <- vect(tile_T2, geom = c("lon_lowestmode", "lat_lowestmode"),
                  crs = "EPSG:4326", keepgeom = FALSE)
  gedi_T2_prj <- project(gedi_T2, "EPSG:6933")
  gcount_T2 <- rasterize(geom(gedi_T2_prj)[, c("x", "y")],
                         GRID.lons.adm.m, fun = "count", background = NA)

  # --- Identify cells with >= min_count points in both T1 & T2 ---
  r_both <- mask(gcount_T1 >= min_count & gcount_T2 >= min_count, GRID.lons.adm.m)

  # --- Check if there are any valid cells before converting to points ---
  valid_cells <- global(!is.na(r_both), "sum", na.rm = TRUE)[1, 1]
  cat(sprintf(" %s valid overlap cells found.\n", valid_cells))
  if (is.na(valid_cells) || valid_cells == 0) {
    cat(sprintf("  → Tile %s: no valid overlap cells.\n", tile_id))
    next
  }

  # --- Conversion to points --- each raster cell becomes a point at its centroid (x, y).
  pts <- try(as.points(r_both, values = FALSE), silent = TRUE)
  if (inherits(pts, "try-error") || nrow(pts) == 0) {
    cat(sprintf("  → Tile %s: empty raster, skipping.\n", tile_id))
    next
  }

  pts_ll <- project(pts, "EPSG:4326")
  coords <- try(terra::crds(pts_ll), silent = TRUE)
  xy_overlap <- as.data.frame(coords[, 1:2, drop = FALSE])
  colnames(xy_overlap) <- c("x", "y")
    
  GRID.coords <- rbind(GRID.coords, xy_overlap)

  cat(sprintf("  → Tile %s: added %d cells.\n", tile_id, nrow(xy_overlap)))
}

# --- Deduplicate grid cells across tiles ---------------------------------
GRID.coords <- unique(GRID.coords)
cat(sprintf("\n✅ Finished: %d grid cells meet min_count = %d in both periods.\n", nrow(GRID.coords), min_count))


# --- Save output to file -------------------------------------------------
GRID.for.matching <- vect(GRID.coords, geom=c("x","y"), crs = "epsg:4326")

filename_out <- paste0("output/", iso3, "_grid.RDS")
#filename_out <- paste0(f.path, "INPUT_grids/", iso3, "_grid.RDS")
print(filename_out)

saveRDS(GRID.for.matching, file = filename_out)

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# STEP2. Clip sampling grid to nonPA areas within country & sample raster layers on nonPA grid
#-------------------------------------------------------------------------------------
#GRID.for.matching <- readRDS(s3_get(paste0(s3.path, "INPUT_grids/", iso3, "_grid.RDS"), force=TRUE))

GRID.pts.nonPA <- project(GRID.for.matching, "epsg:4326")

  for(i in 1:length(allPAs)){
    PA          <- vect(allPAs[i,])
    PA_prj      <- project(PA, "epsg:6933")
    PA_prj_buff <- terra::buffer(PA_prj, width = 10000) ##10km buffer
    PA2         <- project(PA_prj_buff, "epsg:4326")
    overlap     <- GRID.pts.nonPA[PA2]
    if(length(overlap)>0){
      GRID.pts.nonPA0 <- st_difference(sf::st_as_sf(GRID.pts.nonPA), sf::st_as_sf(PA2)) ##remove pts inside poly
      GRID.pts.nonPA <- vect(GRID.pts.nonPA0$geometry)
      GRID.pts.nonPA <- project(GRID.pts.nonPA, "epsg:4326")
    } 
    print(length(GRID.pts.nonPA))
  }

nonPA_xy  <- geom(GRID.pts.nonPA)[,c("x","y")]
  colnames(nonPA_xy)  <- c("x","y")
  nonPA_spdf  <- tryCatch(vect(nonPA_xy, crs="epsg:4326"),      
                          error=function(cond){
                            cat("Country too small - quit processing ", iso3, dim(nonPA_xy),"\n")
                            return(quit(save="no"))})

for (j in 1:length(matching_tifs)){
    ras <- rast(s3_get(paste0(s3.path, "INPUT_covariates_2020/", matching_tifs[j], ".tif")))
    print(matching_tifs[j])
    ras_ex <- extract(ras, nonPA_spdf, method="simple", factors=FALSE)
    nm <- names(ras)
    nonPA_spdf$nm <- ras_ex[, matching_tifs[j]]
    names(nonPA_spdf)[j] <- matching_tifs[j]
  }

nonPA_spdf$x <- geom(nonPA_spdf)[,"x"]
nonPA_spdf$y <- geom(nonPA_spdf)[,"y"]

head(nonPA_spdf)

d_control <- nonPA_spdf
d_control$status <- as.logical("FALSE")
names(d_control) <- make.names(names(d_control), allow_ = FALSE)

d_control <- data.frame(d_control) %>%
    dplyr::rename(
      ### land_cover = lc2000,
      gedi = gedi.l4b,
      land_cover = MapBiomas.brasil.coverage.2020,
      slope = slope,
      elevation = dem,
      popden = pop.den.2020,
      popcnt = pop.cnt.2020,
      min_temp = wc.tmin.2010.2018,
      max_temp = wc.tmax.2010.2018,
      mean_temp = wc.tavg.2010.2018,
      prec = wc.prec.2010.2018,
      tt2city = tt2cities.2015,
      ### wwfbiom = wwf.biomes,
      ### wwfecoreg = wwf.ecoreg,
      d2city = dcities,
      d2road = d2roads,
      lon = x,
      lat = y)
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
  
  d_control$UID <-  seq.int(nrow(d_control))
   
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
    testPA <- project(testPA, "epsg:4326")
    GRID.pts.testPA <- GRID.for.matching[testPA]
    
    #if(length(GRID.pts.testPA)>0){
    if(length(GRID.pts.testPA)>1){
      testPA_xy <- geom(GRID.pts.testPA)[,c("x","y")]
      colnames(testPA_xy) <- c("x","y")
      testPA_spdf  <- vect(testPA_xy, crs="epsg:4326")
                              
        for (j in 1:length(matching_tifs)){
        ras <- rast(s3_get(paste0(s3.path, "INPUT_covariates_2020/", matching_tifs[j], ".tif"), force=TRUE))
        ras <- crop(ras, testPA)
        ras_ex <- extract(ras, testPA_spdf, method="simple", factors=F)
        nm <- names(ras)
        testPA_spdf$nm <- ras_ex[, matching_tifs[j]]
        names(testPA_spdf)[j] <- matching_tifs[j]
      }
    
    testPA_spdf$x <- geom(testPA_spdf)[,"x"]
    testPA_spdf$y <- geom(testPA_spdf)[,"y"]
      
    d_pa <- testPA_spdf
    d_pa$status <- as.logical("TRUE")
    d_pa$DESIG_ENG <- testPA$DESIG_ENG
    d_pa$REP_AREA <- testPA$REP_AREA
    d_pa$PA_STATUS <- testPA$STATUS
    d_pa$PA_STATUSYR <- testPA$STATUS_YR
    d_pa$GOV_TYPE <- testPA$GOV_TYPE
    d_pa$OWN_TYPE <- testPA$OWN_TYPE
    d_pa$MANG_AUTH <- testPA$MANG_AUTH
    names(d_pa) <- make.names(names(d_pa), allow_ = FALSE)
    
    d_pa <- data.frame(d_pa) %>%
            dplyr::rename(
            ### land_cover = lc2000,
            gedi = gedi.l4b,
            land_cover = MapBiomas.brasil.coverage.2020,
            slope = slope,
            elevation = dem,
            popden = pop.den.2020,
            popcnt = pop.cnt.2020,
            min_temp = wc.tmin.2010.2018,
            max_temp = wc.tmax.2010.2018,
            mean_temp = wc.tavg.2010.2018,
            prec = wc.prec.2010.2018,
            tt2city = tt2cities.2015,
            ### wwfbiom = wwf.biomes,
            ### wwfecoreg = wwf.ecoreg,
            d2city = dcities,
            d2road = d2roads,
            lon = x,
            lat = y)
      d_pa$land_cover <- factor(d_pa$land_cover, levels=sequence(10),
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
      
      d_pa$UID <- seq.int(nrow(d_pa))

      filename_out <- paste0("output/", iso3, "_prepped_pa_", testPA$WDPAID, ".RDS")
      #filename_out <- paste0(f.path, "/MATCHING_points/", iso3, "/", iso3, "_prepped_pa_", testPA$WDPAID, ".RDS")
      print(filename_out)
        
      saveRDS(d_pa, file = filename_out)  
    }
  }

png(paste0("output/", iso3, "_matching_points_map.png"), width = 1000, height = 1000, res = 300)
#png(paste0(f.path, iso3, "_matching_points_map.png"), width = 1000, height = 1000, res = 300)
plot(allPAs)
plot(GRID.pts.nonPA, pch=".", col="blue", add=T)
plot(adm, border="red", add=T)
dev.off()

