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


