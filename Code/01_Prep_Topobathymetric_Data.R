#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### LIBRARIES ####
# install required packages (first run only)
# install.packages(c("easypackages", "conflicted", "dplyr", "here", "raster", "terra",
#                    "sp", "sf", "MultiscaleDTM", "fasterize", "ggplot2", "spatialEco"))

# load packages
library(easypackages)
libraries("conflicted", "dplyr", "here", "raster", "terra", "sp", "sf", 
          "MultiscaleDTM", "fasterize", "ggplot2", "spatialEco")

# resolve package conflicts
conflict_prefer("terrain", "terra")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# additional settings
rasterOptions(progress = 'text') # progress info for processing large rasters
options(terra.progress = 1) # progress info for processing large rasters
rasterOptions(tmpdir = "E:/temp_raster_directory") # custom directory for temporary files
terraOptions(tempdir = "E:/temp_raster_directory") # custom directory for temporary files

#### DIRECTORIES ####
# working directory and relative folder path for here()
#setwd("E:/Data/StuartC_DPhil_Ch3/")
#set_here("E:/Data/StuartC_DPhil_Ch3/") # set first-time only
here::i_am(".here")
here::here() # verify

#### COORDINATE SYSTEMS / PROJECTIONS ####
# save PROJ.4 string for WGS 84 / UTM Zone 6S (EPSG WKID 32706) with unit meters
my_crs = CRS("+proj=utm +zone=6 +south +datum=WGS84 +units=m +no_defs +type=crs")

# save PROJ.4 string for the standard geographic coordinate system used by
# Garmin GPS - WGS84 - World Geodetic System 1984 (EPSG WKID 4326) with unit decimal 
# degrees 
gcs = CRS("+proj=longlat +datum=WGS84 +no_defs +type=crs")

#### MOOREA ISLAND ####
# define bounding box coordinates
top = 8066650.744941
bottom = 8050292.944941
left = 188343.404160
right = 208243.844160

# create the polygon coordinates
bbox_coords = matrix(
  c(left, bottom,  # lower-left
    left, top,     # upper-left
    right, top,    # upper-right
    right, bottom, # lower-right
    left, bottom),   # close polygon
  ncol = 2,
  byrow = TRUE)

# create an sf polygon
bbox_poly = 
  st_polygon(list(bbox_coords)) |> 
  st_sfc(crs = 32706) |> 
  st_sf()

# read in the shapefile of all French Polynesia islands from Stanford University,
# available to download at:
# https://geowebservices.stanford.edu:443/geoserver/wfs?outputformat=SHAPE-ZIP&request=GetFeature&service=wfs&srsName=EPSG%3A4326&typeName=druid%3Abh400kc3500&version=2.0.0
fp = st_read(here("Data", "FP_Islands_Shapefile", "bh400kc3500.shp"))
compareCRS(fp, gcs) # should be in GCS WGS 84

# re-project to UTM zone 6S 
fp = st_transform(fp, crs(bbox_poly))
compareCRS(fp, bbox_poly)

# use the bounding box to isolate only Moorea
moorea = st_crop(fp, bbox_poly)

# keep only our desired information
moorea$Island = "Moorea"
moorea$ISO = "FYP"
moorea$English = "French Polynesia"
moorea$French = "Polynésie Française"
moorea$Group = "Society Islands"
moorea = moorea %>% 
  select(Island, ISO, English, French, Group, geometry)

# plot the geometry to make sure that everything looks good
plot(moorea$geometry)

# save the shapefile
st_write(moorea,
         dsn = here("Data", "FP_Islands_Shapefile", "Moorea.shp"),
         append = FALSE)

#### LTER SITES ####
LTER1 = st_polygon(list(matrix(
  c( -149.8455917, -17.48641792,  # lower-left
     -149.8455917, -17.47185366,  # upper-left
     -149.829821,  -17.47185366,  # upper-right
     -149.829821,  -17.48641792,  # lower-right
     -149.8455917, -17.48641792),   # close polygon
  ncol = 2,
  byrow = TRUE))) |>
  st_sfc(crs = 4326) |>
  st_sf() %>%
  st_transform(crs = my_crs)

LTER2 = st_polygon(list(matrix(
  c(-149.8116849, -17.48131958,  # lower-left
    -149.8116849, -17.46576169,  # upper-left
    -149.7961685, -17.46576169,  # upper-right
    -149.7961685, -17.48131958,  # lower-right
    -149.8116849, -17.48131958),   # close polygon
  ncol = 2,
  byrow = TRUE))) |>
  st_sfc(crs = 4326) |>
  st_sf() |>
  st_transform(crs = my_crs)

LTER3 = st_polygon(list(matrix(
  c(-149.7708619, -17.52087158,  # lower-left
    -149.7708619, -17.50382025,  # upper-left
    -149.7519968, -17.50382025,  # upper-right
    -149.7519968, -17.52087158,  # lower-right
    -149.7708619, -17.52087158),   # close polygon
  ncol = 2,
  byrow = TRUE))) |>
  st_sfc(crs = 4326) |>
  st_sf() |>
  st_transform(crs = my_crs)

LTER4 = st_polygon(list(matrix(
  c(-149.7772857, -17.55064263,  # lower-left
    -149.7772857, -17.53305021,  # upper-left
    -149.7566866, -17.53305021,  # upper-right
    -149.7566866, -17.55064263,  # lower-right
    -149.7772857, -17.55064263),   # close polygon
  ncol = 2,
  byrow = TRUE))) |>
  st_sfc(crs = 4326) |>
  st_sf() |>
  st_transform(crs = my_crs)

LTER5 = st_polygon(list(matrix(
  c(-149.8869755, -17.59182383,  # lower-left
    -149.8869755, -17.56818162,  # upper-left
    -149.8561009, -17.56818162,  # upper-right
    -149.8561009, -17.59182383,  # lower-right
    -149.8869755, -17.59182383),   # close polygon
  ncol = 2,
  byrow = TRUE))) |>
  st_sfc(crs = 4326) |>
  st_sf() |>
  st_transform(crs = my_crs)

LTER6 = st_polygon(list(matrix(
  c(-149.934537,  -17.52839766,  # lower-left
    -149.934537,  -17.50735955,  # upper-left
    -149.9115336, -17.50735955,  # upper-right
    -149.9115336, -17.52839766,  # lower-right
    -149.934537,  -17.52839766),   # close polygon
  ncol = 2,
  byrow = TRUE))) |>
  st_sfc(crs = 4326) |>
  st_sf() |>
  st_transform(crs = my_crs)

# rename columns & combine all LTERs into a single, multi-part polygon object
LTER1 = LTER1 %>%
  rename(geometry = st_sfc.st_polygon.list.matrix.c..149.8455917...17.48641792...149.8455917..) %>%
  mutate(Name = "LTER1")
LTER2 = LTER2 %>%
  rename(geometry = st_sfc.st_polygon.list.matrix.c..149.8116849...17.48131958...149.8116849..) %>%
  mutate(Name = "LTER2")
LTER3 = LTER3 %>%
  rename(geometry = st_sfc.st_polygon.list.matrix.c..149.7708619...17.52087158...149.7708619..) %>%
  mutate(Name = "LTER3")
LTER4 = LTER4 %>%
  rename(geometry = st_sfc.st_polygon.list.matrix.c..149.7772857...17.55064263...149.7772857..) %>%
  mutate(Name = "LTER4")
LTER5 = LTER5 %>%
  rename(geometry = st_sfc.st_polygon.list.matrix.c..149.8869755...17.59182383...149.8869755..) %>%
  mutate(Name = "LTER5")
LTER6 = LTER6 %>%
  rename(geometry = st_sfc.st_polygon.list.matrix.c..149.934537...17.52839766...149.934537..) %>%
  mutate(Name = "LTER6")
LTERs = rbind(LTER1, LTER2, LTER3, LTER4, LTER5, LTER6)

# check that everything looks good
ggplot() +
  geom_sf(data = LTERs, aes(fill = Name)) +
  geom_sf(data = moorea, fill = NA, color = "black") +
  theme_minimal() +
  labs(fill = "LTER Site")

# save as a shapefile
st_write(LTERs,
         here("Data", "LTER_Shapefiles", "LTER_Sites.shp"),
         append = FALSE)

#### CALCULATE TOPO-BATHY PREDICTORS ####
##### Bathymetry #####
# read in the Moorea bathymetric datasets - downloaded from SHOM (the Service hydrographique 
# et océanographique de la Marine) - the dataset is called: Lidar - French Polynesia 2015 Moorea 
# https://diffusion.shom.fr/multiproduct/product/configure/id/299

# there are 17 folders that store these data (in a grid system)
base_dir = here("Data", "Topobathy", "SHOM")

# find all .asc files within the sub-folders recursively
asc_files = list.files(
  path = base_dir,
  pattern = "\\.asc$",
  recursive = TRUE,
  full.names = TRUE)

# filter to only include files in MNT1 directories, where MNT1 means surface model
# with 1m resolution
asc_files = asc_files[grepl("/MNT1m/", asc_files)]

cat("Found", length(asc_files), "ASC files in MNT1 folders\n")

# verify the files are from MNT1 folders
print(head(asc_files))

# load all rasters into a list
cat("Loading rasters...\n")

raster_list = lapply(asc_files, function(f) {
  cat("Loading:", basename(f), "\n")
  rast(f)})

# mosaic all tiles together to produce a single topobathy surface
cat("Mosaicking rasters...\n")
bathymetry_mosaic = do.call(mosaic, c(raster_list, fun = "mean"))

# define/assign the CRS based on the SHOM metadata
cat("Defining CRS as WGS 84 / UTM zone 6S...\n")
crs(bathymetry_mosaic) = "EPSG:32706"  # WGS 84 / UTM zone 6S
plot(bathymetry_mosaic)

# this raster actually includes terrestrial coastal areas as well with elevation values! 
# let's save the topobathy data first
# writeRaster(bathymetry_mosaic, 
#             here("Data", "Topobathy", "Topobathymetry.tif"),
#             overwrite = TRUE)

# now, to keep only bathymetry, remove all cells with values > 0.5 m (elevation)
# 0.5 m leaves us a small buffer to account for tidal range
cat("Masking out land portions (elevation > 0.5 m)...\n")
bathymetry_mosaic[bathymetry_mosaic > 0.5] = NA

# in the end, it will be easier to interpret depth coefficients from the models if they
# are positive values, rather than having to mentally flip the negative sign.
bathymetry_mosaic = bathymetry_mosaic*-1

# save the bathymetry raster
writeRaster(bathymetry_mosaic, 
            here("Data", "Topobathy", "Bathymetry.tif"),
            overwrite = TRUE)

res(bathymetry_mosaic) # check resolution, which should be 1 m x 1 m
crs(bathymetry_mosaic) # check the projected coordinate system
minmax(bathymetry_mosaic) # check the range of depth values

# using the bathymetry raster, calculate seascape morphometrics using a 
# quadratic local fit and the Queen's case neighborhood (8 neighboring cells)

##### Slope #####
# the rate of maximum change in depth measured in degrees
slope = Qfit(r = bathymetry_mosaic, 
             w = c(3,3), 
             unit = "degrees",
             metrics = "qslope",
             na.rm = TRUE)

writeRaster(slope, 
            here("Data", "Topobathy", "Slope.tif"),
            overwrite = TRUE)

gc() # free up space if possible

##### Mean Curvature #####
# first derivative of slope and second derivative of depth, measured in m^-1
meanc = Qfit(r = bathymetry_mosaic, 
             w = c(3,3), 
             metrics = "meanc",
             na.rm = TRUE)

writeRaster(meanc, 
            here("Data", "Topobathy", "Curvature.tif"),
            overwrite = TRUE) 

gc() # free up space if possible

##### SAPA Rugosity #####
# surface area to planar area ratio rugosity (slope-corrected)
# surface roughness or complexity, measured as a ratio
sapa = SAPA(r = bathymetry_mosaic, 
            w = c(3,3),
            slope_correction = TRUE,
            na.rm = TRUE)

writeRaster(sapa, 
            here("Data", "Topobathy", "Rugosity.tif"),
            overwrite = TRUE)

gc() # free up space if possible

##### Adjusted Standard Deviation of Bathymetry #####
# standard deviation of bathymetry (a measure of rugosity), adjusted for slope
adjSD = AdjSD(r = bathymetry_mosaic,
              w = c(3, 3),
              na.rm = TRUE)

writeRaster(adjSD, 
            here("Data", "Topobathy", "Adjusted_SD_Bathymetry.tif"),
            overwrite = TRUE)

gc() # free up space if possible

##### Distance from land #####
# convert the Moorea spatial object to a SpatVector
moorea_sv = vect(moorea)

# rasterize the data using bathymetry_mosaic as a template raster
moorea_ras = rasterize(moorea_sv, bathymetry_mosaic, field = 1, background = NA)

# calculate distances from each cell to the nearest point on land
dist_to_moorea = distance(moorea_ras)

# mask to match the bathymetry raster (i.e., keep only cells with a depth value - remove NAs)
dist_to_moorea = mask(dist_to_moorea, bathymetry_mosaic)

# save
writeRaster(dist_to_moorea, 
            here("Data", "Topobathy", "Land_Distance.tif"),
            overwrite = TRUE)

##### Distance from the reef crest #####
# read in the coral reef geomorphic zone data from the Allen Coral Atlas - we will use this layer 
# to define the location of the reef crest. learn more about the geomorphic zones here: 
# https://storage.googleapis.com/coral-atlas-static-files/download-package-materials/Class-Descriptions-Geomorphic-Maps-v3.pdf
geomorphic = st_read(here("Data", "Habitat", "geomorphic.geojson"))
compareCRS(geomorphic, gcs) # should be in GCS WGS 84

# re-project to UTM Zone 6S 
geomorphic = st_transform(geomorphic, my_crs)
compareCRS(geomorphic, my_crs)
plot(geomorphic)

# extract the reef crest from the geomorphic object
reef_crest_sf = geomorphic[geomorphic$class == "Reef Crest", ]

# convert to SpatVector
reef_crest_sv = vect(reef_crest_sf)

# rasterize using bathymetry_mosaic as a template
reef_crest_ras = rasterize(reef_crest_sv, bathymetry_mosaic, field = 1, background = NA)

# calculate distances from each cell to the nearest point along the reef crest
dist_to_reef_crest = distance(reef_crest_ras)

# mask to match the bathymetry raster (i.e., keep only cells with a depth value - remove NAs)
dist_to_reef_crest = mask(dist_to_reef_crest, bathymetry_mosaic)

# save the raster
writeRaster(dist_to_reef_crest, 
            here("Data", "Topobathy", "Crest_Distance.tif"),
            overwrite = TRUE)

# compare the extents, row & column counts, resolutions, and projections of all rasters
# if everything matches, this should return TRUE
compareGeom(bathymetry_mosaic, slope, meanc, sapa, adjSD, dist_to_moorea, dist_to_reef_crest,
            ext = TRUE, res = TRUE, crs = TRUE, rowcol = TRUE)
