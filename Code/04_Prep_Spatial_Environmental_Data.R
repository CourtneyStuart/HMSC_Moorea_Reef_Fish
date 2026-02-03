#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### LIBRARIES ####
# install packages (first run only)
# install.packages(c("easypackages", "raster", "terra", "MultiscaleDTM", "sp",
#                  "sf", "conflicted", "spatialEco", "dplyr", "here", "fasterize",
#                  "nngeo", "ggplot2", "tidyr", "purrr", "stringr"))

# load packages
library(easypackages)
libraries("raster", "terra", "MultiscaleDTM", "sp", "sf", "conflicted", 
          "spatialEco", "dplyr", "here", "fasterize", "nngeo", "ggplot2", 
          "tidyr", "purrr", "stringr")

# resolve package conflicts
conflict_prefer("terrain", "terra")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# additional settings
rasterOptions(progress = 'text') # progress info for processing large rasters
options(terra.progress = 1) # progress info for processing large rasters

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

#### DATA PROCESSING ####
##### SEAFLOOR MORPHOLOGY ####
# we have raster layers of depth, slope, rugosity, etc. across the entire nearshore zone
# of Moorea. we also have fish species presence/absence and abundance data at 
# 6 LTER sites X 3 habitat types/site X 19 years = 342 unique sampling events.
# each habitat type is surveyed on 4 replicate transects (4 reps x 6 sites X 3 habs = 
# 72 surveys annually). for each site-habitat type combo, we know the approx. centroid
# coordinate of the 4 replicate surveys. we will use the approx. coordinates to extract 
# the environmental conditions at each survey location and relate these to the fish data.
# because we know that fish have different home ranges, migration patterns, etc., we will 
# calculate seafloor conditions within 100- and 500-m buffers around these centroids. 

# first read in the rasters prepped in "01_Prep_Topobathymetric_Data.R"
depth = rast(here("Data", "Topobathy", "Bathymetry.tif"))
adjSD_depth = rast(here("Data", "Topobathy", "Adjusted_SD_Bathymetry.tif"))
slope = rast(here("Data", "Topobathy", "Slope.tif"))
curvature = rast(here("Data", "Topobathy", "Curvature.tif"))
rugosity = rast(here("Data", "Topobathy", "Rugosity.tif"))
crest_dist = rast(here("Data", "Topobathy", "Crest_Distance.tif"))
land_dist = rast(here("Data", "Topobathy", "Land_Distance.tif"))

# check the CRS, they should all match
compareGeom(depth, adjSD_depth, slope, curvature, rugosity,
            crest_dist, land_dist, crs = TRUE)
crs(depth)

# read in the transect centroids from MCR LTER
centroids_gcs = st_read(here("Data", "Fish", "DO_NOT_SHARE", "Permanent_Fish_Monitoring_Locations.kml"))

# replace "Fringe" with "Fringing" in the Name column
centroids_gcs = centroids_gcs %>%
  mutate(Name = str_replace(Name, "Fringe", "Fringing"))

# project the data to our desired CRS
centroids = st_transform(centroids_gcs, crs(depth))
st_crs(centroids) == st_crs(depth)

# remove the extra space between "LTER" and the number in Name column
require(stringr)
centroids = centroids %>%
  mutate(Name = str_replace(Name, "LTER ", "LTER"))

# create buffers around the centroids
buff100 = st_buffer(centroids, dist = 100, singleSide = FALSE)
buff500 = st_buffer(centroids, dist = 500, singleSide = FALSE)

# plot the survey locations around the island
moorea = st_read(here("Data", "FP_Islands_Shapefile", "Moorea.shp"))

plot(st_geometry(moorea["Island"]), col = "grey", border = "grey")
plot(st_geometry(buff500), col = "cyan", add = TRUE)
plot(st_geometry(buff100), col = "darkcyan", add = TRUE)

# extract the conditions within the 100 m buffers
rasters = c("depth", "adjSD_depth", "slope", "curvature", "rugosity")
conditions_100m = data.frame()

# loop through rasters
for (i in seq_along(rasters)) {
  raster_name = rasters[i]
  raster_layer = get(raster_name) # load raster by name
  
  # extract values from raster
  extracted = terra::extract(raster_layer, buff100, fun = NULL, cells = FALSE)
  
  # map the name to extracted data
  extracted$Name = buff100$Name[extracted$ID]
  
  # rename raster values column
  colnames(extracted)[2] = "value" # second column contains raster values
  
  # calculate statistics for each raster
  stats = extracted %>%
    group_by(Name) %>%
    summarise(
      mean_100m = mean(value, na.rm = TRUE),
      sd_100m = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(raster = raster_name)
  
  # append to results
  conditions_100m = bind_rows(conditions_100m, stats)
  
  # print progress
  print(paste("Completed raster", i, "of", length(rasters)))
}

# reshape data to wide format, only keeping the required stats
require(stringr)
conditions_100m_wide = conditions_100m %>%
  pivot_wider(
    id_cols = Name,
    names_from = raster,
    values_from = c(mean_100m, sd_100m),
    names_glue = "{raster}_{.value}") %>%
  rename_with(~ str_replace_all(., "_value$", ""), everything()) %>% # remove '_value'
  rename_with(~ str_to_title(str_replace_all(., "_", " ")), everything()) %>% # capitalize each word (with spaces)
  rename_with(~ str_replace_all(., " ", "_"), everything()) %>% # replace spaces with underscores
  rename_with(~ str_replace_all(., "Sd", "SD"), everything()) # replace "Sd" with "SD"

# extract the conditions within the 500 m buffers
conditions_500m = data.frame()

# loop through rasters
for (i in seq_along(rasters)) {
  raster_name = rasters[i]
  raster_layer = get(raster_name) # load raster by name
  
  # extract values from raster
  extracted = terra::extract(raster_layer, buff500, fun = NULL, cells = FALSE)
  
  # map the name to extracted data
  extracted$Name = buff500$Name[extracted$ID]
  
  # rename raster values column
  colnames(extracted)[2] = "value" # second column contains raster values
  
  # calculate statistics for each raster
  stats = extracted %>%
    group_by(Name) %>%
    summarise(
      mean_500m = mean(value, na.rm = TRUE),
      sd_500m = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(raster = raster_name)
  
  # append to results
  conditions_500m = bind_rows(conditions_500m, stats)
  
  # print progress
  print(paste("Completed raster", i, "of", length(rasters)))
}

# reshape data to wide format, only keeping the required stats
conditions_500m_wide = conditions_500m %>%
  pivot_wider(
    id_cols = Name,
    names_from = raster,
    values_from = c(mean_500m, sd_500m),
    names_glue = "{raster}_{.value}") %>%
  rename_with(~ str_replace_all(., "_value$", ""), everything()) %>% # remove '_value'
  rename_with(~ str_to_title(str_replace_all(., "_", " ")), everything()) %>% # capitalize each word (with spaces)
  rename_with(~ str_replace_all(., " ", "_"), everything()) %>% # replace spaces with underscores
  rename_with(~ str_replace_all(., "Sd", "SD"), everything()) # replace "Sd" with "SD"

# save the x and y coordinates of the centroids along with the environmental data
# extract coordinates from centroids
coords = st_coordinates(centroids)

# add coordinates to the centroids data as regular columns
centroids_with_coords = centroids %>%
  mutate(
    Long_UTM6S = coords[, "X"],
    Lat_UTM6S = coords[, "Y"]) %>%
  select(-Description) %>%
  st_drop_geometry()  # remove the geometry column, keeping just the data

# join with the conditions_100m_wide and conditions_500m_wide dataframes
conditions_100m_wide = conditions_100m_wide %>%
  left_join(centroids_with_coords %>% select(Name, Long_UTM6S, Lat_UTM6S), 
            by = "Name") %>%
  relocate(Long_UTM6S, Lat_UTM6S, .after = Name)

conditions_500m_wide = conditions_500m_wide %>%
  left_join(centroids_with_coords %>% select(Name, Long_UTM6S, Lat_UTM6S), 
            by = "Name") %>%
  relocate(Long_UTM6S, Lat_UTM6S, .after = Name)

# join the 100m and 500m conditions dataframes together
environmental_conditions = conditions_100m_wide %>%
  left_join(conditions_500m_wide, by = c("Name", "Long_UTM6S", "Lat_UTM6S"))

# calculate the distance from each fish centroid to the nearest point along the reef 
# crest and along the coastline

# extract raster values at centroid locations
land_dist_values = terra::extract(land_dist, centroids)

# combine into a dataframe
land_dist_df = data.frame(
  Land_Dist = land_dist_values[, 2],  # second column contains the extracted values
  Long_UTM6S = coords[, "X"],
  Lat_UTM6S = coords[, "Y"])

# extract raster values at centroid locations
crest_dist_values = terra::extract(crest_dist, centroids)

# combine into a dataframe
crest_dist_df = data.frame(
  Crest_Dist = crest_dist_values[, 2],  # second column contains the extracted values
  Long_UTM6S = st_coordinates(centroids)[, "X"],
  Lat_UTM6S = st_coordinates(centroids)[, "Y"])

# add these to the environmental conditions dataframe
environmental_conditions = environmental_conditions %>%
  left_join(land_dist_df, by = c("Long_UTM6S", "Lat_UTM6S")) %>%
  left_join(crest_dist_df, by = c("Long_UTM6S", "Lat_UTM6S"))

# remove intermediate data
rm(list = c("land_dist_df", "land_dist_values", "crest_dist_df", "crest_dist_values"))

# we will assume that these conditions are static over time, but we still need to
# include the Year column, though, for subsequent HMSC analyses
year_df = data.frame(Year = 2006:2024)

environmental_expanded = environmental_conditions %>%
  crossing(year_df)

##### HABITAT #####
# read in the habitat data that was previously prepared in script "02_Prep_LTER_Habitat_Data.R"
habitat_data = read.csv(here("Data", "Habitat", "Prepped_MCR_LTER_Habitat_Data.csv"))

# join the habitat data to the environmental_expanded dataframe
environmental_expanded = environmental_expanded %>%
  left_join(habitat_data) %>%
  relocate(Year, .before = Name) %>%
  select(-Site, -Habitat)

##### COTS OUTBREAKS #####
# crown of thorns starfish outbreak data from: 
# MCR LTER: Coral Reef: Long-term Population Dynamics of Acanthaster planci, ongoing since 2005
# led by Andrew Brooks under dataset ID: knb-lter-mcr.1039.
cots = read.csv(here("Data", "COTS", "MCR_LTER_COTS_abundance_2005-2024_20241205.csv"))

# this stores the abundance of Acanthaster planci per transect - using the four
# pseudoreplicate transects conducted within each of the 3 habitat types of each 
# of the 6 LTER sites = 72 transects total.
# we only want data from 2006 through 2024
cots = cots %>%
  filter(Year %in% 2006:2024)

# calculate the total abundance of A. planci seen in each habitat type of each LTER site
cots = cots %>%
  group_by(Year, Site, Habitat) %>%
  summarise(COTS = sum(COTS)) %>%
  select(Year, Site, Habitat, COTS)

plot(cots$Year, cots$COTS)

# add a new Name column to store the sample ID
cots = cots %>%
  mutate(Name = paste0(gsub(" ", "", Site), "_", Habitat)) %>%
  relocate(Name, .after = Year) %>%
  ungroup() %>%
  select(-Site, -Habitat)

##### DEGREE HEATING WEEKS #####
# sea surface temperature, degree heating weeks, and coral bleaching watch data from:
# NOAA Coral Reef Watch (CRW) daily 5km Regional Virtual Stations > Society Archipelago
# https://coralreefwatch.noaa.gov/product/vs/description.php#ascii and 
# https://coralreefwatch.noaa.gov/product/vs/map.php

# NOAA CRW data description: We use the data that fall within the boundary of each Regional
# Virtual Station to create a time series. Rather than provide the value of every 5km satellite
# pixel data point within a Station's boundary, we developed a new Regional Virtual Station 
# algorithm. The algorithm is based on the daily 90th percentile Coral Bleaching HotSpot value
# among a Station's 5km pixels, and the other variables at the pixel where the 90th percentile 
# HotSpot value locates. Daily Regional Virtual Station sea surface temperature (SST),
# SST Anomaly, and Coral Bleaching HotSpot values for a given day are the respective values from 
# the satellite pixel where the 90th percentile Coral Bleaching HotSpot value occurs on that day.
# The Daily Regional Virtual Station Degree Heating Week (DHW) is then calculated (accumulated) 
# from the Daily Regional Virtual Station HotSpots over a consecutive 84 days, following the DHW 
# algorithm. The Daily Regional Virtual Station Bleaching Alert Area single-day value is derived 
# from the Daily Regional Virtual Station Coral Bleaching HotSpot and DHW pair, and a rolling
# Bleaching Alert Area (7-day maximum) composite value is then produced.

crw = read.csv(here("Data", "Heat_Waves", "NOAA_CRW_Society_Archipelago.csv"))

# for each year, calculate the maximum degree heating week (DHW) 
crw = crw %>%
  rename(Year = YYYY) %>%
  group_by(Year) %>%
  summarise(Max_DHW = max(DHW_from_90th_HS.1))

##### CYCLONE OLI #####
# Cyclone Oli hit in 2010 and, for fishes, MCR LTER research staff suspect that there
# were lag effects up to 2 years later
cyclones = data.frame(
  Year = 2006:2024,
  Cyclone = c(
    rep(0, 4),  # 2006-2009: non-cyclone years
    3,          # 2010: cyclone year
    2,          # 2011: 1-year post cyclone
    1,          # 2012: 2-years post cyclone
    rep(0, 12)  # 2013-2024: non-cyclone years & outside of lag time
  ))

##### COMBINE COTS, DHW, CYCLONE OLI #####
# merge all dataframes together
contextual_data = cots %>%
  left_join(crw) %>%
  left_join(cyclones)

# combine the environmental and contextual data
environmental_data = left_join(environmental_expanded, contextual_data) %>%
  relocate(Year, .before = Name) %>%
  relocate(COTS, .after = Lat_UTM6S) %>%
  relocate(Max_DHW, .after = COTS) %>%
  relocate(Cyclone, .after = Max_DHW) %>%
  relocate(Crest_Dist, .after = Cyclone) %>%
  relocate(Land_Dist, .after = Crest_Dist)

# check number of observations per year (should be 18, one per site-habitat combination!)
environmental_data %>%
  group_by(Year) %>%
  summarise(n_obs = n()) %>%
  arrange(Year)

# check number of observations per site-habitat combination (should all be 19, one per year!)
environmental_data %>%
  group_by(Name) %>%
  summarise(n_obs = n()) %>%
  arrange(Name)

# save the dataset
write.csv(environmental_data,
          here("HMSC", "Data", "Intermediate_Datasets", 
               "Full_Spatial_Environmental_Dataset.csv"),
          row.names = FALSE)
