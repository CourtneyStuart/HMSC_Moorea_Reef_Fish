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
# starting with the backreef (lagoon) habitat records from: 
# MCR LTER: Coral Reef: Long-term Community Dynamics: Backreef (Lagoon) Corals Annual Survey, ongoing since 2005. 
# led by Peter Edmunds under dataset ID knb-lter-mcr.1038.
backreef_data = read.csv(here("Data", "Habitat", "knb-lter-mcr.1038",
                              "MCR_LTER_Coral_Cover_Backreef_Long_20250429.csv"))

# ensure percent cover data are in numeric format
backreef_data = backreef_data %>%
  mutate(percent_cover = as.numeric(percent_cover))

# prepare backreef data with standardized names
# exclude the 2005 data because the fish data we will use were only collected following a 
# standardized transect protocol starting in 2006.
backreef_prepared = backreef_data %>%
  mutate(
    Date = paste0(year, "-04"),  # convert year to YYYY-MM format to match later fringing/forereef data 
    Transect = paste("LTER", gsub("LTER_", "", site), "Backreef Coral Transect")) %>%
  filter(!grepl("^2005", Date)) # exclude 2005 observations

# note from Edmunds: The backreef habitat in the lagoons is relatively shallow. There is a 
# minimum depth needed to position the camera frame over corals and still be under water.
# When a quadrat is too shallow, the photo may be tilted off horizontal so that part of the 
# photo is open water, called "blue water", instead of benthic cover. This makes the
# photoquadrat inappropriate for percent cover analysis. For these photoquadrats the code
# "BW" is entered in place of coral cover percent.

# thus, we need to also remove any BW observations.
backreef_prepared = backreef_prepared %>%
  filter(benthic_category!= "BW")

# calculate the mean and SD for observed groups at each date-transect combination.
# no need to specify habitat, because this dataset only includes records from the backreef.
backreef_stats = backreef_prepared %>%
  group_by(Date, Transect, benthic_category) %>%
  summarise(
    Mean_Cover = mean(percent_cover),
    SD_Cover = sd(percent_cover),
    .groups = 'drop')

# get all unique combinations of date, transect, and benthic categories
backreef_dates_transects = backreef_prepared %>%
  distinct(Date, Transect)

backreef_groups = backreef_prepared %>%
  distinct(benthic_category)

# create a complete grid
backreef_grid = backreef_dates_transects %>%
  cross_join(backreef_groups)

# join with summary stats and fill missing with 0s
final_backreef_summary = backreef_grid %>%
  left_join(backreef_stats, 
            by = c("Date", "Transect", "benthic_category")) %>%
  rename(Benthic_Category = benthic_category) %>%
  mutate(
    Mean_Cover = replace_na(Mean_Cover, 0),
    SD_Cover = replace_na(SD_Cover, 0))

# convert to wide format and add a column to clarify that these are all backreef observations
final_backreef_summary_wide = final_backreef_summary %>%
  pivot_wider(
    names_from = Benthic_Category,
    values_from = c(Mean_Cover, SD_Cover),
    names_glue = "{Benthic_Category}_{.value}") %>%
  rename_with(~gsub("Mean_Cover", "Mean_Cover", .x)) %>%
  rename_with(~gsub("SD_cover", "SD_Cover", .x)) %>%
  mutate(Habitat = as.character("Backreef")) %>%
  relocate(Habitat, .after = Date)

# some quick edits to standardize the column names
require(stringr)
final_backreef_summary_wide = final_backreef_summary_wide %>%
  rename_with(~ .x %>%
                # replace spaces and other separators with underscores
                str_replace_all("[ .-]+", "_") %>%
                # split by underscore, capitalize first letter of each word, then rejoin
                str_split("_") %>%
                map_chr(~ .x %>%
                          # preserve "SD" and "CTB" as-is, capitalize others
                          map_chr(function(word) {
                            if (toupper(word) == "SD") {
                              "SD"
                            } else if (toupper(word) == "CTB") {
                              "CTB"
                            } else {
                              str_to_title(word)
                            }
                          }) %>%
                          paste(collapse = "_")))

# now, read in the MCR LTER habitat data from the fringing and forereef zones: 
# MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Corals, ongoing since 2005.
# this is also led by Peter Edmunds, under dataset ID knb-lter-mcr.4.42.
fringefore_data = read.csv(here("Data", "Habitat", "knb-lter-mcr.4.42", "knb-lter-mcr.4_1_20250409.csv"))

# exclude no data records
fringefore_data = fringefore_data %>%
  filter(Percent_Cover != "ND")

# now ensure that percent cover data are in numeric format
fringefore_data = fringefore_data %>%
  mutate(Percent_Cover = as.numeric(Percent_Cover))

# define a new hard coral benthic group
hard_corals = c("Acanthastrea", "Acropora", "Astrea", "Astreopora", "Cyphastrea", 
                 "Danafungia", "Fungiidae unidentified", "Gardineroseris", "Goniastrea", 
                 "Herpolitha", "Leptastrea", "Leptoseris", "Lithophyllon", "Lobactis", 
                 "Lobophyllia", "Montipora", "Pachyseris", "Pavona", "Plesiastrea", 
                 "Pleuractis", "Pocillopora", "Porites", "Porites irregularis", 
                 "Porites rus", "Porites spp. Massive", "Psammocora", "Sandolitha", 
                 "Stylocoeniella", "Tubastrea")

# exclude the 2005 data because the fish data we will use were only collected following a 
# standardized transect protocol starting in 2006.
# exclude the 17 m depth observations because the fish data were recorded at a maximum 
# depth of 10 m along forereef transects. 
# extract transect name from location information by excluding the pole pair and
# quadrat information.
fringefore_data_filtered = fringefore_data %>%
  filter(!grepl("^2005", Date)) %>%
  filter(!Depth == 17) %>%
  mutate(Transect = gsub(" Pole.*$", "", Location))

# create functional group categories and sum within each quadrat
fringefore_data_grouped = fringefore_data_filtered %>%
  mutate(
    Functional_Group = case_when(
      Taxonomy_Substrate_or_Functional_Group %in% hard_corals ~ "Coral",
      TRUE ~ Taxonomy_Substrate_or_Functional_Group)) %>%
  group_by(Date, Habitat, Transect, Section_of_Transect,
           Quadrat_within_section, Functional_Group) %>%
  summarise(Percent_Cover = sum(Percent_Cover), .groups = 'drop')

# calculate the mean and standard deviation (SD) for observed groups at each unique 
# date-habitat-transect combination.
fringefore_stats = fringefore_data_grouped %>%
  group_by(Date, Habitat, Transect, Functional_Group) %>%
  summarise(
    Mean_Cover = mean(Percent_Cover),
    SD_Cover = sd(Percent_Cover),
    .groups = 'drop')

# get all unique combinations of date, transect, and taxonomy groups
fringefore_dates_transects = fringefore_data_grouped %>%
  distinct(Date, Habitat, Transect)

fringefore_groups = fringefore_data_grouped %>%
  distinct(Functional_Group)

# create a complete grid
fringefore_grid = fringefore_dates_transects %>%
  cross_join(fringefore_groups)

# join with the summary stats dataframe and fill missing functional groups with 0s
final_fringefore_summary = fringefore_grid %>%
  left_join(fringefore_stats, 
            by = c("Date", "Habitat", "Transect", "Functional_Group")) %>%
  mutate(
    Mean_Cover = replace_na(Mean_Cover, 0),
    SD_Cover = replace_na(SD_Cover, 0))

# convert to wide format
final_fringefore_summary_wide = final_fringefore_summary %>%
  pivot_wider(
    names_from = Functional_Group,
    values_from = c(Mean_Cover, SD_Cover),
    names_glue = "{Functional_Group}_{.value}") 

# some quick edits to standardize the column names
final_fringefore_summary_wide = final_fringefore_summary_wide %>%
  rename_with(~ .x %>%
                # replace spaces and other separators with underscores
                str_replace_all("[ .-]+", "_") %>%
                # split by underscore, capitalize first letter of each word, then rejoin
                str_split("_") %>%
                map_chr(~ .x %>%
                          # preserve "SD" and "CTB" as-is, capitalize others
                          map_chr(function(word) {
                            if (toupper(word) == "SD") {
                              "SD"
                            } else if (toupper(word) == "CTB") {
                              "CTB"
                            } else {
                              str_to_title(word)
                            }
                          }) %>%
                          paste(collapse = "_")))

# combine the backreef, fringing reef, and forereef data
habitat_data = rbind(final_backreef_summary_wide,
                     (final_fringefore_summary_wide %>% 
                        select(-Non_Coralline_Crustose_Algae_Mean_Cover, 
                               -Non_Coralline_Crustose_Algae_SD_Cover,
                               -Soft_Coral_Mean_Cover,
                               -Soft_Coral_SD_Cover,
                               -Unknown_Or_Other_Mean_Cover,
                               -Unknown_Or_Other_SD_Cover)))
# check for NAs
anyNA(habitat_data)

# some quick tidying (we will use the same naming scheme across habitat and fish data)
habitat_data = habitat_data %>%
  mutate(
    Date = substr(Date, 1, 4),
    Transect = str_extract(Transect, "LTER \\d+") %>%
      str_replace(" ", "_")) %>%
  rename(Year = Date, Site = Transect) %>%
  mutate(Site = str_replace(Site, "LTER_", "LTER"),
         Name = paste(Site, Habitat, sep = "_")) %>% 
  relocate(Name, .after = Year) %>%
  relocate(Site, .before = Habitat)

# save the prepped habitat data
write.csv(habitat_data,
          here("Data", "Habitat", "Prepped_MCR_LTER_Habitat_Data.csv"),
          row.names = FALSE)
