#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### LIBRARIES ####
# install packages (first run only)
# install.packages(c("here", "easypackages", "raster", "terra", "sp", "sf",
#                    "conflicted", "spatialEco", "dplyr", "ggplot2", "tidyr",
#                    "purrr", "stringr", "corrplot", "Cairo", "usdm",
#                    "PNWColors", "ape", "Hmsc"))

# load packages
library(easypackages)
libraries("here", "raster", "terra", "sp", "sf", "conflicted", "spatialEco",
          "dplyr", "ggplot2", "tidyr", "purrr", "stringr", "corrplot", "Cairo", 
          "usdm", "PNWColors", "ape", "Hmsc")

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

#### COORDINATE SYSTEMS ####
# save PROJ.4 string for WGS 84 / UTM Zone 6S (EPSG WKID 32706) with unit meters
my_crs = CRS("+proj=utm +zone=6 +south +datum=WGS84 +units=m +no_defs +type=crs")

# save PROJ.4 string for the standard geographic coordinate system used by
# Garmin GPS - WGS84 - World Geodetic System 1984 (EPSG WKID 4326) with unit decimal
# degrees 
gcs = CRS("+proj=longlat +datum=WGS84 +no_defs +type=crs")

#### DATA PROCESSING ####
# read in our datasets
fish_PA = read.csv(here("HMSC", "Data", "Intermediate_Datasets", 
                        "Fish_Presence_Absence_Dataset.csv"))

fish_ABU = read.csv(here("HMSC", "Data", "Intermediate_Datasets", 
                        "Fish_Abundance_Dataset.csv"))

spatenv = read.csv(here("HMSC", "Data", "Intermediate_Datasets", 
                        "Full_Spatial_Environmental_Dataset.csv"))

traits = read.csv(here("HMSC", "Data", "Intermediate_Datasets",
                       "Fish_Traits_Dataset.csv")) %>%
  mutate(Species = gsub(" ", ".", Species)) # format Genus.species

tree = read.tree(here("HMSC", "Data", "Intermediate_Datasets", 
                      "Fish_Taxonomic_Tree.tre"))

# check how species names are stored in our taxonomic tree
tree$tip.label

# use . instead of _ as separators to match species names in other datasets
tree$tip.label = gsub("_", ".", tree$tip.label) 

# create master datasets that have ALL data
master_PA = spatenv %>%
  left_join(fish_PA)

master_ABU = spatenv %>%
  left_join(fish_ABU)

# check for collinearity among the numeric spatial-environmental predictors
full_preds = master_PA %>%
  select(5:6, 8:39) %>%
  mutate(COTS = as.integer(COTS)) %>%
  mutate(Max_DHW = as.integer(Max_DHW))

# create a correlation matrix
cormat = cor(full_preds, 
             use = "complete.obs")

# save a color palette for plotting
palette = pnw_palette("Shuksan2", 200, type = "continuous")

# set plotting margins
par(mar = c(0,0,0,0))

# create and save the full correlation plot
corrplot(cormat, method = "color", col = palette, type = "upper",
         order = "original", addCoef.col = "black", number.cex = 0.50, 
         number.digits = 2, tl.col = "black", tl.srt = 40, tl.cex = 0.70)
Cairo(file = here("Figures", "Full_Correlation_Matrix.png"), 
      bg = "white", type = "png", units = "in", width = 10, height = 10, 
      pointsize = 12,  dpi = 600)
par(mar = c(0,0,0,0))
corrplot(cormat, method = "color", col = palette, type = "upper",
         order = "original", addCoef.col = "black", number.cex = 0.50, 
         number.digits = 2, tl.col = "black", tl.srt = 40, tl.cex = 0.70)
dev.off()

# the "usdm" package's vifcor function runs through all pairs of variables and checks
# whether any exceed the r threshold we define
print(vifcor(full_preds, th = 0.7))

# the vifstep function runs through all pairs of variables and checks whether any exceed 
# the VIF threshold we define
print(vifstep(full_preds, th = 10))

# the AdjSD of depth variables show a strong positive correlation (unsurprisingly) with 
# the depth, slope, and curvature variables at both scales. because many fish species have
# specific depth preferences, and because depth influences temperature/salinity/DO/light/etc., 
# I want to keep depth itself. 

# I also want to capture seafloor structure, and it seems the best way to do this is using
# curvature. AdjSD of depth and slope have collineariy issues, and SAPA rugosity has a very
# small range of values (at both scales), so I expect it to have limited explanatory power.
# for these reasons, keep curvature rather than AdjSD of depth, slope, and rugosity. 

# While it was flagged for collinearity issues in the full predictor suite, I want to keep
# the land_dist due to its influence on exposure to anthropogenic pollutants, soils, 
# tourism, etc. try keeping land_dist instead of crest_dist (which is itself in some ways
# redundant with our other predictors).

selected_preds = full_preds %>%
  select(COTS, Max_DHW, Land_Dist,
         Depth_Mean_100m, Depth_Mean_500m,
         Curvature_Mean_100m, Curvature_Mean_500m,
         Coral_Mean_Cover, Macroalgae_Mean_Cover, CTB_Mean_Cover)

cormat2 = cor(selected_preds, 
             use = "complete.obs")

# create and save the restricted correlation plot
corrplot(cormat2, method = "color", col = palette, type = "upper",
         order = "original", addCoef.col = "black", number.cex = 0.8, 
         number.digits = 2, tl.col = "black", tl.srt = 60, tl.cex = 1)
Cairo(file = here("Figures", "Restricted_Correlation_Matrix.png"),
      bg = "white", type = "png", units = "in", width = 10, height = 10,
      pointsize = 12,  dpi = 600)
par(mar = c(0,0,0,0))
corrplot(cormat2, method = "color", col = palette, type = "upper",
         order = "original", addCoef.col = "black", number.cex = 0.8, 
         number.digits = 2, tl.col = "black", tl.srt = 60, tl.cex = 1)
dev.off()

# checking correlation and VIF scores
print(vifcor(selected_preds, th = 0.7))
print(vifstep(selected_preds, th = 10))
print(vifstep(selected_preds, th = 5)) # even stricter VIF assessment

# while curvature (500 m) and land_dist are correlated (-0.74), I will keep both 
# considering they're conceptually distinct and that all variables have VIF < 5. 

# creating the desired dataframes
# presence-absence fish data
data_PA = master_PA %>%
  select(Year, Name, Long_UTM6S, Lat_UTM6S,
         COTS, Max_DHW, Cyclone, Land_Dist, 
         Depth_Mean_100m, Depth_Mean_500m,
         Curvature_Mean_100m, Curvature_Mean_500m,
         Coral_Mean_Cover, Macroalgae_Mean_Cover, CTB_Mean_Cover,
         Abudefduf.septemfasciatus:Zebrasoma.velifer) # select all species columns
anyNA(data_PA) # there should be no NAs

# fish abundance data
data_ABU = master_ABU %>%
  select(Year, Name, Long_UTM6S, Lat_UTM6S,
         COTS, Max_DHW, Cyclone, Land_Dist, 
         Depth_Mean_100m, Depth_Mean_500m,
         Curvature_Mean_100m, Curvature_Mean_500m,
         Coral_Mean_Cover, Macroalgae_Mean_Cover, CTB_Mean_Cover,
         Abudefduf.septemfasciatus:Zebrasoma.velifer) # select all species columns
anyNA(data_ABU) # there should be no NAs

# separate out the site and habitat information from the current Name column so that 
# we can estimate the effects of both independently
# also ensure that the columns are formatted properly at this stage
data_PA = data_PA %>%
  separate(Name, into = c("Site", "Habitat"), sep = "_", remove = TRUE) %>%
  relocate(Site, .after = Year) %>%
  relocate(Habitat, .after = Site) %>%
  mutate(
    Year = as.factor(Year),
    Site = as.factor(Site),
    Habitat = as.factor(Habitat),
    Cyclone = as.factor(Cyclone))

data_ABU = data_ABU %>%
  separate(Name, into = c("Site", "Habitat"), sep = "_", remove = TRUE) %>%
  relocate(Site, .after = Year) %>%
  relocate(Habitat, .after = Site) %>%
  mutate(
    Year = as.factor(Year),
    Site = as.factor(Site),
    Habitat = as.factor(Habitat),
    Cyclone = as.factor(Cyclone))

#### CREATING REQUIRED HMSC INPUTS ####
##### Y species data #####
species_cols = 17:159 # species column numbers from data_PA

# presence-absence matrix
Y_PA = as.matrix(data_PA[, species_cols])

# abundance matrix
Y_ABU = as.matrix(data_ABU[, species_cols])

# dimension check
dim(Y_PA) 
dim(Y_ABU)

##### X spatial-environmental data #####
X = model.matrix(~ Habitat + COTS + Max_DHW + Cyclone + 
                   Land_Dist + Depth_Mean_100m + Depth_Mean_500m +
                   Curvature_Mean_100m + Curvature_Mean_500m +
                   Coral_Mean_Cover + Macroalgae_Mean_Cover + CTB_Mean_Cover,
                 data = data_PA)

# scale numeric and integer predictors (NOT factors)
continuous_vars = c("COTS", "Max_DHW", "Land_Dist", 
                    "Depth_Mean_100m", "Depth_Mean_500m",
                    "Curvature_Mean_100m", "Curvature_Mean_500m",
                    "Coral_Mean_Cover", "Macroalgae_Mean_Cover", "CTB_Mean_Cover")  
X_scaled = X
X_scaled[, continuous_vars] = scale(X[, continuous_vars])

##### Study design #####
study_design = data.frame(
  site = factor(data_PA$Site),
  year = factor(data_PA$Year))

##### Spatial coordinates for sites #####
# get unique coordinates for each site
unique_sites = data_PA %>%
  group_by(Site) %>%
  slice(1) %>%
  select(Site, Long_UTM6S, Lat_UTM6S)

xy_sites = as.matrix(unique_sites[, c("Long_UTM6S", "Lat_UTM6S")])
rownames(xy_sites) = unique_sites$Site

##### Random levels #####
# site random effect with spatial structure
rL.site = HmscRandomLevel(sData = xy_sites)

# Year random effect  
rL.year = HmscRandomLevel(units = study_design$year)

##### CHECKS #####
# check alignment
nrow(Y_PA) == nrow(X_scaled)  # should be TRUE
nrow(Y_PA) == nrow(study_design)  # should be TRUE
nrow(Y_ABU) == nrow(X_scaled)  # should be TRUE
nrow(Y_ABU) == nrow(study_design)  # should be TRUE

# check for NAs
sum(is.na(Y_PA))  # preferably no NAs
sum(is.na(Y_ABU))  # preferably no NAs
sum(is.na(X_scaled))  # should be 0 (no NAs in predictors!)
sum(is.na(study_design))  # should be 0 (no NAs in study design!)

# check spatial data
nrow(xy_sites)  # should be 6
rownames(xy_sites)  # should match site names

# check species names match between Y and trait data
colnames(Y_PA) == colnames(Y_ABU) # should be the same (all TRUE)
species_Y = colnames(Y_PA) # should be the same (all TRUE)
species_trait = traits$Species

# check for matches
length(species_Y)  # there should be 143 species
length(species_trait)  # how many species have trait data?

# species that match across the trait and response data
matching_species = dplyr::intersect(species_Y, species_trait)
length(matching_species)

# HMSC requires:
# trait data as a dataframe
# rownames must exactly match column names in Y
# only include species that are in Y

# to be sure, filter trait data to only species in Y
trait_hmsc = traits[traits$Species %in% species_Y, ]

# set species names as rownames
rownames(trait_hmsc) = trait_hmsc$Species

# remove the Species column (now it's in rownames)
trait_hmsc$Species = NULL

# reorder to match Y column order
trait_hmsc = trait_hmsc[colnames(Y_PA), ]

# check alignment
all(rownames(trait_hmsc) == colnames(Y_PA))  # must be TRUE!
all(rownames(trait_hmsc) == colnames(Y_ABU))  # must be TRUE!

# check trait data
head(trait_hmsc)
str(trait_hmsc)

# scale continuous traits
continuous_traits = c("Trophic_Level", "Max_TL_cm")
trait_hmsc[, continuous_traits] = scale(trait_hmsc[, continuous_traits])

# make sure categorical traits are factors
categorical_traits = c("Body_Shape", "Column_Position", "Schooling",
                       "Shoaling", "Solitary", "Reproductive_Mode",
                       "Spawn_Agg")
trait_hmsc[, categorical_traits] = lapply(trait_hmsc[, categorical_traits], factor)

# whoops, all species in our data are reef-associated!
levels(trait_hmsc$Column_Position)

# so remove the Column_Position trait!
trait_hmsc = trait_hmsc %>%
  select(-Column_Position)

# define the models just to ensure that we've formatted the data properly -
# we will not run any MCMC sampling yet!!!! if everything is formatted properly,
# these should run with no errors!!!!

# define the presence-absence model
PA_model = Hmsc(Y = Y_PA,
          XData = as.data.frame(X_scaled),
          # specify which predictors to use
          XFormula = ~ HabitatForereef + HabitatFringing + COTS +
            Max_DHW + Cyclone1 + Cyclone2 + Cyclone3 + Land_Dist + 
            Depth_Mean_100m + Depth_Mean_500m +
            Curvature_Mean_100m + Curvature_Mean_500m +
            Coral_Mean_Cover + Macroalgae_Mean_Cover + CTB_Mean_Cover,
          TrData = trait_hmsc,  # add trait data
          # specify which traits to use
          TrFormula = ~ Body_Shape + Max_TL_cm + Trophic_Level +
            Reproductive_Mode + Spawn_Agg, 
          phyloTree = tree, # add taxonomic tree
          studyDesign = study_design, # add study design
          ranLevels = list(site = rL.site, 
                           year = rL.year),
          distr = "probit")

# define the abundance model
ABU_model = Hmsc(Y = Y_ABU,
                 XData = as.data.frame(X_scaled),
                 # specify which predictors to use
                 XFormula = ~ HabitatForereef + HabitatFringing + COTS +
                   Max_DHW + Cyclone1 + Cyclone2 + Cyclone3 + Land_Dist + 
                   Depth_Mean_100m + Depth_Mean_500m +
                   Curvature_Mean_100m + Curvature_Mean_500m +
                   Coral_Mean_Cover + Macroalgae_Mean_Cover + CTB_Mean_Cover,
                 TrData = trait_hmsc,  # add trait data
                 # specify which traits to use
                 TrFormula = ~ Body_Shape + Max_TL_cm + Trophic_Level +
                   Reproductive_Mode + Spawn_Agg, 
                 phyloTree = tree, # add taxonomic tree
                 studyDesign = study_design, # add study design
                 ranLevels = list(site = rL.site, 
                                  year = rL.year),
                distr = "lognormal poisson")

# save the files ready for Hmsc MCMC sampling
save(Y_PA, Y_ABU, X_scaled, trait_hmsc, tree, 
     study_design, rL.site, rL.year, xy_sites,
     file = here("HMSC", "Data", "Hmsc_Data_Ready.RData"))

write.csv(Y_PA,
          here("HMSC", "Data", "Y_Presence_Absence.csv"),
          row.names = FALSE)

write.csv(Y_ABU,
          here("HMSC", "Data", "Y_Abundance.csv"),
          row.names = FALSE)

write.csv(X_scaled,
          here("HMSC", "Data", "X_Scaled.csv"),
          row.names = FALSE)

write.csv(study_design,
          here("HMSC", "Data", "Study_Design.csv"),
          row.names = FALSE)

write.csv(trait_hmsc,
          here("HMSC", "Data", "Traits.csv"),
          row.names = TRUE)

write.tree(tree,
           here("HMSC", "Data", "Tree.tre"))

write.csv(master_PA,
          here("HMSC", "Data", "Master_Presence_Absence_Dataset.csv"),
          row.names = FALSE)

write.csv(master_ABU,
          here("HMSC", "Data", "Master_Abundance_Dataset.csv"),
          row.names = FALSE)
