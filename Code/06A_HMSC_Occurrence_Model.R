#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### SYSTEM SETTINGS ####
# set environment variables BEFORE any heavy linear-algebra is used and before loading 
# packages that call BLAS.
Sys.setenv(OPENBLAS_NUM_THREADS = "1",
           OMP_NUM_THREADS      = "1",
           MKL_NUM_THREADS      = "1")
           
#### LIBRARIES ####
# install packages (first run only)
#install.packages(c("here", "easypackages", "Hmsc", "tidyr", "dplyr", "ggplot2", "coda"))

# load packages
library(easypackages)
libraries("here", "Hmsc", "tidyr", "dplyr", "ggplot2", "coda")

#### DIRECTORIES ####
# working directory and relative folder path for here()
#setwd("E:/Data/StuartC_DPhil_Ch3/")
#set_here("E:/Data/StuartC_DPhil_Ch3/") # set first-time only
here::i_am(".here")
here::here() # verify
model_directory = as.character(here("HMSC", "Models"))

# load prepared data from script "05_Check_And_Prep_Model_Inputs.R"
load(here("HMSC", "Data", "Hmsc_Data_Ready.RData"))

# set random seed
set.seed(1)

# check the trait data really quickly to confirm that they'll all be helpful
summary(trait_hmsc)

# the schooling, shoaling, and solitary trait columns are very zero-inflated
# with such few observations of each trait, I don't think it would be helpful
# or wise to include them as the estimates would have extremely high uncertainty. 
# exclude these three traits from the models...

#### DEFINE THE MODELS ####
# define the presence-absence model
PA_model = Hmsc(Y = Y_PA,
                XData = as.data.frame(X_scaled),
                # specify which predictors to use
                XFormula = ~ HabitatForereef + HabitatFringing + COTS +
                  Max_DHW + Cyclone1 + Cyclone2 + Cyclone3 + Land_Dist + 
                  Depth_Mean_100m + Depth_Mean_500m +
                  Curvature_Mean_100m + Curvature_Mean_500m +
                  Coral_Mean_Cover + Macroalgae_Mean_Cover + CTB_Mean_Cover,
                TrData = trait_hmsc,  
                # specify which traits to use
                TrFormula = ~ Body_Shape + Max_TL_cm + Trophic_Level +
                  Reproductive_Mode + Spawn_Agg, 
                phyloTree = tree, 
                studyDesign = study_design, 
                ranLevels = list(site = rL.site, 
                                 year = rL.year),
                distr = "probit")

#### MCMC SETTINGS ####
nParallel = 2
nChains = 2
samples = 300 # 300 samples per chain = 600 total
thin = 100
transient = 5000
# total iterations per chain = transient + (samples * thin) = 
# 5000 + 30000 = 35000

##### PA Model #####
cat("PA Model - thin =", thin, ", transient =", transient, "\n")
cat("Total iterations per chain:", transient + (samples * thin), "\n")
cat("Total posterior samples:", samples * nChains, "\n\n")

start_time = Sys.time()
PA_model = sampleMcmc(PA_model, 
                      thin = thin, 
                      samples = samples, 
                      transient = transient,
                      nChains = nChains, 
                      nParallel = nParallel,
                      initPar = "fixed effects")
end_time = Sys.time()

cat("Completed in:", difftime(end_time, start_time, units = "mins"), "minutes\n")

# save outputs
filename = file.path(model_directory, 
                     paste0("PA_model_chains_", nChains, 
                            "_samples_", samples * nChains,  # total samples
                            "_thin_", thin, ".rda"))
save(PA_model, file = filename)
cat("Saved:", filename, "\n\n")