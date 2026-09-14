#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### SYSTEM SETTINGS ####
# Set environment variables BEFORE any heavy linear-algebra is used / before loading
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
setwd("E:/Data/StuartC_DPhil_Ch3/")
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

# #### RAW ABUNDANCE MODEL ####
# ##### DEFINE THE MODEL #####
# ABU_model = Hmsc(
#   Y = Y_ABU,
#   XData = data_ABU,
#   XFormula = ~ Habitat + COTS + Max_DHW + Cyclone +
#     Land_Dist + Depth_Mean_100m + Depth_Mean_500m +
#     Curvature_Mean_100m + Curvature_Mean_500m +
#     Coral_Mean_Cover + Macroalgae_Mean_Cover + CTB_Mean_Cover,
#   TrData = trait_hmsc,
#   TrFormula = ~ Body_Shape + Max_TL_cm + Trophic_Level +
#     Reproductive_Mode + Spawn_Agg,
#   phyloTree = tree,
#   studyDesign = study_design,
#   ranLevels = list(site = rL.site,
#                    year = rL.year),
#   distr = "lognormal poisson",
#   XScale = TRUE,
#   TrScale = TRUE)
# 
# ##### MCMC SETTINGS #####
# nParallel = 4
# nChains = 4
# samples = 1000 # 1000 samples per chain = 4000 total
# thin = 50
# transient = 10000
# # total iterations per chain = transient + (samples * thin) = 
# # 60000 iterations = 10000 transient + (1000 samples * 50 thin)
# 
# ##### RUN THE MODEL #####
# cat("ABU Model - thin =", thin, ", transient =", transient, "\n")
# cat("Total iterations per chain:", transient + (samples * thin), "\n")
# cat("Total posterior samples:", samples * nChains, "\n\n")
# 
# start_time = Sys.time()
# ABU_model = sampleMcmc(ABU_model, 
#                        thin = thin, 
#                        samples = samples, 
#                        transient = transient,
#                        nChains = nChains, 
#                        nParallel = nParallel,
#                        initPar = "fixed effects")
# end_time = Sys.time()
# 
# cat("Completed in:", difftime(end_time, start_time, units = "mins"), "minutes\n")
# 
# ##### SAVE THE OUTPUTS #####
# filename = file.path(model_directory, 
#                      paste0("ABU_model_chains_", nChains, 
#                             "_total_samples_", samples * nChains,  # total samples
#                             "_thin_", thin, ".rda"))
# save(ABU_model, file = filename)
# cat("Saved:", filename, "\n\n")

#### CONDITIONAL ABUNDANCE MODEL ####
##### CREATE CONDITIONAL ABUNDANCE DATA #####
# Convert zeros to NA (i.e., abundance conditional on presence)
Y_ABU_conditional = Y_ABU
Y_ABU_conditional[Y_ABU_conditional == 0] = NA

##### DEFINE THE MODEL #####
ABU_conditional_model = Hmsc(
  Y = Y_ABU_conditional,
  XData = data_ABU,
  XFormula = ~ Habitat + COTS + Max_DHW + Cyclone +
    Land_Dist + Depth_Mean_100m + Depth_Mean_500m +
    Curvature_Mean_100m + Curvature_Mean_500m +
    Coral_Mean_Cover + Macroalgae_Mean_Cover + CTB_Mean_Cover,
  TrData = trait_hmsc,
  TrFormula = ~ Body_Shape + Max_TL_cm + Trophic_Level +
    Reproductive_Mode + Spawn_Agg,
  phyloTree = tree,
  studyDesign = study_design,
  ranLevels = list(site = rL.site,
                   year = rL.year),
  distr = "lognormal poisson",
  XScale = TRUE,
  TrScale = TRUE)

##### MCMC SETTINGS #####
nParallel = 4
nChains = 4
samples = 1000 # 1000 samples per chain = 4000 total
thin = 50
transient = 10000
# total iterations per chain = transient + (samples * thin) = 
# 60000 iterations = 10000 transient + (1000 samples * 50 thin)

##### RUN THE MODEL #####
cat("ABU conditional model - thin =", thin, ", transient =", transient, "\n")
cat("Total iterations per chain:", transient + (samples * thin), "\n")
cat("Total posterior samples:", samples * nChains, "\n\n")

start_time = Sys.time()

ABU_conditional_model = sampleMcmc(
  ABU_conditional_model, 
  thin = thin, 
  samples = samples, 
  transient = transient,
  nChains = nChains, 
  nParallel = nParallel,
  initPar = "fixed effects")

end_time = Sys.time()

cat("Completed in:", difftime(end_time, start_time, units = "mins"), "minutes\n")

##### SAVE THE OUTPUTS #####
filename = file.path(
  model_directory, 
  paste0("ABU_conditional_model_chains_", nChains, 
         "_total_samples_", samples * nChains,
         "_thin_", thin, ".rda")
)

save(ABU_conditional_model, file = filename)

cat("Saved:", filename, "\n\n")