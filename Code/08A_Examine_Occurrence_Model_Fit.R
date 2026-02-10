#### CONTACT ####
# Courtney Stuart (courtney.e.stuart@gmail.com; courtney.stuart@mansfield.ox.ac.uk)

#### LIBRARIES ####
# install required packages (first run only)
# install.packages(c("easypackages", "conflicted", "tidyr", "dplyr", "here", "ggplot2",
#                    "Hmsc", "coda", "corrplot", "Cairo", "stringr", "PNWColors",
#                    "colorspace", "grDevices", "colorRamps"))

# load packages
library(easypackages)
libraries("conflicted", "tidyr", "dplyr", "here", "ggplot2",
          "Hmsc", "coda", "corrplot", "Cairo", "stringr", 
          "PNWColors", "colorspace", "grDevices", "colorRamps")

# resolve package conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

set.seed(1)

# save some potential color palettes for plotting later
starfish = pnw_palette(name = "Starfish", n = 8, type = "continuous")
bay = pnw_palette(name = "Bay", n = 8, type = "continuous")
anemone = pnw_palette(name = "Anemone", n = 8, type = "continuous")
sailboat = pnw_palette(name = "Sailboat", n = 8, type = "continuous")

#### DIRECTORIES ####
# working directory and relative folder path for here
setwd("E:/Data/StuartC_DPhil_Ch3/")
#set_here("E:/Data/StuartC_DPhil_Ch3/") # set first-time only
here::i_am(".here")
here::here() # verify

# paths to HMSC data and results
data.directory = here("HMSC", "Data")
model.directory = here("HMSC", "Models")

list.files(model.directory)

# read in the results file
nChains = 2
samples = 600
thin = 100
filename = file.path(model.directory, 
                     paste0("PA_model_chains_", as.character(nChains),
                            "_samples_",as.character(samples),
                            "_thin_",as.character(thin),".rda"))
load(filename)

# convert to coda object
mpost = convertToCodaObject(PA_model, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

#### EXPLANATORY POWER ####
# compute fitted values (expected values)
fitted = computePredictedValues(PA_model, expected = TRUE)

# evaluate model fit
MF_fit = evaluateModelFit(hM = PA_model, predY = fitted)

# check model fit metrics
round((mean(MF_fit$TjurR2)), digits = 2)  # mean Tjur R2 across species
round((summary(MF_fit$TjurR2)), digits = 2)  # distribution of R2 values
round((mean(MF_fit$AUC)), digits = 2) # mean AUC across species
round((summary(MF_fit$AUC)), digits = 2) # distribution of AUC values

# create simple histogram plots of TjurR2 and AUC
jpeg(filename = here("Figures", "PA_Model", "AUC_and_TjurR2_Explanatory_Power.jpg"),
     width = 16,
     height = 8,
     units = "in",
     res = 450)

par(mfrow = c(1,2))
hist(MF_fit$AUC, xlab = expression("AUC"~"(explanatory)"), main = NULL, xlim = c(0,1),
     col = bay[1])
abline(v = 0.5, col = "black", lty = "dashed")
hist(MF_fit$TjurR2, xlab = expression("Tjur R"^2~"(explanatory)"), main = NULL, 
     col = bay[1], xlim = c(0,1))
abline(v = 0, col = "black", lty = "dashed")
dev.off()

tjur_fit = data.frame(
  species = colnames(PA_model$Y),
  TjurR2 = MF_fit$TjurR2)

auc_fit = data.frame(
  species = colnames(PA_model$Y),
  AUC = MF_fit$AUC)

#### PREDICTIVE POWER ####
# generates new predictions by sampling from the posterior predictive distribution
# compute predicted values
predicted = computePredictedValues(PA_model, expected = FALSE)

# evaluate predictive power
MF_predictive = evaluateModelFit(hM = PA_model, predY = predicted)

# check metrics
round((mean(MF_predictive$TjurR2)), digits = 2)  # mean Tjur R2 across species
round((summary(MF_predictive$TjurR2)), digits = 2)  # distribution of R2 values
round((mean(MF_predictive$AUC)), digits = 2) # mean AUC across species
round((summary(MF_predictive$AUC)), digits = 2) # distribution of AUC values

# create simple histogram plots of TjurR2 and AUC
jpeg(filename = here("Figures", "PA_Model", "AUC_and_TjurR2_Predictive.jpg"),
     width = 16,
     height = 8,
     units = "in",
     res = 450)
par(mfrow = c(1,2))
hist(MF_predictive$AUC, xlab = expression("AUC"~"(predictive)"), main = NULL,
     xlim = c(0,1), col = bay[2])
abline(v = 0.5, col = "black", lty = "dashed")
hist(MF_predictive$TjurR2, xlab = expression("Tjur R"^2~"(predictive)"), main = NULL,
     xlim = c(0,1), col = bay[2])
abline(v = 0, col = "black", lty = "dashed")
dev.off()

tjur_predictive = data.frame(
  species = colnames(PA_model$Y),
  TjurR2 = MF_predictive$TjurR2)

auc_predictive = data.frame(
  species = colnames(PA_model$Y),
  AUC = MF_predictive$AUC)

#### MULTI-PANEL PERFORMANCE PLOT ####
# four panel plot to show both explanatory and predictive power
jpeg(filename = here("Figures", "PA_Model", "AUC_and_TjurR2_Explanatory_and_Predictive.jpg"),
     width = 8,
     height = 8,
     units = "in",
     res = 450)
par(mar = c(5.1, 4.1, 1.5, 1.5),
    mfrow = c(2,2))
hist(MF_fit$AUC, xlab = expression("AUC"~"(explanatory)"), main = NULL, 
     xlim = c(0,1), col = bay[1])
abline(v = 0.5, col = "black", lty = "dashed")
hist(MF_fit$TjurR2, xlab = expression("Tjur R"^2~"(explanatory)"), main = NULL,
     xlim = c(0,1), col = bay[1])
abline(v = 0, col = "black", lty = "dashed")
hist(MF_predictive$AUC, xlab = expression("AUC"~"(predictive)"), main = NULL,
     xlim = c(0,1), col = bay[3])
abline(v = 0.5, col = "black", lty = "dashed")
hist(MF_predictive$TjurR2, xlab = expression("Tjur R"^2~"(predictive)"), main = NULL,
     xlim = c(0,1), col = bay[3])
abline(v = 0, col = "black", lty = "dashed")
dev.off()

# create a table that includes both the explanatory and predictive metrics
model_fit = (auc_fit %>%
               rename(AUC_Explanatory = AUC) %>%
               left_join(auc_predictive %>%
                           rename(AUC_Predictive = AUC)) %>%
               left_join(tjur_fit %>%
                           rename(TjurR2_Explanatory = TjurR2)) %>%
               left_join(tjur_predictive %>%
                           rename(TjurR2_Predictive = TjurR2)) %>%
               mutate(across(c(AUC_Explanatory, AUC_Predictive, 
                               TjurR2_Explanatory, TjurR2_Predictive), 
                             ~round(., 3))) %>%
               rename(Species = species) %>%
               mutate(Species = str_replace(Species, "\\.", " ")))

# save the data for inclusion in supplementary materials
write.csv(model_fit,
          here("HMSC", "Data", "PA_Model_Fit_Metrics.csv"),
          row.names = FALSE)
