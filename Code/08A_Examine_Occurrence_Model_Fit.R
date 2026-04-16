#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

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
# working directory and relative folder path for here()
setwd("E:/Data/StuartC_DPhil_Ch3/")
#set_here("E:/Data/StuartC_DPhil_Ch3/") # set first-time only
here::i_am(".here")
here::here() # verify

# paths to HMSC data and results
data.directory = here("HMSC", "Data")
model.directory = here("HMSC", "Models")

list.files(model.directory)

# read in the results file
nChains = 4
samples = 1000
thin = 100
filename = file.path(model.directory, 
                     paste0("PA_model_chains_", as.character(nChains),
                            "_samples_",as.character(samples),
                            "_thin_",as.character(thin),".rda"))
load(filename)

# convert to coda object
mpost = convertToCodaObject(PA_model, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

#### EXPLANATORY PERFORMANCE ####
# compute fitted values (expected values)
fitted = computePredictedValues(PA_model, expected = TRUE)

# evaluate model fit
MF_fit = evaluateModelFit(hM = PA_model, predY = fitted)

# check model fit metrics
round((mean(MF_fit$TjurR2)), digits = 2)  # mean Tjur R2 across species
round((summary(MF_fit$TjurR2)), digits = 2)  # distribution of Tjur R2 values
round((mean(MF_fit$AUC)), digits = 2) # mean AUC across species
round((summary(MF_fit$AUC)), digits = 2) # distribution of AUC values

# create simple histogram plots of Tjur R2 and AUC
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

#### IN-SAMPLE POSTERIOR PREDICTIVE PERFORMANCE ####
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

# create simple histogram plots of AUC and Tjur R2
jpeg(filename = here("Figures", "PA_Model", "AUC_and_TjurR2_In_Sample_Predictive.jpg"),
     width = 16,
     height = 8,
     units = "in",
     res = 450)
par(mfrow = c(1,2))
hist(MF_predictive$AUC, xlab = expression("AUC"~"(in-sample posterior predictive)"), main = NULL,
     xlim = c(0,1), col = bay[2])
abline(v = 0.5, col = "black", lty = "dashed")
hist(MF_predictive$TjurR2, xlab = expression("Tjur R"^2~"(in-sample posterior predictive)"), main = NULL,
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
# four panel plot to show both explanatory and in-sample poserior predictive performance
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
hist(MF_predictive$AUC, xlab = expression("AUC"~"(in-sample posterior predictive)"), main = NULL,
     xlim = c(0,1), col = bay[3])
abline(v = 0.5, col = "black", lty = "dashed")
hist(MF_predictive$TjurR2, xlab = expression("Tjur R"^2~"(in-sample posterior predictive)"), main = NULL,
     xlim = c(0,1), col = bay[3])
abline(v = 0, col = "black", lty = "dashed")
dev.off()

# create a table that includes both the explanatory and in-sample posterior predictive 
# performance metrics
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

# how many species had AUC >= 0.99?
model_fit %>%
  summarise(
    n_auc_expl = sum(AUC_Explanatory >= 0.99, na.rm = TRUE),
    n_auc_pred = sum(AUC_Predictive >= 0.99, na.rm = TRUE))

# which species had the lowest explanatory AUC?
model_fit %>%
  slice_min(AUC_Explanatory, n = 1)

# which species had the highest explanatory AUC?
model_fit %>%
  slice_max(AUC_Explanatory, n = 1)

# which species had the lowest in-sample posterior predictive AUC?
model_fit %>%
  slice_min(AUC_Predictive, n = 1)

# which species had the highest in-sample posterior predictive AUC?
model_fit %>%
  slice_max(AUC_Predictive, n = 1)

# which species had the lowest explanatory Tjur R2?
model_fit %>%
  slice_min(TjurR2_Explanatory, n = 1)

# which species had the highest explanatory Tjur R2?
model_fit %>%
  slice_max(TjurR2_Explanatory, n = 1)

# which species had the lowest in-sample posterior predictive Tjur R2?
model_fit %>%
  slice_min(TjurR2_Predictive, n = 1)

# which species had the highest in-sample posterior predictive Tjur R2?
model_fit %>%
  slice_max(TjurR2_Predictive, n = 1)

# what percentage of species had AUCs >= 0.8 and Tjur R2s >= 0.3?
model_fit %>%
  summarise(
    pct_auc_expl = mean(AUC_Explanatory >= 0.8) * 100,
    pct_auc_pred = mean(AUC_Predictive >= 0.8) * 100,
    pct_tjur_expl = mean(TjurR2_Explanatory >= 0.3) * 100,
    pct_tjur_pred = mean(TjurR2_Predictive >= 0.3) * 100)
