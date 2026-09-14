#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### LIBRARIES ####
# install required packages (first run only)
# install.packages(c("easypackages", "conflicted", "tidyr", "dplyr", "here", "ggplot2",
#                    "Hmsc", "coda", "stringr", "patchwork", "PNWColors"))

# load packages
library(easypackages)
libraries("conflicted", "tidyr", "dplyr", "here", "ggplot2",
          "Hmsc", "coda", "stringr", "patchwork", "PNWColors")

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

#### EFFECTIVE SAMPLE SIZES & PSRFs ####
# examine MCMC convergence using the potential scale reduction factor (PSRF) of the beta, 
# omega, and gamma parameters.

# read in the results file
nChains = 4
samples = 4000
thin = 50
filename = file.path(model.directory, 
                     paste0("PA_model_chains_", as.character(nChains),
                            "_total_samples_",as.character(samples),
                            "_thin_",as.character(thin),".rda"))
load(filename)

# convert to coda object
mpost = convertToCodaObject(PA_model, 
                            spNamesNumbers = c(T,F),
                            covNamesNumbers = c(T,F))

##### BETA #####
# compute the effective sample sizes for beta
es.beta = effectiveSize(mpost$Beta)
summary(es.beta) # look at the spread of effective sample sizes

# which parameters had an effective sample size for beta <= 100?
low.es.beta = es.beta[es.beta <= 100]
print(low.es.beta)

# is this a problem relating to species rarity? look at some examples...
Y_data = as.data.frame(PA_model$Y)
sum(Y_data$Centropyge.flavissima)
sum(Y_data$Neocirrhites.armatus)
sum(Y_data$Zebrasoma.scopas)
# rarity doesn't seem to be the issue, let's move forward but keep this in mind...

# calculate the PSRF values for beta - ideally, we want all PSRF <= 1.1
psrf.beta = gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
summary(psrf.beta) # look at the spread of values

# identify the beta parameters that did not converge
unconverged_beta = as.data.frame(which(psrf.beta[, "Point est."] > 1.1))

# what percentage of the beta point estimates are <= 1.1?
round((sum(psrf.beta[, "Point est."] <= 1.1) / 
         length(psrf.beta[, "Point est."]) * 100),
      digits = 2)

# what percentage of the beta upper CI estimates are <= 1.1? (stricter assessment)
round((sum(psrf.beta[, "Upper C.I."] <= 1.1) / 
         length(psrf.beta[, "Upper C.I."]) * 100),
      digits = 2)

# PSRF values for beta indicate excellent overall convergence, with 98.53% of point
# estimates and 95.93% of upper CI estimates <= 1.1. this suggests reliable estimation
# of species’ responses to environmental covariates (fixed effects).

# save data for species–covariate pairs for which β parameters did not achieve
# satisfactory convergence
i = which(psrf.beta[, "Point est."] > 1.1)
x = strsplit(gsub("^B\\[|\\]$", "", rownames(psrf.beta)[i]), ", ")
supp_table = data.frame(
  Species = sapply(x, `[`, 2), Covariate = sapply(x, `[`, 1),
  ESS = es.beta[rownames(psrf.beta)[i]],
  PSRF_point_est = psrf.beta[i, 1], PSRF_upper_CI = psrf.beta[i, 2])
supp_table[sapply(supp_table, is.numeric)] = round(supp_table[sapply(supp_table, is.numeric)], 2)
write.csv(supp_table, here("HMSC", "Data", "Unconverged_Betas_PA_Model.csv"),
          row.names = FALSE)

##### OMEGA #####
# to look at all omega PSRFs we run the line below
# WARNING, we have many species pairs (with 149 unique species) so this takes
# a lot of time and computational effort!!!
psrf.omega = gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf

# if we're short on time, we can instead take a sub-sample of 5000 randomly
# selected species pairs to avoid excessive computations and get a rough look
# at things
# tmp = mpost$Omega[[1]]
# z = ncol(tmp[[1]])
# sel = sample(z, size = 5000)
# 
# # here we take the subset of species pairs + loop over the 4 MCMC chains
# for(i in 1:length(tmp)){
#    tmp[[i]] = tmp[[i]][,sel]}
# 
# psrf.omega = gelman.diag(tmp, multivariate = FALSE)$psrf
# summary(psrf.omega) # look at the spread of values

# keep only the upper triangle of the 149 x 149 omega matrix, excluding the
# diagonal, so that each unique species pair is counted only once
keep = upper.tri(matrix(FALSE, 149, 149), diag = FALSE)

# select the corresponding PSRF rows
psrf.omega.unique = psrf.omega[as.vector(keep), , drop = FALSE]
summary(psrf.omega.unique) # look at the spread of values

# check number of unique species pairs
nrow(psrf.omega.unique) # this should be 11026

# identify omega estimates that did not converge (unique species pairs with point
# estimate PSRF > 1.1)
unconverged_omega = as.data.frame(
  which(psrf.omega.unique[, "Point est."] > 1.1))
nrow(unconverged_omega)

# what percentage of the omega point estimates are <= 1.1?
round((sum(psrf.omega.unique[, "Point est."] <= 1.1) /
         length(psrf.omega.unique[, "Point est."]) * 100),
      digits = 2)

# what percentage of the omega upper CI estimates are <= 1.1? (stricter assessment)
round((sum(psrf.omega.unique[, "Upper C.I."] <= 1.1) /
         length(psrf.omega.unique[, "Upper C.I."]) * 100),
      digits = 2)

# PSRF values for omega indicate good overall convergence, with 94.6% of point
# estimates and 87.23% of upper CI estimates <= 1.1.

##### GAMMA #####
# now check the gamma parameters
gamma_params = as.mcmc.list(mpost$Gamma)

# look at the spread of ESS for the gamma parameters
es.gamma = effectiveSize(gamma_params)
print(summary(es.gamma))

# which parameters had an effective sample size for gamma < 100?
low.es.gamma = es.gamma[es.gamma <= 100]
print(low.es.gamma)

# calculate the PSRFs
psrf.gamma = gelman.diag(gamma_params, multivariate = FALSE)$psrf

# look at the spread of gamma PSRFs
summary(psrf.gamma)

# what percentage of the gamma point estimates are <= 1.1?
round((sum(psrf.gamma[, "Point est."] <= 1.1) /
         length(psrf.gamma[, "Point est."]) * 100),
      digits = 2)

# what percentage of the gamma upper CI estimates are <= 1.1? (stricter assessment)
round((sum(psrf.gamma[, "Upper C.I."] <= 1.1) /
         length(psrf.gamma[, "Upper C.I."]) * 100),
      digits = 2)

# identify the gamma parameters that potentially did not converge
# (point estimate PSRF > 1.1)
unconverged_gamma = as.data.frame(
  which(psrf.gamma[, "Point est."] > 1.1))

# PSRF values for gamma indicate excellent convergence, with 100% of point
# estimates and upper CI estimates <= 1.1.

#### MULTI-PANEL PSRF PLOT ####
# make a three panel plot to show the PSRF spread for each of the parameters

# prepare the data (ordered here to match the order of the manuscript results section)
df_beta = data.frame(psrf = psrf.beta[, "Point est."], 
                     parameter = "Beta")

df_gamma = data.frame(psrf = psrf.gamma[, "Point est."], 
                     parameter = "Gamma")

df_omega = data.frame(psrf = psrf.omega.unique[, "Point est."], 
                      parameter = "Omega")

# create the individual plots
p1 = ggplot(df_beta, aes(x = psrf)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), 
                     labels = scales::label_number(accuracy = 0.1)) +
  labs(x = "PSRF (beta)", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

p2 = ggplot(df_gamma, aes(x = psrf)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), 
                     labels = scales::label_number(accuracy = 0.1), 
                     breaks = scales::breaks_width(0.1)) +
  labs(x = "PSRF (gamma)", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

p3 = ggplot(df_omega, aes(x = psrf)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), 
                     labels = scales::label_number(accuracy = 0.1)) +
  labs(x = "PSRF (omega)", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# combine plots with automatic labels
combined_plot = p1 + p2 + p3 + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_suffix = ')')
combined_plot

# save the plot as a figure
ggsave(here("Figures", "PA_model", "Convergence_PSRFs.jpeg"),
       plot = combined_plot,
       width = 8,
       height = 4,
       dpi = 450)

#### TRACE PLOTS ####
# extract beta parameters from all four MCMC chains
beta_chain1 = mpost$Beta[[1]]
beta_chain2 = mpost$Beta[[2]]
beta_chain3 = mpost$Beta[[3]]
beta_chain4 = mpost$Beta[[4]]

# create a data frame for each chain
df_beta1 = data.frame(
  iteration = 1:nrow(beta_chain1),
  chain = "Chain 1",
  beta_chain1)

df_beta2 = data.frame(
  iteration = 1:nrow(beta_chain2),
  chain = "Chain 2",
  beta_chain2)

df_beta3 = data.frame(
  iteration = 1:nrow(beta_chain3),
  chain = "Chain 3",
  beta_chain3)

df_beta4 = data.frame(
  iteration = 1:nrow(beta_chain4),
  chain = "Chain 4",
  beta_chain4)

# combine chains
df_beta = bind_rows(df_beta1, df_beta2, df_beta3, df_beta4)

# reshape to long format
df_beta_long = df_beta %>%
  pivot_longer(cols = -c(iteration, chain),
               names_to = "parameter",
               values_to = "value")

# get all unique parameter names
all_params = unique(df_beta_long$parameter)

# randomly sample 10 parameters
sampled_params = sample(all_params, 10)

# filter to just the sampled parameters
df_beta_sample = df_beta_long %>%
  dplyr::filter(parameter %in% sampled_params)

# plot the sampled parameters
mar = c(5.1, 4.1, 4.1, 2.1) 
ggplot(df_beta_sample, 
       aes(x = iteration, 
           y = value, 
           color = chain)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) +
  scale_fill_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) + 
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_Random_Sample_4Chains_4000Samples_50Thin.jpg"),
       width = 8, height = 10, units = "in", dpi = 300)

# look at a few random species from the community as examples...

# look at all parameters for Abudefduf.septemfasciatus
df_beta_sp1 = df_beta_long %>%
  dplyr::filter(grepl("Abudefduf\\.septemfasciatus", parameter))

ggplot(df_beta_sp1, 
       aes(x = iteration, 
           y = value, 
           color = chain)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) +
  scale_fill_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) + 
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_A_septemfasciatus_4Chains_4000Samples_50Thin.jpg"),
       width = 8, height = 10, units = "in", dpi = 300)

# look at all parameters for Acanthurus.triostegus
df_beta_sp2 = df_beta_long %>%
  dplyr::filter(grepl("Acanthurus\\.triostegus", parameter))

ggplot(df_beta_sp2, 
       aes(x = iteration, 
           y = value, 
           color = chain)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) +
  scale_fill_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) + 
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_A_triostegus_4Chains_4000Samples_50Thin.jpg"),
       width = 8, height = 10, units = "in", dpi = 300)

# look at all parameters for Scarus.altipinnis
df_beta_sp3 = df_beta_long %>%
  dplyr::filter(grepl("Scarus\\.altipinnis", parameter))

ggplot(df_beta_sp3, 
       aes(x = iteration, 
           y = value, 
           color = chain)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) +
  scale_fill_manual(values = c(starfish[5], anemone[1], starfish[6], anemone[3])) + 
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_S_altipinnis_4Chains_4000Samples_50Thin.jpg"),
       width = 8, height = 10, units = "in", dpi = 300)
