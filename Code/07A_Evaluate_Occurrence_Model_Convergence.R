#### CONTACT ####
# Courtney Stuart (courtney.e.stuart@gmail.com; courtney.stuart@mansfield.ox.ac.uk)

#### LIBRARIES ####
# install required packages (first run only)
# install.packages(c("easypackages", "conflicted", "tidyr", "dplyr", "here", "ggplot2",
#                    "Hmsc", "coda", "stringr", "patchwork", "PNWColors))

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
# working directory and relative folder path for here
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
nChains = 2
samples = 600
thin = 100
filename = file.path(model.directory, 
                     paste0("PA_model_chains_", as.character(nChains),
                            "_samples_",as.character(samples),
                            "_thin_",as.character(thin),".rda"))
load(filename)

# convert to coda object
mpost = convertToCodaObject(PA_model, 
                            spNamesNumbers = c(T,F),
                            covNamesNumbers = c(T,F))

# compute the effective sample sizes for beta
es.beta = effectiveSize(mpost$Beta)
summary(es.beta) # look at the spread of effective sample sizes

# which parameters had an effective sample size for beta < 100?
low.es.beta = es.beta[es.beta < 100]
print(low.es.beta)

# interesting, there seems to be a localized issue with 9 parameters and
# most of these appear to be for Neocirrhites armatus! 

# is this a problem relating to species rarity?
Y_data = as.data.frame(PA_model$Y)
sum(Y_data$Neocirrhites.armatus)
# there were 131  Neocirrhites.armatus presences, so they're not rare...
sum(Y_data$Chromis.xanthura)
# there were 56 Chromis.xanthura presences, so they are a less prevalent species
# but still met the requirement for inclusion

# let's move forward but keep this in mind...

# calculate the PSRF values for beta - ideally, we want all PSRF <= 1.1
psrf.beta = gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
summary(psrf.beta) # look at the spread of values

# check the problematic (PSRF > 1.1) parameters for Neocirrhites.armatus
# and Chromis.xanthura
N.armatus.psrf = (as.data.frame(psrf.beta))[grepl("Neocirrhites\\.armatus", 
                                                  rownames(psrf.beta)), ]
N.armatus.psrf[N.armatus.psrf$`Point est.` > 1.1, ]

C.xanthura.psrf = (as.data.frame(psrf.beta))[grepl("Chromis\\.xanthura", 
                                                  rownames(psrf.beta)), ]
C.xanthura.psrf[C.xanthura.psrf$`Point est.` > 1.1, ]

# as a reminder, here are the ESS issues:
low.es.beta

# the parameters with ESS < 100 above are also those with PSRF > 1.1! this suggests
# highly localized, species-specific issues. the community results should still be 
# solid for interpretation thanks to the strong convergence results overall, but
# use caution if interpreting environmental relationships for these two species 
# specifically (as we have limited confidence in these beta parameter estimates).

# to look at all omega PSRFs we would run the line below; however, we have many
# species pairs (with 143 unique species) so this would take a lot of time and 
# computational effort!!!
# psrf = gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf

# instead, for omega, we will take a sub-sample of 5000 randomly selected species 
# pairs to avoid excessive computations.
tmp = mpost$Omega[[1]]
z = ncol(tmp[[1]])
sel = sample(z, size = 5000)

# here we take the subset of species pairs + loop over the 2 MCMC chains
for(i in 1:length(tmp)){ 
  tmp[[i]] = tmp[[i]][,sel]}

psrf.omega = gelman.diag(tmp, multivariate = FALSE)$psrf
summary(psrf.omega) # look at the spread of values

# what percentage of the beta point estimates are <= 1.1?
round((sum(psrf.beta[, "Point est."] <= 1.1) / 
         length(psrf.beta[, "Point est."]) * 100),
      digits = 2)

# what percentage of the beta upper CI estimates are <= 1.1? (stricter assessment)
round((sum(psrf.beta[, "Upper C.I."] <= 1.1) / 
         length(psrf.beta[, "Upper C.I."]) * 100),
      digits = 2)

# what percentage of the omega point estimates are <= 1.1?
round((sum(psrf.omega[, "Point est."] <= 1.1) / 
         length(psrf.omega[, "Point est."]) * 100),
      digits = 2)

# what percentage of the omega upper CI estimates are <= 1.1? (stricter assessment)
round((sum(psrf.omega[, "Upper C.I."] <= 1.1) / 
         length(psrf.omega[, "Upper C.I."]) * 100),
      digits = 2)

# identify the beta parameters that did not converge
unconverged_beta = as.data.frame(which(psrf.beta[, "Point est."] > 1.1))

# identify the omega parameters that did not converge (of the 5000 subsampled)
unconverged_omega = as.data.frame(which(psrf.omega[, "Point est."] > 1.1))

# the PSRF values for beta look great - over 99% convergence for point estimates
# and approx. 95% convergence for the stricter upper CI assessment!

# the PSRF values for omega look good as well - over 99% convergence for point
# estimates and approx. 93% convergence for the stricter upper CI assessment!

# now check the gamma parameters
if("Gamma" %in% names(mpost)) {
  gamma_params = mpost$Gamma
  
  # convert to mcmc.list if needed
  gamma_params = as.mcmc.list(gamma_params)
  
  # check convergence
  library(coda)
  psrf.gamma = gelman.diag(gamma_params, multivariate = FALSE)
}

# look at the spread of ESS for the gamma parameters
es.gamma = effectiveSize(gamma_params)
print(summary(es.gamma))

# calculate the PSRFs
gamma_params = mpost$Gamma
psrf.gamma = gelman.diag(gamma_params, multivariate = FALSE)

# look at the spread of gamma PSRFs
summary(psrf.gamma$psrf)
max((psrf.gamma$psrf)[, "Point est."])

# what percentage of the gamma point estimates are <= 1.1?
round((sum((psrf.gamma$psrf)[, "Point est."] <= 1.1) /
         length((psrf.gamma$psrf)[, "Point est."]) * 100),
      digits = 2)

# what percentage of the gamma upper CI estimates are <= 1.1? (stricter assessment)
round((sum((psrf.gamma$psrf)[, "Upper C.I."] <= 1.1) / 
         length((psrf.gamma$psrf)[, "Upper C.I."]) * 100),
      digits = 2)

#### MULTI-PANEL PSRF PLOT ####
# make a three panel plot to show the PSRF spread for each of the paramaters

# prepare the data
df_beta = data.frame(psrf = psrf.beta[, "Point est."], 
                     parameter = "Beta")
df_omega = data.frame(psrf = psrf.omega[, "Point est."], 
                      parameter = "Omega")
df_gamma = data.frame(psrf = (psrf.gamma$psrf)[, "Point est."], 
                      parameter = "Gamma")

# create the individual plots
p1 = ggplot(df_beta, aes(x = psrf)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "PSRF (beta)", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

p2 = ggplot(df_omega, aes(x = psrf)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "PSRF (omega)", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

p3 = ggplot(df_gamma, aes(x = psrf)) +
  geom_histogram(bins = 30, fill = "gray70", color = "black") +
  geom_vline(xintercept = 1.1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "PSRF (gamma)", y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# combine plots with automatic labels
combined_plot = p1 + p2 + p3 + 
  plot_annotation(tag_levels = 'A')

# save the plot as a figure
ggsave(here("Figures", "PA_model", "Convergence_PSRFs.jpeg"),
       plot = combined_plot,
       width = 8,
       height = 4,
       dpi = 450)

#### TRACE PLOTS ####
# extract beta parameters from all chains
beta_chain1 = mpost$Beta[[1]]
beta_chain2 = mpost$Beta[[2]]

# create data frames for each chain
df_beta1 = data.frame(
  iteration = 1:nrow(beta_chain1),
  chain = "Chain 1",
  beta_chain1)

df_beta2 = data.frame(
  iteration = 1:nrow(beta_chain2),
  chain = "Chain 2",
  beta_chain2)

# combine chains
df_beta = bind_rows(df_beta1, df_beta2)

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

# plot the sampled beta parameters
mar = c(5.1, 4.1, 4.1, 2.1) 
ggplot(df_beta_sample, 
       aes(x = iteration, 
           y = value, 
           color = chain)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(bay[1], anemone[3])) +
  scale_fill_manual(values = c(bay[1], anemone[3])) + 
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_Random_Sample_2Chains_600Samples_100Thin.jpg"),
       width = 8, height = 10, units = "in", dpi = 300)

# look at a few random species from the community as examples
# look at all parameters for Abudefduf.septemfasciatus
df_beta_sp1 = df_beta_long %>%
  dplyr::filter(grepl("Abudefduf\\.septemfasciatus", parameter))

ggplot(df_beta_sp1, 
       aes(x = iteration, 
           y = value, 
           color = chain)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(bay[1], anemone[3])) +
  scale_fill_manual(values = c(bay[1], anemone[3])) +
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_A_septemfasciatus_2Chains_600Samples_100Thin.jpg"),
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
  scale_color_manual(values = c(bay[1], anemone[3])) +
  scale_fill_manual(values = c(bay[1], anemone[3])) +
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_A_triostegus_2Chains_600Samples_100Thin.jpg"),
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
  scale_color_manual(values = c(bay[1], anemone[3])) +
  scale_fill_manual(values = c(bay[1], anemone[3])) +
  theme_bw() +
  labs(x = "Posterior sample",
       y = "Parameter value") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       filename = here("Figures", "PA_Model",
                       "Beta_Trace_S_altipinnis_2Chains_600Samples_100Thin.jpg"),
       width = 8, height = 10, units = "in", dpi = 300)
