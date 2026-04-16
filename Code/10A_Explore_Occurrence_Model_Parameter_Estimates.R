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

#### SAVE ABBREVIATED SPECIES NAMES ####
make_unique_abbrev = function(names) {
  # expects "Genus.species" (first dot used as separator)
  genus  = sub("\\..*$", "", names)
  species = sub("^[^.]*\\.", "", names)
  n = length(genus)
  len = rep(1L, n)                        # start with 1-letter abbrev
  build = function(l) paste0(substr(genus, 1, l), ". ", species)
  
  labels = build(len)
  while(any(duplicated(labels))) {
    # indices involved in duplicates
    dup = labels[duplicated(labels) | duplicated(labels, fromLast = TRUE)]
    idx = which(labels %in% dup)
    # increase only those indices (cap at full genus length)
    len[idx] = pmin(nchar(genus[idx]), len[idx] + 1L)
    labels = build(len)
    # if we've used full genus and still have duplicates, force uniqueness
    if(all(len == nchar(genus)) && any(duplicated(labels))) {
      labels = make.unique(labels, sep = " ")
      break}}
  labels}

abbrv_species = make_unique_abbrev(PA_model$spNames)
full_species = PA_model$spNames

#### SPECIES NICHES ####
# I want to use abbreviated species names for plotting, set these as these column names 
colnames(PA_model$Y) = abbrv_species
PA_model$spNames = abbrv_species

# construct a beta plot showing the estimates of species' niche parameters
postBeta = getPostEstimate(PA_model, parName = "Beta")

# create the plot
jpeg(filename = here("Figures", "PA_Model", "Beta_Plot.jpg"), 
     width = 10, 
     height = 12, 
     units = "in", 
     res = 450)

par(mar = c(9, 11, 0, 0), 
    mgp = c(0, 0, 0), 
    oma = c(0, 0, 0, 0),
    omi = c(0, 0, 0, 0),
    xpd = NA)

plotBeta(PA_model, 
         post = postBeta, 
         plotTree = FALSE, 
         spNamesNumbers = c(TRUE, FALSE),
         covNamesNumbers = c(TRUE, FALSE), 
         colors = colorRampPalette(c("darkred","white","darkblue")),
         mar = c(9, 11, 0, 0),
         mgp = c(0, 2, 0),
         supportLevel = 0.95)

dev.off()

# the intercept here refers to the reference level - HabitatBackreef in a non-cyclone
# year (0).

# the full beta plot may be hard to read because of the large number of species (n = 143)!
# instead, split the species in half and make two beta plots.

# use numeric indices instead of species names to define groups
group1 = rev(1:72)
group2 = rev(73:142)

# create and save the first beta plot
jpeg(filename = here("Figures", "PA_Model", "Beta_Plot_1.jpg"), 
     width = 6, 
     height = 12, 
     units = "in", 
     res = 600)

par(mar = c(10, 6.75, 0, 0), 
    mgp = c(0, 0, 0), 
    oma = c(0, 0, 0, 0),
    omi = c(0, 0, 0, 0),
    xpd = NA)

plotBeta(PA_model,
         post = postBeta, 
         param = "Support", 
         plotTree = FALSE, 
         SpeciesOrder = "Vector",
         SpVector = group1,  # numeric indices, not names
         spNamesNumbers = c(TRUE, FALSE),
         covNamesNumbers = c(TRUE, FALSE), 
         colors = colorRampPalette(c("darkred","white","darkblue")),
         mar = c(10, 6.75, 0, 0),
         mgp = c(0, 2, 0),
         supportLevel = 0.95)

dev.off()

# and now do the same for the second beta plot
jpeg(filename = here("Figures", "PA_Model", "Beta_Plot_2.jpg"), 
     width = 6, 
     height = 12, 
     units = "in", 
     res = 600)

par(mar = c(10, 6.75, 0, 0), 
    mgp = c(0, 0, 0), 
    oma = c(0, 0, 0, 0),
    omi = c(0, 0, 0, 0),
    xpd = NA)

plotBeta(PA_model, 
         post = postBeta, 
         param = "Support",
         plotTree = FALSE, 
         SpeciesOrder = "Vector",
         SpVector = group2,
         spNamesNumbers = c(TRUE, FALSE),
         covNamesNumbers = c(TRUE, FALSE), 
         colors = colorRampPalette(c("darkred","white","darkblue")),
         mar = c(10, 6.75, 0, 0),
         mgp = c(0, 2, 0),
         supportLevel = 0.95)

dev.off()

# closer look
thresh = 0.95

# tidy the matrices (postBeta: rows = variables, cols = species)
rownames(postBeta$mean) = PA_model$covNames
beta_df = as.data.frame(as.table(postBeta$mean))
names(beta_df) = c("Variable", "Species", "Mean")
beta_df$Support_Pos = as.vector(postBeta$support)
beta_df$Support_Neg = as.vector(postBeta$supportNeg)

# classify by 95% posterior sign
beta_df = beta_df %>%
  mutate(Sign95 = case_when(
    Support_Pos >= thresh ~ "Positive",
    Support_Neg >= thresh ~ "Negative",
    TRUE                 ~ "No_Effect"))

# we don't want to plot beta parameters for the species-variable combinations with 
# ESS < 100, so identify and remove these
es.beta = effectiveSize(mpost$Beta)
low.es.beta = es.beta[es.beta < 100]
print(low.es.beta)

bad_pairs = do.call(rbind, lapply(names(low.es.beta), function(x) {
  s = sub("^B\\[", "", x)
  s = sub("\\]$", "", s)
  parts = strsplit(s, ",\\s*")[[1]]
  data.frame(
    Variable = parts[1],
    Species  = parts[2],
    stringsAsFactors = FALSE)}))

bad_pairs = bad_pairs %>%
  mutate(Species = sub("^(.)([^.]*)\\.([^.]*)$", "\\1. \\3", Species))

# remove these species-environment pairs
beta_df_filt = beta_df %>%
  anti_join(bad_pairs, by = c("Variable", "Species"))

# counts per predictor
counts = beta_df_filt %>% group_by(Variable) %>%
  summarise(N_Positive = sum(Sign95 == "Positive"),
            N_Negative = sum(Sign95 == "Negative"),
            N_No_Effect = sum(Sign95 == "No_Effect")) %>%
  mutate(label = paste0("+:", N_Positive, " / -:", N_Negative, " / 0:", N_No_Effect))

print(counts)

# quick faceted plot (means only; color by sign)
ann = beta_df %>% 
  group_by(Variable) %>% 
  summarise(ypos = max(Mean, na.rm = TRUE) * 1.05) %>% 
  left_join(counts, by = "Variable")

ggplot(beta_df, 
       aes(x = Species, y = Mean, colour = Sign95)) +
  geom_point(size = 1) +
  facet_wrap(~Variable, scales = "free_y", ncol = 4) +
  coord_flip() +
  theme_minimal() +
  geom_text(data = ann, 
            aes(x = Inf,
                y = ypos, 
                label = label), 
            inherit.aes = FALSE, hjust = 1.05, size = 3)

# one figure per variable:
# now, for each variable, plot the posterior mean beta for each species. shade the dots
# based on whether the coefficient is negative, positive, or negligible at the 95% support
# level. because there are so many species, also add a tally in the legend at the bottom.
# save a separate jpeg file for each variable. 
vars = setdiff(PA_model$covNames, "(Intercept)")
out_dir = here("Figures", "PA_Model")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (v in vars) {
  
  d = filter(beta_df_filt, Variable == v)
  
  # counts for legend labels
  n_pos = sum(d$Sign95 == "Positive")
  n_neg = sum(d$Sign95 == "Negative")
  n_none = sum(d$Sign95 == "No_Effect")
  
  legend_labels = c(
    Positive  = paste0("Positive (n = ", n_pos, ")"),
    Negative  = paste0("Negative (n = ", n_neg, ")"),
    No_Effect = paste0("None (n = ", n_none, ")"))
  
  p = ggplot(d, aes(Species, Mean, colour = Sign95)) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 6, face = "italic"),
      plot.title = element_text(hjust = 0.4),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box.margin = margin(t = -7.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 10)) +
    scale_colour_manual(
      name   = "Effect (95% posterior support)",
      values = c(
        Positive  = bay[1],
        Negative  = bay[8],
        No_Effect = bay[5]),
      labels = legend_labels,
      guide  = guide_legend(nrow = 1)) +
    labs(title = v, x = NULL, y = "Posterior mean (Beta)")

  ggsave(
    filename = file.path(out_dir, paste0(v, ".jpeg")),
    plot = p,
    device = "jpeg",
    width = 8, height = 10, units = "in",
    dpi = 450,
    limitsize = FALSE)
}

#### TRAIT SIGNALS ####
# examine if the species niches are linked to their traits with a gamma plot
postGamma = getPostEstimate(PA_model, parName = "Gamma")
plotGamma(PA_model, 
          post = postGamma, 
          trNamesNumbers = c(TRUE, FALSE),
          covNamesNumbers = c(TRUE, FALSE),
          colors = colorRampPalette(c("darkred","white","darkblue")),
          supportLevel = 0.95)

# save the figure
jpeg(filename = here("Figures", "PA_Model", "Gamma_Plot.jpg"), 
     width = 8, 
     height = 5, 
     units = "in", 
     res = 450)

par(mar = c(9, 11, 0, 0),
    mgp = c(0, 0, 0))

plotGamma(PA_model, 
          post = postGamma, 
          trNamesNumbers = c(TRUE, FALSE),
          covNamesNumbers = c(TRUE, FALSE),
          colors = colorRampPalette(c("darkred","white","darkblue")),
          supportLevel = 0.95,
          mar = c(9, 11, 0, 0))

dev.off()

# another way of examining the influence of traits is to see how much of the variation
# they explain among the responses of the species to the environmental covariates:
vp_ungrouped = computeVariancePartitioning(PA_model)
vp_ungrouped$R2T$Beta
round(((vp_ungrouped$R2T$Beta)*100), 2) # as a percentage, rather than proportion

# these results are consistent with the above figures: the traits explain only a very 
# minor part of the variation. the same negative result is obtained also in the sense
# that traits explain only a negligible proportion of variation in species' occurrence:
vp_ungrouped$R2T$Y
round(((vp_ungrouped$R2T$Y)*100), 2) # as a percentage, rather than proportion

#### PHYLOGENETIC SIGNAL ####
# next evaluate the posterior distribution of the phylogenetic signal in species niches.
# a rho of 0 would mean phylogeny explains nothing (species niches are independent of
# evolutionary history), while rho of 1 would mean phylogeny perfectly predicts niches 
# (all variation is explained by evolutionary relationships). 
summary(mpost$Rho)

# visualize the posterior distribution of Rho
# extract posterior samples of Rho
rho_samples_matrix = as.matrix(mpost$Rho)

# flatten the matrix into a single vector
rho_samples_flat = as.vector(rho_samples_matrix)

# plot the density of the posterior samples
ggplot(data.frame(Rho = rho_samples_flat), aes(x = Rho)) +
  geom_density(fill = "darkblue", alpha = 0.5) +
  labs(title = "Posterior Distribution of Phylogenetic Signal (Rho)",
       x = "Rho (Phylogenetic Signal)",
       y = "Density") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"))

# these results indicate a moderate phylogenetic signal - closely related
# species tend to respond similarly to environmental conditions. about 
# 51% of the variation in species' niches can be attributed to shared evolutionary 
# history (BUT NOTE that this assumes we are not missing any important environmental
# covariates or functional traits).

#### RESIDUAL SPECIES ASSOCIATIONS ####
##### co-occurrence matrix for only significant pairwise associations (+/-1) ##### 
require(corrplot)
OmegaCor = computeAssociations(PA_model)
supportLevel = 0.95

# build sign matrix (1 = positive assoc, -1 = negative assoc, 0 = no strong support)
toPlot = ((OmegaCor[[1]]$support > supportLevel) |
            (OmegaCor[[1]]$support < (1 - supportLevel))) * OmegaCor[[1]]$mean
toPlot = sign(toPlot)
diag(toPlot) = 0 # remove the diagonal (species self pairs automatically have +1 association)

# species names
sp_names = rownames(OmegaCor[[1]]$mean)
if (is.null(sp_names)) sp_names = as.character(seq_len(nrow(OmegaCor[[1]]$mean)))

# find species with at least one non-zero association
keep_species = sp_names[rowSums(abs(toPlot)) > 0]

if (length(keep_species) == 0) {
  stop("No species have >= 95% posterior support for positive OR negative residual correlations.")
} else if (length(keep_species) == 1) {
  stop("Only one species has significant associations — nothing to plot as a matrix.")
}

# preserve ordering
plotOrder = corrMatOrder(OmegaCor[[1]]$mean, order = "AOE")
plotOrder_names = sp_names[plotOrder]
keep_order_names = plotOrder_names[plotOrder_names %in% keep_species]

# subset matrix
toPlot_sub = toPlot[keep_order_names, keep_order_names]

# plot
corrplot(toPlot_sub, 
         method = "color",
         type = "lower",
         font = 3,
         outline = TRUE,
         addgrid.col = "gray40",
         tl.cex = 0.625,
         tl.col = "black",
         tl.srt = 45,
         cl.cex = 1.5,
         col = c("darkred", "white", "darkblue"),
         is.corr = FALSE,
         cl.length = 3)

# save the figure
jpeg(filename = here("Figures", "PA_Model", "Omega_Plot_Significant_Only.jpg"), 
     width = 10, height = 10, units = "in", res = 450)
par(mar = c(0, 0, 0, 0), xpd = TRUE)
corrplot(toPlot_sub, 
         method = "color",
         type = "lower",
         font = 3,
         outline = TRUE,
         addgrid.col = "gray40",
         tl.cex = 0.625,
         tl.col = "black",
         tl.srt = 45,
         cl.cex = 1.5,
         col = c("darkred", "white", "darkblue"),
         is.corr = FALSE,
         cl.length = 3)
dev.off()

gc()

##### full co-occurrence matrix #####
# this is the correlation matrix for the full community (all omegas), including 
# non-significant pairwise associations
OmegaCor = computeAssociations(PA_model)
supportLevel = 0.95

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot = sign(toPlot)
plotOrder = corrMatOrder(OmegaCor[[1]]$mean, order = "AOE")

corrplot(toPlot[plotOrder, plotOrder], 
         method = "color",
         type = "lower",
         font = 3,
         outline = TRUE,
         addgrid.col = "gray40",
         tl.cex = 0.625,
         tl.col = "black",
         tl.srt = 45,
         cl.cex = 1.5,
         col = c("darkred", "white", "darkblue"),
         is.corr = FALSE,
         cl.length = 3)

# save the figure
jpeg(filename = here("Figures", "PA_Model", "Omega_Plot_Full.jpg"), 
     width = 16, 
     height = 16, 
     units = "in", 
     res = 450)
par(mar = c(0, 0, 0, 0),
    xpd = TRUE)
corrplot(toPlot[plotOrder, plotOrder], 
         method = "color",
         type = "lower",
         font = 3,
         outline = TRUE,
         addgrid.col = "gray40",
         tl.cex = 0.625,
         tl.col = "black",
         tl.srt = 45,
         cl.cex = 1.5,
         col = c("darkred", "white", "darkblue"),
         is.corr = FALSE,
         cl.length = 3)
dev.off()

gc()

# identify those species with the greatest positive and negative associations
# restore the full species names for these calculations
colnames(PA_model$Y) = full_species
PA_model$spNames = full_species

# get the residual species co-occurrence matrix
postOmega = getPostEstimate(PA_model, parName = "Omega")
support = postOmega$support

# use the upper triangle only (exclude diagonal self-pairs)
upper = upper.tri(support)

# find positive associations with strong support (posterior probability >= 95%)
idx_pos = which(support >= 0.95 & upper, arr.ind = TRUE)

# associations with <=0.05 posterior probability of a positive association
# can be interpreted as having >=0.95 support for a negative association
idx_neg = which(support <= 0.05 & upper, arr.ind = TRUE)

# identify the positive pairs
positive_pairs = data.frame(
  Species_1 = colnames(PA_model$Y)[idx_pos[, 1]],
  Species_2 = colnames(PA_model$Y)[idx_pos[, 2]],
  Support   = support[idx_pos])

# quick double check to confirm there are no self pairs (A-B B-A)
positive_pairs = positive_pairs %>%
  mutate(sp1 = pmin(Species_1, Species_2),
         sp2 = pmax(Species_1, Species_2)) %>%
  distinct(sp1, sp2, .keep_all = TRUE) %>%
  select(-sp1, -sp2)

# identify the negative pairs
negative_pairs = data.frame(
  Species_1 = colnames(PA_model$Y)[idx_neg[, 1]],
  Species_2 = colnames(PA_model$Y)[idx_neg[, 2]],
  Support   = support[idx_neg])

# quick double check to confirm there are no self pairs (A-B B-A)
negative_pairs = negative_pairs %>%
  mutate(sp1 = pmin(Species_1, Species_2),
         sp2 = pmax(Species_1, Species_2)) %>%
  distinct(sp1, sp2, .keep_all = TRUE) %>%
  select(-sp1, -sp2)

# count how many strong residual positive associations there are
n_positive = sum(support[upper] >= 0.95)

#  count how many strong residual negative associations there are
n_negative = sum(support[upper] <= 0.05)

# who had the most positive and negative associations?
positive_counts = sort(
  table(c(positive_pairs$Species_1, positive_pairs$Species_2)),
  decreasing = TRUE) |>
  as.data.frame()
colnames(positive_counts) = c("Species", "N_positive_associations")

# who had the most positive and negative associations?
negative_counts = sort(
  table(c(negative_pairs$Species_1, negative_pairs$Species_2)),
  decreasing = TRUE) |>
  as.data.frame()
colnames(negative_counts) = c("Species", "N_negative_associations")

##### butterflyfish interactions #####
# identify positive residual species-to-species associations involving Chaetodon 
# butterflyfishes keeping only unique pairs (no A-B B-A duplicates)
positive_chaetodon_pairs = positive_pairs %>% 
  filter(str_detect(Species_1, "Chaetodon\\.") | 
           str_detect(Species_2, "Chaetodon\\.")) %>%
  rowwise() %>%
  mutate(sp1 = min(c(Species_1, Species_2)),
         sp2 = max(c(Species_1, Species_2))) %>%
  ungroup() %>%
  distinct(sp1, sp2, .keep_all = TRUE) %>%
  select(-sp1, -sp2)

# identify negative residual species-to-species associations involving Chaetodon
# butterflyfishes keeping only unique pairs (no A-B B-A duplicates)
negative_chaetodon_pairs = negative_pairs %>% 
  filter(str_detect(Species_1, "Chaetodon\\.") | 
           str_detect(Species_2, "Chaetodon\\.")) %>%
  rowwise() %>%
  mutate(sp1 = min(c(Species_1, Species_2)),
         sp2 = max(c(Species_1, Species_2))) %>%
  ungroup() %>%
  distinct(sp1, sp2, .keep_all = TRUE) %>%
  select(-sp1, -sp2)

# what percentage of ALL positive and negative residual species-to-species associations
# involved butterflyfish?
100 * (nrow(positive_chaetodon_pairs) / nrow(positive_pairs))
100 * (nrow(negative_chaetodon_pairs) / nrow(negative_pairs))

# how many positive pairs were Chaetodon-Chaetodon?
positive_chaetodon_chaetodon = positive_chaetodon_pairs %>%
  filter(str_detect(Species_1, "Chaetodon\\.") &
           str_detect(Species_2, "Chaetodon\\.")) %>%
  nrow()

# how many negative pairs were Chaetodon-Chaetodon?
negative_chaetodon_chaetodon = negative_chaetodon_pairs %>%
  filter(str_detect(Species_1, "Chaetodon\\.") &
           str_detect(Species_2, "Chaetodon\\.")) %>%
  nrow()

# do chaetodontids co-occur with acanthurids and pomacentrids more (positive) or 
# less (negative) than expected based on their responses to environmental covariates?
acanthuridae = c("Acanthurus", "Ctenochaetus", "Naso")
pomacentridae = c("Chromis", "Dascyllus", "Pycnochromis", "Stegastes")
target_genera = c(acanthuridae, pomacentridae)

# positive residual associations
positive_acanthurid_pomacentrid = positive_chaetodon_pairs %>%
  filter(
    (str_detect(Species_1, "Chaetodon\\.") & 
       str_detect(Species_2, paste0("^(", paste(target_genera, collapse="|"), ")\\."))) |
      (str_detect(Species_2, "Chaetodon\\.") & 
         str_detect(Species_1, paste0("^(", paste(target_genera, collapse="|"), ")\\."))))
nrow(positive_acanthurid_pomacentrid)

# negative residual associations
negative_acanthurid_pomacentrid = negative_chaetodon_pairs %>%
  filter(
    (str_detect(Species_1, "Chaetodon\\.") & 
       str_detect(Species_2, paste0("^(", paste(target_genera, collapse="|"), ")\\."))) |
      (str_detect(Species_2, "Chaetodon\\.") & 
         str_detect(Species_1, paste0("^(", paste(target_genera, collapse="|"), ")\\."))))
nrow(negative_acanthurid_pomacentrid)
