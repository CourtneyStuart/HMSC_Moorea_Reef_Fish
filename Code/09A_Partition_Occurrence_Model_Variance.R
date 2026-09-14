#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### LIBRARIES ####
# install required packages (first run only)
# install.packages(c("easypackages", "conflicted", "tidyr", "dplyr", "here", "ggplot2",
#                    "Hmsc", "coda", "corrplot", "Cairo", "stringr", "PNWColors",
#                    "colorspace", "grDevices", "colorRamps", "tibble", "ape"))

# load packages
library(easypackages)
libraries("conflicted", "tidyr", "dplyr", "here", "ggplot2",
          "Hmsc", "coda", "corrplot", "Cairo", "stringr", 
          "PNWColors", "colorspace", "grDevices", "colorRamps", "tibble", "ape")

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
samples = 4000
thin = 50
filename = file.path(model.directory, 
                     paste0("PA_model_chains_", as.character(nChains),
                            "_total_samples_",as.character(samples),
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
  len = rep(1L, n) # start with 1-letter abbrev
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

species_labels = make_unique_abbrev(PA_model$spNames)

#### FIXED VS. RANDOM ####
# group all fixed effects together
group_fixed = rep(1, ncol(PA_model$X))

# partition the variance explained between the fixed vs. random effects
vp_simple = computeVariancePartitioning(
  PA_model,
  group = group_fixed,
  groupnames = "Fixed effects")

# save the values
vp_simple_vals = vp_simple$vals

# extract the fixed and random effects
fixed  = vp_simple_vals["Fixed effects", ]
random = colSums(vp_simple_vals[c("Random: site", "Random: year"), ])

# combine into a matrix where rows = components and columns = species
FR_vals = rbind(
  "Fixed effects" = fixed,
  "Random effects" = random)

# calculate mean variance explained for each component
mean_fixed = round(mean(fixed), digits = 2)
round(summary(fixed), digits = 2)
mean_random = round(mean(random), digits =  2)
round(summary(random), digits = 2)

# calculate summary stats over all species
# convert to data frame
FR_df = as.data.frame(FR_vals)

# add component column
FR_df$Component = rownames(FR_df)

# convert to long format
FR_long = reshape(
  FR_df,
  varying = colnames(FR_df)[colnames(FR_df) != "Component"],
  v.names = "Variance",
  timevar = "Species",
  times = colnames(FR_df)[colnames(FR_df) != "Component"],
  direction = "long")

# clean up
rownames(FR_long) = NULL
FR_long = FR_long[, c("Species", "Component", "Variance")]

# how many species had greater than one-third of variance captured by random effects?
random_sp_third = FR_long %>%
  filter(Component == "Random effects" &
           Variance >= (1/3))

# how many species had greater than one-half of variance captured by random effects?
random_sp_half = FR_long %>%
  filter(Component == "Random effects" &
           Variance >= (1/2))

# choose colors for plotting
cols = c("blue4", "darkorange")  # blue = fixed, orange = random

# set the file name, save location, and settings
jpeg(filename = here("Figures", "PA_Model",
                     "Variance_Partitioning_Fixed_vs_Random.jpg"),
     width = 14,
     height = 7,
     units = "in",
     res = 450,
     quality = 100)

par(mar = c(5.5, 3.5, 0.0, 0.0) + 0.1,
    mgp = c(2.15, 0.5, 0),
    xpd = TRUE)

# create the bar plot
bp = barplot(
  FR_vals,
  las = 1,
  cex.axis = 1,
  names.arg = rep("", ncol(FR_vals)),
  col = cols,
  ylab = "Proportion of variance in occurrence",
  border = "black",
  legend = FALSE,
  ylim = c(0, max(colSums(FR_vals)) * 1.1),
  space = c(0, rep(0, ncol(FR_vals) - 1)))

# add species names
text(
  x = bp,
  y = -0.01,
  labels = species_labels, #use colnames(FR_vals) if you want full names instead
  srt = 90,
  adj = 1,
  xpd = TRUE,
  cex = 0.65,
  font = 3)

# create legend labels with mean variance (2 decimal places)
legend_labels = c(
  sprintf("Fixed effects (mean = %.2f)", mean_fixed),
  sprintf("Random effects (mean = %.2f)", mean_random))

# legend centered at top
legend(
  x = "top",
  legend = rev(legend_labels),
  fill = rev(cols),
  cex = 1,
  xpd = TRUE,
  bty = "n",
  horiz = TRUE,
  inset = c(0, 0.015))

dev.off()

# calculate summary stats over all species
# mean variance explained across species
mean_FR = aggregate(
  Variance ~ Component,
  data = FR_long,
  FUN = mean)
print(mean_FR)

# max variance explained across species
max_FR = aggregate(
  Variance ~ Component,
  data = FR_long,
  FUN = max)
print(max_FR)

# min variance explained across species
min_FR = aggregate(
  Variance ~ Component,
  data = FR_long,
  FUN = min)
print(min_FR)

# show both proportion (0 to 1) and percentage (0 to 100)
mean_FR$Percent = 100 * mean_FR$Variance
max_FR$Percent  = 100 * max_FR$Variance
min_FR$Percent  = 100 * min_FR$Variance
print(mean_FR)
print(max_FR)
print(min_FR)

#### GROUPED FIXED EFFECTS ####
# check model structure before grouping
head(PA_model$X)  # fixed effects design matrix
names(PA_model$rL)  # random effects: should be "site" and "year"

# set groups
group = c(
  0,     # (Intercept) - typically coded as 0
  1, 1,  # HabitatForereef, HabitatFringing (Habitat type)
  2,     # COTS abundance
  3,     # Max_DHW
  4, 4, 4,  # Cyclone1, Cyclone2, Cyclone3 (Cyclone exposure)
  5,     # Land_Dist 
  6, 6,  # Depth_Mean_100m, Depth_Mean_500m (Depth)
  7, 7,  # Curvature_Mean_100m, Curvature_Mean_500m (Curvature)
  8, 8, 8   # Coral_Mean_Cover, Macroalgae_Mean_Cover, CTB_Mean_Cover (Benthic cover)
)

# set group names
group_names = c("Habitat type", "COTS abundance", "Maximum DHW", "Cyclone exposure", 
                "Distance to land", "Depth", "Curvature", "Benthic cover")

# verify the length matches
length(group)  # should be 16
ncol(PA_model$X)  # should be 16

# compute variance partitioning using grouped covariates
vp_grouped = computeVariancePartitioning(PA_model, 
                                 group = group, 
                                 groupnames = group_names)

# examine results
vp_grouped$vals  # variance proportion for each group and species
vp_grouped$R2T   # R2T gives the variance among species explained by traits, measured for 
# species' responses to covariates ($R2T$Beta) and species occurrences ($R2T$Y)
round(mean(vp_grouped$R2T$Beta), digits = 2)
round(mean(vp_grouped$R2T$Y), digits = 2)

# plot variance partitioning
par(mfrow = c(1, 1))

# extract variance partitioning values
vp_grouped_vals = vp_grouped$vals

# calculate mean variance explained by each component
vp_grouped_means = rowMeans(vp_grouped_vals) # convert to percentage

# create legend labels with means
legend_labels = paste0(rownames(vp_grouped_vals), 
                       " (mean = ", round(vp_grouped_means, 2), ")")

##### landscape orientation #####
# create the plot
jpeg(filename = here("Figures", "PA_Model", "Variance_Partitioning_Grouped.jpg"), 
     width = 16, 
     height = 8, 
     units = "in", 
     res = 450,
     quality = 100)

par(mar = c(5.5, 3.5, 0.25, 9.25) + 0.1,   # increase bottom and right margins
    mgp = c(2.15, 0.5, 0),  # moves axis title and labels closer
    xpd = TRUE)  # allow plotting outside plot region

bp = barplot(vp_grouped_vals, 
             las = 1,
             cex.axis = 0.8,
             names.arg = rep("", ncol(vp_grouped_vals)),
             col = colorRamps::matlab.like(nrow(vp_grouped_vals)),
             ylab = "Proportion of variance in occurrence",
             border = "black",
             legend = FALSE,
             ylim = c(0, max(colSums(vp_grouped_vals)) * 1.02),
             space = c(0, rep(0, ncol(vp_grouped_vals) - 1)))

text(x = bp, 
     y = -0.01,
     labels = species_labels, #use colnames(vp_grouped_vals) if you want full names instead
     srt = 90,
     adj = 1, 
     xpd = TRUE,
     cex = 0.6,
     font = 3)

# add legend outside plot area (to the top-right)
legend(x = par("usr")[2] - 6,
       y = par("usr")[4],
       legend = rev(legend_labels),
       fill = rev(colorRamps::matlab.like(nrow(vp_grouped_vals))),
       cex = 0.8,
       xpd = TRUE,
       bty = "n")
dev.off()

##### portrait orientation #####
jpeg(
  filename = here("Figures", "PA_Model", "Variance_Partitioning_Grouped_Portrait.jpg"),
  width = 9,
  height = 13,
  units = "in",
  res = 450,
  quality = 100)

par(
  mar = c(5.5, 5.5, 0, 0.5),
  mgp = c(-1.75, 0.75, 0),
  xpd = NA)

bp = barplot(
  vp_grouped_vals,
  horiz = TRUE,
  las = 1,
  cex.axis = 0.8,
  cex.lab = 0.8,
  names.arg = rep("", ncol(vp_grouped_vals)),
  col = colorRamps::matlab.like(nrow(vp_grouped_vals)),
  xlab = "Proportion of variance in occurrence",
  border = "black",
  legend = FALSE,
  xlim = c(0, max(colSums(vp_grouped_vals)) * 1.02),
  space = 0)

# add italicized species names on y-axis
axis(
  side = 2,
  at = bp,
  labels = species_labels, #use colnames(vp_ungrouped_vals) if you want full names instead
  las = 2,
  font = 3,
  cex.axis = 0.4)

# # bottom legend with extra separation
legend(
  "bottom",
  legend = legend_labels,
  fill = colorRamps::matlab.like(nrow(vp_grouped_vals)),
  ncol = 3,
  bty = "n",
  cex = 0.7,
  inset = c(0, -0.0875))

dev.off()

# save the grouped variance partitioning results as a csv file
require(tibble)
vp_grouped_transposed = vp_grouped_vals %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Species") %>%
  mutate(across(-Species, ~ round(.x, 2)))

write.csv(vp_grouped_transposed,
          here("HMSC", "Data", "Variance_Partitioning_Grouped.csv"),
          row.names = FALSE)

#### UNGROUPED FIXED EFFECTS ####
# variance partitioning without grouping fixed-effect covariates
# the intercept is assigned to the first group because it has no variation
group = c(1, 1:(ncol(PA_model$X) - 1))
groupnames = colnames(PA_model$X)[-1]

vp_ungrouped = computeVariancePartitioning(
  PA_model,
  group = group,
  groupnames = groupnames)

# examine results
vp_ungrouped$vals # variance proportion for each variable and species

# check that the group vector has one value for every column of X
length(group) # should be 16
ncol(PA_model$X) # should be 16

# check group assignments
group
groupnames

# these will be the same as for vp_grouped above
# vp_ungrouped$R2T
# round(mean(vp_grouped$R2T$Beta), digits = 2) 
# round(mean(vp_grouped$R2T$Y), digits = 2)

plotVariancePartitioning(PA_model, VP = vp_ungrouped)

# extract variance partitioning values
vp_ungrouped_vals = vp_ungrouped$vals

# save as a data frame
vp_ungrouped_vals_df = as.data.frame(vp_ungrouped_vals)

# find species where any single variable explains >= 50% of variance
high_influence_species = data.frame(
  species = colnames(vp_ungrouped_vals_df),
  max_var_explained = apply(vp_ungrouped_vals_df, 2, max),
  top_variable = rownames(vp_ungrouped_vals_df)[apply(vp_ungrouped_vals_df, 2, which.max)]) %>%
  filter(max_var_explained >= 0.5) %>%
  arrange(desc(max_var_explained))

# count how many species have HabitatForereef as the top variable
sum(high_influence_species$top_variable == "HabitatForereef")

# count how many species have depth_mean_100m as the top variable
sum(high_influence_species$top_variable == "Depth_Mean_100m")

# calculate mean variance explained by each component
vp_ungrouped_means = rowMeans(vp_ungrouped_vals)

# create legend labels with means
legend_labels = paste0(rownames(vp_ungrouped_vals),
                       " (mean = ", round(vp_ungrouped_means, 2), ")")

##### landscape orientation #####
# create the plot
jpeg(filename = here("Figures", "PA_Model", "Variance_Partitioning_Not_Grouped.jpg"), 
     width = 16, 
     height = 8, 
     units = "in", 
     res = 450,
     quality = 100)

par(mar = c(5.5, 3, 0.25, 11.25) + 0.1,
    mgp = c(2, 0.5, 0),
    xpd = TRUE)

bp = barplot(vp_ungrouped_vals, 
             las = 1,
             cex.axis = 0.8,
             names.arg = rep("", ncol(vp_ungrouped_vals)),
             col = colorRamps::matlab.like(nrow(vp_ungrouped_vals)),
             ylab = "Proportion of variance in occurrence",
             border = "black",
             legend = FALSE,
             ylim = c(0, max(colSums(vp_ungrouped_vals)) * 1.02),
             space = c(-0.5, rep(0, ncol(vp_ungrouped_vals) - 1)))

text(x = bp, 
     y = -0.01,
     labels = species_labels, #use colnames(vp_ungrouped_vals) if you want full names instead
     srt = 90,
     adj = 1,
     xpd = TRUE,
     cex = 0.6,
     font = 3)

legend(x = par("usr")[2] - 6,
       y = par("usr")[4],
       legend = rev(legend_labels),
       fill = rev(colorRamps::matlab.like(nrow(vp_ungrouped_vals))),
       cex = 0.8,
       xpd = TRUE,
       bty = "n")

dev.off()

##### portrait orientation #####
par(mfrow = c(1, 1))

# calculate mean variance explained by each component
vp_means = rowMeans(vp_ungrouped_vals)

legend_labels = paste0(
  rownames(vp_ungrouped_vals), " (mean = ", round(vp_ungrouped_means, 2), ")")

jpeg(
  filename = here("Figures", "PA_Model", "Variance_Partitioning_Not_Grouped_Portrait.jpg"),
  width = 8,
  height = 14,
  units = "in",
  res = 450,
  quality = 100)

par(
  mar = c(9.75, 5.5, 0, 0.5),
  mgp = c(-1.5, 0.75, 0),
  xpd = NA)

bp = barplot(
  vp_ungrouped_vals,
  horiz = TRUE,
  las = 1,
  cex.axis = 0.8,
  names.arg = rep("", ncol(vp_ungrouped_vals)),
  col = colorRamps::matlab.like(nrow(vp_ungrouped_vals)),
  xlab = "Proportion of variance in occurrence",
  border = "black",
  legend = FALSE,
  xlim = c(0, max(colSums(vp_ungrouped_vals)) * 1.02),
  space = 0)

# add italicized species names on y-axis
axis(
  side = 2,
  at = bp,
  labels = species_labels, #use colnames(vp_ungrouped_vals) if you want full names instead
  las = 2,
  font = 3,
  cex.axis = 0.4)

# # bottom legend with extra separation
legend(
  "bottom",
  legend = legend_labels,
  fill = colorRamps::matlab.like(nrow(vp_ungrouped_vals)),
  ncol = 2,
  bty = "n",
  cex = 0.8,
  inset = c(0, -0.16))

dev.off()

# save the ungrouped variance partitioning results as a csv file
vp_ungrouped_transposed = vp_ungrouped_vals %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Species") %>%
  mutate(across(-Species, ~ round(.x, 2)))

write.csv(vp_ungrouped_transposed,
          here("HMSC", "Data", "Variance_Partitioning_Ungrouped.csv"),
          row.names = FALSE)

#### GROUPED - PHYLOGENETIC ORDERING ####
# extract the phylogenetic tree directly from the Hmsc model
phylo_tree = PA_model$phyloTree

# create a copy of the tree for plotting
phylo_tree_plot = phylo_tree

# create abbreviated species labels for plotting
phylo_tree_plot$tip.label = sapply(
  strsplit(phylo_tree$tip.label, "\\."),
  function(x) paste0(
    substr(x[1], 1, 1),
    ". ",
    paste(x[-1], collapse = " ")))

# extract the phylogenetic tip order
phylo_order = match(
  phylo_tree$tip.label,
  colnames(vp_grouped_vals))

# check that all species are present in the variance partitioning results
if (any(is.na(phylo_order))) {
  stop("some species in the phylogenetic tree are missing from vp_grouped_vals.")
}

# reorder variance partitioning values according to phylogenetic tip order
vp_grouped_vals_phylo = vp_grouped_vals[, phylo_order]

# calculate mean variance explained by each component
vp_grouped_means_phylo = rowMeans(vp_grouped_vals_phylo)

# create legend labels with means
legend_labels_phylo = paste0(
  rownames(vp_grouped_vals_phylo),
  " (mean = ",
  round(vp_grouped_means_phylo, 2),
  ")")

##### identify species with poor beta convergence #####
# species for which variance partitioning should not be interpreted
problem_species = c(
  "Neocirrhites.armatus",
  "Zebrasoma.scopas",
  "Centropyge.flavissima")

# identify their positions in the phylogenetically ordered results
problem_idx = phylo_tree$tip.label %in% problem_species

# check that all specified species occur in the tree
if (!all(problem_species %in% phylo_tree$tip.label)) {
  stop("one or more problem_species are missing from the phylogenetic tree.")
}

##### determine how much of the phylogeny to show #####
# calculate the depth of each node from the root
node_depth = ape::node.depth.edgelength(phylo_tree)

# find the maximum root-to-tip distance
max_tree_depth = max(node_depth)

# proportion of the tree to show
recent_tree_fraction = 0.10

# calculate the left-hand limit of the displayed tree
tree_xmin = max_tree_depth - (
  max_tree_depth * recent_tree_fraction)

##### create portrait figure #####
jpeg(
  filename = here(
    "Figures","PA_Model",
    "Variance_Partitioning_Grouped_Portrait_Phylogenetic_Recent.jpg"),
  width = 11,
  height = 14,
  units = "in",
  res = 600,
  quality = 100)

##### colors #####
vp_cols = colorRamps::matlab.like(
  nrow(vp_grouped_vals_phylo))

# color used to flag species with poor beta convergence
poor_convergence_col = "grey40"

##### main plot area #####
# tree, labels and variance plot use exactly the same vertical coordinates
main_bottom = 0.07
main_top = 1

##### phylogenetic tree #####
par(
  fig = c(0.00, 0.1, main_bottom, main_top),
  mar = c(2.2, 0, 0, 0),
  xpd = FALSE)

# plot the terminal portion of the phylogeny without tip labels
ape::plot.phylo(
  phylo_tree_plot,
  type = "phylogram",
  direction = "rightwards",
  show.tip.label = FALSE,
  edge.width = 0.75,
  x.lim = c(
    tree_xmin,
    max_tree_depth),
  y.lim = c(
    0.5,
    length(phylo_tree$tip.label) + 0.5),
  no.margin = FALSE)

# retrieve the coordinates used by plot.phylo
phylo_coords = get(
  "last_plot.phylo",
  envir = .PlotPhyloEnv)

# extract the y-position of each species tip
tip_y = phylo_coords$yy[
  seq_len(length(phylo_tree$tip.label))]

##### species labels #####
par(
  fig = c(0.1, 0.21, main_bottom, main_top),
  mar = c(2.2, 0, 0, 0),
  xpd = FALSE,
  new = TRUE)

# create empty plotting region for species labels
plot(
  NA,
  xlim = c(0, 1),
  ylim = c(
    0.5,
    length(phylo_tree$tip.label) + 0.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  bty = "n")

# add abbreviated species names
text(
  x = 0,
  y = tip_y,
  labels = phylo_tree_plot$tip.label,
  pos = 4,
  offset = 0,
  font = 3,
  cex = 0.6)

##### variance partitioning #####
par(
  fig = c(0.17, 1.00, main_bottom, main_top),
  mar = c(2.2, 0, 0, 0),
  mgp = c(1.7, 0.5, 0),
  xpd = FALSE,
  new = TRUE)

# set up plotting region
plot(
  NA,
  xlim = c(
    0,
    max(colSums(vp_grouped_vals_phylo)) * 1.02),
  ylim = c(
    0.5,
    length(phylo_tree$tip.label) + 0.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  bty = "n")

# height of each horizontal bar
bar_height = 1

# draw stacked variance partitioning bars
for (i in seq_len(ncol(vp_grouped_vals_phylo))) {
  
  # NEW: plot a single grey bar for species with poor beta convergence
  if (problem_idx[i]) {
    
    rect(
      xleft = 0,
      ybottom = tip_y[i] - bar_height / 2,
      xright = sum(vp_grouped_vals_phylo[, i]),
      ytop = tip_y[i] + bar_height / 2,
      col = poor_convergence_col,
      border = "black"
    )
    
  } else {
    
    # plot the usual stacked variance partition for all other species
    x_left = 0
    
    for (j in seq_len(nrow(vp_grouped_vals_phylo))) {
      
      x_right = x_left + vp_grouped_vals_phylo[j, i]
      
      rect(
        xleft = x_left,
        ybottom = tip_y[i] - bar_height / 2,
        xright = x_right,
        ytop = tip_y[i] + bar_height / 2,
        col = vp_cols[j],
        border = "black"
      )
      
      x_left = x_right
    }
  }
}

# add x-axis
axis(
  side = 1,
  cex.axis = 0.9,
  tck = -0.02)

# add x-axis title
mtext(
  "Proportion of variance in occurrence",
  side = 1,
  line = -1.75,
  cex = 0.9)

##### legend #####
par(
  fig = c(0.00, 1.00, 0.00, main_bottom),
  mar = c(0, 0, 0, 0),
  xpd = FALSE,
  new = TRUE)

# create empty plotting region for legend
plot.new()

# add poor-convergence category to legend
legend(
  "top",
  legend = c(
    legend_labels_phylo,
    "Poor β convergence – not interpreted"
  ),
  fill = c(
    vp_cols,
    "grey40"
  ),
  ncol = 3,
  bty = "n",
  cex = 0.9,
  x.intersp = 0.5,
  y.intersp = 0.9)

dev.off()

#### UNGROUPED - PHYLOGENETIC ORDERING ####
# extract the phylogenetic tree directly from the Hmsc model
phylo_tree = PA_model$phyloTree

# create a copy of the tree for plotting
phylo_tree_plot = phylo_tree

# create abbreviated species labels for plotting
phylo_tree_plot$tip.label = sapply(
  strsplit(phylo_tree$tip.label, "\\."),
  function(x) paste0(
    substr(x[1], 1, 1),
    ". ",
    paste(x[-1], collapse = " ")))

# extract the phylogenetic tip order
phylo_order = match(
  phylo_tree$tip.label,
  colnames(vp_ungrouped_vals))

# check that all species are present in the variance partitioning results
if (any(is.na(phylo_order))) {
  stop("some species in the phylogenetic tree are missing from vp_ungrouped_vals.")
}

# reorder variance partitioning values according to phylogenetic tip order
vp_ungrouped_vals_phylo = vp_ungrouped_vals[, phylo_order]

# calculate mean variance explained by each component
vp_ungrouped_means_phylo = rowMeans(vp_ungrouped_vals_phylo)

# create legend labels with means
legend_labels_phylo = paste0(
  rownames(vp_ungrouped_vals_phylo),
  " (mean = ",
  round(vp_ungrouped_means_phylo, 2),
  ")")

##### identify species with poor beta convergence #####
# species for which variance partitioning should not be interpreted
problem_species = c(
  "Neocirrhites.armatus",
  "Zebrasoma.scopas",
  "Centropyge.flavissima")

# identify their positions in the phylogenetically ordered results
problem_idx = phylo_tree$tip.label %in% problem_species

# check that all specified species occur in the tree
if (!all(problem_species %in% phylo_tree$tip.label)) {
  stop("one or more problem_species are missing from the phylogenetic tree.")
}

##### determine how much of the phylogeny to show #####
# calculate the depth of each node from the root
node_depth = ape::node.depth.edgelength(phylo_tree)

# find the maximum root-to-tip distance
max_tree_depth = max(node_depth)

# proportion of the tree to show
recent_tree_fraction = 0.10

# calculate the left-hand limit of the displayed tree
tree_xmin = max_tree_depth - (
  max_tree_depth * recent_tree_fraction)

##### create portrait figure #####
jpeg(
  filename = here(
    "Figures","PA_Model",
    "Variance_Partitioning_Ungrouped_Portrait_Phylogenetic_Recent.jpg"),
  width = 11,
  height = 14,
  units = "in",
  res = 600,
  quality = 100)

##### colors #####
vp_cols = colorRamps::matlab.like(
  nrow(vp_ungrouped_vals_phylo))

# color used to flag species with poor beta convergence
poor_convergence_col = "grey40"

##### main plot area #####
# tree, labels and variance plot use exactly the same vertical coordinates
main_bottom = 0.07
main_top = 1

##### phylogenetic tree #####
par(
  fig = c(0.00, 0.1, main_bottom, main_top),
  mar = c(2.2, 0, 0, 0),
  xpd = FALSE)

# plot the terminal portion of the phylogeny without tip labels
ape::plot.phylo(
  phylo_tree_plot,
  type = "phylogram",
  direction = "rightwards",
  show.tip.label = FALSE,
  edge.width = 0.75,
  x.lim = c(
    tree_xmin,
    max_tree_depth),
  y.lim = c(
    0.5,
    length(phylo_tree$tip.label) + 0.5),
  no.margin = FALSE)

# retrieve the coordinates used by plot.phylo
phylo_coords = get(
  "last_plot.phylo",
  envir = .PlotPhyloEnv)

# extract the y-position of each species tip
tip_y = phylo_coords$yy[
  seq_len(length(phylo_tree$tip.label))]

##### species labels #####
par(
  fig = c(0.1, 0.21, main_bottom, main_top),
  mar = c(2.2, 0, 0, 0),
  xpd = FALSE,
  new = TRUE)

# create empty plotting region for species labels
plot(
  NA,
  xlim = c(0, 1),
  ylim = c(
    0.5,
    length(phylo_tree$tip.label) + 0.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  bty = "n")

# add abbreviated species names
text(
  x = 0,
  y = tip_y,
  labels = phylo_tree_plot$tip.label,
  pos = 4,
  offset = 0,
  font = 3,
  cex = 0.6)

##### variance partitioning #####
par(
  fig = c(0.17, 1.00, main_bottom, main_top),
  mar = c(2.2, 0, 0, 0),
  mgp = c(1.7, 0.5, 0),
  xpd = FALSE,
  new = TRUE)

# set up plotting region
plot(
  NA,
  xlim = c(
    0,
    max(colSums(vp_ungrouped_vals_phylo)) * 1.02),
  ylim = c(
    0.5,
    length(phylo_tree$tip.label) + 0.5),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  bty = "n")

# height of each horizontal bar
bar_height = 1

# draw stacked variance partitioning bars
for (i in seq_len(ncol(vp_ungrouped_vals_phylo))) {
  
  # plot a single grey bar for species with poor beta convergence
  if (problem_idx[i]) {
    
    rect(
      xleft = 0,
      ybottom = tip_y[i] - bar_height / 2,
      xright = sum(vp_ungrouped_vals_phylo[, i]),
      ytop = tip_y[i] + bar_height / 2,
      col = poor_convergence_col,
      border = "black"
    )
    
  } else {
    
    # plot the usual stacked variance partition for all other species
    x_left = 0
    
    for (j in seq_len(nrow(vp_ungrouped_vals_phylo))) {
      
      x_right = x_left + vp_ungrouped_vals_phylo[j, i]
      
      rect(
        xleft = x_left,
        ybottom = tip_y[i] - bar_height / 2,
        xright = x_right,
        ytop = tip_y[i] + bar_height / 2,
        col = vp_cols[j],
        border = "black"
      )
      
      x_left = x_right
    }
  }
}

# add x-axis
axis(
  side = 1,
  cex.axis = 0.8,
  tck = -0.02)

# add x-axis title
mtext(
  "Proportion of variance in occurrence",
  side = 1,
  line = -1.75,
  cex = 0.8)

##### legend #####
par(
  fig = c(0.00, 1.00, 0.00, main_bottom),
  mar = c(0, 0, 0, 0),
  xpd = FALSE,
  new = TRUE)

# create empty plotting region for legend
plot.new()

# add poor-convergence category to legend
legend(
  "top",
  legend = c(
    legend_labels_phylo,
    "Poor β convergence – not interpreted"
  ),
  fill = c(
    vp_cols,
    "grey40"
  ),
  ncol = 3,
  bty = "n",
  cex = 0.8,
  x.intersp = 0.5,
  y.intersp = 0.9)

dev.off()
