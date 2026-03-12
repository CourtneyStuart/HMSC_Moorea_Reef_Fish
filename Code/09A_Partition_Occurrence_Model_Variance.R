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
#setwd("E:/Data/StuartC_DPhil_Ch3/")
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

species_labels = make_unique_abbrev(PA_model$spNames)

#### FIXED VS. RANDOM SIMPLE ####
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
random_sp = FR_long %>%
  filter(Component == "Random effects" &
           Variance >= (1/3))

# choose colors for plotting
cols = c("blue4", "darkorange")  # blue = fixed, orange = random

# set the file name, save location, and settings
jpeg(filename = here("Figures", "PA_Model", "Variance_Partitioning_Fixed_vs_Random.jpg"),
     width = 14,
     height = 7,
     units = "in",
     res = 450)

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
# check model structure
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

# compute variance partitioning
vp_grouped = computeVariancePartitioning(PA_model, 
                                 group = group, 
                                 groupnames = group_names)

# examine results
vp_grouped$vals  # variance proportion for each group and species
vp_grouped$R2T   # R2T gives the variance among species explained by traits, measured for 
# species' responses to covariates ($R2T$Beta) and species occurrences ($R2T$Y)

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
     res = 450)

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
     cex = 0.65,
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
  res = 450)

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
  cex.axis = 0.475)

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
vp_grouped_transposed = vp_grouped_vals %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Species") %>%
  mutate(across(-Species, ~ round(.x, 2)))

write.csv(vp_grouped_transposed,
          here("HMSC", "Data", "Variance_Partitioning_Grouped.csv"),
          row.names = FALSE)

#### UNGROUPED FIXED EFFECTS ####
# variance partitioning without groupings 
vp_ungrouped = computeVariancePartitioning(PA_model)

# examine results
vp_ungrouped$vals # variance proportion for each variable and species
vp_ungrouped$R2T # R2T gives the variance among species explained by traits, measured for 
# species' responses to covariates ($R2T$Beta) and species occurrences ($R2T$Y)

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
     res = 450)

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
     cex = 0.625,
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
  res = 450)

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
  cex.axis = 0.5)

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