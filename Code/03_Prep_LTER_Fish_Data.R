#### CONTACT ####
# Courtney Stuart (courtney.seascape@gmail.com)

#### LIBRARIES ####
# install packages (first run only)
# install.packages(c("easypackages", "conflicted", "dplyr", "here", "ggplot2", 
#                    "ape", "data.tree", "tidyverse", "rfishbase", "stringdist", 
#                    "purrr", "tictoc"))

# load packages
library(easypackages)
libraries("conflicted", "dplyr", "here", "ggplot2", "ape", "data.tree", 
          "tidyverse", "rfishbase", "stringdist", "purrr", "tictoc")

# resolve package conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

#### DIRECTORIES ####
# working directory and relative folder path for here()
#setwd("E:/Data/StuartC_DPhil_Ch3/")
#set_here("E:/Data/StuartC_DPhil_Ch3/") # set first-time only
here::i_am(".here")
here::here() # verify

#### FULL TRANSECT SURVEY DATA #### 
# read in the full Moorea LTER annual fish data downloaded from 
# https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-mcr&identifier=6&revision=64
data = read.csv(here("Data", "Fish", "MCR_LTER_Annual_Fish_Survey_20250324.csv"))

# keep only observations of bony, ray-finned fishes (i.e., drop the cartilaginous fishes)
# drop transects where there were no fish observed at all
# drop observations where only family could be identified (not genus or species)
data = data %>%
  filter(! Family %in% c("Carcharhinidae", "Ginglymostomatidae", 
                         "Myliobatidae", "Dasyatidae")) %>%
  filter(! Taxonomy == "No fish observed") %>%
  filter(!str_detect(Taxonomy, "unidentified"))

# remove any observations where the species binomial name is NULL (if there are any) or
# where multiple unidentified species (spp.) were listed together (if there are any)
data = data %>% 
  filter(!Taxonomy == "NULL") %>%
  filter(!str_detect(Taxonomy, "spp\\."))

# there may be multiple sp. per genus, some of which may already be numbered 
# (e.g., Acanthurus sp. 1) --> extract genus + any existing number from 
# "Genus sp" or "Genus sp. X"
data = data %>%
  mutate(
    match = str_match(Taxonomy, "^([A-Z][a-z]+) (sp\\.?)(?: ([0-9]+))?$"),
    genus_extracted = match[,2],
    sp_marker = match[,3],
    existing_num = match[,4])

# get the current max assigned number for each genus (if any)
existing_nums = data %>%
  filter(!is.na(existing_num)) %>%
  mutate(existing_num = as.integer(existing_num)) %>%
  group_by(genus_extracted) %>%
  summarise(max_num = max(existing_num), .groups = "drop")

# assign new numbers to unnumbered "Genus sp" in increasing numeric order
unnumbered_sp = data %>%
  mutate(row_id = row_number()) %>%
  filter(!is.na(sp_marker) & is.na(existing_num)) %>%
  group_by(genus_extracted) %>%
  mutate(new_num = row_number()) %>%
  ungroup() %>%
  left_join(existing_nums, by = "genus_extracted") %>%
  mutate(
    final_num = new_num + coalesce(max_num, 0),
    Taxonomy_clean = paste0(genus_extracted, " sp. ", final_num)) %>%
  select(row_id, Taxonomy_clean)

# merge these records back into the full dataset and add a final Taxonomy_clean column
data = data %>%
  mutate(row_id = row_number()) %>%
  left_join(unnumbered_sp, by = "row_id") %>%
  mutate(
    Taxonomy_clean = coalesce(
      Taxonomy_clean,  # from unnumbered_sp if applicable
      ifelse(
        !is.na(sp_marker) & !is.na(existing_num),
        paste0(genus_extracted, " sp. ", existing_num),
        Taxonomy))) %>%
  select(-match, -sp_marker, -existing_num, -row_id, -genus_extracted)

# (cf) likely means confer or compare to, suggesting a confident but uncertain ID -> 
# drop (cf) and assume that the IDs are correct because (cf) will cause problems when
# later constructing a taxonomic tree and extracting traits
data = data %>%
  mutate(
    # add "1" to any (cf) sp. without a number
    Taxonomy_clean = ifelse(
      str_detect(Taxonomy_clean, "\\(cf\\) sp\\.(?!\\s*\\d)"),
      str_replace(Taxonomy_clean, "\\(cf\\) sp\\.", "(cf) sp. 1"),
      Taxonomy_clean),
    # remove ALL (cf), regardless of position or surrounding spaces
    Taxonomy_clean = str_replace_all(Taxonomy_clean, "\\(?cf\\)?", "") %>% 
      str_squish())  # remove any extra spaces left behind

# confirm that it worked - this returns TRUE if any row still contains "cf", FALSE otherwise
any(str_detect(data$Taxonomy_clean, "cf"))

# now check for near duplicates in taxonomic names, these may be typos/digitization errors
species = unique(data$Taxonomy_clean)

# calculate string distance matrix
dist_matrix = stringdistmatrix(species, species, method = "lv")

# convert to a matrix with row and column names
dist_matrix_named = as.matrix(dist_matrix)
rownames(dist_matrix_named) = species
colnames(dist_matrix_named) = species

# convert to a data frame of pairs
similar_names = as.data.frame(as.table(dist_matrix_named))

# filter out self-pairs and keep only small distances (e.g., ≤3)
near_duplicates = similar_names %>%
  filter(Var1 != Var2, Freq <= 3) %>%
  rowwise() %>%
  mutate(
    Name1 = min(as.character(Var1), as.character(Var2)),
    Name2 = max(as.character(Var1), as.character(Var2))
  ) %>%
  ungroup() %>%
  distinct(Name1, Name2, .keep_all = TRUE) %>%
  arrange(Freq)

# clean up column names
near_duplicates = near_duplicates %>%
  rename(species1 = Var1, species2 = Var2, distance = Freq) %>%
  filter(
    !grepl("\\bsp\\.", species1) &
      !grepl("\\bsp\\.", species2))

# look for the entries that have both genus and species names, check if any are
# duplicates with little typos.
# I don't see any issues - these are truly different species
rm(list = c("near_duplicates", "similar_names", "dist_matrix", "dist_matrix_named"))

# how many species were recorded over the full 2006-2024 period? these must be ID'd
# to the species level - not containing "sp."
species_IDd = as.data.frame(species) %>%
  filter(!str_detect(species, "sp\\.")) %>%
  distinct(species)

# how many unique transects are there (ignoring the width of the transect belt)
data = data %>%
  mutate(
    Transect_no_width = str_remove_all(Location, "\\s*Swath\\s*\\d+\\s*m?") %>%
      str_trim())

# there should be 72 unique transects:
# 6 LTER sites X 3 habitats per site X 4 permanent transects per habitat
unique_transects = data %>%
  distinct(Transect_no_width) 

# still ignoring the width of the transect belt, record whether each species was observed 
# or not along the transect line in each year
transect_data_PA = data %>%
  group_by(Year, Date, Site, Transect_no_width, Taxonomy_clean) %>% 
  summarise(Count = as.integer(sum(Count) > 0), .groups = "drop") %>%
  rename(PA = Count) # 1 if seen, 0 if not

# do the same, but for abundance instead of presence-absence
transect_data_ABU = data %>%
  group_by(Year, Date, Site, Transect_no_width, Taxonomy_clean) %>% 
  summarise(Count = sum(Count), .groups = "drop") %>%
  rename(Abundance = Count) # sum of abundance per species per transect per year

# at this point, we're missing the transects where each species was entirely absent 
# (PA = 0 and Abundance = 0). we need to add these absence locations back in

# get all unique year × transect combinations (ignoring transect width)
unique_transects = data %>%
  distinct(Year, Date, Site, Transect_no_width)

# get all unique species
unique_species = data %>%
  distinct(Taxonomy_clean)

# create all possible year × transect × species combos
all_transect_species = expand_grid(unique_transects, unique_species)

# now complete the presence–absence data
transect_data_PA = data %>%
  group_by(Year, Date, Site, Transect_no_width, Taxonomy_clean) %>%
  summarise(PA = as.integer(sum(Count) > 0), .groups = "drop") %>%
  right_join(all_transect_species, 
             by = c("Year", "Date", "Site", "Transect_no_width", "Taxonomy_clean")) %>%
  mutate(PA = replace_na(PA, 0))

# now complete the abundance data
transect_data_ABU = data %>%
  group_by(Year, Date, Site, Transect_no_width, Taxonomy_clean) %>%
  summarise(Abundance = sum(Count), .groups = "drop") %>%
  right_join(all_transect_species, 
             by = c("Year", "Date", "Site", "Transect_no_width", "Taxonomy_clean")) %>%
  mutate(Abundance = replace_na(Abundance, 0))

# the data are still at the level of the four pseudoreplicate transects within
# each year-site-habitat combination - we want to know presence-absence aggregated
# across those four transects. 

# extract habitat type and aggregate presence-absence
aggregate_data_PA = transect_data_PA %>%
  mutate(Habitat = case_when(
    str_detect(Transect_no_width, "Backreef") ~ "Backreef",
    str_detect(Transect_no_width, "Forereef") ~ "Forereef",
    str_detect(Transect_no_width, "Fringing") ~ "Fringing",
    TRUE ~ NA_character_)) %>%
  # if present in ANY of the 4 transects, mark as present (1)
  group_by(Year, Site, Habitat, Taxonomy_clean) %>%
  summarise(PA = as.integer(max(PA) > 0),
            .groups = "drop")

# we should have 6 sites X 3 habitats X 19 years = 342 observations
aggregate_data_PA %>%
  distinct(Year, Site, Habitat) %>%
  nrow()

# get the full species list for each Year-Site-Habitat combination
species_per_combination = aggregate_data_PA %>%
  group_by(Year, Site, Habitat) %>%
  summarise(species_list = list(sort(unique(Taxonomy_clean))),
            n_species = n_distinct(Taxonomy_clean),
            .groups = "drop")

# check if all combinations have the same number of species
species_per_combination %>%
  count(n_species)
length(unique(species_per_combination$n_species)) == 1

# now repeat the aggregation process for the abundance data
# the data are still at the level of the four pseudoreplicate transects within
# each year-site-habitat combination - we want to sum abundance across those four transects. 
aggregate_data_ABU = transect_data_ABU %>%
  mutate(Habitat = case_when(
    str_detect(Transect_no_width, "Backreef") ~ "Backreef",
    str_detect(Transect_no_width, "Forereef") ~ "Forereef",
    str_detect(Transect_no_width, "Fringing") ~ "Fringing",
    TRUE ~ NA_character_)) %>%
  # sum abundances across the 4 transects within each habitat
  group_by(Year, Site, Habitat, Taxonomy_clean) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE),
            .groups = "drop")

# we should have 6 sites X 3 habitats X 19 years = 342 observations
aggregate_data_ABU %>%
  distinct(Year, Site, Habitat) %>%
  nrow()

# get the complete species list for each Year-Site-Habitat combination
species_per_combination = aggregate_data_ABU %>%
  group_by(Year, Site, Habitat) %>%
  summarise(species_list = list(sort(unique(Taxonomy_clean))),
            n_species = n_distinct(Taxonomy_clean),
            .groups = "drop")

# check if all combinations have the same number of species
species_per_combination %>%
  count(n_species)
length(unique(species_per_combination$n_species)) == 1

# finally, pivot the data to wide format
aggregate_data_PA_wide = aggregate_data_PA %>%
  pivot_wider(
    names_from = Taxonomy_clean,
    values_from = PA,
    values_fill = 0)
  
aggregate_data_ABU_wide = aggregate_data_ABU %>%
  pivot_wider(
    names_from = Taxonomy_clean,
    values_from = Abundance,
    values_fill = 0)

# remove the long format tables 
rm(list = c("aggregate_data_PA", "aggregate_data_ABU"))

# triple check that the PA and Abundance data align. 
# species should only be present (PA = 1) when abundance is >0 and vice versa. 
# get the species columns (exclude Year, Site, Habitat)
species_cols = setdiff(names(aggregate_data_ABU_wide), c("Year", "Site", "Habitat"))

inconsistencies = tibble()

for (sp in species_cols) {
  # where PA = 1 but Abundance = 0
  pa_but_no_abu = sum(aggregate_data_PA_wide[[sp]] == 1 & aggregate_data_ABU_wide[[sp]] == 0)
  
  # where Abundance > 0 but PA = 0
  abu_but_no_pa = sum(aggregate_data_ABU_wide[[sp]] > 0 & aggregate_data_PA_wide[[sp]] == 0)
  
  if (pa_but_no_abu > 0 | abu_but_no_pa > 0) {
    inconsistencies = bind_rows(inconsistencies, 
                                 tibble(species = sp,
                                        pa_but_no_abundance = pa_but_no_abu,
                                        abundance_but_no_pa = abu_but_no_pa))
  }
}

# check results
if (nrow(inconsistencies) == 0) {
  print("TRUE: All PA and Abundance data are consistent!")
} else {
  print("FALSE: Inconsistencies found!:")
  print(inconsistencies)
}

# the aggregate PA and aggregate ABU wide dataframes now include all possible species, 
# before accounting for prevalence/rarity or traits.

#### FILTERING BASED ON # OF OBSERVATIONS/RARITY ####

# filter the data further to remove species that were observed less than 10 times
aggregate_data_PA_wide = aggregate_data_PA_wide %>%
  select(1:3, 
         all_of(names(.)[4:418][colSums(.[4:418], na.rm = TRUE) >= 10]))

# filter species based on occurrence frequency --> they should present at at least 10% of 
# sampling units, but no more than 90%. this removes too rare and too ubiquitous species. 
aggregate_data_PA_wide = aggregate_data_PA_wide %>%
  select(1:3, 
         where(~ if(is.numeric(.x)) {
           occurrence_rate = sum(.x, na.rm = TRUE) / nrow(aggregate_data_PA_wide)
           occurrence_rate >= 0.10 & occurrence_rate <= 0.90} else {
           TRUE  # keep columns Year, Site, Habitat
         }))

# remove those species from the abundance data, too
aggregate_data_ABU_wide = aggregate_data_ABU_wide %>%
  select(1:3, all_of(names(aggregate_data_PA_wide)[4:ncol(aggregate_data_PA_wide)]))

# remove temporary data we no longer need 
rm(list = c("all_transect_species", "existing_nums", "inconsistencies", "transect_data_ABU",
            "transect_data_PA", "unique_species", "unique_transects", "unnumbered_sp", "sp", 
            "species", "species_cols", "abu_but_no_pa", "pa_but_no_abu"))

# here are the species once accounting for rarity
species_names = names(aggregate_data_PA_wide)[4:ncol(aggregate_data_PA_wide)]

#### TAXONOMIC DATA ####
# read in the taxonomic information provided by the MCR LTER program
taxo = read.csv(here("Data", "Fish", "MCR_LTER_Fish_Species_List_20250325.csv"))

# keep only bony fishes in class Actinopterygii
taxo = taxo %>% filter(class == "Actinopterygii") %>%
  filter(! speciesbinomial == "NULL")

# check that we know at least down to the genus level
any(taxo$speciesbinomial == "NULL")
any(taxo$genus == "NULL")
any(taxo$family == "NULL")

# drop any observations where multiple unidentified species were listed together (spp.)
taxo = taxo %>% filter(!str_detect(speciesbinomial, "spp\\."))

# drop the unidentified records
taxo = taxo %>% filter(!str_detect(speciesbinomial, "unidentified"))

# remove common names
taxo = taxo %>% select(-commonname)
taxo[is.na(taxo)] = "Unknown"

# extract genus + any existing number from "Genus sp" or "Genus sp. X"
taxo = taxo %>%
  mutate(
    match = str_match(speciesbinomial, "^([A-Z][a-z]+) (sp\\.?)(?: ([0-9]+))?$"),
    genus_extracted = match[,2],
    sp_marker = match[,3],
    existing_num = match[,4])

# get current max assigned number for each genus (if any)
existing_nums = taxo %>%
  filter(!is.na(existing_num)) %>%
  mutate(existing_num = as.integer(existing_num)) %>%
  group_by(genus_extracted) %>%
  summarise(max_num = max(existing_num), .groups = "drop")

# assign new numbers to unnumbered "Genus sp"
unnumbered_sp = taxo %>%
  mutate(row_id = row_number()) %>%
  filter(!is.na(sp_marker) & is.na(existing_num)) %>%
  group_by(genus_extracted) %>%
  mutate(new_num = row_number()) %>%
  ungroup() %>%
  left_join(existing_nums, by = "genus_extracted") %>%
  mutate(
    final_num = new_num + coalesce(max_num, 0),
    speciesbinomial_clean = paste0(genus_extracted, " sp. ", final_num)) %>%
  select(row_id, speciesbinomial_clean)

# merge back into full dataset and construct final speciesbinomial_clean column
taxo = taxo %>%
  mutate(row_id = row_number()) %>%
  left_join(unnumbered_sp, by = "row_id") %>%
  mutate(
    speciesbinomial_clean = coalesce(
      speciesbinomial_clean,  # from unnumbered_sp if applicable
      ifelse(
        !is.na(sp_marker) & !is.na(existing_num),
        paste0(genus_extracted, " sp. ", existing_num),
        speciesbinomial))) %>%
  select(-match, -sp_marker, -existing_num, -row_id, -genus_extracted)

taxo = taxo %>%
  mutate(speciesbinomial_clean = ifelse(
    str_detect(speciesbinomial_clean, "\\(cf\\) sp\\.(?!\\s*\\d)"),
    str_replace(speciesbinomial_clean, "\\(cf\\) sp\\.", "(cf) sp. 1"),
    speciesbinomial_clean))

# check for near duplicates (these may be typos)
# get unique species names
species = unique(taxo$speciesbinomial_clean)

# calculate string distance matrix
dist_matrix = stringdistmatrix(species, species, method = "lv")

# convert to a matrix with row and column names
dist_matrix_named = as.matrix(dist_matrix)
rownames(dist_matrix_named) = species
colnames(dist_matrix_named) = species

# convert to a data frame of pairs
similar_names = as.data.frame(as.table(dist_matrix_named))

# filter out self-pairs and keep only small distances (e.g., ≤3)
near_duplicates = similar_names %>%
  filter(Var1 != Var2, Freq <= 3) %>%
  arrange(Freq)

# clean up column names
near_duplicates = near_duplicates %>%
  rename(species1 = Var1, species2 = Var2, distance = Freq) %>%
  filter(
    !grepl("\\bsp\\.", species1) &
      !grepl("\\bsp\\.", species2))

# there are three typos!
# fix the typos based on the accepted names on FishBase
# Anampses twisti (incorrect) --> Anampses twistii (correct)
# Antennarius tuberosus (incorrect) --> Antennatus tuberosus  (correct)
# Pseudogramma polyacanthum (incorrect) --> Pseudogramma polyacantha (correct)
taxo = taxo %>%
  mutate(
    speciesbinomial = case_when(
      speciesbinomial == "Anampses twisti" ~ "Anampses twistii",
      speciesbinomial == "Antennarius tuberosus" ~ "Antennatus tuberosus",
      speciesbinomial == "Pseudogramma polyacanthum" ~ "Pseudogramma polyacantha",
      TRUE ~ speciesbinomial),
    speciesbinomial_clean = case_when(
      speciesbinomial_clean == "Anampses twisti" ~ "Anampses twistii",
      speciesbinomial_clean == "Antennarius tuberosus" ~ "Antennatus tuberosus",
      speciesbinomial_clean == "Pseudogramma polyacanthum" ~ "Pseudogramma polyacantha",
      TRUE ~ speciesbinomial_clean))

# replace the original speciesbinomial with the new, corrected speciesbinomial_clean
taxo = taxo %>%
  select(-speciesbinomial) %>%
  rename(speciesbinomial = speciesbinomial_clean)

# (cf) likely means confer or compare to, suggesting a confident but uncertain ID
# drop (cf) and assume that the IDs are correct because (cf) will cause problems when
# constructing the taxonomic tree
taxo = taxo %>%
  mutate(speciesbinomial = str_replace_all(speciesbinomial, "\\s*\\(cf\\)", ""))
any(str_detect(taxo$speciesbinomial, "cf"))

# finally, keep only the species that we have in the community/transect survey data 
taxo = taxo %>%
  filter(speciesbinomial %in% species_names)
setdiff(taxo$species_binomial, species_names)

# remove the full list of species in the original taxonomy spreadsheet (no longer needed)
rm(species) 

# are we missing taxonomic information for any species?
lost_taxa = as.data.frame(setdiff(species_names, taxo$speciesbinomial))

# yes, we're missing information for blennies -> add this in manually
blennies = tibble(
  phylum = rep("Chordata", 2),
  subphylum = rep("Craniata", 2),
  infraphylum = rep("Vertebrate", 2),
  class = rep("Actinopterygii", 2),
  subclass = rep("Neopterygii", 2),
  division = rep("Teleostei", 2),
  subdivision = rep("Euteleostei", 2),
  order = rep("Blenniiformes", 2),
  family = rep("Blenniidae", 2),
  genus = c("Plagiotremus", "Cirripectes"),
  speciesbinomial = c(
    "Plagiotremus tapeinosoma",
    "Cirripectes variolosus"))

# add the missing species info back to our taxonomic data
taxo = rbind(taxo, blennies) %>%
  distinct(speciesbinomial)

# are there any other missing species?
lost_taxa = as.data.frame(setdiff(species_names, taxo$speciesbinomial))
# 0 - no more problems

# remove some of the temp data we no longer need
rm(list = c("lost_taxa", "blennies", "existing_nums", "near_duplicates",
            "similar_names", "dist_matrix", "dist_matrix_named", "unnumbered_sp"))

# now build a taxonomic tree from the taxonomic paths
taxo$path = apply(taxo, 1, function(x) paste(na.omit(x), collapse = ";"))
taxo$pathString = paste("Life", taxo$phylum, taxo$subphylum, taxo$infraphylum,
                        taxo$class, taxo$subclass, taxo$division, taxo$subdivision,
                        taxo$order, taxo$family, taxo$genus,
                        taxo$speciesbinomial, sep = "/")
taxo_tree = as.Node(taxo, pathName = "pathString")
phylo_tree = as.phylo.Node(taxo_tree)

# save the taxonomic tree for the full LTER dataset (before considering traits)
svg(here("Figures", "Full_LTER_Taxonomic_Tree.svg"), width = 20, height = 40)
plot.phylo(phylo_tree, cex = 0.35)
dev.off()

# save the tree as a .tre file (Newick format)
write.tree(phylo_tree, file = here("Data", "Fish", "FULL_LTER_Taxonomic_Tree.tre"))

#### TRAITS ####
# download different trait tables, using our species list and the rfishbase package
traits_species = species(species_names, 
                         fields = c("Species", "BodyShapeI", "DemersPelag"))

traits_ecology = ecology(species_names, 
                         fields = c("Species", "Schooling", "Shoaling",
                                    "Solitary", "DietTroph", "FoodTroph"))

traits_reproduction = reproduction(species_names,
                                   fields = c("Species", "ReproMode", "SpawnAgg"))

# ! NOTE here that we will add species' maximum lengths later in the script (this will
# require a larger code chunk, which I'd like to keep separate for clarity)

# from traits_species, do some quick renaming of columns/levels
unique(traits_species$BodyShapeI)
traits_species = traits_species %>%
  mutate(BodyShapeI = recode(BodyShapeI,
                             "elongated" = "Elongated",
                             "short and / or deep" = "Short/Deep",
                             "fusiform / normal" = "Fusiform",
                             "eel-like" = "Eel-Like"))

unique(traits_species$DemersPelag)
traits_species = traits_species %>%
  mutate(DemersPelag = recode(DemersPelag,
                              "reef-associated" = "Reef-Associated",
                              "benthopelagic" = "Benthopelagic"))

# from traits_ecology, create a unified "Troph" column where quantitative stomach content data 
# (in FoodTroph) are prioritized over qualitative estimates from the literature (in DietTroph)
traits_ecology = traits_ecology %>%
  mutate(Troph = coalesce(FoodTroph, DietTroph)) %>%
  select(-FoodTroph, -DietTroph)

# now some quick cleaning 
unique(traits_reproduction$ReproMode)
traits_reproduction = traits_reproduction %>%
  mutate(ReproMode = recode(ReproMode,
                            "dioecism" = "Dioecism",
                            "protandry" = "Protandry",
                            "protogyny" = "Protogyny"))

# use +1 for yes, 0 for no in schooling, solitary, and shoaling
traits_ecology = traits_ecology %>%
  mutate(Schooling = Schooling*-1) %>%
  mutate(Shoaling = Shoaling*-1) %>%
  mutate(Solitary = Solitary*-1)

# use +1 for yes, 0 for no in SpawnAgg
traits_reproduction = traits_reproduction %>%
  mutate(SpawnAgg = SpawnAgg*-1) %>%
  rename(Spawn_Agg = SpawnAgg)

# join all three trait datasets by "Species"
traits_all = traits_species %>%
  inner_join(traits_ecology, by = "Species") %>%
  inner_join(traits_reproduction, by = "Species")

# a few small changes to the column names in the trait dataframe
traits_all = traits_all %>%
  rename(Body_Shape = BodyShapeI,
         Column_Position = DemersPelag,
         Trophic_Level = Troph,
         Reproductive_Mode = ReproMode)

# any missing trait data?
any(is.na(traits_all))

# it's likely that we're missing some trophic info, check this
summary(traits_all$Trophic_Level)

# yes, we are! what's missing?
missing_troph_data = traits_all %>%
  filter(is.na(Trophic_Level)) %>%
  pull(Species) %>%
  as.data.frame()

# add in the missing Trophic_Level values using information from the FishBase website 
# directly, where values are means based on size and trophs of closest relatives
traits_all = traits_all %>%
  mutate(Trophic_Level = case_when(
    Species == "Amblygobius nocturnus"      ~ 3.0,
    Species == "Caracanthus maculatus"      ~ 3.3,
    Species == "Cheilinus oxycephalus"     ~ 3.4,
    Species == "Chlorurus microrhinos"     ~ 2.0,
    Species == "Dascyllus flavicaudus"     ~ 3.0,
    Species == "Gnatholepis anjerensis"    ~ 3.4,
    Species == "Pomachromis fuscidorsalis" ~ 3.5,
    Species == "Pseudocheilinus tetrataenia" ~ 3.3,
    Species == "Pseudojuloides atavai"     ~ 3.4,
    Species == "Pycnochromis acares"       ~ 3.0,
    Species == "Synodus binotatus"         ~ 4.0,
    TRUE ~ Trophic_Level  # keep existing Troph values for all other species
  ))

# any other missing trait data?
any(is.na(traits_all))

# where are those NAs?
colSums(is.na(traits_all)) > 0

# we're missing reproductive mode info for a few species...this is also unavailable
# on the FishBase website - so drop these species
traits_all = traits_all %>% tidyr::drop_na()

# add the maximum total length, in centimeters, for each species
###### start maximum TL calculations ######

# canonicalize species - try to use rfishbase validate_names if available
canonicalize_species = function(species) {
  if ("rfishbase" %in% .packages()) {
    # try-catch in case validate_names not available or errors
    tryCatch({
      vn = rfishbase::validate_names(species)
      if (!is.null(vn) && length(vn) > 0 && !is.na(vn[1])) return(vn[1])
    }, error = function(e) species)
  }
  species
}

# robust convert_to_TL function that uses all conversions and inverts reversed equations
convert_to_TL = function(species, length, from_type, conversion_df) {
  # basic guards
  if (is.na(length) || is.na(from_type) || is.null(species) || species == "") return(NA_real_)
  if (toupper(from_type) == "TL") return(as.numeric(length))
  
  species_can = canonicalize_species(species)
  
  # try to find candidate rows for that species (be flexible)
  conv_rows = conversion_df %>%
    filter(!is.na(Species)) %>%
    filter(Species == species_can | Species == species |
             str_detect(tolower(Species), fixed(tolower(species), ignore_case = TRUE)))
  # if none found at all, try fuzzy partial match to be diagnostic
  if (nrow(conv_rows) == 0) {
    conv_rows = conversion_df %>%
      filter(!is.na(Species)) %>%
      filter(str_detect(tolower(Species), fixed(tolower(str_replace_all(species, "\\s+.*$", "")), ignore_case = TRUE)))
  }
  if (nrow(conv_rows) == 0) return(NA_real_)
  
  # normalize column names for detection
  nm = names(conv_rows) %>% tolower()
  
  # ensure numeric a, b and coerce when possible
  conv_rows = conv_rows %>%
    mutate(
      a_num = suppressWarnings(as.numeric(.data$a)),
      b_num = suppressWarnings(as.numeric(.data$b)),
      # treat missing intercept as 0 when slope exists
      a_num = ifelse(is.na(a_num) & !is.na(b_num), 0, a_num)
    )
  
  # loop over rows and attempt to compute predicted TL from each row,
  # using heuristics to determine orientation/direction
  predicted_list = vector("numeric", length = 0)
  used_rows = 0L
  
  for (i in seq_len(nrow(conv_rows))) {
    row = conv_rows[i, ]
    a_i = row$a_num
    b_i = row$b_num
    
    # skip if b is missing or zero (can't invert or compute)
    if (is.na(b_i) || b_i == 0) next
    
    # determine type info if present
    type_field = if ("type" %in% names(row)) as.character(row$Type) else NA_character_
    
    # look for explicit length-type columns if they exist - common naming variants
    # try multiple probable column name variants
    length_type1 = NA_character_; length_type2 = NA_character_
    possible_names = names(row)
    # common variants:
    possible_t1 = c("lengthtype1", "length_type1", "length1_type", "lengthtype_1", "len1_type")
    possible_t2 = c("lengthtype2", "length_type2", "length2_type", "lengthtype_2", "len2_type")
    for (pn in possible_t1) if (pn %in% tolower(possible_names)) length_type1 = as.character(row[[which(tolower(possible_names) == pn)]])
    for (pn in possible_t2) if (pn %in% tolower(possible_names)) length_type2 = as.character(row[[which(tolower(possible_names) == pn)]])
    
    # normalize from_type and TL
    from_up = toupper(from_type)
    TL_up = "TL"
    
    # heuristic checks to infer mapping:
    # 1) if length_type2 == TL and length_type1 == from_type => equation likely: Length2 (TL) = a + b * Length1 (from)
    # 2) if length_type1 == TL and length_type2 == from_type => equation likely: Length1 (TL) = a + b * Length2 (from)
    # 3) if Type string contains both TL and from_type, try to detect order "TL.*FL" etc.
    # 4) fallback: try both orientations (direct and inverted), but only keep valid results.
    predicted_here = c()
    
    # helper to safely compute direct TL = a + b * length
    direct_predict = function(a, b, len) {
      as.numeric(a + b * len)
    }
    # helper to invert equation when equation is: from = a + b * TL  -> TL = (from - a) / b
    invert_predict = function(a, b, len) {
      as.numeric((len - a) / b)
    }
    
    used_this_row = FALSE
    
    if (!is.na(length_type2) && !is.na(length_type1)) {
      lt1 = toupper(trimws(length_type1)); lt2 = toupper(trimws(length_type2))
      # If lt2 is TL and lt1 is from_type => LT2 = a + b * LT1  -> compute
      if (lt2 == TL_up && lt1 == from_up) {
        predicted_here = c(predicted_here, direct_predict(a_i, b_i, length))
        used_this_row = TRUE
      } else if (lt1 == TL_up && lt2 == from_up) {
        # equation: LT1 = a + b * LT2  -> LT2 is from -> invert
        predicted_here = c(predicted_here, invert_predict(a_i, b_i, length))
        used_this_row = TRUE
      }
      # if explicit columns exist but don't match, fall through to other heuristics
    }
    
    # check Type string heuristics if not used yet
    if (!used_this_row && !is.na(type_field)) {
      tf = toupper(type_field)
      # if pattern looks like "TL-FL" or "TL vs FL" etc, try to interpret order
      # check for "TL" before "FL" meaning likely TL = a + b * FL
      if (str_detect(tf, TL_up) && str_detect(tf, from_up)) {
        # check which appears earlier in string (order)
        pos_TL = str_locate(tf, TL_up)[1]
        pos_from = str_locate(tf, from_up)[1]
        if (!is.na(pos_TL) && !is.na(pos_from)) {
          if (pos_TL < pos_from) {
            # TL appears before FROM -> commonly means "TL-FL" -> TL = a + b * FL
            predicted_here = c(predicted_here, direct_predict(a_i, b_i, length))
            used_this_row = TRUE
          } else {
            # FROM appears before TL -> e.g. "FL-TL" -> FL = a + b * TL  -> invert
            predicted_here = c(predicted_here, invert_predict(a_i, b_i, length))
            used_this_row = TRUE
          }
        } else {
          # ambiguous Type string; try both direct and inverted as fallback
          predicted_here = c(predicted_here, direct_predict(a_i, b_i, length))
          predicted_here = c(predicted_here, invert_predict(a_i, b_i, length))
          used_this_row = TRUE
        }
      }
    }
    
    # if still not used, fallback: attempt both orientations but check for reasonable numeric outputs
    if (!used_this_row) {
      # direct TL = a + b * from_length
      direct_val = tryCatch(direct_predict(a_i, b_i, length), error = function(e) NA_real_)
      invert_val = tryCatch(invert_predict(a_i, b_i, length), error = function(e) NA_real_)
      
      # decide whether direct or inverted seems plausible:
      # heuristic: predicted TL should be positive and not absurd (e.g., not orders of magnitude off)
      plausible = function(x) { !is.na(x) && is.finite(x) && x > 0 && x < 10000 } # coarse bounds
      if (plausible(direct_val)) predicted_here = c(predicted_here, direct_val)
      if (plausible(invert_val)) predicted_here = c(predicted_here, invert_val)
      
      # if both plausible, keep both (we will average). If neither plausible, skip this row.
    }
    
    # accept only finite positive predictions
    predicted_here = predicted_here[is.finite(predicted_here) & !is.na(predicted_here) & predicted_here > 0]
    
    if (length(predicted_here) > 0) {
      predicted_list = c(predicted_list, predicted_here)
      used_rows = used_rows + 1L
    }
  } # end loop over conversion rows
  
  # if we got no predictions, return NA
  if (length(predicted_list) == 0) return(NA_real_)
  # here we return mean
  mean_pred = mean(predicted_list, na.rm = TRUE)
  # and here we attach attributes useful for debugging
  attr(mean_pred, "n_equations_used") = used_rows
  attr(mean_pred, "n_predictions") = length(predicted_list)
  attr(mean_pred, "sd_predicted_TL") = ifelse(length(predicted_list) > 1, sd(predicted_list, na.rm = TRUE), NA_real_)
  mean_pred
}

# now run the functions
species_list = traits_all$Species

length_data = rfishbase::species(species_list) %>%
  select(Species, Length, LTypeMaxM)

length_conversions = rfishbase::length_length(species_list)

max_length_data = length_data %>%
  rowwise() %>%
  mutate(
    Max_TL_cm = convert_to_TL(Species, Length, LTypeMaxM, length_conversions),
    # diagnostics: how many equations were used (from attribute)
    n_equations_used = ifelse(!is.na(Max_TL_cm), as.integer(attr(Max_TL_cm, "n_equations_used")), 0L),
    n_predictions = ifelse(!is.na(Max_TL_cm), as.integer(attr(Max_TL_cm, "n_predictions")), 0L),
    sd_predicted_TL = ifelse(!is.na(Max_TL_cm), as.numeric(attr(Max_TL_cm, "sd_predicted_TL")), NA_real_)
  ) %>%
  ungroup() %>%
  select(Species,
         Original_Length = Length,
         Original_Type = LTypeMaxM,
         Max_TL_cm,
         n_equations_used,
         n_predictions,
         sd_predicted_TL)

###### end maximum TL calculations ###### 

# add the maximum total length data to our other traits
traits_all = traits_all %>%
  left_join((max_length_data %>%
          select(Species, Max_TL_cm)),
       by = "Species")

# one final check to make sure there's no missing trait data
anyNA(traits_all)

# how many species were retained from the full LTER list?
length(unique(traits_all$Species))

# how many were lost along the way as we added more trait data?
length(unique(traits_species$Species)) - length(unique(traits_all$Species))
length(unique(traits_ecology$Species)) - length(unique(traits_all$Species))
length(unique(traits_reproduction$Species)) - length(unique(traits_all$Species))

# filter the taxo dataframe to keep only species for which trait data are available
taxo_filtered = taxo %>%
  filter(speciesbinomial %in% traits_all$Species)

# make a new taxonomic tree
taxo_tree_filtered = as.Node(taxo_filtered, pathName = "pathString")
phylo_tree_filtered = as.phylo.Node(taxo_tree_filtered)

# save the graphic
svg(here("Figures", "Filtered_LTER_Taxonomic_Tree.svg"), width = 20, height = 40)
plot.phylo(phylo_tree_filtered, cex = 0.75)
dev.off()

# save the tree as a .tre file (Newick format)
write.tree(phylo_tree_filtered, file = here("Data", "Fish", "Filtered_LTER_Taxonomic_Tree.tre"))

# final step! we need to now go back to the MCR LTER survey data and keep only 
# those species that we have taxonomic and trait data for! 
# get species names from columns 4 to end
species_from_transects = names(aggregate_data_PA_wide)[4:ncol(aggregate_data_PA_wide)]

# keep only species that are in traits_all$Species
species_to_keep = species_from_transects[species_from_transects %in% traits_all$Species]

# filter the PA and Abundance dataframes
aggregate_data_PA_wide = aggregate_data_PA_wide %>%
  select(1:3, all_of(species_to_keep))

aggregate_data_ABU_wide = aggregate_data_ABU_wide %>%
  select(1:3, all_of(species_to_keep))

# create a new Name Id column by combining the Site and Habitat columns
aggregate_data_PA_wide = aggregate_data_PA_wide %>%
  mutate(Name = paste0(gsub("_", "", Site), "_", Habitat)) %>%
  relocate(Name, .after = Year) %>%
  select(-Site, -Habitat)

aggregate_data_ABU_wide = aggregate_data_ABU_wide %>%
  mutate(Name = paste0(gsub("_", "", Site), "_", Habitat)) %>%
  relocate(Name, .after = Year) %>%
  select(-Site, -Habitat)

#### SAVING DATASETS FOR MODELLING ####
write.csv(aggregate_data_PA_wide,
          here("HMSC", "Data", "Intermediate_Datasets", "Fish_Presence_Absence_Dataset.csv"),
          row.names = FALSE)
write.csv(aggregate_data_ABU_wide,
          here("HMSC", "Data", "Intermediate_Datasets", "Fish_Abundance_Dataset.csv"),
          row.names = FALSE)
write.csv(traits_all,
          here("HMSC", "Data", "Intermediate_Datasets", "Fish_Traits_Dataset.csv"),
          row.names = FALSE)
write.tree(phylo_tree_filtered, 
           file = here("HMSC", "Data", "Intermediate_Datasets", "Fish_Taxonomic_Tree.tre"))

# save R data with these final datasets
save(list = c("aggregate_data_PA_wide", "aggregate_data_ABU_wide",
              "traits_all", "phylo_tree_filtered", "species_to_keep"),
     file = here("Code", "Prepped_Community_Taxonomy_Trait_Data.RData"))
