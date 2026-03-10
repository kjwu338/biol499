###
#BIOL 499 Script 2: Clean out the seed traits and species you want
#Last updated: 2/17/2026
#Author: Kelly Wu
#Description: Only include species from germination trait
#Included traits for analysis of covariance
###

#load library
library(dplyr)
library(tidyr)

#load data
seed_traits <-read.csv("5_outputtables/st_merged.csv")
species_list <- read.csv("2_data_raw/gc_seedlingdata.csv")

#PART 1. Clean data-----
length(unique(species_list$Species))

#select species you want
species_filter <- seed_traits %>%
  semi_join(species_list, by = c("species" = "Species"))
length(unique(species_filter$species))

unique_species_list <- unique(species_filter$species)

#select traits
#seed dry mass, seed length, seed longevity (33), seed width, seed shape
traits <- c(26, 27, 239, 349)

trait_filter <- species_filter %>%
  filter(TraitID %in% traits)

#select units
units <- c("g / 1000 seeds", "g/1000", "mm")

unit_filter <- trait_filter %>%
  filter(OrigUnitStr %in% units)


#Part 2. Restructure dataframe----
cleanspeciesdf <- unit_filter %>%
  # a) normal seed dry mass units
  # make "g / 1000 seeds" and "g/1000" the same unit for averaging
  mutate(OrigUnitStr = ifelse(OrigUnitStr %in% c("g / 1000 seeds", "g/1000"), "g/1000", OrigUnitStr)) %>%
  
  # b) group by species and trait name 
  group_by(AccSpeciesID, species, traitname) %>%
  
  # c) calculate weighted average based on sample size (n_samples)
  summarize(
    weighted_mean = sum(mean_value * n_samples, na.rm = TRUE) / sum(n_samples, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # d) pivot the data so each trait becomes a column
  pivot_wider(
    names_from = traitname, 
    values_from = weighted_mean
  ) %>%
  
  # e) replace missing values with "n/a"
  # Note: This converts numeric values to characters to allow the "n/a" string
  mutate(across(-c(AccSpeciesID, species), ~ as.character(.))) %>%
  mutate(across(-c(AccSpeciesID, species), ~ replace_na(., "n/a")))

# View result
print(cleanspeciesdf)

#Part 3: Merge with species list --------
# code for creating a base dataframe of all unique species from the original seedling data
# a) Left join the cleaned trait data so all original species are preserved
final_species_df <- species_list %>%
  select(Species) %>%
  distinct() %>%
  left_join(cleanspeciesdf, by = c("Species" = "species")) %>%
  
  # b) replace n/a with actual"n/a"
  # This ensures consistency with the "n/a" format in your cleanspeciesdf
  mutate(across(everything(), ~ replace_na(as.character(.), "n/a")))

# View the final dataframe containing all original species and their available traits
print(final_species_df)

write.csv(final_species_df, "5_outputtables/clean_gctrytraits.csv")

###READ ME BEFORE CLOSING SCRIPT: The output was highly modified using excel to incorporate values from various sources 
# before being input into script 3 and 4 for data analysis 
# modified df can be found in 6_manualediting => "gc_speciestraitsfinal.csv"