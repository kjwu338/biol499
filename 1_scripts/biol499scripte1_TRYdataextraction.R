###
#TRY data extraction code
#Last updated: 2/17/2026
#Author: Kelly Wu
#Description: Script to extract averages for traits per species from TRY data sets
#Problem: mean is not calculating properly, seed dry mass is huge
###

#Load packages
library(tidyverse)
library(readr)

#Input data
df <-read.delim("2_data_raw/traitsubmission_jan_katarinna.txt")

#When the cell in column "TraitID" is not blank, extract the whole row and create new table
cleandf <- df %>% filter(!is.na(TraitID))

#write.csv(cleandf,"5_outputtables/achifol_clean.csv")

#For every unique AccSpeciesID and TraitID + (UNIT) in the columns:
###
#Kelly's notes: 
#I was going to extract for most common unit but I realized there are some traits where there are lot of common units 
#So I just filtered for species, trait and unit. Can calculated means for all units
###

results <- cleandf %>% 
  mutate(
    is_numeric = OrigValueStr == "" | !is.na(as.numeric(as.character(OrigValueStr))),
    OrigValueStr_numeric = ifelse(is_numeric, as.numeric(as.character(OrigValueStr)), NA) #if else
  ) %>%  
  group_by(AccSpeciesID, TraitID, OrigUnitStr) %>% #group BY SPECIES AND TRAIT AND UNIT
  
  #Finding most common unit
  #add_count(OrigUnitStr, name = "unit_count") %>%  # Count occurrences of each unit
  #filter(unit_count == max(unit_count)) %>%     # Keep only rows with max count
  #slice(1:n()) %>%  # If there's a tie, keep first unit alphabetically
  
  #Calculation
  summarise(
    # Check if ALL values in the group are numeric
    all_numeric = all(is_numeric),  # Add this check first
    
    # Calculate mean ONLY if ALL values are numeric
    mean_value = ifelse(all_numeric,
                        mean(OrigValueStr_numeric, na.rm = TRUE), #mean calculation 
                        NA),
    
    #OUTPUT DATA FRAME
    # Concatenate unique words for text values
    text_values = ifelse(!all_numeric,
                         paste(unique(OrigValueStr), collapse = "; "),
                         NA),
    
    # Metadata/other info
    species = first(AccSpeciesName),
    traitname = first(TraitName), 
    reference = paste(unique(Reference), collapse = "; "),
    
    # Summary stats (only for numeric)
    n_samples = n(), 
    sd_value = ifelse(all_numeric,
                      sd(OrigValueStr_numeric, na.rm = TRUE), 
                      NA),
    min_value = ifelse(all_numeric,
                       min(OrigValueStr_numeric, na.rm = TRUE), 
                       NA),
    max_value = ifelse(all_numeric,
                       max(OrigValueStr_numeric, na.rm = TRUE), 
                       NA),
    
    .groups = 'drop'
  )

results_sorted <- results %>%
  arrange(species, traitname) %>%
  select(AccSpeciesID, species, TraitID, traitname, OrigUnitStr, n_samples, mean_value, everything())

####MAKE SURE YOU WANT THE OUTPUT ON YOUR COMPUTER: 
#Output table eg. possible csv
write.csv(results_sorted, "5_outputtables/kattrydataoutput.csv", row.names = FALSE)

