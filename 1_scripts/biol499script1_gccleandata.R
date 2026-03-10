##
#BIOL 499 GC Script 1: Clean Data & Descriptive Analysis
#Date: 2/11/2026
#Description:
#Author: 
##

# load packages
library(dplyr)
library(tidyverse)

# read data
df <- read.csv("2_data_raw/gc_seedlingdata.csv") #load data

# PART 1: clean data ----
# keep only the first 3 columns 
df <- df %>% 
  dplyr::select(1:3)

# Split the sample column into three columns and deal with double control 
df_cleaned <- df %>%
  separate("sample", 
           into = c("site", "microsite_condition", "surface_depth"), 
           sep = "-", 
           remove = FALSE) %>% #keep original
  mutate(
    double_control = is.na(microsite_condition),
    microsite_condition = if_else(double_control, "double_control", microsite_condition),
    surface_depth = if_else(double_control, "double_control", surface_depth)
  )
    
df_cleannames <- df_cleaned %>%
    mutate(
    #recode microsite_condition with full names
    microsite_condition = case_when(
      microsite_condition == "D" ~ "dung",
      microsite_condition == "B" ~ "biocrust",
      microsite_condition == "C" ~ "control",
      TRUE ~ microsite_condition
    ),
    surface_depth = case_when(
      surface_depth == "S" ~ "surface",
      surface_depth == "0" ~ "0-5",
      surface_depth == "5" ~ "5-10",
      TRUE ~ surface_depth
    )
  ) 

write.csv(df_cleannames, "5_outputtables/gc_wcategories.csv")

# PART 2: Descriptive analysis ----

#1) Number of species + count 
speciescount <- df_cleannames %>%
  group_by(Species) %>%
  summarise (
    count = n()
  )
print(speciescount)

#2) Count per site
sitecount <- df_cleannames %>%
  group_by(site) %>%
  summarise(
    count = n()
  )
print(sitecount)

#3) gc sample count
gctraycount <- df_cleannames %>%
  group_by(Sample) %>%
  summarize(
    count = n()
  )

print(gctraycount)

#4) Group By condition: Tray count, Seedling Count and Species Count
intx_trays <- df_cleannames %>%
  group_by(microsite_condition,surface_depth) %>%
  summarize(
    tray_num = n_distinct(Sample), 
    seedling_count = n(), 
    species_count = n_distinct(Species)) %>%
  arrange(desc(tray_num))
print(intx_trays)

write.csv(intx_trays, "5_outputtables/intx_trays.csv")

#5) Group by just microsites
microsites <- df_cleannames %>%
  group_by(microsite_condition) %>%
  summarize(
    seedling_count = n(), 
    species_count = n_distinct(Species), 
    tray_count = n_distinct(Sample)
  )
print(microsites)

#Richness per surface depth


#output another dataframe