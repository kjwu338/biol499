#biol499 merging files
#date: 2/23/2026

#load data
st_kat <- read.csv("5_outputtables/kattrydataoutput.csv")
st_aly <- read.csv("5_outputtables/alytrydataoutput.csv")

# unique combinations
unique_keys <- c("AccSpeciesID", "TraitID", "OrigUnitStr", "n_samples")

merged_unique <- bind_rows(st_aly, st_kat) %>%
  distinct(across(all_of(unique_keys)), .keep_all = TRUE)

#save
write.csv(merged_unique, "5_outputtables/st_merged.csv", row.names = FALSE)
