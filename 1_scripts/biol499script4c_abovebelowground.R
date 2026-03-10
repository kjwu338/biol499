####
# BIOL 499 Script 4c: GLMMs/Species metrics for above/belowground
# Last Updated: March 4, 2026
# Description: running GLMMs based on depth groupings
####
#questions I have: 
# 1. Should I compare above/belowground community compositions? Yes. How do I do that?

library(glmmTMB) # zero inflated data
library(lme4) # mixed-effects models
library(lmerTest) # provide p values for mixed-effect models
library(dplyr) # summarizing and tranform tables
library(DHARMa) # testing for zero inflation
library(tidyr)
library(car)


# load data
gc <- read.csv("5_outputtables/gc_wcategories.csv")

#clean data (same between script a & c)
gc_clean <- gc %>%
  filter(microsite_condition != "double_control") %>%
  mutate(
    # Create the Above vs Belowground category
    depth_simple = case_match(surface_depth, 
                              "surface" ~ "aboveground",
                              c("0-5", "5-10", "10-May") ~ "belowground"),
    # Set factors and order
    depth_simple = factor(depth_simple, levels = c("aboveground", "belowground")),
    microsite_condition = factor(microsite_condition, 
                                 levels = c("biocrust", "dung", "control")),
    site = as.factor(site)
  )

gc_collapsed <- gc_clean %>%
  group_by(site, microsite_condition, depth_simple, species) %>%
  summarise(total_count = n(), .groups = "drop")

  # create master 120 list
all_trays_120 <- expand_grid(
  site = factor(1:20),
  microsite_condition = factor(c("biocrust", "dung", "control"), 
                               levels = c("biocrust", "dung", "control")),
  depth_simple = factor(c("aboveground", "belowground"), 
                        levels = c("aboveground", "belowground"))
)


#1. ABUNDANCE/COUNT ------

# STEP 1: normalize abundance by volume + clean data ----
# aggregate seedling counts by treatment
abundance_120 <- all_trays_120 %>%
  left_join(
    gc_collapsed %>% 
      group_by(site, microsite_condition, depth_simple) %>% 
      summarise(abundance = sum(total_count), .groups = "drop"),
    by = c("site", "microsite_condition", "depth_simple")
  ) %>%
  mutate(abundance = replace_na(abundance, 0))

abundance120_normalized <- abundance_120 %>%
  mutate(
    tray = row_number(),
    # 1. Determine volume based on the new depth category
    # Surface = 15.71, Belowground (0-5 + 5-10) = 31.42
    sample_volume = if_else(depth_simple == "aboveground", 15.71, 31.42),
    
    # 2. Normalize abundance by the specific volume
    seedpercm3 = abundance / sample_volume
  )

# set factors
tray_normal <- tray_normal %>%
  mutate(
    surface_depth = factor(surface_depth, levels = c("surface", "0-5", "5-10")),
    microsite_condition = factor(microsite_condition, levels = c("control", "biocrust", "dung"))
  )

# STEP 2: run glmm-----
#attempt 1: failed assumptions : glmer-----
#model_abundance <- glmer(seedpercm3 ~ microsite_condition * surface_depth + (1 | site), 
#data = tray_normal)

#anova(model_abundance) # anova table
#attempt 2: glmmTMB-----this model meets all assumptions----
model_abundance_120 <- glmmTMB(seedpercm3 ~ microsite_condition * depth_simple + (1|site), 
                               data = abundance120_normalized,
                               family = tweedie())

Anova(model_abundance, type = "II") # is this line of code okay? random effect model

# STEP 3: test assumptions ----
# check normality = right skewed
hist(abundance120_normalized$seedpercm3)

#tweedie model passed all assumptions
simulationOutput_120 <- simulateResiduals(fittedModel = model_abundance_120)
plot(simulationOutput_120) 
testOutliers(simulationOutput_120) 
testDispersion(simulationOutput_120)


#2. RICHNESS

# STEP 1: clean data and create species count
richness_counts_120 <- gc_clean %>%
  group_by(site, microsite_condition, depth_simple) %>%
  summarise(
    richness = n_distinct(species), 
    .groups = "drop"
  )

# STEP 2: join with master list to include empty trays ----
tray_richness_120 <- all_trays_120 %>%
  left_join(richness_counts_120, by = c("site", "microsite_condition", "depth_simple")) %>%
  mutate(
    tray = row_number(),
    richness = replace_na(richness, 0)
  )

# STEP 3: run glmmTMB for R]richness ----
# Richness is count data -> poisson + zero Inflation is standard
model_richness_120 <- glmmTMB(richness ~ microsite_condition * depth_simple + (1|site), 
                              data = tray_richness_120,
                              ziformula = ~1, 
                              family = poisson)

Anova(model_richness_120, type = "II")

# STEP 4: test assumptions ----
sim_richness_120 <- simulateResiduals(model_richness_120)
plot(sim_richness_120)

# 3. EVENNESS

# STEP 1: calculate shannon h for pielou’s j (evenness) ----
evenness_calc_120 <- gc_clean %>%
  # Count each species per microsite/depth combination
  group_by(site, microsite_condition, depth_simple, species) %>%
  summarise(sp_count = n(), .groups = "drop_last") %>%
  # Calculate proportions based on the NEW collapsed totals
  mutate(
    p_i = sp_count / sum(sp_count),
    ln_p_i = p_i * log(p_i)
  ) %>%
  summarise(
    shannon_h = -sum(ln_p_i),
    richness_val = n(), 
    .groups = "drop"
  ) %>%
  mutate(
    evenness = if_else(richness_val > 1, shannon_h / log(richness_val), 0)
  )

# STEP 2: join with master list to include empty trays ----
tray_evenness_120 <- all_trays_120 %>%
  left_join(evenness_calc_120, by = c("site", "microsite_condition", "depth_simple")) %>%
  mutate(
    tray = row_number(),
    # Beta family cannot handle exactly 0 or 1
    evenness_raw = case_when(
      evenness <= 0 ~ 0.001,
      evenness >= 1 ~ 0.999,
      TRUE ~ evenness
    )
  )

# STEP 3: run glmmTMB for evenness ----
model_evenness_120 <- glmmTMB(evenness_raw ~ microsite_condition * depth_simple + (1|site), 
                              dispformula = ~microsite_condition, 
                              data = tray_evenness_120,
                              family = beta_family())

Anova(model_evenness_120, type = "II")

# STEP 4: check assumptions -------
sim_evenness_120 <- simulateResiduals(model_evenness_120)
plot(sim_evenness_120)
