####
# BIOL 499 Script 4a: GLMMs/Species metrics by microsite & depth
# **univariate analysis
# Last Updated: March 5, 2026
# Description: determine if microsites (biocrusts and cow pats) act as selective 
# filters based on functional traits like seed mass, length, and shape.
####
#questions I have: 
# 1. Should I compare above/belowground community compositions? Yes. How do I do that?
# Notes: Since there is an unbalanced design, we will not cross microsite and depth, 
# rather we just do a one way comparison

####

library(glmmTMB) # zero inflated data
library(lme4) # mixed-effects models
library(lmerTest) # provide p values for mixed-effect models
library(dplyr) # summarizing and tranform tables
library(DHARMa) # testing for zero inflation
library(tidyr)
library(car)

# STEP 0: prep data
# load data
gc <- read.csv("5_outputtables/gc_wcategories.csv")

#clean data (same between script a & c)
gc_clean <- gc %>%
  filter(microsite_condition != "double_control") %>%
  mutate(
    surface_depth = case_match(surface_depth, 
                               "10-May" ~ "5-10", 
                               .default = surface_depth),
    surface_depth = factor(surface_depth, levels = c("surface", "0-5", "5-10")),
    microsite_condition = factor(microsite_condition, 
                                 levels = c("biocrust", "dung", "control")),
    site = as.factor(site)
  )

# create master 160 list
all_trays <- expand_grid(
  site = factor(1:20),
  combo = c("biocrust_surface", "biocrust_0-5", "biocrust_5-10",
            "dung_surface", "dung_0-5", "dung_5-10",
            "control_0-5", "control_5-10")
) %>%
  separate(combo, into = c("microsite_condition", "surface_depth"), sep = "_")
######

#1. ABUNDANCE/COUNT ------

# STEP 1: CLEAN DATA
# not normalizing abundance by volume (15.71) bc surface treatments are not standarized ----
# aggregate seedling counts by trays
tray_counts <- gc_clean %>%
  group_by(Tray.., site, microsite_condition, surface_depth) %>%
  summarise(total_seeds = n(), .groups = "drop") #drop grouping structure

tray_normal <- all_trays %>%
  left_join(tray_counts, by = c("site", "microsite_condition", "surface_depth")) %>%
  mutate(
    tray_id = row_number(),
    total_seeds = replace_na(total_seeds, 0), # assign 0 to trays with no seedlings
    treatment = paste(microsite_condition, surface_depth, sep = "_")
)

# set factors
treatment_levels <- c(
  "control_surface", "biocrust_surface", "dung_surface",
  "control_0-5",     "biocrust_0-5",     "dung_0-5",
  "control_5-10",    "biocrust_5-10",    "dung_5-10"
)

tray_normal <- tray_normal %>%
  mutate(
    treatment = factor(treatment, levels = treatment_levels),
    site = factor(site)
  )

# STEP 2: run glmm-----
#attempt 1: failed assumptions : glmer-----
#model_abundance <- glmer(seedpercm3 ~ microsite_condition * surface_depth + (1 | site), 
                        #data = tray_normal)

#anova(model_abundance) # anova table
#attempt 2: glmmTMB-----this model meets all assumptions----
#model_abundance <- glmmTMB(seedpercm3 ~ microsite_condition * surface_depth + (1|site), 
                                #data = tray_normal,
                                #family = tweedie())
#attempt 3: 
# count data; non-normal data, treatment = predictor, site as random effect
# glmmTMB for zero-inflated; Should i use?
model_abundance <- glmer(total_seeds ~ treatment + (1 | site), 
                         data = tray_normal,
                         family = poisson())

Anova(model_abundance, type = "III") # is this line of code okay? random effect model

# STEP 3: test assumptions ----
# check normality = right skewed
hist(tray_normal$total_seeds)

simulationOutput <- simulateResiduals(fittedModel = model_abundance)
plot(simulationOutput) # ks test(poor goodness of fit) and levene test(non-equal variance) is significant 
testOutliers(simulationOutput) #outlier test - not significant
testDispersion(simulationOutput)


#2. RICHNESS

richness_counts <- gc_clean %>%
  group_by(site, microsite_condition, surface_depth) %>%
  summarise(richness = n_distinct(Species), .groups = "drop")

tray_richness <- all_trays %>%
  left_join(richness_counts, by = c("site", "microsite_condition", "surface_depth")) %>%
  mutate(
    richness = replace_na(richness, 0),
    treatment = factor(paste(microsite_condition, surface_depth, sep = "_"), 
                       levels = treatment_levels)
  )

# Run model with combined treatment
model_richness <- glmer(richness ~ treatment + (1|site), 
                          data = tray_richness,
                          family = poisson)

Anova(model_richness, type = "III")

# STEP 4: test assumptions ----
sim_richness <- simulateResiduals(model_richness) # non significant/pass assumptions
plot(sim_richness)

# 3. EVENNESS

# STEP 1: calculate shannon h for pielou’s j (evenness) ----
evenness_calc <- gc_clean %>%
  # grouping 
  group_by(site, microsite_condition, surface_depth, Species) %>% 
  # identifying species 
  summarise(sp_count = n(), .groups = "drop_last") %>%
  # calculates relative abundance
  mutate(
    p_i = sp_count / sum(sp_count), 
    # natural log
    ln_p_i = p_i * log(p_i)
  ) %>%
  summarise(
    # shannon diversity
    shannon_h = -sum(ln_p_i),
    # richness = tray_richness$richness
    richness_val = n(), 
    .groups = "drop"
  ) %>%
  mutate(
    # Pielou's evenness = H′/ln(S)
    # evenness requires at least 2 species to be mathematically defined
    evenness = if_else(richness_val > 1, shannon_h / log(richness_val), NA_real_)
  )

# STEP 2: join with master list to include empty trays ----
tray_evenness <- all_trays %>%
  left_join(evenness_calc, by = c("site", "microsite_condition", "surface_depth")) %>%
  mutate(
    # Create the combined treatment factor
    treatment = factor(paste(microsite_condition, surface_depth, sep = "_"), 
                       levels = treatment_levels),
    # Beta distribution cannot handle exact 0 or 1.
    # Trays with 0 or 1 species result in NA
    evenness_raw = case_when(
      evenness <= 0 ~ 0.001,
      evenness >= 1 ~ 0.999,
      TRUE ~ evenness
    )
  )

# STEP 3: run glmmTMB for evenness ----
# We use dispformula to fix the Levene test by modeling variance per microsite
model_evenness <- glmmTMB(evenness_raw ~ treatment + (1 | site), 
                          dispformula = ~microsite_condition, 
                          data = tray_evenness, 
                          family = beta_family())

Anova(model_evenness, type = "III")

# STEP 4: check assumptions ----
sim_evenness <- simulateResiduals(model_evenness)
plot(sim_evenness) # good enough

