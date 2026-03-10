####
# BIOL 499 Script 3: Just Trait Analysis 
# Last Updated: March 5, 2026
# Description: determine if microsites (biocrusts and cow pats) act as selective 
# filters based on functional traits like seed mass, length, and shape.
# Notes: find better data for seed mass; maybe drop seed length

####
# PART 1: Data Cleaning & Merging
# PART 2: Descriptive Statistics (Mean, SD, CWM)
# PART 3: testing filtering (CV & KS Test)
# PART 4: Linear Mixed Models (LMM)
####

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(DHARMa)
library(glmmTMB)

# 1. DATA PREP ----
seedlingdf <- read.csv("6_manualediting/gc_seedtraitsfinal.csv")
gcseedling <- read.csv("5_outputtables/gc_wcategories.csv")

seedlingdf_clean <- seedlingdf %>%
  mutate(across(everything(), ~ na_if(as.character(.), "n/a"))) %>%
  mutate(final_mass = as.numeric(final_mass),
         final_length = as.numeric(final_length))

merged_seedling_traits <- gcseedling %>%
  left_join(seedlingdf_clean, by = c("species" = "Species")) %>% 
  filter(microsite_condition != "double_control") %>% 
  mutate(
    surface_depth = case_match(surface_depth, 
                               "10-May" ~ "5-10", 
                               .default = surface_depth),
    surface_depth = factor(surface_depth, levels = c("surface", "0-5", "5-10")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

# 2. DESCRIPTIVE STATS ----
# group level stats
trait_descriptive_stats <- merged_seedling_traits %>%
  group_by(microsite_condition, surface_depth) %>%
  summarise(
    mean_mass = mean(final_mass, na.rm = TRUE),
    sd_mass = sd(final_mass, na.rm = TRUE),
    mean_length = mean(final_length, na.rm = TRUE),
    sd_length = sd(final_length, na.rm = TRUE),
    n_seedlings = n(), .groups = "drop"
  )

# Community Weighted Means (CWM)
cwm_tray <- merged_seedling_traits %>%
  group_by(tray, site, microsite_condition, surface_depth) %>%
  summarise(cwm_mass = mean(final_mass, na.rm = TRUE),
            cwm_length = mean(final_length, na.rm = TRUE), .groups = "drop")

# 3. testing filtering  ----
# Coefficient of Variation (CV)
cv_summary <- trait_descriptive_stats %>%
  mutate(cv_mass = (sd_mass / mean_mass) * 100,
         cv_length = (sd_length / mean_length) * 100)

# KS Test for mass (Variety Clumping)
dung_mass_dist <- merged_seedling_traits %>% filter(microsite_condition == "dung") %>% pull(final_mass) %>% na.omit()
biocrust_mass_dist <- merged_seedling_traits %>% filter(microsite_condition == "biocrust") %>% pull(final_mass) %>% na.omit()
ks_result_mass <- ks.test(dung_mass_dist, biocrust_mass_dist)

# 4. glmm/lmm ----
# Full Models
#model_mass <- lmer(final_mass ~ microsite_condition * surface_depth + (1 | site/tray), data = merged_seedling_traits)
#model_length <- lmer(final_length ~ microsite_condition * surface_depth + (1 | site/tray), data = merged_seedling_traits)
model_mass <- glmmTMB(final_mass ~ microsite_condition * surface_depth + (1 | site/tray),
        data = merged_seedling_traits,
        family = Gamma(link = "log"))

# Simple Models (Microsite only)
#model_mass_micro <- lmer(final_mass ~ microsite_condition + (1 | site/tray), data = merged_seedling_traits)
#model_length_micro <- lmer(final_length ~ microsite_condition + (1 | site/tray), data = merged_seedling_traits)

# Post-hoc Significance Letters
emm_mass <- emmeans(model_mass_micro, ~ microsite_condition)
emm_letters_mass <- cld(emm_mass, Letters = letters, adjust = "tukey")

#5. Testing assumptions
hist(merged_seedling_traits$final_mass)
sim_mass <- simulateResiduals(model_mass)
plot(sim_mass)
testResiduals(sim_mass) #KS test (Normality), Dispersion test (Variance), Outlier test

#all lmer models failed assumptions
# for mass glmm loglink still failed, should I log transform?
