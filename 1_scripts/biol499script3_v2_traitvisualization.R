####
# BIOL 499 Script 3: Just Trait Visualization
# Last Updated: March 5, 2026
# Description: the graphs for seed trait: 
###

####
# PART 1: CV Bar Charts (Filtering)
# PART 2: Trait Boxplots (Raw Data)
# PART 3: Density Plots (KS Distribution)
# PART 4: CWM Bar Charts
####

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. CV BAR CHARTS ----
cv_long <- cv_summary %>%
  pivot_longer(cols = c(cv_mass, cv_length), names_to = "trait", values_to = "cv") %>%
  mutate(trait = recode(trait, "cv_mass" = "Mass CV (%)", "cv_length" = "Length CV (%)"))

ggplot(cv_long, aes(x = microsite_condition, y = cv, fill = microsite_condition)) +
  geom_col(width = 0.6, alpha = 0.8) +
  facet_grid(trait ~ surface_depth, scales = "free_y") +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A")) +
  theme_classic() +
  labs(y = "Coefficient of Variation (%)", fill = "Microsite", x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

# 2. RAW TRAIT BOXPLOTS ----
ggplot(merged_seedling_traits, aes(x = microsite_condition, y = final_mass, fill = microsite_condition)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  facet_wrap(~ surface_depth) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A")) +
  labs(y = "Seed Mass (g)", x = NULL, fill = "Microsite") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

# 3. KS DENSITY PLOT ----
ggplot(merged_seedling_traits %>% filter(microsite_condition %in% c("dung", "biocrust")), 
       aes(x = final_mass, fill = microsite_condition)) +
  geom_density(alpha = 0.4) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("KS D = ", round(ks_result_mass$statistic, 3), "\np = ", round(ks_result_mass$p.value, 3))) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A")) +
  theme_classic() + labs(x = "Seed Mass (g)", y = "Density")

# 4. CWM PLOT ----
cwm_plot_data <- cwm_tray %>%
  group_by(microsite_condition, surface_depth) %>%
  summarise(mean_cwm = mean(cwm_mass, na.rm = TRUE),
            se_cwm = sd(cwm_mass, na.rm = TRUE) / sqrt(n()), .groups = "drop")

ggplot(cwm_plot_data, aes(x = microsite_condition, y = mean_cwm, fill = microsite_condition)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_cwm - se_cwm, ymax = mean_cwm + se_cwm), width = 0.2) +
  facet_wrap(~ surface_depth) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A")) +
  labs(y = "CWM Seed Mass (g)", x = NULL) +
  theme_classic() + theme(axis.text.x = element_blank(), legend.position = "bottom")

