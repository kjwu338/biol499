####
# BIOL 499 Script 4a: visualization script
# Last Updated: March 5, 2026
# Description: making graphs for abundance, richness and evenness
# Notice: MUST RUN script4a_v2 first 
####

library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(dplyr)

# Abundance visualization -----

# STEP 1: get letters from the poisson count model
emm_counts <- emmeans(model_abundance, ~ treatment, type = "response")
emm_letters <- as.data.frame(cld(emm_counts, Letters = letters, adjust = "tukey")) %>%
  mutate(.group = trimws(.group)) %>%
  # IMPORTANT: Turn "biocrust_surface" back into two columns for the plot
  separate(treatment, into = c("microsite_condition", "surface_depth"), sep = "_") %>%
  mutate(
    surface_depth = factor(surface_depth, levels = c("surface", "0-5", "5-10")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

# STEP 2: calculate label positions
letter_pos <- tray_normal %>%
  group_by(microsite_condition, surface_depth) %>%
  summarise(y_pos = max(total_seeds, na.rm = TRUE) + 2, .groups = "drop") 

# Join letters to their positions
emm_letters <- left_join(emm_letters, letter_pos, by = c("microsite_condition", "surface_depth"))

# STEP 3: create the plot 
ggplot(tray_normal, aes(x = microsite_condition, y = total_seeds, fill = microsite_condition)) +
  # Whiskers
  stat_boxplot(geom = "errorbar", width = 0.25) +
  # Main boxplot (outliers shown as open circles)
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  # Significance Letters
  geom_text(data = emm_letters, 
            aes(x = microsite_condition, y = y_pos, label = .group),
            size = 4.5, fontface = "bold", vjust = 0) +
  # Faceting by Depth
  facet_wrap(~ surface_depth, nrow = 1,
             labeller = labeller(surface_depth = c("surface" = "Surface", 
                                                   "0-5" = "0-5 cm", 
                                                   "5-10" = "5-10 cm"))) +
  # Professional Colors
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", 
                               "dung" = "#C4944A", 
                               "control" = "#D4C09A"),
                    labels = c("Biocrust", "Dung", "Control")) +
  # Labels
  labs(x = NULL, 
       y = "Seedling Count (per tray)", 
       fill = "Microsite") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

# Richness visualization -----

# STEP 1: get sig letters for treatment model
# = combined 'treatment' variable as the predictor
emm_rich <- emmeans(model_richness, ~ treatment, type = "response")
emm_rich_letters <- as.data.frame(cld(emm_rich, Letters = letters, adjust = "tukey")) %>%
  mutate(.group = trimws(.group)) %>%
  # SPLIT the treatment back into two columns for the plot facets and fill
  separate(treatment, into = c("microsite_condition", "surface_depth"), sep = "_") %>%
  mutate(
    surface_depth = factor(surface_depth, levels = c("surface", "0-5", "5-10")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

# STEP 2: calculate letter positions
letter_pos_rich <- tray_richness %>%
  group_by(microsite_condition, surface_depth) %>%
  summarise(y_pos = max(richness, na.rm = TRUE) + 0.8, .groups = "drop") # 0.8 to clear whiskers

emm_rich_letters <- left_join(emm_rich_letters, letter_pos_rich, 
                              by = c("microsite_condition", "surface_depth"))

# STEP 3: plot richness
ggplot(tray_richness, aes(x = microsite_condition, y = richness, fill = microsite_condition)) +
  # whisker bars
  stat_boxplot(geom = "errorbar", width = 0.25) +
  # boxplot
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  # letters
  geom_text(data = emm_rich_letters, 
            aes(x = microsite_condition, y = y_pos, label = .group),
            size = 4.5, fontface = "bold", vjust = 0) +
  # split data into subplots
  facet_wrap(~ surface_depth, nrow = 1,
             labeller = labeller(surface_depth = c(
               "surface" = "Surface",
               "0-5"     = "0-5 cm",
               "5-10"    = "5-10 cm"))) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B",
                               "dung"      = "#C4944A",
                               "control"   = "#D4C09A"),
                    labels = c("Biocrust", "Dung", "Control")) +
  labs(x = NULL, y = "Number of Species", fill = "Microsite") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

# EVENNESS------

# STEP 1: Get significance letters from the treatment-based model
emm_even <- emmeans(model_evenness, ~ treatment, type = "response")
emm_even_letters <- as.data.frame(cld(emm_even, Letters = letters, adjust = "tukey")) %>%
  mutate(.group = trimws(.group)) %>%
  # split the treatment back into microsite and depth
  separate(treatment, into = c("microsite_condition", "surface_depth"), sep = "_") %>%
  mutate(
    surface_depth = factor(surface_depth, levels = c("surface", "0-5", "5-10")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

# STEP 2: calculate letter positions
letter_pos_even <- tray_evenness %>%
  filter(!is.na(evenness)) %>%
  group_by(microsite_condition, surface_depth) %>%
  # position letters above max value
  summarise(y_pos = max(evenness, na.rm = TRUE) + 0.05, .groups = "drop") 

emm_even_letters <- left_join(emm_even_letters, letter_pos_even, 
                              by = c("microsite_condition", "surface_depth"))

# STEP 3: plot evenness
ggplot(tray_evenness %>% filter(!is.na(evenness)), 
       aes(x = microsite_condition, y = evenness, fill = microsite_condition)) +
  # Whiskers
  stat_boxplot(geom = "errorbar", width = 0.25) +
  # Main Boxplot
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  # Letters
  geom_text(data = emm_even_letters, 
            aes(x = microsite_condition, y = y_pos, label = .group),
            size = 4.5, fontface = "bold", vjust = 0) +
  # Faceting
  facet_wrap(~ surface_depth, nrow = 1,
             labeller = labeller(surface_depth = c(
               "surface" = "Surface", 
               "0-5" = "0-5 cm", 
               "5-10" = "5-10 cm"))) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", 
                               "dung" = "#C4944A", 
                               "control" = "#D4C09A"),
                    labels = c("Biocrust", "Dung", "Control")) +
  # keeps scale clean
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Pielou's Evenness (J)", fill = "Microsite") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

# Pielou's values go from 0-1. 
# Close to 1 means high evenness, ie there is no one dominant species
# Close to 0 means dominance by a single species