####
# BIOL 499 Script 4c: visualization script
# Last Updated: March 6, 2026
# Description: making graphs for abundance, richness and evenness with above and belowground
# Notice: MUST RUN script4c_abovebelowground first 
####

library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(dplyr)

# 1. ABUNDANCE VISUALIZATION -----

# STEP 1: get letters from the updated model
emm_norm <- emmeans(model_abundance_120, ~ microsite_condition | depth_simple, type = "response")
emm_letters <- as.data.frame(cld(emm_norm, Letters = letters, adjust = "tukey")) %>%
  mutate(
    .group = trimws(.group),
    depth_simple = factor(depth_simple, levels = c("aboveground", "belowground")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

# STEP 2: calculate label positions
letter_pos <- abundance120_normalized %>%
  group_by(microsite_condition, depth_simple) %>%
  summarise(y_pos = max(seedpercm3, na.rm = TRUE) + 0.02, .groups = "drop")

emm_letters <- left_join(emm_letters, letter_pos, by = c("microsite_condition", "depth_simple"))

# STEP 3: create the plot
ggplot(abundance120_normalized, aes(x = microsite_condition, y = seedpercm3, fill = microsite_condition)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  geom_text(data = emm_letters, aes(x = microsite_condition, y = y_pos, label = .group),
            size = 4, fontface = "bold", vjust = 0) +
  facet_wrap(~ depth_simple, nrow = 1,
             labeller = labeller(depth_simple = c("aboveground" = "Aboveground (Surface)", 
                                                  "belowground" = "Belowground (0-10 cm)"))) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A"),
                    labels = c("Biocrust", "Dung", "Control")) +
  labs(x = NULL, 
       y = expression(bold(paste("Seedling Density (counts / ", cm^3, ")"))), 
       fill = "Microsite") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

# 2. RICHNESS VISUALIZATION -----

emm_rich <- emmeans(model_richness_120, ~ microsite_condition | depth_simple, type = "response")
emm_rich_letters <- as.data.frame(cld(emm_rich, Letters = letters, adjust = "tukey")) %>%
  mutate(
    .group = trimws(.group),
    depth_simple = factor(depth_simple, levels = c("aboveground", "belowground")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

letter_pos_rich <- tray_richness_120 %>%
  group_by(microsite_condition, depth_simple) %>%
  summarise(y_pos = max(richness, na.rm = TRUE) + 0.5, .groups = "drop")

emm_rich_letters <- left_join(emm_rich_letters, letter_pos_rich, by = c("microsite_condition", "depth_simple"))

ggplot(tray_richness_120, aes(x = microsite_condition, y = richness, fill = microsite_condition)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  geom_text(data = emm_rich_letters, aes(x = microsite_condition, y = y_pos, label = .group),
            size = 4, fontface = "bold", vjust = 0) +
  facet_wrap(~ depth_simple, nrow = 1,
             labeller = labeller(depth_simple = c("aboveground" = "Aboveground", "belowground" = "Belowground"))) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A")) +
  labs(x = NULL, y = "Species Richness", fill = "Microsite") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

# 3. EVENNESS VISUALIZATION -----

emm_even <- emmeans(model_evenness_120, ~ microsite_condition | depth_simple, type = "response")
emm_even_letters <- as.data.frame(cld(emm_even, Letters = letters, adjust = "tukey")) %>%
  mutate(
    .group = trimws(.group),
    depth_simple = factor(depth_simple, levels = c("aboveground", "belowground")),
    microsite_condition = factor(microsite_condition, levels = c("biocrust", "dung", "control"))
  )

letter_pos_even <- tray_evenness_120 %>%
  filter(!is.na(evenness)) %>%
  group_by(microsite_condition, depth_simple) %>%
  summarise(y_pos = max(evenness, na.rm = TRUE) + 0.05, .groups = "drop")

emm_even_letters <- left_join(emm_even_letters, letter_pos_even, by = c("microsite_condition", "depth_simple"))

ggplot(tray_evenness_120 %>% filter(!is.na(evenness)), 
       aes(x = microsite_condition, y = evenness, fill = microsite_condition)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  geom_boxplot(outlier.shape = 21, width = 0.6, alpha = 0.8) +
  geom_text(data = emm_even_letters, aes(x = microsite_condition, y = y_pos, label = .group),
            size = 4, fontface = "bold", vjust = 0) +
  facet_wrap(~ depth_simple, nrow = 1,
             labeller = labeller(depth_simple = c("aboveground" = "Aboveground", "belowground" = "Belowground"))) +
  scale_fill_manual(values = c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Pielou's Evenness (J)", fill = "Microsite") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
