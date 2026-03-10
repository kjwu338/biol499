####
# BIOL 499 Script 4b: Multivariate analysis eg. Beta Diversity, PERMANOVA, and NMDS
# Description: NMDS, PERMANOVA, Beta Dispersion
# last updated: 3/8/2026
# author: Kelly Wu
####

#NMDS VS PERMANOVA
# NMDS is a visualization tool; shows shape of data
# PERMANOVA: calculates actual distance between groups

library(dplyr)
library(ggplot2)
library(vegan)
library(tidyr)

# STEP 0: Load and clean data
gc <- read.csv("5_outputtables/gc_wcategories.csv")

#clean data (same between script a, b, c)
# seedings from germination experiment (n = 372)
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

# create species matrix wide format
sp_wide <- gc_clean %>%
  # Note: Join by site/cond/depth
  group_by(tray, site, microsite_condition, surface_depth, species) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(!is.na(species)) %>%  
  pivot_wider(names_from = species, values_from = count, values_fill = 0) 

# Metadata and Matrix
meta_perm <- sp_wide %>% 
  dplyr::select(site, microsite_condition, surface_depth)
sp_mat <- sp_wide %>% 
  dplyr::select(-site, -microsite_condition, -surface_depth) %>% 
  as.matrix()

# remove empty trays (Bray-Curtis cannot handle rows with 0 sums)
valid_rows <- rowSums(sp_mat) > 0
sp_mat     <- sp_mat[valid_rows, ]
meta_perm  <- meta_perm[valid_rows, ]

# STEP 1. MULTIVARIATE ANALYSIS -----

# a. ORDINATION CALCs
sp_mat_sqrt <- sqrt(sp_mat) # sqrt to deemphasis abundant species
bray_dist   <- vegdist(sp_mat_sqrt, method = "bray") # dissimilarity matrix

#moved from 2 dimensions to 3 dimensions to reduce stress for better match
# 3D map of relationships
nmds_gc <- metaMDS(sp_mat_sqrt, distance = "bray", k = 3, trymax = 100)
nmds_gc$stress

# b. plot with ggplot----
# NMDS scores (coordinates)
nmds_scores <- as.data.frame(scores(nmds_gc, display = "sites"))

# add 101 row metadata
plot_data <- cbind(nmds_scores, meta_perm)

# colors
cols <- c("biocrust" = "#8B9E6B", "dung" = "#C4944A", "control" = "#D4C09A")

ggplot(plot_data, aes(x = NMDS1, y = NMDS2)) +
  # add 95% Confidence Ellipses
  stat_ellipse(aes(color = microsite_condition), type = "t", level = 0.95, linewidth = 1) +
  
  # add points + alpha makes them slightly transparent
  geom_point(aes(color = microsite_condition, shape = surface_depth), 
             size = 3, alpha = 0.7) +
  
  # color and shapes
  scale_color_manual(values = cols, name = "Microsite") +
  scale_shape_manual(values = c(16, 17, 15), name = "Depth") + # 16=circle, 17=tri, 15=square
  
  # legend
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1 # NMDS should usually be square
  ) +
  
  # stress value
  labs(title = "NMDS of Seed Bank Communities",
       subtitle = paste0("Stress = ", round(nmds_gc$stress, 3)))

# STEP 2. MUltivariate stats ----

# a. Beta dispersion: selective filtering 
# is variation within dung different than biocrust?
disp_micro <- betadisper(bray_dist, meta_perm$microsite_condition)

anova(disp_micro) # p < 0.05 means one microsite is more selective/constrained
#notes: p value = 0.2353, no one microsite is more selective for species
permutest(disp_micro, pairwise = TRUE) # indicates which pair differs

# b. PERMANOVA 

#make 8 treatment groups
meta_perm$combo <- interaction(meta_perm$microsite_condition, meta_perm$surface_depth)

# are species composition different between groups?
perm_results <- adonis2(sp_mat_sqrt ~ combo, 
                        data = meta_perm, 
                        strata = meta_perm$site, # Block by site x20
                        method = "bray")
print(perm_results)

