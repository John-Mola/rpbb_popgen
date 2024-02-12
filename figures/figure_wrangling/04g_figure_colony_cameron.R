##%######################################################%##
#                                                          #
####       FIGURE 5 COMPARISON WITH CAMERON DATA        ####
#                                                          #
##%######################################################%##


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(ggsignif)



# DATA --------------------------------------------------------------------

df_affinis_east <- readRDS("./analyses/analyses_output/03d_cameron_comparison.Rdata")
df_cld_colony_compare <- readRDS("./analyses/analyses_output/03d_cameron_comparison_cld.Rdata")

# WRANGLING ---------------------------------------------------------------

df_wrangled <- df_affinis_east %>% 
  group_by(species) %>% 
  summarise(mean_detect = mean(n_detected), se_detect = std.error(n_detected)) %>% 
  full_join(., df_cld_colony_compare)

level_order <- c("affinis", "pensylvanicus", "bimaculatus", "impatiens")

p_cam <- df_wrangled %>% 
  ggplot(., aes(x = factor(species, levels = level_order), y = mean_detect)) +
  geom_pointrange(aes(ymin = mean_detect - se_detect, ymax = mean_detect + se_detect, fill = species), size = 2, linewidth = 1, shape = 21) +
  geom_text(aes(y = mean_detect + se_detect + 0.01, label = letters)) +
  
  # geom_signif(y_position = c(0.98, 0.94, 0.87), xmin = c(1, 1, 1), xmax = c(4, 3, 2), annotations = c("P = 0.009", "P = 0.076", "P = 0.992"), textsize = 6, vjust = 0) +
  theme_classic(base_size = 20) +
  labs(x = "Species", y = "Mean Colonies Detected (± SE)") +
  scale_x_discrete(labels = c("B. affinis", "B. pensylvanicus", "B. bimaculatus", "B. impatiens")) +
  theme(axis.text.x = element_text(face = "italic"), legend.position = "none") +
  # coord_cartesian(ylim = c(0,1)) +
  scale_fill_manual(values = c("#F4D75D", "black", "black", "black")) 

# save it

ggsave(p_cam, filename =  "./figures/figures_output/04g_figure_compare_to_cameron.svg", width = 8, height = 5)

