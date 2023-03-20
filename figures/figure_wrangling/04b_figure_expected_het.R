##%######################################################%##
#                                                          #
####                  MAKING FIGURE OF                  ####
####       EXPECTED HETEROZYGOSITY BY POPULATION        ####
#                                                          #
##%######################################################%##

# PURPOSE - this script is used to create a barplot of expected heterozygosity across 100km "population" clusters; alternatively, we may want to make plots at the region or k-level.


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(plotrix)

#affinis colors
v_affinis_colors <- c("#060200", "#F4D75D", "#D17E23")

# DATA --------------------------------------------------------------------

# df_het_rpbb <- readRDS("./analyses/analyses_output/03a02_output_df_he_ho_by_site.Rdata") %>% 
#   mutate(region = if_else(Site == "Appalachian", "Appalachian", if_else(Site == "Twin Cities", "Twin Cities", "Central")))

df_basic_rpbb <- readRDS("./analyses/analyses_output/03a02_output_basic_rpbb.Rdata")
df_sum_stats <- readRDS("./analyses/analyses_output/03a02_output_df_sum_stats_by_site.Rdata")



# WRANGLE POP LEVEL HE ----------------------------------------------------

df_he_sum <- df_basic_rpbb$Hs %>% 
  as.tibble(rownames = "locus") %>% 
  pivot_longer(!locus, names_to = "site", values_to = "he") %>% 
  group_by(site) %>% 
  summarise(mean_he = mean(he), se_he = std.error(he, na.rm = TRUE)) %>% 
  inner_join(., df_sum_stats, by = c("site" = "name")) %>% 
  mutate(region = if_else(site == "Appalachian", "Appalachian", if_else(site == "Twin Cities", "Twin Cities", "Central")))


# df_het_sum <- df_basic_rpbb$Ho %>% 
#   as.tibble(rownames = "locus") %>% 
#   pivot_longer(!locus, names_to = "site", values_to = "ho") %>% 
#   group_by(site) %>% 
#   summarise(mean_ho = mean(ho), se_ho = std.error(ho, na.rm = TRUE)) %>% 
#   inner_join(., df_he_sum, by = c("site" = "site"))

df_hs_long <- df_basic_rpbb$Hs %>% 
  as.tibble(rownames = "locus") %>% 
  pivot_longer(!locus, names_to = "site", values_to = "value") %>% 
  mutate(measure = "he")

df_ho_long <- df_basic_rpbb$Ho %>% 
  as.tibble(rownames = "locus") %>% 
  pivot_longer(!locus, names_to = "site", values_to = "value") %>% 
  mutate(measure = "ho")

df_het_sum <- df_hs_long %>% 
  bind_rows(., df_ho_long) %>% 
  group_by(site, measure) %>% 
  summarise(mean_value = mean(value), se_value = std.error(value, na.rm = TRUE)) %>% 
  inner_join(., df_sum_stats, by = c("site" = "name")) %>% 
  mutate(region = if_else(site == "Appalachian", "Appalachian", if_else(site == "Twin Cities", "Twin Cities", "Central")))

# FIGURE OF HETEROZYGOSITY BY POPULATION ----------------------------------

# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])


# He only
(p_he_site <- df_he_sum %>%
    ggplot(., aes(x = reorder(site, N), y = mean_he, fill = region))+
      geom_pointrange(aes(ymin = mean_he - se_he, ymax = mean_he + se_he, size = N), pch = 21, linewidth = 1, fatten = 1)+
      geom_text(aes(y = 1, label = N)) +
      # geom_hline(yintercept = 0.678, linetype = 2) +
      # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
      scale_fill_manual(values = v_affinis_colors[1:3])+
  scale_size(guide = "none") +
      labs(x = "Region", y = "Heterozygosity") +
      theme_classic(base_size = 15) +
      # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip(ylim = c(0,1))
  )

# He and Ho

(p_het_site <- df_het_sum %>%
    ggplot(., aes(x = reorder(site, -N), y = mean_value, fill = measure))+
    geom_pointrange(aes(ymin = mean_value - se_value, ymax = mean_value + se_value, size = N, shape = measure), linewidth = 1, fatten = 0.8, position = position_dodge(width = 0.5))+
    # geom_text(aes(y = 1, label = N)) +
    # geom_hline(yintercept = 0.678, linetype = 2) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
    scale_fill_manual(values = v_affinis_colors[1:3])+
    scale_shape_manual(values = c(21, 22)) +
    scale_size(guide = "none") +
    labs(x = "Region", y = "Heterozygosity") +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0,1))
    # coord_flip(ylim = c(0,1))
)


# SAVE OUTPUTS ------------------------------------------------------------

# For now, just saving this single figure. Might want to change things or have different display items later. Want feedback from co-authors so moving on for now.

ggsave(p_he_site, filename = "./figures/figures_output/04b_figre_he_by_site.svg")

ggsave(p_het_site, filename = "./figures/figures_output/04b_figre_heho_by_site.svg")
