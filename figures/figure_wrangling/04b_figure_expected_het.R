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
library(patchwork)
library(ggsignif)
library(adegenet)
library(pegas)
library(poppr)
library(hierfstat)
#affinis colors
v_affinis_colors <- c("#060200", "#F4D75D", "#D17E23")

# DATA --------------------------------------------------------------------

# df_het_rpbb <- readRDS("./analyses/analyses_output/03a02_output_df_he_ho_by_site.Rdata") %>% 
#   mutate(region = if_else(Site == "Appalachian", "Appalachian", if_else(Site == "Twin Cities", "Twin Cities", "Central")))

df_basic_rpbb <- readRDS("./analyses/analyses_output/03a02_output_basic_rpbb.Rdata")
df_sum_stats <- readRDS("./analyses/analyses_output/03a02_output_df_sum_stats_by_site.Rdata")


df_centroids <- readRDS("./data/data_output/output_01c_df_cluster100_centroids.Rdata") %>% 
  arrange(longitude_center)

gen_rpbb_flt <- readRDS("./analyses/analyses_output/03a01_gen_rpbb_flt.Rdata")


# loading cld outputs

df_cld_pop_ar <- read_rds("./analyses/analyses_output/03a02_output_df_cld_pop_ar.Rdata")
df_cld_pop_he <- read_rds("./analyses/analyses_output/03a02_output_df_cld_pop_he.Rdata")
df_cld_pop_fis <- read_rds("./analyses/analyses_output/03a02_output_df_cld_pop_fis.Rdata")

df_cld_reg_ar <- read_rds("./analyses/analyses_output/03a02_output_df_cld_reg_ar.Rdata")
df_cld_reg_he <- read_rds("./analyses/analyses_output/03a02_output_df_cld_reg_he.Rdata")
df_cld_reg_fis <- read_rds("./analyses/analyses_output/03a02_output_df_cld_reg_fis.Rdata")


# WRANGLE POP LEVEL HE ----------------------------------------------------

df_he_sum <- df_basic_rpbb$Hs %>% 
  as.tibble(rownames = "locus") %>% 
  pivot_longer(!locus, names_to = "site", values_to = "he") %>% 
  group_by(site) %>% 
  summarise(mean_he = mean(he), se_he = std.error(he, na.rm = TRUE)) %>% 
  inner_join(., df_sum_stats, by = c("site" = "name")) %>% 
  mutate(region = if_else(site == "Appalachian", "Appalachian", if_else(site == "Twin Cities", "Twin Cities", "Central"))) %>% 
  full_join(., df_centroids, by = c("site" = "named_cluster100")) %>% 
  full_join(., df_cld_pop_he)



# WRANGLE FIS -----------------------------------------------------------


df_fis_sum <- df_basic_rpbb$Fis %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "site", values_to = "fis") %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  group_by(site) %>% 
  summarise(mean_fis = mean(fis, na.rm = TRUE), se_fis = std.error(fis, na.rm = TRUE)) %>% 
  inner_join(., df_sum_stats, by = c("site" = "name")) %>% 
  mutate(region = if_else(site == "Appalachian", "Appalachian", if_else(site == "Twin Cities", "Twin Cities", "Central"))) %>% 
  full_join(., df_centroids, by = c("site" = "named_cluster100"))  %>% 
  full_join(., df_cld_pop_fis)



# WRANGLE AR ------------------------------------------------------------


df_ar_sum <- allelic.richness(genind2hierfstat(gen_rpbb_flt))$Ar %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "site", values_to = "richness") %>% 
  group_by(site) %>% 
  summarise(mean_ar = mean(richness, na.rm = TRUE), se_ar = std.error(richness, na.rm = TRUE)) %>% 
  inner_join(., df_sum_stats, by = c("site" = "name")) %>% 
  mutate(region = if_else(site == "Appalachian", "Appalachian", if_else(site == "Twin Cities", "Twin Cities", "Central"))) %>% 
  full_join(., df_centroids, by = c("site" = "named_cluster100"))  %>% 
  full_join(., df_cld_pop_ar)

# df_het_sum <- df_basic_rpbb$Ho %>% 
#   as.tibble(rownames = "locus") %>% 
#   pivot_longer(!locus, names_to = "site", values_to = "ho") %>% 
#   group_by(site) %>% 
#   summarise(mean_ho = mean(ho), se_ho = std.error(ho, na.rm = TRUE)) %>% 
#   inner_join(., df_he_sum, by = c("site" = "site"))

# df_hs_long <- df_basic_rpbb$Hs %>% 
#   as.tibble(rownames = "locus") %>% 
#   pivot_longer(!locus, names_to = "site", values_to = "value") %>% 
#   mutate(measure = "he")
# 
# df_ho_long <- df_basic_rpbb$Ho %>% 
#   as.tibble(rownames = "locus") %>% 
#   pivot_longer(!locus, names_to = "site", values_to = "value") %>% 
#   mutate(measure = "ho")
# 
# df_het_sum <- df_hs_long %>% 
#   bind_rows(., df_ho_long) %>% 
#   group_by(site, measure) %>% 
#   summarise(mean_value = mean(value), se_value = std.error(value, na.rm = TRUE)) %>% 
#   inner_join(., df_sum_stats, by = c("site" = "name")) %>% 
#   mutate(region = if_else(site == "Appalachian", "Appalachian", if_else(site == "Twin Cities", "Twin Cities", "Central")))

# FIGURE OF HETEROZYGOSITY BY POPULATION ----------------------------------

# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])



# POP LEVEL HE --------------------------------------------------------



(p_he_site <- df_he_sum %>%
    ggplot(., aes(x = reorder(site, longitude_center), y = mean_he, fill = region))+
      geom_pointrange(aes(ymin = mean_he - se_he, ymax = mean_he + se_he, size = N), pch = 21, linewidth = 1, fatten = 1)+
   geom_text(aes(y = mean_he + se_he + 0.1, label = letters)) +
      #geom_text(aes(y = 1, label = N)) +
      # geom_hline(yintercept = 0.678, linetype = 2) +
      # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
      scale_fill_manual(values = c( "#D17E23", "#F4D75D", "#060200"))+
  scale_size(guide = "none") +
      labs(x = "Region", y = "Heterozygosity") +
      theme_classic(base_size = 15) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
   coord_cartesian(ylim = c(0.1, 1))
    # coord_flip(ylim = c(0,1))
  )



# POP LEVEL FIS ------------------------------------------------------------




# fis 
(p_fis_site <- df_fis_sum %>%
    ggplot(., aes(x = reorder(site, longitude_center), y = mean_fis, fill = region))+
    geom_pointrange(aes(ymin = mean_fis - se_fis, ymax = mean_fis + se_fis, size = N), pch = 21, linewidth = 1, fatten = 1)+
   geom_text(aes(y = mean_fis + se_fis + 0.07, label = letters)) +
   
    #geom_text(aes(y = 1, label = N)) +
    # geom_hline(yintercept = 0.678, linetype = 2) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
    scale_fill_manual(values = c( "#D17E23", "#F4D75D", "#060200"))+
    scale_size(guide = "none") +
    labs(x = "Region", y ="") +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
   coord_cartesian(ylim = c(-0.2, 0.7))
 
  # coord_flip(ylim = c(0,1))
)





# POP LEVEL AR ------------------------------------------------------------





# ar 
(p_ar_site <- df_ar_sum %>%
    ggplot(., aes(x = reorder(site, longitude_center), y = mean_ar, fill = region))+
    geom_pointrange(aes(ymin = mean_ar - se_ar, ymax = mean_ar + se_ar, size = N), pch = 21, linewidth = 1, fatten = 1)+
   geom_text(aes(y = mean_ar + se_ar + 0.1, label = letters)) +
   
    #geom_text(aes(y = 1, label = N)) +
    # geom_hline(yintercept = 0.678, linetype = 2) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
    scale_fill_manual(values = c( "#D17E23", "#F4D75D", "#060200"))+
    scale_size(guide = "none") +
    labs(x = "Region", y = "AR") +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # coord_flip(ylim = c(0,1))
)


# He and Ho

# (p_het_site <- df_het_sum %>%
#     ggplot(., aes(x = reorder(site, -N), y = mean_value, fill = measure))+
#     geom_pointrange(aes(ymin = mean_value - se_value, ymax = mean_value + se_value, size = N, shape = measure), linewidth = 1, fatten = 0.8, position = position_dodge(width = 0.5))+
#     # geom_text(aes(y = 1, label = N)) +
#     # geom_hline(yintercept = 0.678, linetype = 2) +
#     # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
#     scale_fill_manual(values = v_affinis_colors[1:3])+
#     scale_shape_manual(values = c(21, 22)) +
#     scale_size(guide = "none") +
#     labs(x = "Region", y = "Heterozygosity") +
#     theme_classic(base_size = 15) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     coord_cartesian(ylim = c(0,1))
#     # coord_flip(ylim = c(0,1))
# )


# SAVE OUTPUTS ------------------------------------------------------------

# For now, just saving this single figure. Might want to change things or have different display items later. Want feedback from co-authors so moving on for now.

# ggsave(p_he_site, filename = "./figures/figures_output/04b_figre_he_by_site.svg")

# ggsave(p_het_site, filename = "./figures/figures_output/04b_figre_heho_by_site.svg")




# RENEWING ----------------------------------------------------------------

#AR, He, Fis

df_ar_reg <- read_rds("./analyses/analyses_output/03a02_output_df_ar_reg.Rdata")
df_hs_reg <- read_rds("./analyses/analyses_output/03a02_output_df_hs_reg.Rdata")
df_fis_reg <- read_rds("./analyses/analyses_output/03a02_output_df_fis_reg.Rdata")

level_order <- c("north west", "central", "appalachian")

# summarizing AR

df_reg_ar_sum <- df_ar_reg %>% 
  group_by(region) %>% 
  summarise(mean_ar = mean(ar), se_ar = std.error(ar, na.rm = TRUE))  %>% 
  full_join(., df_cld_reg_ar)

# summarizing He

df_reg_he_sum <- df_hs_reg %>% 
  group_by(region) %>% 
  summarise(mean_he = mean(he), se_he = std.error(he, na.rm = TRUE))   %>% 
  full_join(., df_cld_reg_he)

# summarizing fis

df_reg_fis_sum <- df_fis_reg %>% 
  group_by(region) %>% 
  summarise(mean_fis = mean(fis), se_fis = std.error(fis, na.rm = TRUE))   %>% 
  full_join(., df_cld_reg_fis)

# # plot df_ar
# 
# p_ar <- df_ar_reg %>% 
#   ggplot(., aes(x = factor(region, level = level_order), y = ar)) +
#   geom_boxplot(aes(fill = region), width = 0.25, alpha = 0.8) +
#   geom_point(aes(fill = region), shape = 21, size = 3, alpha = 0.5) +
#   geom_signif(y_position = c(23, 20), xmin = c(1, 2), xmax = c(3, 3), annotations = c("·", "*"), textsize = 8, vjust = 0.5) +
#   theme_classic(base_size = 15) +
#   scale_x_discrete(labels = c("Twin Cities", "Central", "Appalachian")) +
#   labs(x = "Genetic Cluster", y = "Rarefied \n Allelic Richness") +
#   scale_fill_manual(values = v_affinis_colors[3:1]) +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()) +
#   coord_cartesian(ylim = c(0, 23))
 

# ar 
(p_ar <- df_reg_ar_sum %>%
    ggplot(., aes(x = factor(region, level = level_order), y = mean_ar, fill = region))+
    geom_pointrange(aes(ymin = mean_ar - se_ar, ymax = mean_ar + se_ar), pch = 21, linewidth = 1, fatten = 1, size = 4)+
    geom_text(aes(y = mean_ar + se_ar + 0.25, label = letters)) +
    
    #geom_text(aes(y = 1, label = N)) +
    # geom_hline(yintercept = 0.678, linetype = 2) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
    scale_fill_manual(values = c( "#D17E23", "#F4D75D", "#060200"))+
    scale_size(guide = "none") +
    labs(x = "Cluster", y = "Rarefied Allelic \n Richness (AR)") +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank())
  # coord_flip(ylim = c(0,1))
)



 
# # plot df_hs
# 
# p_he <- df_hs_reg %>% 
#   ggplot(., aes(x = factor(region, level = level_order), y = he)) +
#   geom_boxplot(aes(fill = region), width = 0.25, alpha = 0.8, outlier.colour = "white") +
#   geom_point(aes(fill = region), shape = 21, size = 3, alpha = 0.5) +
#   theme_classic(base_size = 15) +
#   scale_x_discrete(labels = c("Twin Cities", "Central", "Appalachian")) +
#   labs(x = "Genetic Cluster", y = "Expected \n Heterozygosity") +
#   scale_fill_manual(values = v_affinis_colors[3:1]) +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()) +
#   coord_cartesian(ylim = c(0.1, 1))

(p_he <- df_reg_he_sum %>%
    ggplot(., aes(x = factor(region, level = level_order), y = mean_he, fill = region))+
    geom_pointrange(aes(ymin = mean_he - se_he, ymax = mean_he + se_he), pch = 21, linewidth = 1, fatten = 1, size = 4)+
    geom_text(aes(y = mean_he + se_he + 0.025, label = letters)) +
    
    #geom_text(aes(y = 1, label = N)) +
    # geom_hline(yintercept = 0.678, linetype = 2) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
    scale_fill_manual(values = c( "#D17E23", "#F4D75D", "#060200"))+
    scale_size(guide = "none") +
    labs(x = "Cluster", y = expression("Expected \n Heterozygosity (H"[E]*")")) +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank()) +
      coord_cartesian(ylim = c(0.1, 1))
    
  # coord_flip(ylim = c(0,1))
)



# # plot df_fis
# 
# 
# p_fis <- df_fis_reg %>% 
#   ggplot(., aes(x = factor(region, level = level_order), y = fis)) +
#   geom_boxplot(aes(fill = region), width = 0.25, alpha = 0.8, outlier.colour = "white") +
#   geom_point(aes(fill = region), shape = 21, size = 3, alpha = 0.5) +
#   geom_signif(y_position = c(0.55, 0.45), xmin = c(1, 2), xmax = c(3, 3), annotations = c("***", "***"), textsize = 8, vjust = 0.5) +
#   theme_classic(base_size = 15) +
#   scale_x_discrete(labels = c("Twin Cities", "Central", "Appalachian")) +
#   labs(x = "Genetic Cluster", y = expression(F[IS])) +
#   scale_fill_manual(values = v_affinis_colors[3:1]) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   coord_cartesian(ylim = c(-0.2, 0.7))

(p_fis <- df_reg_fis_sum %>%
    ggplot(., aes(x = factor(region, level = level_order), y = mean_fis, fill = region))+
    geom_pointrange(aes(ymin = mean_fis - se_fis, ymax = mean_fis + se_fis), pch = 21, linewidth = 1, fatten = 1, size = 4)+
    geom_text(aes(y = mean_fis + se_fis + 0.025, label = letters)) +
    
    #geom_text(aes(y = 1, label = N)) +
    # geom_hline(yintercept = 0.678, linetype = 2) +
    # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
    scale_fill_manual(values = c( "#D17E23", "#F4D75D", "#060200"))+
    scale_size(guide = "none") +
    labs(x = "Genetic Cluster", y = expression("Coefficient of \n Inbreeding (F"[IS]*")")) +
    theme_classic(base_size = 15) +
    scale_x_discrete(labels = c("Twin Cities", "Central", "Appalachian")) +
    
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_cartesian(ylim = c(-0.2, 0.7))
  
  # coord_flip(ylim = c(0,1))
)


# combine them

(p_combined_region <-  (p_ar / p_he / p_fis) + plot_annotation(tag_levels = 'A'))

# super combo

(p_combined_site <- ((p_ar_site + theme(axis.text.x = element_blank(), axis.title = element_blank())) / (p_he_site + theme(axis.text.x = element_blank(), axis.title = element_blank())) / (p_fis_site + labs(x = "100-km Population"))) + plot_layout(guides = "collect"))


p_combined_site_region <- (p_combined_region | p_combined_site) + plot_layout(widths = c(1,3)) + plot_annotation(tag_levels = 'A') & theme(legend.position = "none")



# SAVE IT -----------------------------------------------------------------

ggsave(p_combined_site_region, filename = "./figures/figures_output/04b_figure_regional_diversity.svg")

