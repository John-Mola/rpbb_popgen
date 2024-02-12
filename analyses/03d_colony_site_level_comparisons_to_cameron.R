##%######################################################%##
#                                                          #
#### COMPARING B AFFINIS UNIQUE COLONIES PER INDIVIDUAL ####
####         TO RESULTS FROM CAMERON ET AL 2011         ####
#                                                          #
##%######################################################%##

library(tidyverse)
library(report)
library(plotrix)
library(multcomp)




# DATA --------------------------------------------------------------------


df_rpbb_colonizer <- read_csv("./analyses/outputs_colony/r_colonizer/rpbb_batch_colonizeR.csv")


df_cameron <- read_csv("./data/data_raw/cameron2011_colony_counts.csv") %>% 
  dplyr::select(-subpopulation, -latitude, -Longitude, -n_colony_detected) %>% 
  rename(site = locality) %>% 
  mutate(n_detected = n_colonies/n_individuals) %>% 
  filter(n_individuals > 1)

# WRANGLING ---------------------------------------------------------------


df_ncol_site <- df_rpbb_colonizer %>% 
  group_by(cluster10, year) %>% 
  mutate(n_individuals = n()) %>% 
  group_by(cluster10, year, family_index, n_individuals) %>% 
  tally() %>% 
  group_by(cluster10, year, n_individuals) %>% 
  tally(name = "n_colonies") %>% 
  filter(n_individuals > 1) %>% 
  mutate(species = "affinis") %>% 
  ungroup() %>% 
  dplyr::select(species, cluster10, n_individuals, n_colonies, year) %>% 
  mutate(n_detected = n_colonies/n_individuals,
         site = paste0(year, "_", cluster10)) %>% 
  dplyr::select(-year)



df_affinis_cameron <- bind_rows(df_ncol_site, df_cameron)


# only eastern species from cameron and affinis

df_affinis_east <- df_affinis_cameron %>% 
  filter(species %in% c("impatiens", "pensylvanicus", "bimaculatus", "affinis"))


glm_combined <- df_affinis_east %>% 
  mutate(species = factor(species, levels = c("pensylvanicus", "bimaculatus", "impatiens", "affinis"))) %>% 
  glm(data =., n_detected ~ species, family = "quasibinomial")


summary(glm_combined)       
posthoc_glm_combined <- glht(glm_combined, linfct = mcp(species = "Tukey"))
summary(posthoc_glm_combined)
cld(posthoc_glm_combined)
df_cld_colony_compare <- broom::tidy(cld(posthoc_glm_combined))



saveRDS(df_affinis_east, file = "./analyses/analyses_output/03d_cameron_comparison.Rdata")
saveRDS(df_cld_colony_compare, file = "./analyses/analyses_output/03d_cameron_comparison_cld.Rdata")




