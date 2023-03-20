#RESAMPLING VOSNESENSKII TO SEE FAMILY DISTRIBUTION

library(tidyverse)

# datasets

df_sra <- read_csv("data/data_raw/misc/sra_captures_post_COLONY.csv")

# saveRDS(df_capwire_results, "./analyses/analyses_output/03c_bigsites_capwireoutput.Rdata")

df_capwire_results <- read_rds("./analyses/analyses_output/03c_bigsites_capwireoutput.Rdata")

#Do this for 10, 20, and 40 specimens. Instead of plotting as before, calculate the range of proportion detected (the mean, 5% and 95% range)


# SHUFFLING FUNCTION ------------------------------------------------------


type_shuffle_f_m = function(df, v_sites, ssize) {

  df %>% 
    filter(site %in% v_sites) %>% 
    sample_n(size = ssize) %>% 
    distinct(ClusterIndex) %>% 
    summarise(nrows = n())

}


# WRANGLING SIERRA DATASETS -----------------------------------------------



#Target sites (to keep comparison fair)
sra_sites <- c("JM01", "JM02", "JM03", "JM04", "JM05", "JM06", "JM07", "JM10", "JMX", "JMY", "forest_north", "forest_south", "mdw_north", "mdw_south", "north_east")

#datasets
df_sra15_vos <- df_sra %>% filter(species == "vosnesenskii", year == 2015)
df_sra15_bif <- df_sra %>% filter(species == "bifarius", year == 2015)
df_sra18_bif <- df_sra %>% filter(species == "bifarius", year == 2018)



# SHUFFLE SAMPLE SIZE OF 9 ------------------------------------------------

#shuffle and match for sample size of 9 individuals
df_shuf_vos15_n9 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 9))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n9 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 9))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n9 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 9))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n9 <- bind_rows(df_shuf_bif15_n9, df_shuf_bif18_n9, df_shuf_vos15_n9) %>% mutate(raw_detected = nrows/9)

#shuffle and match for sample size of 10 individuals
df_shuf_vos15_n10 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 10))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n10 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 10))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n10 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 10))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n10 <- bind_rows(df_shuf_bif15_n10, df_shuf_bif18_n10, df_shuf_vos15_n10) %>% mutate(raw_detected = nrows/10)

#shuffle and match for sample size of 11 individuals
df_shuf_vos15_n11 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 11))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n11 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 11))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n11 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 11))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n11 <- bind_rows(df_shuf_bif15_n11, df_shuf_bif18_n11, df_shuf_vos15_n11) %>% mutate(raw_detected = nrows/11)

#shuffle and match for sample size of 12 individuals
df_shuf_vos15_n12 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 12))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n12 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 12))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n12 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 12))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n12 <- bind_rows(df_shuf_bif15_n12, df_shuf_bif18_n12, df_shuf_vos15_n12) %>% mutate(raw_detected = nrows/12)

#shuffle and match for sample size of 13 individuals
df_shuf_vos15_n13 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 13))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n13 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 13))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n13 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 13))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n13 <- bind_rows(df_shuf_bif15_n13, df_shuf_bif18_n13, df_shuf_vos15_n13) %>% mutate(raw_detected = nrows/13)

#shuffle and match for sample size of 15 individuals
df_shuf_vos15_n15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 15))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 15))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 15))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n15 <- bind_rows(df_shuf_bif15_n15, df_shuf_bif18_n15, df_shuf_vos15_n15) %>% mutate(raw_detected = nrows/15)

#shuffle and match for sample size of 18 individuals
df_shuf_vos15_n18 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 18))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n18 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 18))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n18 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 18))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n18 <- bind_rows(df_shuf_bif15_n18, df_shuf_bif18_n18, df_shuf_vos15_n18) %>% mutate(raw_detected = nrows/18)

#shuffle and match for sample size of 19 individuals
df_shuf_vos15_n19 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 19))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n19 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 19))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n19 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 19))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n19 <- bind_rows(df_shuf_bif15_n19, df_shuf_bif18_n19, df_shuf_vos15_n19) %>% mutate(raw_detected = nrows/19)

#shuffle and match for sample size of 42 individuals
df_shuf_vos15_n42 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites, 42))) %>% bind_rows() %>% mutate(set = "vos15")
df_shuf_bif15_n42 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites, 42))) %>% bind_rows() %>% mutate(set = "bif15")
df_shuf_bif18_n42 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites, 42))) %>% bind_rows() %>% mutate(set = "bif18")
df_all_shuf_n42 <- bind_rows(df_shuf_bif15_n42, df_shuf_bif18_n42, df_shuf_vos15_n42) %>% mutate(raw_detected = nrows/42)



# COMPARISONS TO OBSERVED VALUES ------------------------------------------

#my brain hurts and I cannot think of how to do this even remotely elegantly (just like above...lol)

#Turtle21, 42 individuals

df_compare_Turtle21 <- df_all_shuf_n42 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Turtle21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Turtle21")

#MNZoo21, 19 individuals

df_compare_MNZoo21 <- df_all_shuf_n19 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "MNZoo21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "MNZoo21")

#MNZoo20, 18 individuals

df_compare_MNZoo20 <- df_all_shuf_n18 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "MNZoo20"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "MNZoo20")

#Nachusa21, 15 individuals

df_compare_Nachusa21 <- df_all_shuf_n15 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Nachusa21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Nachusa21")

#Foxglove21, 13 individuals

df_compare_Foxglove21 <- df_all_shuf_n13 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Foxglove21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Foxglove21")

#Hamline21, 12 individuals

df_compare_Hamline21 <- df_all_shuf_n12 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Hamline21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Hamline21")

#Cherokee20, 11 individuals

df_compare_Cherokee20 <- df_all_shuf_n11 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Cherokee20"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Cherokee20")

#Illiniwek21, 11 individuals

df_compare_Illiniwek21 <- df_all_shuf_n11 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Illiniwek21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Illiniwek21")

#Elaine20, 10 individuals

df_compare_Elaine20 <- df_all_shuf_n10 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Elaine20"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Elaine20")

#Appalachian21, 9 individuals

df_compare_Appalachian21 <- df_all_shuf_n9 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "Appalachian21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "Appalachian21")

#CapSprings21, 9 individuals

df_compare_CapSprings21 <- df_all_shuf_n9 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "CapSprings21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "CapSprings21")

#SeedSavers21, 9 individuals

df_compare_SeedSavers21 <- df_all_shuf_n9 %>% filter(nrows <= pull(filter(df_capwire_results, custom_name == "SeedSavers21"), n_colonies)) %>% group_by(set) %>% tally() %>% pivot_wider(names_from = set, values_from = n) %>% mutate(custom_name = "SeedSavers21")

## combine them all - a few get "lost" because they are all 0s (meaning the observed value was less than all reshuffled values)
df_compare_all_reshuffles <- bind_rows(df_compare_Appalachian21, df_compare_CapSprings21, df_compare_Cherokee20, df_compare_Elaine20, df_compare_Foxglove21, df_compare_Hamline21, df_compare_Illiniwek21, df_compare_MNZoo20, df_compare_MNZoo21, df_compare_Nachusa21, df_compare_SeedSavers21, df_compare_Turtle21)

# merge with the capwire results table

df_joined_capwire_compare <- full_join(df_capwire_results, df_compare_all_reshuffles) %>% replace_na(replace = list(bif15 = 0, bif18 = 0, vos15 = 0)) %>% 
  mutate(prob_bif15 = bif15/1000,
         prob_bif18 = bif18/1000,
         prob_vos15 = vos15/1000) %>% 
  dplyr::select(-bif15, -vos15, -bif18) %>% arrange(prop_detected)

df_joined_capwire_compare %>% arrange(prop_detected)


# SAVE OUTPUT -------------------------------------------------------------

saveRDS(df_joined_capwire_compare, "./analyses/analyses_output/03d_joined_comparison_to_sierras.Rdata")

# saving examples from 9, 11, 19, and 42 datasets as examples for plotting. Later, may want to use all of them. 

saveRDS(df_all_shuf_n9, "./analyses/analyses_output/03d_shuffle_results9.Rdata")
saveRDS(df_all_shuf_n11, "./analyses/analyses_output/03d_shuffle_results11.Rdata")
saveRDS(df_all_shuf_n19, "./analyses/analyses_output/03d_shuffle_results19.Rdata")
saveRDS(df_all_shuf_n42, "./analyses/analyses_output/03d_shuffle_results42.Rdata")








# readRDS("./analyses/analyses_output/03d_joined_comparison_to_sierras.Rdata")
# 
# ## look at output
# # 
# df_all_shuf_n10 %>%
#   group_by(set) %>%
#   summarise(q01 = quantile(raw_detected, 0.01),
#             q05 = quantile(raw_detected, 0.05),
#             q95 = quantile(raw_detected, 0.95),
#             q99 = quantile(raw_detected, 0.99))
# # 
# df_all_shuf_n19 %>%
#   group_by(set) %>%
#   summarise(q01 = quantile(raw_detected, 0.01),
#             q05 = quantile(raw_detected, 0.05),
#             q95 = quantile(raw_detected, 0.95),
#             q99 = quantile(raw_detected, 0.99))
# 
# df_all_shuf_n30 %>%
#   group_by(set) %>% 
#   summarise(q01 = quantile(raw_detected, 0.01),
#             q05 = quantile(raw_detected, 0.05),
#             q95 = quantile(raw_detected, 0.95),
#             q99 = quantile(raw_detected, 0.99))
