##%######################################################%##
#                                                          #
####            DIPLOID MALES COUNTING USING            ####
####            NOMURA AND TANIGUCHI METHOD             ####
#                                                          #
##%######################################################%##

#Nomura, T., and Y. Taniguchi. 2022. Simple method for combining multiple-loci marker genotypes to estimate diploid male proportion, with an application to a threatened bumble bee population in Japan. Insectes Sociaux.

# This script is intended to yield an estimate of the number of diploid males in the population, corrected for level of homozygosity across the population. The intention is that by just simply counting the number of males with at least one (or two, or three) heterozygous locus we may be underestimating the true frequency of diploid males in the population. This would occur because if the population exhibits low heterozygosity in general, we may have some number of females that are homozygous anyway...so we might expect that some diploid males are homozygous, but that is unobservable because of the way ploidy is assessed using MSATs


# PACKAGES ----------------------------------------------------------------

library(tidyverse)


# DATA --------------------------------------------------------------------

#dataset of genotypes
# this is set up for source in the markdown document, it will not run if you don't update the directory depth
df_merged <- read_rds("./data/data_output/output_01d_merged_genotypes.Rdata")



# FINDING HET COUNT FOR ALL -----------------------------------------------

df_genos <- df_merged %>% 
  dplyr::select(internal_barcode, 28:53)

# map(c("B124", "BTERN01", "BT28", "BT10", "B96", "BT30", "BTMS0081", "BTMS0066", "BTMS0083", "B126", "BTMS0062", "BTERN02", "BTMS0086", "BL13", "BTMS0059"),

#This works but I don't really understand how

#From: https://stackoverflow.com/questions/49838191/r-iterate-over-pairs-of-columns-in-dataframe-by-regex-match


df_matched <- purrr::map(c("btern01", "bt28", "b96", "bt30", "btms0081", "btms0066", "btms0083", "b126", "btms0062", "btern02", "btms0086", "bl13", "btms0059"), ~ df_genos %>%
                    mutate(
                      !!(paste0(.x, "_ans", collapse = "")) :=
                        UQ(rlang::sym(paste0(.x, "_1", collapse = ""))) == UQ(rlang::sym(paste0(.x, "_2", collapse = ""))) )) %>%
  reduce(., left_join) %>%
  dplyr::select(internal_barcode, contains("_ans")) %>%
  rowwise() %>%
  mutate(het_count = sum(c_across(all_of(contains("_ans"))) == "FALSE", na.rm = TRUE))

df_matched %>% group_by(het_count) %>% tally()

df_het_matched <- inner_join(df_merged, dplyr::select(df_matched, internal_barcode, het_count)) %>% 
  filter(loci_w_data >= 10)


# CALCULATING FREQUENCIES FOR ALL SPECIMENS AND THEN SPECIFIC SITES -------

# Only Appalachians, Milwaukee, and Twin Cities have more than 6 males with 26, 14, and 50. 


# #freq at 1 more more loci
# 
# (df_1het <- df_het_matched %>% 
#   mutate(is_het = if_else(het_count >= 1, "yes", "no")) %>% 
#   count(sex, is_het))
# 
# (v_1het <- df_1het %>% 
#   filter(is_het == "yes") %>% 
#   mutate(prop_diploid = n/sum(n)) %>% 
#   filter(sex == "male") %>% 
#   pull(prop_diploid))
# 
# #freq at 2 more more loci
# 
# 
# (df_2het <- df_het_matched %>% 
#   mutate(is_het = if_else(het_count >= 2, "yes", "no")) %>% 
#   count(sex, is_het))
# 
# (v_2het <- df_2het %>% 
#   filter(is_het == "yes") %>% 
#   mutate(prop_diploid = n/sum(n)) %>% 
#   filter(sex == "male") %>% 
#   pull(prop_diploid))
# 
# 
# # REPEATED AT SITE LEVEL FOR MN ZOO ---------------------------------------
# 
# #for mnzoo at level 1
# (df_1het_mnzoo <- df_het_matched %>% 
#    filter(str_detect(site, "Minnesota Zoo")) %>% 
#    mutate(is_het = if_else(het_count >= 1, "yes", "no")) %>% 
#    count(sex, is_het))
#   
# (v_1het_mnzoo <- df_1het_mnzoo %>% 
#   filter(is_het == "yes") %>% 
#   mutate(prop_diploid = n/sum(n)) %>% 
#   filter(sex == "male") %>% 
#   pull(prop_diploid))
# 
# #for mnzoo at level 2
# (df_2het_mnzoo <- df_het_matched %>% 
#     filter(str_detect(site, "Minnesota Zoo")) %>% 
#     mutate(is_het = if_else(het_count >= 2, "yes", "no")) %>% 
#     count(sex, is_het))
# 
# (v_2het_mnzoo <- df_2het_mnzoo %>% 
#     filter(is_het == "yes") %>% 
#     mutate(prop_diploid = n/sum(n)) %>% 
#     filter(sex == "male") %>% 
#     pull(prop_diploid))
# 
# 
# # for Appalachian
# 
# #for appalachian at level 1
# (df_1het_app <- df_het_matched %>% 
#     filter(named_cluster100 == "Appalachian") %>% 
#     mutate(is_het = if_else(het_count >= 1, "yes", "no")) %>% 
#     count(sex, is_het))
# 
# (v_1het_app <- df_1het_app %>% 
#     filter(is_het == "yes") %>% 
#     mutate(prop_diploid = n/sum(n)) %>% 
#     filter(sex == "male") %>% 
#     pull(prop_diploid))
# 
# #for appalachian at level 2
# (df_2het_app <- df_het_matched %>% 
#     filter(named_cluster100 == "Appalachian") %>% 
#     mutate(is_het = if_else(het_count >= 2, "yes", "no")) %>% 
#     count(sex, is_het))
# 
# (v_2het_app <- df_2het_app %>% 
#     filter(is_het == "yes") %>% 
#     mutate(prop_diploid = n/sum(n)) %>% 
#     filter(sex == "male") %>% 
#     pull(prop_diploid))
# 




# for all sites -----------------------------------------------------------

# actually no point in doing this unless merging various site names because there is splitting like "MN Zoo - site 1" "MN Zoo - Site 2" etc that should be merged; probably merge with same clusters used by colony count dataset

# df_2het_allsites <- df_het_matched %>%
#   mutate(is_het = if_else(het_count >= 2, "yes", "no")) %>%
#   count(named_cluster100, sex, is_het) %>%
#   complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>% 
#   # filter(is_het == "yes") %>% 
#   group_by(named_cluster100, is_het) %>% 
#   mutate(prop_diploid = n/sum(n),
#          n_mf = sum(n)) %>% 
#   # filter(sex == "male") %>% 
#   pivot_wider(names_from = is_het, values_from = n)



# giant summary table of all sites ----------------------------------------

# almost assuredly there are several better ways of doing this but somehow this is what I ended up doing...it's really kinda terrible but it works!
#n_male is the number of males from the site, n_female number of females, n_total is total, n_dip_total is the number of diploid/heterozygous individuals from the site, n_dip_male and n_dip_female are the numbers of heterozygous individuals by sex, cor_prop_dip_male is the "corrected" proportion of diploid males as in the proportion of males that are diploid compared to the total number of diploid individuals in the dataset, and total_prop_dip_male is just the number of diploid males/number of males

#there are apparently only 3 clusters with >10 males so maybe it only really makes sense to do the total and then the few focused looks (as already done above) anyway!

df_het_2locus <- df_het_matched %>%
  mutate(is_het = if_else(het_count >= 2, "yes", "no")) %>%
  count(named_cluster100, sex, is_het) %>%
  complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>%
  pivot_wider(names_from = is_het, values_from = n) %>%
  rowwise() %>%
  mutate(n_sex = sum(no, yes), n_dips = sum(yes)) %>%
  ungroup() %>%
  group_by(named_cluster100) %>%
  mutate(n_tot_dips = sum(n_dips),
         n_total = sum(n_sex)) %>%
  filter(sex == "male") %>%
  rename(n_male = n_sex,
         n_dip_male = yes) %>%
  mutate(n_female = n_total-n_male,
         n_dip_female = n_tot_dips - n_dip_male,
         n_dip_total = n_tot_dips) %>%
  select(-no, -sex, -n_dips) %>%
  select(n_male, n_female, n_total, n_dip_total, n_dip_male, n_dip_female) %>%
  mutate(cor_prop_dip_male = n_dip_male/n_dip_total,
         tot_prop_dip_male = n_dip_male/n_male)



#same as above but only 1 locus required

df_het_1locus <- df_het_matched %>%
  mutate(is_het = if_else(het_count >= 1, "yes", "no")) %>%
  count(named_cluster100, sex, is_het) %>%
  complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>%
  pivot_wider(names_from = is_het, values_from = n) %>%
  rowwise() %>%
  mutate(n_sex = sum(no, yes), n_dips = sum(yes)) %>%
  ungroup() %>%
  group_by(named_cluster100) %>%
  mutate(n_tot_dips = sum(n_dips),
         n_total = sum(n_sex)) %>%
  filter(sex == "male") %>%
  rename(n_male = n_sex,
         n_dip_male = yes) %>%
  mutate(n_female = n_total-n_male,
         n_dip_female = n_tot_dips - n_dip_male,
         n_dip_total = n_tot_dips) %>%
  select(-no, -sex, -n_dips) %>%
  select(n_male, n_female, n_total, n_dip_total, n_dip_male, n_dip_female) %>%
  mutate(cor_prop_dip_male = n_dip_male/n_dip_total,
         tot_prop_dip_male = n_dip_male/n_male)



# SAVE OUTPUTS ------------------------------------------------------------

saveRDS(df_het_1locus, "./analyses/analyses_output/03b_output_het1locus.Rdata")
saveRDS(df_het_2locus, "./analyses/analyses_output/03b_output_het2locus.Rdata")
