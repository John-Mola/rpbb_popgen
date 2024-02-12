##%######################################################%##
#                                                          #
####             DIPLOID MALE CALCULATIONS              ####
#                                                          #
##%######################################################%##


# PACKAGES ----------------------------------------------------------------

library(tidyverse)


# DATA --------------------------------------------------------------------

#dataset of genotypes
# this is set up for source in the markdown document, it will not run if you don't update the directory depth
df_merged <- read_rds("./data/data_output/output_01d_merged_genotypes.Rdata")



# FINDING HET COUNT FOR ALL -----------------------------------------------

df_genos <- df_merged %>% 
  dplyr::select(internal_barcode, 28:53)

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



# MALE DIPLOID FREQUENCY AT 3 LOCI OR MORE --------------------------------

df_het_3locus <- df_het_matched %>%
  mutate(is_het = if_else(het_count >= 3, "yes", "no")) %>%
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


saveRDS(df_het_3locus, "./analyses/analyses_output/03b_output_het3locus.Rdata")







# CHECKING SOME STUFF FOLLOWING REVIEWS -----------------------------------


# removing btern02 to see if that changes stuff

# df_matched_nobtern <- purrr::map(c("btern01", "bt28", "b96", "bt30", "btms0081", "btms0066", "btms0083", "b126", "btms0062", "btern02", "btms0086", "bl13", "btms0059"), ~ df_genos %>%
#                            mutate(
#                              !!(paste0(.x, "_ans", collapse = "")) :=
#                                UQ(rlang::sym(paste0(.x, "_1", collapse = ""))) == UQ(rlang::sym(paste0(.x, "_2", collapse = ""))) )) %>%
#   reduce(., left_join) %>%
#   dplyr::select(internal_barcode, contains("_ans")) %>%
#   rowwise() %>%
#   dplyr::select(-btern02_ans) %>% 
#   mutate(het_count = sum(c_across(all_of(contains("_ans"))) == "FALSE", na.rm = TRUE))
# 
# # df_matched %>% group_by(het_count) %>% tally()
# 
# df_het_matched_nobtern02 <- inner_join(df_merged, dplyr::select(df_matched_nobtern, internal_barcode, het_count)) %>% 
#   filter(loci_w_data >= 10)
# 
# 
# df_het3loci_nobtern02 <- df_het_matched_nobtern02 %>%
#   mutate(is_het = if_else(het_count >= 3, "yes", "no")) %>%
#   count(named_cluster100, sex, is_het) %>%
#   complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>%
#   pivot_wider(names_from = is_het, values_from = n) %>%
#   rowwise() %>%
#   mutate(n_sex = sum(no, yes), n_dips = sum(yes)) %>%
#   ungroup() %>%
#   group_by(named_cluster100) %>%
#   mutate(n_tot_dips = sum(n_dips),
#          n_total = sum(n_sex)) %>%
#   filter(sex == "male") %>%
#   rename(n_male = n_sex,
#          n_dip_male = yes) %>%
#   mutate(n_female = n_total-n_male,
#          n_dip_female = n_tot_dips - n_dip_male,
#          n_dip_total = n_tot_dips) %>%
#   select(-no, -sex, -n_dips) %>%
#   select(n_male, n_female, n_total, n_dip_total, n_dip_male, n_dip_female) %>%
#   mutate(cor_prop_dip_male = n_dip_male/n_dip_total,
#          tot_prop_dip_male = n_dip_male/n_male)
# 
# 
# df_het1loci <- df_het_matched %>%
#   mutate(is_het = if_else(het_count >= 1, "yes", "no")) %>%
#   count(named_cluster100, sex, is_het) %>%
#   complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>%
#   pivot_wider(names_from = is_het, values_from = n) %>%
#   rowwise() %>%
#   mutate(n_sex = sum(no, yes), n_dips = sum(yes)) %>%
#   ungroup() %>%
#   group_by(named_cluster100) %>%
#   mutate(n_tot_dips = sum(n_dips),
#          n_total = sum(n_sex)) %>%
#   filter(sex == "male") %>%
#   rename(n_male = n_sex,
#          n_dip_male = yes) %>%
#   mutate(n_female = n_total-n_male,
#          n_dip_female = n_tot_dips - n_dip_male,
#          n_dip_total = n_tot_dips) %>%
#   select(-no, -sex, -n_dips) %>%
#   select(n_male, n_female, n_total, n_dip_total, n_dip_male, n_dip_female) %>%
#   mutate(cor_prop_dip_male = n_dip_male/n_dip_total,
#          tot_prop_dip_male = n_dip_male/n_male)
# 
# 
# df_het2loci <- df_het_matched %>%
#   mutate(is_het = if_else(het_count >= 2, "yes", "no")) %>%
#   count(named_cluster100, sex, is_het) %>%
#   complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>%
#   pivot_wider(names_from = is_het, values_from = n) %>%
#   rowwise() %>%
#   mutate(n_sex = sum(no, yes), n_dips = sum(yes)) %>%
#   ungroup() %>%
#   group_by(named_cluster100) %>%
#   mutate(n_tot_dips = sum(n_dips),
#          n_total = sum(n_sex)) %>%
#   filter(sex == "male") %>%
#   rename(n_male = n_sex,
#          n_dip_male = yes) %>%
#   mutate(n_female = n_total-n_male,
#          n_dip_female = n_tot_dips - n_dip_male,
#          n_dip_total = n_tot_dips) %>%
#   select(-no, -sex, -n_dips) %>%
#   select(n_male, n_female, n_total, n_dip_total, n_dip_male, n_dip_female) %>%
#   mutate(cor_prop_dip_male = n_dip_male/n_dip_total,
#          tot_prop_dip_male = n_dip_male/n_male)
# 
# 
# df_het3loci <- df_het_matched %>%
#   mutate(is_het = if_else(het_count >= 3, "yes", "no")) %>%
#   count(named_cluster100, sex, is_het) %>%
#   complete(named_cluster100, nesting(sex, is_het), fill = list(n = 0)) %>%
#   pivot_wider(names_from = is_het, values_from = n) %>%
#   rowwise() %>%
#   mutate(n_sex = sum(no, yes), n_dips = sum(yes)) %>%
#   ungroup() %>%
#   group_by(named_cluster100) %>%
#   mutate(n_tot_dips = sum(n_dips),
#          n_total = sum(n_sex)) %>%
#   filter(sex == "male") %>%
#   rename(n_male = n_sex,
#          n_dip_male = yes) %>%
#   mutate(n_female = n_total-n_male,
#          n_dip_female = n_tot_dips - n_dip_male,
#          n_dip_total = n_tot_dips) %>%
#   select(-no, -sex, -n_dips) %>%
#   select(n_male, n_female, n_total, n_dip_total, n_dip_male, n_dip_female) %>%
#   mutate(cor_prop_dip_male = n_dip_male/n_dip_total,
#          tot_prop_dip_male = n_dip_male/n_male)
