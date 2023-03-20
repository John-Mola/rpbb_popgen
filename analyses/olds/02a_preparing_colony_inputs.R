##%######################################################%##
#                                                          #
####             PREPARING FILES FOR COLONY             ####
#                                                          #
##%######################################################%##

# Purpose: This script is used to take the output of 01d (merged genotype file) and prepare various inputs for COLONY. 

# We need an ERR file and a GENO file. We then paste those together to make the DAT file. 

#FUTURE DATS TO MAKES
#TODO - only females, known colonies excluded (done here!; I think the known colonies could throw off COLONY's assumptions)
#TODO - all females
#TODO - only known colonies
#TODO - males

# PACKAGES ----------------------------------------------------------------

library(tidyverse)


# DATA ---------------------------------------------------------------------

#Genotype data
df_rpbb_fulldata <- readRDS("./data/data_output/output_01d_merged_genotypes.Rdata")

#Error rate data of kept loci
df_error_rates <- readRDS("./data/data_output/output_01d_error_rates.Rdata")

#Helper vector with names of loci
v_loci_kept <- readRDS("./data/data_output/output_01d_vector_good_loci.Rdata")

# WRANGLING GENO FILE -----------------------------------------------------


# filter out individuals with less than 10 markers, select only the shortname and loci, convert NAs to 0, drop column names 

df_rpbb_female_geno <- df_rpbb_fulldata %>% 
  filter(loci_w_data >= 10,
         #NOTE that this is ONLY females in this dataset now!!
         sex == "female",
         #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
         is.na(which_nest)) %>% 
  dplyr::select(internal_barcode, all_of(v_loci_kept)) %>% 
  mutate(across(v_loci_kept, ~replace_na(.x, 0)))


write_tsv(df_rpbb_female_geno, file = "./analyses/inputs_colony/geno_files/output_02a_rpbb_female_geno.txt", col_names = FALSE)


# MN ZOO 2020 ONLY --------------------------------------------------------

# Filter to only MN Zoo specimens in 2020, then do rest

# filter out individuals with less than 10 markers, select only the shortname and loci, convert NAs to 0, drop column names

df_rpbb_mnz_2020 <- df_rpbb_fulldata %>% 
  # site specific filters
  filter(str_detect(site, pattern = "Zoo"),
         year == 2020) %>% 
  # general filters (copied from above)
  filter(loci_w_data >= 10,
         #NOTE that this is ONLY females in this dataset now!!
         sex == "female",
         #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
         is.na(which_nest)) %>% 
  dplyr::select(internal_barcode, all_of(v_loci_kept)) %>% 
  mutate(across(v_loci_kept, ~replace_na(.x, 0)))


write_tsv(df_rpbb_mnz_2020, file = "./analyses/inputs_colony/geno_files/output_02a_rpbb_mnz_2020.txt", col_names = FALSE)


# MN ZOO 2021 ONLY --------------------------------------------------------

# Filter to only MN Zoo specimens in 2021, then do rest

# filter out individuals with less than 10 markers, select only the shortname and loci, convert NAs to 0, drop column names

df_rpbb_mnz_2021 <- df_rpbb_fulldata %>% 
  # site specific filters
  filter(str_detect(site, pattern = "Zoo"),
         year == 2021) %>% 
  # general filters (copied from above)
  filter(loci_w_data >= 10,
         #NOTE that this is ONLY females in this dataset now!!
         sex == "female",
         #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
         is.na(which_nest)) %>% 
  dplyr::select(internal_barcode, all_of(v_loci_kept)) %>% 
  mutate(across(v_loci_kept, ~replace_na(.x, 0)))


write_tsv(df_rpbb_mnz_2021, file = "./analyses/inputs_colony/geno_files/output_02a_rpbb_mnz_2021.txt", col_names = FALSE)


# TURTLE VALLEY ONLY ------------------------------------------------------

# Filter to only turtle valley, then do rest

# filter out individuals with less than 10 markers, select only the shortname and loci, convert NAs to 0, drop column names

df_rpbb_turtle_2021 <- df_rpbb_fulldata %>% 
  # site specific filters
  filter(str_detect(site, pattern = "Turtle"),
         year == 2021) %>% 
  # general filters (copied from above)
  filter(loci_w_data >= 10,
         #NOTE that this is ONLY females in this dataset now!!
         sex == "female",
         #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
         is.na(which_nest)) %>% 
  dplyr::select(internal_barcode, all_of(v_loci_kept)) %>% 
  mutate(across(v_loci_kept, ~replace_na(.x, 0)))


write_tsv(df_rpbb_turtle_2021, file = "./analyses/inputs_colony/geno_files/output_02a_rpbb_turtle_2021.txt", col_names = FALSE)



# WRANGLING ERROR RATE FILE -----------------------------------------------

#this seems really dumb not sure why i couldn't get pivot_wider to behave as expected AFTER adding blank mutate columns but but oh well it works

df_rpbb_err <- df_error_rates %>% 
  #pivot existing error rate data to wide format so all loci are the column names with their error rate as only value
  pivot_wider(names_from = locus, values_from = error_rate) %>% 
  # add two dummy rows of NAs
  add_row(.before = 1) %>% 
  add_row(.before = 1) %>% 
  # replace NAs with 0s to satisfy COLONY criteria
  replace(is.na(.), 0)

write_tsv(df_rpbb_err, file = "./analyses/inputs_colony/err_files/output_02a_rpbb_err.txt", col_names = TRUE)



# NEXT STEPS ARE DONE OUTSIDE OF R ----------------------------------------

# In command line, generate a header file (steal my old one), then paste the geno and error rate files onto that to make the .DAT file. 
# Run .DAT file using colony_runner.opt or whatever it's named 




