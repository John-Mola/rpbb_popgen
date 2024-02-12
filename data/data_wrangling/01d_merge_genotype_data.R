##%######################################################%##
#                                                          #
####           MERGING METADATA WITH GENOTYPE           ####
####              DATA FROM USDA DATABASE               ####
#                                                          #
##%######################################################%##



#PURPOSE - to take the metadata output of 01c (where individuals are assigned to "clusters") and merge it with the genotype data provided by the USDA


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(janitor)


# DATA --------------------------------------------------------------------

# Output of 01c, the wrangled metadata assigned to clusters
df_clustered <- readRDS("./data/data_output/output_01c_df_rpbb_clustered.Rdata")

# USDA genotype data

df_geno <- readxl::read_xlsx("./data/data_raw/Bombus_affinis_repository__msatdata_ver_22September2022.xlsx", sheet = 2) %>% clean_names()


# USDA error rate data
df_err_rates <- readxl::read_xlsx("./data/data_raw/Bombus_affinis_repository__msatdata_ver_22September2022.xlsx", sheet = 3) %>% clean_names() %>% dplyr::select(locus, error_rate)



# REMOVE BAD LOCI ---------------------------------------------------------


# Using a threshold of 0.3 (matches USDA sheet recommend)
accepted_error = 0.3


df_good_loci <- df_err_rates %>%
  filter(error_rate < accepted_error) %>%
  mutate(locus = tolower(locus))

v_good_loci <- df_good_loci %>% pull(locus)

v_good_loci_cols <- df_good_loci %>%
  mutate(
    loci_1 = paste0(locus, "_1"),
    loci_2 = paste0(locus, "_2")
  ) %>%
  select(loci_1, loci_2) %>%
  pivot_longer(cols = c(loci_1, loci_2)) %>%
  pull(value)
  

# saving the bad loci column names as that seems cleaner to filter out below
v_bad_loci_cols <- df_err_rates %>%
  mutate(locus = tolower(locus)) %>% 
  filter(error_rate > accepted_error) %>% 
  mutate(
    loci_1 = paste0(locus, "_1"),
    loci_2 = paste0(locus, "_2")
  ) %>%
  select(loci_1, loci_2) %>%
  pivot_longer(cols = c(loci_1, loci_2)) %>%
  pull(value)


# WRANGLING OF USDA SHEET -------------------------------------------------

df_flt_geno <- df_geno %>% 
  dplyr::select(-sample_tube_id_with_description, -name_3, -name_18)


# MERGE METADATA WITH GENOTYPES BY BARCODE --------------------------------

# # Then adding a loci count. Below are two vectors of loci names. 
# v_loci <- dplyr::select(df_flt_geno, -internal_barcode) %>% names()
# #Not sure if we'll need this, but hey. 
# v_loci_raw <- str_remove(v_loci, "\\_[^.]*$") %>% unique()

df_merged <- inner_join(df_clustered, df_flt_geno, by = c("internal_barcode")) %>% 
  #just some rearranging of column order for prettiness
  dplyr::select(internal_barcode, longname, sex, from_nest, which_nest, state, county, year, date, site, fl_genus, fl_species, cluster05, cluster10, cluster25, cluster50, cluster100, n_cluster05, n_cluster10, n_cluster25, n_cluster50, n_cluster100, named_cluster100, latitude, longitude, everything()) %>% 
  # REMOVE the bad loci
  dplyr::select(-v_bad_loci_cols) %>% 
  # COUNT the number of loci each individual has
    rowwise() %>%
    # this is kinda janky because I'm counting values above zero instead of just non-NA values but whatever it works
    mutate(loci_w_data = sum(c_across(all_of(v_good_loci_cols)) > 0, na.rm = TRUE)/2) %>%
    ungroup()

  

# SAVE RESULTS ------------------------------------------------------------

# Genotype file
saveRDS(df_merged, "./data/data_output/output_01d_merged_genotypes.Rdata")

# Error rates df
saveRDS(df_good_loci, "./data/data_output/output_01d_error_rates.Rdata")

# Error rates vector
saveRDS(v_good_loci_cols, "./data/data_output/output_01d_vector_good_loci.Rdata")
