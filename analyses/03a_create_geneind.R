##%######################################################%##
#                                                          #
####               CREATE GENEIND OBJECT                ####
#                                                          #
##%######################################################%##

#PURPOSE - this script is used to create a geneind object from our genotype file (really from the output of COLONY combined back with our genotype file!)


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(adegenet)
library(poppr)
library(hierfstat)
library(reshape2)

# DATA --------------------------------------------------------------------

#THIS CURRENT DATASET is females only and excluding specimens from known colonies; to match the 02a01 script where this filtering is done, I simply just copied below (as that script outputs only a txt file of a later step)

# # pre-COLONY siblings data
# 
# df_rpbb_fulldata <- readRDS("./data/data_output/output_01d_merged_genotypes.Rdata") %>% 
#   filter(loci_w_data >= 10,
#          #NOTE that this is ONLY females in this dataset now!!
#          sex == "female",
#          #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
#          is.na(which_nest)) #%>% mutate(cluster = paste0(state, " (",cluster,")"))


# post-COLONY siblings data

df_rpbb_colonizer <- read_csv("./analyses/outputs_colony/r_colonizer/rpbb_batch_colonizeR.csv")

v_rpbb_keepers <- readRDS("./analyses/outputs_colony/r_colonizer/02c_v_rpbb_keepers.Rdata")

# Loci used in analysis

v_loci_kept <- readRDS("./data/data_output/output_01d_vector_good_loci.Rdata")


# PREPARE DATA FOR CREATING GENIND OBJECT ---------------------------------

# Need to have the genotypes in dataframe with no other information, then a vector of individual IDs and a vector of site IDs. We then provide that to df2genind() to create the genind object


# select only the individual shortnames, sites, and genotypes from the dataframe
df_rpbb_simple <- df_rpbb_colonizer %>% 
  dplyr::select(internal_barcode, named_cluster100, all_of(v_loci_kept)) %>%
  # !!!!!!!! THIS IS WHERE THE REMOVAL OF SIBLINGS HAPPENS!!!!!!!!
  filter(internal_barcode %in% v_rpbb_keepers)

# create vector of individual IDs
v_rpbb_shortnames <- pull(df_rpbb_simple, internal_barcode)

# create vector of site names
v_rpbb_sites <- pull(df_rpbb_simple, named_cluster100)

# remove shortname, site names, then combine adjacent columns to get each allele into a single cell, last a little ditty to remove the "_1" in the merged columns
df_rpbb_onlyGeno <- df_rpbb_simple %>% 
  dplyr::select(-internal_barcode, -named_cluster100)

df_rpbb_mergeGenos <- mapply(function(x, y) {
  paste(x, y, sep = ",")},
  df_rpbb_onlyGeno[ ,seq(1, ncol(df_rpbb_onlyGeno), by = 2)],
  df_rpbb_onlyGeno[ ,seq(2, ncol(df_rpbb_onlyGeno), by = 2)])

colnames(df_rpbb_mergeGenos) <- gsub(x = colnames(df_rpbb_mergeGenos), pattern = "_1", replacement = "", fixed = TRUE)  

df_rpbb_alleles <- as.data.frame(df_rpbb_mergeGenos) %>% 
  mutate_all(funs(str_replace_all(., "NA,NA", NA_character_)))



# CREATE GENIND OBJECT ----------------------------------------------------

gen_rpbb = df2genind(df_rpbb_alleles, ploidy = 2, ind.names = v_rpbb_shortnames, pop = v_rpbb_sites, sep = ",")



# SAVE GENEIND OBJECT -----------------------------------------------------

saveRDS(gen_rpbb, "./analyses/analyses_output/03a_rpbb_femaleNOknown_NOSibs_genind.Rdata")




