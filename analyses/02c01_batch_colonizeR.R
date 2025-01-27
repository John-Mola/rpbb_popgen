##%######################################################%##
#                                                          #
####             WRANGLING OUTPUT OF COLONY             ####
#                                                          #
##%######################################################%##

#PURPOSE - this script is to take the output of COLONY, filter out "bad" families, do some merging with the original metadata file, and otherwise "wrangle" the information to be usable for downstream analysis; filename is a bit weird because of carryover from an old version


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(janitor)


# DATA --------------------------------------------------------------------

# group keys helper, needed for easy matching back to the long df of best clusters
df_keys_match <- readRDS("./analyses/inputs_colony/02a01_rpbb_batch_group_keys.Rdata")

# pull all over the bestclusters from COLONY outputs...make a list of file paths, then load all of them, bind them together with a key, join that key to the group keys, janitor
list_path_bcs = list.files(path = "analyses/outputs_colony/batch_rpbb/", pattern="*.BestCluster", recursive = TRUE, full.names = TRUE)
df_bcs = lapply(list_path_bcs, read_table, col_names = TRUE) %>% bind_rows(.id = "rowid") %>% inner_join(., df_keys_match) %>% clean_names()


# Metadata from the 01d step; we may need the genotype data when loading into future steps, so keep it all connected
df_md <- readRDS("data/data_output/output_01d_merged_genotypes.Rdata")

# Name the output file for below
pn_batch = "rpbb_batch_colonizeR"

# Output directory for the final csv
od = "analyses/outputs_colony/r_colonizer/"

# BESTCLUSTER FILTERING ---------------------------------------------------

#ONLY KEEP families with a probability equal to/greater than 80% -- might want to make this more stringent later!!
fam_threshold <- 0.8

# filtering poor families. Adds a value to make the number way over anything generated by COLONY to ensure no dupes. Then I bind the cluster and year information to that family name to avoid duplicate numbers between colony batches

df_flt_batch <- df_bcs %>%  mutate(family_index = if_else(probability < fam_threshold, 1000+row_number(), cluster_index),
                                   family_index = paste0(family_index, "_", clst_yr))



# COLONIZER FUNCTION ------------------------------------------------------

#Copied from my old github repository https://github.com/John-Mola/SNPs_to_COLONY/blob/master/Scripts_MasterCopies/02a_colonizeR.R

# Not sure I care that this is in a function still. Kinda seems unnecessary. Maybe it's good though. 

# Provide function with a .BestCluster output from COLONY, a metadata file with a column named "unique.ID" that matches column "OffspringID" from BestCluster (ostensibly, these things already exist since that's needed for the input to DAT_maker), a name for the output file, and an output directory.

# Run the function
colonizeR<-function(bestCluster, metadat, projectName, outdir)
{
  require(tidyverse)

  # Trim BestCluster
  bestC = bestCluster %>%
    select(family_index, offspring_id)

  # Join, remove individuals absent in BestCluster, create FamilyID column, remove extra variables
  dat_full <- inner_join(metadat, bestC, by=c("internal_barcode"="offspring_id"))  %>%
    arrange(family_index) %>%
    mutate(occur = !family_index %in% family_index[duplicated(family_index)]
    ) %>%
    mutate(family_id = if_else(occur == TRUE, "s", as.character(family_index))) #%>%
  #select(-occur,-gps.ID,-Library.ID,-unique.ID.old,-plate,-occur,-plate.ID)

  ## Save the dataframe
  write_csv(dat_full, file = paste0(outdir,projectName,".csv"), col_names = TRUE)

}


# Run it with objects defined above (change names if desired)
colonizeR(df_flt_batch,df_md,pn_batch,od)



# # FILTERING 1 REPRESENTATIVE PER FAMILY -----------------------------------

# If we want to have a file where only 1 sibling per family is retained for downstream analysis. This does that.Not sure it's really necessary though...

# Load the file we just created for viewing and further filtering

df_batch_colonizer <- read_csv("./analyses/outputs_colony/r_colonizer/rpbb_batch_colonizeR.csv")



df_rpbb_keepers <- df_batch_colonizer %>%
  group_by(family_index) %>%
  sample_n(1)

v_rpbb_keepers <- df_rpbb_keepers %>% pull(internal_barcode)

saveRDS(v_rpbb_keepers, file = "./analyses/outputs_colony/r_colonizer/02c_v_rpbb_keepers.Rdata")



