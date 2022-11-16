##%######################################################%##
#                                                          #
####             WRANGLING OUTPUT OF COLONY             ####
#                                                          #
##%######################################################%##

#PURPOSE - this script is to take the output of COLONY, filter out "bad" families, do some merging with the original metadata file, and otherwise "wrangle" the information to be usable for downstream analysis

# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(janitor)


# DATA --------------------------------------------------------------------

# BestCluster from COLONY; for the run using ONLY females EXCLUDING known colonies
df_bc <- read_table(file = "./analyses/outputs_colony/rpbb_femalesNOknown/rpbb_femalesNOknown.BestCluster", col_names = TRUE) %>% clean_names()

# Metadata from the 01d step; we may need the genotype data when loading into future steps, so keep it all connected
df_md <- readRDS("data/data_output/output_01d_merged_genotypes.Rdata")

# Name the output file for below
pn = "rpbb_femalesNOknown_2022-11-14-colonizeR"

# Output directory for the final csv
od = "analyses/outputs_colony/r_colonizer/"

# BESTCLUSTER FILTERING ---------------------------------------------------

#ONLY KEEP families with a probability equal to/greater than 80%
fam_threshold <- 0.8


df_flt_bc <- df_bc %>%
  #1000+row_number just adds a high number so it doesn't get confused with any real families. Made it 1000 because that means anything with a four digit number was reassigned. This may prove to be a bad decision later if something required 3 digits for some dumb reason. 
  #changing name to "family_index" because I use "cluster" for assignment to population clusters. To mimimize confusion, let's change that...
  mutate(family_index = if_else(probability < fam_threshold, 1000+row_number(), cluster_index))
  

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
colonizeR(df_flt_bc,df_md,pn,od)

# Load the file we just created for viewing and further filtering

df_rpbb_colonizer <- read_csv("./analyses/outputs_colony/r_colonizer/rpbb_femalesNOknown_2022-11-14-colonizeR.csv")



# SANITY CHECK OF LOCATIONS -----------------------------------------------

#TODO - figure out why this is happening...cause it seems wrong. Re-ran COLONY turning off polyandry but that doesn't seem to fix it completely...it's certainly improved though. 
#TODO - re-running AGAIN allowing inbreeding in the calculation. The manual says this may slow down computation a lot...but does not say it would "hurt" per se

# this takes the output, removes singletons, groups by family, and then determines which clusters each family (colony) belongs to. Because clusters are grouped by proximity...families with individuals across clusters seem kind of suspicious...There are several colonies with representatives from multiple clusters. 
df_rpbb_colonizer %>% filter(family_id != "s") %>% group_by(family_id) %>% summarise(which_clusters = toString(unique(cluster))) %>% View()


# FILTERING 1 REPRESENTATIVE PER FAMILY -----------------------------------

#We want to have a file where only 1 sibling per family is retained for downstream analysis. This does that. 

#TODO revisit this - I am not satisfied with COLONY output but for now, let's move forward to get a "complete" analysis 


#oh yeah this is completely unnecessary because as long as you select *at least* one per family then all singletons will be selected too...
# df_rpbb_singletons <- df_rpbb_colonizer %>% filter(family_id == "s")
# 
# df_rpbb_lonesibling <- df_rpbb_colonizer %>%
#   filter(family_id != "s") %>%
#   group_by(family_id) %>%
#   sample_n(1)
# 
# v_singletons <- df_rpbb_singletons %>% pull(internal_barcode)
# v_familyreps <- df_rpbb_lonesibling %>% pull(internal_barcode)

df_rpbb_keepers <- df_rpbb_colonizer %>% 
  group_by(family_index) %>% 
  sample_n(1)

v_rpbb_keepers <- df_rpbb_keepers %>% pull(internal_barcode)

saveRDS(v_rpbb_keepers, file = "./analyses/outputs_colony/r_colonizer/02c_v_rpbb_keepers.Rdata")
