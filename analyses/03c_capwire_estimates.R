##%######################################################%##
#                                                          #
#### ESTIMATING ML ABUNDANCE OF COLONIES ACROSS SITES   ####
#                                                          #
##%######################################################%##

# This script is intended to take the output of colonizer, create family size tables for sites, and then run capwire to get the ML estimate of family size. 


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(capwire)


# DATA --------------------------------------------------------------------

df_rpbb_colonizer <- read_csv("./analyses/outputs_colony/r_colonizer/rpbb_batch_colonizeR.csv")



# FINDING SAMPLE SIZES AND COLONY COUNTS ALL SITES ------------------------


# find the sample size per 10km

df_N_10kmYear <- df_rpbb_colonizer %>% 
  count(cluster10, year, name = "N") 

# Find the number of colonies (including singletons!) by 10km and year

df_Ncol_10kmYear <- df_rpbb_colonizer %>% 
  group_by(cluster10, year, family_index) %>% 
  tally() %>% 
  group_by(cluster10, year) %>% 
  tally(name = "Ncol")
  
# Join those two together

df_N_Ncol_10kmYear <- inner_join(df_N_10kmYear, df_Ncol_10kmYear, by = c("cluster10", "year"))


# WRANGLING FOR "BIG SITES" ONLY STUFF ------------------------------------

## save a df of all sites with >=10 samples. 

df_sites_10kmN9 <- df_rpbb_colonizer %>% 
  group_by(cluster10, year, site, named_cluster100) %>% 
  summarise(n=n()) %>% 
  group_by(cluster10, year, named_cluster100) %>% 
  summarise(n = sum(n), sites = toString(site)) %>% 
  filter(n >= 9) %>% 
  arrange(desc(n)) %>% 
  # giving stuff meaningful names manually
  mutate(custom_name = case_when(
    cluster10 == 50 & year == 2021 ~ "Turtle21",
    cluster10 == 14 & year == 2021 ~ "MNZoo21",
    cluster10 == 14 & year == 2020 ~ "MNZoo20",
    cluster10 == 20 & year == 2021 ~ "Nachusa21",
    cluster10 == 53 & year == 2021 ~ "Foxglove21",
    cluster10 == 8 & year == 2021 ~ "Hamline21",
    cluster10 == 6 & year == 2020 ~ "Cherokee20",
    cluster10 == 19 & year == 2021 ~ "Illiniwek21",
    cluster10 == 8 & year == 2020 ~ "Elaine20",
    cluster10 == 5 & year == 2021 ~ "SeedSavers21",
    cluster10 == 16 & year == 2021 ~ "CapSprings21",
    cluster10 == 29 & year == 2021 ~ "Appalachian21",
  )) %>% 
  dplyr::select(cluster10, year, custom_name)

## filter the rpbb_colonizer df to include only the target sites

df_bigsite_colonizer <- df_rpbb_colonizer %>% 
  inner_join(., df_sites_10kmN9, by = c("cluster10", "year"))

v_bigsites <- df_sites_10kmN9 %>% pull(custom_name)



# GENERATE FAMILY DISTRIBUTIONS -------------------------------------------

#make a df with each family count for pulling site-level vectors from

df_bigsite_counts <- df_bigsite_colonizer %>% 
  group_by(custom_name, family_index) %>% 
  tally() %>% 
  arrange(desc(n)) 


# Get the family size vector for each site
## could automate this but probably too much of a PITA so....

## "Cherokee20"    "Elaine20"      "Foxglove21"    "Hamline21"     "Illiniwek21" "MNZoo20"       "MNZoo21"       "Nachusa21"     "SeedSavers21"  "Turtle21"  
#NOTE - occassionally the bootstrap function is failing with "Error in max(lik) : invalid 'type' (list) of argument"...not sure why, because then it will run fine on another attempt

# Appalachian21----

v_Appalachian21 <- df_bigsite_counts %>% 
  filter(custom_name == "Appalachian21") %>% 
  pull(n)

captable_Appalachian21 <- buildClassTable(as.numeric(v_Appalachian21))
rs_Appalachian21 <- fitTirm(data=captable_Appalachian21, max.pop=1000)
#throws an error because there are only singletons...
#bs_Appalachian21 <- bootstrapCapwire(rs_Appalachian21, bootstraps = 1000)
output_Appalachian21 <- tibble(n_genotyped = sum(v_Appalachian21), n_colonies = length(v_Appalachian21), ml.colony.num =  NA_integer_, CI.lower= NA_integer_, CI.upper = NA_integer_)


# CapSprings21----

v_CapSprings21 <- df_bigsite_counts %>% 
  filter(custom_name == "CapSprings21") %>% 
  pull(n)

captable_CapSprings21 <- buildClassTable(as.numeric(v_CapSprings21))
rs_CapSprings21 <- fitTirm(data=captable_CapSprings21, max.pop=1000)
bs_CapSprings21 <- bootstrapCapwire(rs_CapSprings21, bootstraps = 1000)
output_CapSprings21 <- tibble(n_genotyped = sum(v_CapSprings21), n_colonies = length(v_CapSprings21), ml.colony.num =  rs_CapSprings21[[3]], CI.lower= bs_CapSprings21[[2]][[1]], CI.upper = bs_CapSprings21[[2]][[2]])

# Cherokee20---- 

v_Cherokee20 <- df_bigsite_counts %>% 
  filter(custom_name == "Cherokee20") %>% 
  pull(n)

captable_Cherokee20 <- buildClassTable(as.numeric(v_Cherokee20))
rs_Cherokee20 <- fitTirm(data=captable_Cherokee20, max.pop=1000)
bs_Cherokee20 <- bootstrapCapwire(rs_Cherokee20, bootstraps = 1000)
output_Cherokee20 <- tibble(n_genotyped = sum(v_Cherokee20), n_colonies = length(v_Cherokee20), ml.colony.num =  rs_Cherokee20[[3]], CI.lower= bs_Cherokee20[[2]][[1]], CI.upper = bs_Cherokee20[[2]][[2]])

# Elaine20----    

v_Elaine20 <- df_bigsite_counts %>% 
  filter(custom_name == "Elaine20") %>% 
  pull(n)

captable_Elaine20 <- buildClassTable(as.numeric(v_Elaine20))
rs_Elaine20 <- fitTirm(data=captable_Elaine20, max.pop=1000)
bs_Elaine20 <- bootstrapCapwire(rs_Elaine20, bootstraps = 1000)
output_Elaine20 <- tibble(n_genotyped = sum(v_Elaine20), n_colonies = length(v_Elaine20), ml.colony.num =  rs_Elaine20[[3]], CI.lower= bs_Elaine20[[2]][[1]], CI.upper = bs_Elaine20[[2]][[2]])

# Foxglove21----

v_Foxglove21 <- df_bigsite_counts %>% 
  filter(custom_name == "Foxglove21") %>% 
  pull(n)

captable_Foxglove21 <- buildClassTable(as.numeric(v_Foxglove21))
rs_Foxglove21 <- fitTirm(data=captable_Foxglove21, max.pop=1000)
bs_Foxglove21 <- bootstrapCapwire(rs_Foxglove21, bootstraps = 1000)
output_Foxglove21 <- tibble(n_genotyped = sum(v_Foxglove21), n_colonies = length(v_Foxglove21), ml.colony.num =  rs_Foxglove21[[3]], CI.lower= bs_Foxglove21[[2]][[1]], CI.upper = bs_Foxglove21[[2]][[2]])

# Hamline21  ----

v_Hamline21 <- df_bigsite_counts %>% 
  filter(custom_name == "Hamline21") %>% 
  pull(n)

captable_Hamline21 <- buildClassTable(as.numeric(v_Hamline21))
rs_Hamline21 <- fitTirm(data=captable_Hamline21, max.pop=1000)
bs_Hamline21 <- bootstrapCapwire(rs_Hamline21, bootstraps = 1000)
output_Hamline21 <- tibble(n_genotyped = sum(v_Hamline21), n_colonies = length(v_Hamline21), ml.colony.num =  rs_Hamline21[[3]], CI.lower= bs_Hamline21[[2]][[1]], CI.upper = bs_Hamline21[[2]][[2]])

# Illiniwek21----

v_Illiniwek21 <- df_bigsite_counts %>% 
  filter(custom_name == "Illiniwek21") %>% 
  pull(n)

captable_Illiniwek21 <- buildClassTable(as.numeric(v_Illiniwek21))
rs_Illiniwek21 <- fitTirm(data=captable_Illiniwek21, max.pop=1000)
bs_Illiniwek21 <- bootstrapCapwire(rs_Illiniwek21, bootstraps = 1000)
output_Illiniwek21 <- tibble(n_genotyped = sum(v_Illiniwek21), n_colonies = length(v_Illiniwek21), ml.colony.num =  rs_Illiniwek21[[3]], CI.lower= bs_Illiniwek21[[2]][[1]], CI.upper = bs_Illiniwek21[[2]][[2]])

# MNZoo20   ----

v_MNZoo20 <- df_bigsite_counts %>% 
  filter(custom_name == "MNZoo20") %>% 
  pull(n)

captable_MNZoo20 <- buildClassTable(as.numeric(v_MNZoo20))
rs_MNZoo20 <- fitTirm(data=captable_MNZoo20, max.pop=1000)
bs_MNZoo20 <- bootstrapCapwire(rs_MNZoo20, bootstraps = 1000)
output_MNZoo20 <- tibble(n_genotyped = sum(v_MNZoo20), n_colonies = length(v_MNZoo20), ml.colony.num =  rs_MNZoo20[[3]], CI.lower= bs_MNZoo20[[2]][[1]], CI.upper = bs_MNZoo20[[2]][[2]])

# MNZoo21   ----

v_MNZoo21 <- df_bigsite_counts %>% 
  filter(custom_name == "MNZoo21") %>% 
  pull(n)

captable_MNZoo21 <- buildClassTable(as.numeric(v_MNZoo21))
rs_MNZoo21 <- fitTirm(data=captable_MNZoo21, max.pop=1000)
bs_MNZoo21 <- bootstrapCapwire(rs_MNZoo21, bootstraps = 1000)
output_MNZoo21 <- tibble(n_genotyped = sum(v_MNZoo21), n_colonies = length(v_MNZoo21), ml.colony.num =  rs_MNZoo21[[3]], CI.lower= bs_MNZoo21[[2]][[1]], CI.upper = bs_MNZoo21[[2]][[2]])

# Nachusa21  ----

v_Nachusa21 <- df_bigsite_counts %>% 
  filter(custom_name == "Nachusa21") %>% 
  pull(n)

captable_Nachusa21 <- buildClassTable(as.numeric(v_Nachusa21))
rs_Nachusa21 <- fitTirm(data=captable_Nachusa21, max.pop=1000)
bs_Nachusa21 <- bootstrapCapwire(rs_Nachusa21, bootstraps = 1000)
output_Nachusa21 <- tibble(n_genotyped = sum(v_Nachusa21), n_colonies = length(v_Nachusa21), ml.colony.num =  rs_Nachusa21[[3]], CI.lower= bs_Nachusa21[[2]][[1]], CI.upper = bs_Nachusa21[[2]][[2]])

# SeedSavers21----

v_SeedSavers21 <- df_bigsite_counts %>% 
  filter(custom_name == "SeedSavers21") %>% 
  pull(n)

captable_SeedSavers21 <- buildClassTable(as.numeric(v_SeedSavers21))
rs_SeedSavers21 <- fitTirm(data=captable_SeedSavers21, max.pop=1000)
bs_SeedSavers21 <- bootstrapCapwire(rs_SeedSavers21, bootstraps = 1000)
output_SeedSavers21 <- tibble(n_genotyped = sum(v_SeedSavers21), n_colonies = length(v_SeedSavers21), ml.colony.num =  rs_SeedSavers21[[3]], CI.lower= bs_SeedSavers21[[2]][[1]], CI.upper = bs_SeedSavers21[[2]][[2]])

# Turtle21 ----

v_Turtle21 <- df_bigsite_counts %>% 
  filter(custom_name == "Turtle21") %>% 
  pull(n)

captable_Turtle21 <- buildClassTable(as.numeric(v_Turtle21))
rs_Turtle21 <- fitTirm(data=captable_Turtle21, max.pop=1000)
bs_Turtle21 <- bootstrapCapwire(rs_Turtle21, bootstraps = 1000)
output_Turtle21 <- tibble(n_genotyped = sum(v_Turtle21), n_colonies = length(v_Turtle21), ml.colony.num =  rs_Turtle21[[3]], CI.lower= bs_Turtle21[[2]][[1]], CI.upper = bs_Turtle21[[2]][[2]])



# COMBINE ALL OUTPUTS INTO COMMON DATAFRAME -------------------------------

df_capwire_results <- bind_rows(output_Appalachian21, output_CapSprings21, output_Cherokee20, output_Elaine20, output_Foxglove21, output_Hamline21, output_Illiniwek21, output_MNZoo20, output_MNZoo21, output_Nachusa21, output_SeedSavers21, output_Turtle21) %>% mutate(custom_name = sort(v_bigsites)) %>% dplyr::select(custom_name, everything()) %>% 
  mutate(prop_detected = n_colonies/ml.colony.num,
         raw_detected = n_colonies/n_genotyped)

df_capwire_results %>% arrange(desc(n_genotyped))



# SAVE VARIOUS OUTPUTS ----------------------------------------------------

# Full table of all sites

saveRDS(df_N_Ncol_10kmYear, "./analyses/analyses_output/03c_allsite_colonycounts.Rdata")

# Capwrite results output

saveRDS(df_capwire_results, "./analyses/analyses_output/03c_bigsites_capwireoutput.Rdata")


