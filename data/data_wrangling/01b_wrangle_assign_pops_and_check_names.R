##%######################################################%##
#                                                          #
####           SECOND DATA WRANGLING FILE FOR           ####
####     2020/2021 RPBB POPULATION GENETICS PROJECT     ####
#                                                          #
##%######################################################%##


# The purpose of this document is to continue the wrangling of the 2020/2021 RPBB genetic data. Within this document, we ensure consistent naming of sites, states, etc. Originally I planned to do the spatial clustering here, too...but that might be a fraught step so I'll do it in the next script. 


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(naniar)

`%ni%` <- Negate(`%in%`)


# DATA --------------------------------------------------------------------

df_meta <- readRDS("./data/data_output/output_01a_df_rpbb_metadata.Rdata")


# VARIOUS NAME STANDARDIZATION --------------------------------------------

# Below, I want to wrangle to ensure consistent naming in sex, state, county, site, and floral_host (i.e. no redundant names or different names for the same thing)

# Consistent SEX labels ---

# For females there is "f", "female", and "worker"
# there are NO gynes in the dataset, so caste is unnecessary
# For males there is "m" "male" and "Male"
# There are 5 NAs for some reason...
#TODO - !!!!!!!!!!!!FOR NOW WE EXCLUDE THE NA SEX INDIVIDUALS!!!!!!!!!!!!!!!

v_males <- c("m", "male", "Male")
v_females <- c("f", "female", "Worker")

# Consistent STATE labels ---

v_wisconsin <- c("WI", "wisconsin")
v_illinois <- c("IL", "Illinois")
v_iowa <- c("IA", "IO", "Iowa")
v_minnesota <- c("minnesota", "Minnesota", "MN")
v_virginia <- c("Virginia")
v_west_virginia <- c("West Virginia", "WV")

# Consistent COUNTY labels ---

# the only issues with counties seem to be capitalization, so I just add that to the wrangling below

# Consistent SITE labels ---

# actually...I'm not sure this really matters and it's kinda a pain. If it matters down the line, I'll do it. But for now, I should probably focus on doing smarter things like measuring the distance between everyone and then clustering them

# Consistent FLORAL_HOST labels

# done below as it's a lot...

# The actual wrangling step; doing it kinda in a redundant way for better readability IMO

df_wrg_meta <- df_meta %>% 
  # SEX names and filter
  mutate(sex = case_when(
    sex %in% v_males ~ "male",
    sex %in% v_females ~ "female",
    sex == "NA" ~ NA_character_
  )) %>% 
  filter(!is.na(sex)) %>% 
  # STATE names
  mutate(state = case_when(
    state %in% v_wisconsin ~ "wisconsin",
    state %in% v_illinois ~ "illinois",
    state %in% v_iowa ~ "iowa",
    state %in% v_minnesota ~ "minnesota",
    state %in% v_virginia ~ "virginia",
    state %in% v_west_virginia ~ "west virginia",
    TRUE ~ NA_character_
  )) %>% 
  # COUNTY names
  mutate(county = tolower(county)) %>% 
  mutate(floral_host = case_when(
    floral_host == "agastache foeniculum" ~ "Agastache foeniculum",
    floral_host == "astilbe" ~ "Astilbe sp.",
    floral_host == "Cirsium species" ~ "Cirsium sp.",
    floral_host == "Daucus carota" ~ "Daucus carrota",
    floral_host == "eutrochium purpureum" ~ "Eurtrochium purpureum",
    floral_host == "linaria vulgaris" ~ "Linaria vulgaris",
    floral_host == "lotus corniculatus" ~ "Lotus corniculatus",
    floral_host == "melilotus albus" ~ "Melilotus albus",
    floral_host == "monarda fistulosa" ~ "Monarda fistulosa",
    floral_host %in% c("pycanthemum", "Pycnanthemum (virginiana or tenuiflora)", "pycnanthemum virginianum", "Pycnanthemum virginianum") ~ "Pycnanthemum sp.",
    floral_host == "Rudbeckia laciniata 'Goldkugel'" ~ "Rudbeckia laciniata",
    floral_host == "securigera varia" ~ "Securigera varia",
    floral_host == "sorbaria sorbifolia" ~ "Sorbaria sorbifolia",
    floral_host == "teucrium canadense" ~ "Teucrium canadense",
    floral_host == "thalictrum" ~ "Thalictrum sp.",
    floral_host == "Unknown" ~ NA_character_,
    floral_host == "vicia americana" ~ "Vicia americana",
    floral_host == "in flight" ~ NA_character_,
    TRUE ~ floral_host
  )) %>% 
  separate(col = floral_host, into = c("fl_genus", "fl_species"), sep = " ", extra = "merge")


saveRDS(df_wrg_meta, "./data/data_output/output_01b_df_rpbb_metadata.Rdata")
