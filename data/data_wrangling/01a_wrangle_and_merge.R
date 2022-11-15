## %######################################################%##
#                                                          #
####              WRANGLING AND MERGING OF              ####
####           2020 AND 2021 RPBB COLLECTIONS           ####
#                                                          #
## %######################################################%##

# This document is used to wrange and merge the 2020 and 2021 RPBB genetic collections. Field data and USDA-called genotype data are both loaded in and combined. Some additional helper columns are added and other changes to prepare the data for use in COLONY and downstream analyses.
# Some metadata from 2020 were cleaned using a manual step due to confusion arising from collecting data from several collaborators during a global pandemic...
# Field data from my own collections were also cleaned in a separate step in a prior database. Should the need arise, this file can be provided (but is somewhat unnecessary as we don't have that same cleaning file for other collaborators who submitted data and used a mix of hand and automated cleaning methods)

# Initial tasks

# DONE - which metadata do we want? From USDA or from collectors? What is missing (if anything relevant) from USDA?
#### --> basically only need caste, whether or not it's from a nest, and floral associations. Everything else is redundant. So let's go with whatever the USDA sheet says for all of that.

# Quality checks

# TODO - confirm missing specimens with Jay Watson
# TODO - clean up this document

# Things to do in the next step 01b

# TODO - ensure specimens have real lat-long and not obscured, if obscured get real and merge
# TODO merge with genotype data
# TODO create a "number of loci" column
# TODO - assign all specimens to consistent site names
# TODO group nearby specimens into "populations"...


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(naniar)

`%ni%` <- Negate(`%in%`)

# READ RAW DATA, CLEAN NAMES ----------------------------------------------

# cleaned metadata from 2020 specimens

df_2020_meta <- read_csv("./data/data_output/df_wrangled_2020rpbbMSATs_2021_07_09.csv")

# metadata and genotypes from USDA database

# this was the old call - the new data has a different format with metadata and genotype data on different pages
# df_raw_2021_meta <- read_excel("./metadata/meta_raw/AffinsRepUSDA_Koch_21Dec2021.xlsx", sheet = 2) %>% clean_names()
# maybe for now we just merge metadata with USDA database and make sure all specimens are accounted for...then, we merge the genotype data in later

df_raw_USDA_meta <- read_excel("./data/data_raw/Bombus_affinis_repository__msatdata_ver_22September2022.xlsx", sheet = 1) %>% clean_names()

# USDA data sheet two to help remove metadata that doesn't have any genotype data

df_raw_USDA_2meta <- read_excel("./data/data_raw/Bombus_affinis_repository__msatdata_ver_22September2022.xlsx", sheet = 2) %>% clean_names()

# additional metadata from collectors ----

df_boone <- read_csv("./data/data_raw/from_collectors/boone_raw_data_10-26-2021.csv") %>% clean_names()

df_hepner <- read_excel("./data/data_raw/from_collectors/hepner_raw_data_2021.xlsx") %>% clean_names()

df_watson <- read_excel("./data/data_raw/from_collectors/watson_raw_data_2021.xlsx") %>% clean_names()

df_jean <- read_excel("./data/data_raw/from_collectors/jean_raw_data_2021.xlsx") %>% clean_names()

df_kochanski <- read_excel("./data/data_raw/from_collectors/kochanski_raw_data_2021.xlsx", sheet = 1) %>% clean_names()

df_runquist <- read_excel("./data/data_raw/from_collectors/runquist_raw_data_2021.xlsx", sheet = 1) %>% clean_names()

df_mola <- read_csv("./data/data_output/2021-08-07-rpbb-ISP-JM-tarsi-matched.csv") %>% clean_names()



# REMOVING METADATA WITHOUT GENOTYPES -------------------------------------

# in the USDA sheet2 the specimens that do not have genotypes have a blank "Name" column (some specimens that have a name lack sufficient genotype data but we'll filter that at a later step)

df_has_geno <- df_raw_USDA_2meta %>%
  filter(!is.na(name_3))

# need to merge sheet 2 and sheet one by the internal_barcode because I need dates in order to filter out old specimens from the expected metadata

df_modern <- inner_join(df_has_geno, df_raw_USDA_meta, by = c("internal_barcode", "sample_tube_id_with_description")) %>%
  mutate(year = year(caught_dd_mmm_yyyy)) %>%
  filter(year > 2019)



# CLEANING FIELD COLLECTOR METADATA ---------------------------------------

df_flt_boone <- df_boone %>%
  filter(
    bombus_species == "affinis",
    leg_or_swab == "leg"
  ) %>%
  dplyr::select(unique_id, floral_host, notes) %>%
  mutate(
    caste = NA_character_,
    nest = NA_character_,
    id2 = NA_character_
  ) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)

df_flt_hepner <- df_hepner %>%
  dplyr::select(unique_id, caste, nest, floral_host, notes) %>%
  mutate(id2 = NA_character_) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)

df_flt_jean <- df_jean %>%
  dplyr::select(unique_id, caste, nest, floral_host, notes) %>%
  mutate(id2 = NA_character_) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)

df_flt_kochanski <- df_kochanski %>%
  dplyr::select(unique_id, caste, nest, floral_host, notes) %>%
  mutate(id2 = NA_character_) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)

df_flt_runquist <- df_runquist %>%
  dplyr::select(unique_id, caste, nest, floral_host, notes) %>%
  mutate(id2 = NA_character_) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)

df_flt_watson <- df_watson %>%
  dplyr::select(usda_tag, notes) %>%
  rename(unique_id = usda_tag) %>%
  mutate(
    caste = NA_character_,
    nest = NA_character_,
    id2 = NA_character_,
    floral_host = NA_character_
  ) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)

df_flt_mola <- df_mola %>%
  dplyr::select(un_id, vial_top, caste, nectar_plant, notes) %>%
  rename(unique_id = vial_top, id2 = un_id, floral_host = nectar_plant) %>%
  mutate(
    unique_id = as.character(unique_id),
    nest = NA_character_,
    id2 = NA_character_
  ) %>%
  dplyr::select(unique_id, id2, caste, nest, floral_host, notes)


df_collectors_mrg <- bind_rows(df_flt_boone, df_flt_hepner, df_flt_jean, df_flt_kochanski, df_flt_mola, df_flt_watson, df_flt_runquist) %>%
  mutate(
    nest = if_else(floral_host %in% c("nest", "entering nest", "leaving nest"), "nest", nest),
    floral_host = if_else(floral_host %in% c("nest", "entering nest", "leaving nest"), NA_character_, floral_host)
  ) %>%
  replace_with_na_all(condition = ~ .x %in% common_na_strings) %>%
  # why did I add this "attached" column? Great documentation, John.
  mutate(attached = "attached")



# SOME WRANGLING TO 2021 DATA  --------------------------------------------

# TODO - Some values are duplicated. Search for dupes and missing values
# CLEAN - these are the duplicated values, for now I have simply removed all instances as we cannot be sure what is the "true" one or if they came from separate specimens or what:

v_dupes <- c("jmk-12jul21-01", "jmk-12jul21-02")

df_flt_2021_meta <- df_raw_USDA_meta %>%
  mutate(year = year(caught_dd_mmm_yyyy)) %>%
  inner_join(., df_collectors_mrg, by = c("sample_tube_id_with_description" = "unique_id")) %>%
  filter(sample_tube_id_with_description %ni% v_dupes) %>%
  rename(
    longname = sample_tube_id_with_description,
    sex = male_m_female_f_unk,
    date = caught_dd_mmm_yyyy,
    from_nest = nest,
    site = site_description_name
  )


# MERGE 2021 WITH 2020 ----------------------------------------------------

# CHECK - this is weird because the 2020 data is entirely from its own metadata file and the 2021 data is from the USDA mostly and then combined with a bit of collector data. Potentially for 2020 we want to pull only the relevant columns (nest, floral associations, notes) and an ID column and then join them....
## --> I think this ends up being fine as this was previously sorted. Just join the barcode for reference.

df_USDA_codes <- df_raw_USDA_meta %>%
  select(internal_barcode, sample_tube_id_with_description)

df_partial_2020_meta <- dplyr::select(df_2020_meta, longname, sex, from_nest, which_nest, latitude, longitude, state, county, year, date, site_longname) %>%
  rename(site = site_longname) %>%
  mutate(floral_host = NA_character_)

df_partial_2021_meta <- df_flt_2021_meta %>%
  dplyr::select(longname, sex, from_nest, latitude, longitude, state, county, year, date, site, floral_host) %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )

# This merge ends up containing missing 10 specimens that are in the USDA database but do not have appropriate metadata. 4 are the duplicates from Jade's collections (with no real way of figuring out exactly where they came from), 1 is an untracked specimen from Elaine (similarly with no real way of knowing whats up), but the remaining five are from Jay. There's a chance some of them are from the known nest...so for now, this seems like an appropriate list.

df_meta_merged <- bind_rows(df_partial_2020_meta, df_partial_2021_meta) %>%
  inner_join(., df_USDA_codes, by = c("longname" = "sample_tube_id_with_description")) %>%
  filter(year > 2019)

# There were some records with missing lat-long coordinates. I merge them here for simplicity sake.

df_missing_ll <- read_csv("./data/data_raw/from_collectors/df_mlb_missing_coord_fixed.csv") %>%
  dplyr::select(sample_tube_id_with_description, latitude, longitude)

df_meta_merged_corrected <- full_join(df_meta_merged, df_missing_ll, by = c("longname" = "sample_tube_id_with_description")) %>%
  # this is a dumb way to do it but I'm lazy rn and it works. Not sure how to do the join above in the way I want and not end up with missing/redundant coordinate values
  mutate(
    latitude = case_when(
      !is.na(latitude.x) ~ latitude.x,
      is.na(latitude.x) ~ latitude.y,
    ),
    longitude = case_when(
      !is.na(longitude.x) ~ longitude.x,
      is.na(longitude.x) ~ longitude.y
    )
  ) %>%
  dplyr::select(-latitude.x, -longitude.x, -latitude.y, -longitude.y) %>%
  # there's one coordinate MLB provided that is missing a USDA barcode (presumably it was not genotyped) so filter that out
  filter(!is.na(internal_barcode))


saveRDS(df_meta_merged_corrected, file = "./data/data_output/output_01a_df_rpbb_metadata.Rdata")

# SOME QUALITY CHECKS -----------------------------------------------------

# Who is missing?
#
# # df_raw_USDA_2meta has all of the specimens that actually have genetic data
# nrow(df_has_geno) #498 rows - so there's this many missing:
# nrow(df_has_geno) - nrow(df_meta_merged) # 21 missing!
#
# df_missing_meta <- anti_join(df_has_geno, df_meta_merged, by = c("sample_tube_id_with_description" = "longname"))
#
# #ok - some of these missing are old specimens! Easy peasy. Made a new df called df_modern which keeps only specimens from 2020 onwards
#
# df_missing_meta2 <- anti_join(df_modern, df_meta_merged, by = c("sample_tube_id_with_description" = "longname"))
# nrow(df_missing_meta2)

# ok ...now we "only" have 10 missing...4 missing ones are the duplicates that was intentionally removed. So there are 6 that I dunno why they're missing. All from Elaine (1) and Jay (5)

# can simply pull the Elaine and Jay ones from here over into the merged metadata. The Kochanski ones that are duplicated will need to be dropped completely unless there's a way of separating them out/determining true information for each. Might also need the non-obscured locations.



# SOME APPALACHIA STUFF ---------------------------------------------------

# df_meta_merged %>%
#   filter(state %in% c("West Virginia", "WV", "Virginia")) %>%
#   group_by(county) %>%
#   tally()
#
# df_meta_merged %>%
#   filter(state %in% c("West Virginia", "WV", "Virginia")) %>%
#   group_by(site) %>%
#   tally()
#
# df_meta_merged %>%
#   filter(state %in% c("West Virginia", "WV", "Virginia")) %>%
#   group_by(sex) %>%
#   tally()
#
# df_meta_merged %>%
#   filter(state %in% c("West Virginia", "WV", "Virginia"),
#          sex == "Male") %>%
#   mutate(dayno = lubridate::yday(date)) %>%
#   group_by(dayno) %>%
#   tally() %>%
#   ggplot(., aes(x = dayno, y = n)) +
#   geom_point()
#
# #####
####
###
##
#
##
###
####
##### --- from previous

# REMOVE OLD SPECIMENS ----------------------------------------------------

# NOTE - a few specimens had incorrect year (2021 instead of 2020), I corrected these in Excel for simplicity-sake
#
# df_flt <- df %>%
#   mutate(year = year(caught_dd_mmm_yyyy),
#          longitude = as.numeric(longitude),
#          latitude = as.numeric(latitude)) %>%
#   # this includes specimens from 2017 onwards
#   filter(year > 2015)
#
#
# # CORRECTING SOME ERRORS AND ASSIGNING TO COLONIES ------------------------
#
# #NOTE - the specimens from "Redwing" in Hennepin are actually from a site just called "Minneapolis Nest" and were entered wrong at some point
#
# #NOTE - 828 East Ave is the redwing nest...
#
# nest_names <- c("Redwing", "Minneapolis Nest", "828 East Ave Redwing")
#
# df_corrected <- df_flt %>%
#   mutate(site_description_name = if_else(site_description_name == "Redwing" & county == "Hennepin", "Minneapolis Nest", site_description_name),
#          site_description_name = if_else(site_description_name == "828 East Ave Redwing", "Redwing", site_description_name)) %>%
#   # add a column that specifies if a bee comes from a known nest
#   mutate(from_nest = if_else(site_description_name %in% nest_names, "yes", "no"),
#          which_nest = if_else(from_nest == "yes", site_description_name, NA_character_),
#          site_shortname = abbreviate(site_description_name))
#
#
# # ASSIGNING TO SIMPLIFIED REGIONAL SITES ----------------------------------
#
# v_iowa_south <- c("Hawk Eye - Iowa City", "Hawk Eye Male - Iowa City")
# v_iowa_east <- c("East of Andrew")
# v_iowa_north <- c("Seed Saver", "Hutchinson Farm GG")
# v_mn_metro <- c("Pickerell Lake - found dead", "MSP Airport on Tarmac - found dead", "Minneapolis - 39th Ave South", "Ramsey MN", "Carver Park Reserve - Found Dead", "Minneapolis - Mt. Fistulosa", "Across the stree from MN nest", "18069 Rawlings St - Elkriver City", "Eden Prarie", "Redwing", "Minneapolis Nest", "Ames Lake Park", "929 Central Ave Red Wing", "Cherokee Heights", "Elaine's Garden", "828 East Ave Redwing", "Minnesota Zoo")
# # for now, lumping MN Zoo (and all Twin Cities region sites) into one category
# #v_mn_zoo <- c("Minnesota Zoo")
# v_wi_jeff <- c("Natural Resources Conservation Service Wetland Reserve Easement--private property")
# v_wi_dane <- c("Behind 7637 Farmington Way", "Near Kwik Trip", "Lamp Park", "High Point Church", "Madison Country Day School", "Merry Street")
# v_wv <- c("USFS Marlinton Ranger Station Pollinator Garden", "USFS Monongahela Nat'l Forest - Neola")
# v_ill <- c("unk")
#
#
# df_broadsites <- df_corrected %>%
#   mutate(reg_site = case_when(
#     site_description_name %in% v_iowa_south ~ "iowa south",
#     site_description_name %in% v_iowa_east ~ "iowa east",
#     site_description_name %in% v_iowa_north ~ "iowa north",
#     site_description_name %in% v_mn_metro ~ "minneapolis metro",
#     #site_description_name %in% v_mn_zoo ~ "minneapolis zoo",
#     site_description_name %in% v_wi_jeff ~ "wisconsin east",
#     site_description_name %in% v_wi_dane ~ "wisconsin dane",
#     site_description_name %in% v_wv ~ "west virginia",
#     site_description_name %in% v_ill ~ "northern illinois",
#     site_description_name == TRUE ~ NA_character_
#   ))
#
# # ASSIGN SHORTNAMES -------------------------------------------------------
#
# # shortname should include sex, some sort of location, unique number for duplicates; needs to be less than 20 characters for COLONY
#
# df_shortnamed <- df_broadsites %>%
#   mutate(region_shortname = abbreviate(reg_site)) %>%
#   group_by(male_m_female_f_unk, region_shortname) %>%
#   mutate(group_count = str_pad(row_number(), width = 2, side = "left", pad = 0),
#          shortname = paste0(region_shortname, group_count, male_m_female_f_unk))
#
# # LOCI COUNT COLUMN -------------------------------------------------------
#
# # copy-paste loci names then add _1 and _2 to them to match the df
# v_loci_raw <- c("b124", "btern01", "bt28", "bt10", "b96", "bt30", "btms0081", "btms0066", "btms0083", "b126", "btms0062", "btern02", "btms0086", "bl13", "btms0059")
# v_loci_1 <- paste0(v_loci_raw, "_1")
# v_loci_2 <- paste0(v_loci_raw, "_2")
# v_loci <- c(v_loci_1, v_loci_2)
#
#
# # Do a rowwise count of the number of alleles scored (divide by 2 to get number of loci with data)
#
# df_loci_count <- df_shortnamed %>%
#   rowwise() %>%
#   # this is kinda janky because I'm counting values above zero instead of just non-NA values but whatever it works
#   mutate(loci_w_data = sum(c_across(v_loci) > 0, na.rm = TRUE)/2) %>%
#   ungroup()
#
#
# # FINAL CHECKS AND SAVE ---------------------------------------------------
#
#
# # yincredibly inefficient but just want it the way i want it
#
# df_wrangled <- df_loci_count %>%
#   dplyr::select(1:12, year, from_nest, which_nest, site_shortname, reg_site, region_shortname, group_count, loci_w_data, everything()) %>%
#   mutate(specimen_type = if_else(leg_1_or_zero == 1, "leg", if_else(whole_body_1_or_zero == 1, "body", "unknown"))) %>%
#   dplyr::select(1, specimen_type, everything(), -whole_body_1_or_zero, -leg_1_or_zero) %>%
#   rename(longname = sample_description_a_l, sex = male_m_female_f_unk, date = caught_dd_mmm_yyyy, site_longname = site_description_name) %>%
#   dplyr::select(longname, shortname, specimen_type, sex, from_nest, which_nest, latitude, longitude, reg_site, region_shortname, group_count, site_longname, site_shortname, county, state, date, time_hrs, year, loci_w_data, everything(), -time_hrs)
#
#
#
#
# # SAVE IT!
#
# write_csv(df_wrangled, "./data/data_output/df_wrangled_2020rpbbMSATs_2021_07_09.csv")
#
#











# test area ---------------------------------------------------------------

# # this finds all of the site names that are in the same state and county. I could then measure the distances between these to ensure they are >10km (or whatever) from any other site. If they are, then they are a unique "population". If they are not, then the sites should be lumped as a "population". Also need to check possibility of sites in neighbor counties being <10km apart. Some specimens collected by Elaine in Ramsey, Hennepin, and Dakota county might be like this. The two WV specimens as well.
# df_corrected %>% group_by(state, county) %>% summarise(sites = toString(unique(site_description_name)))
#
# ## MN Zoo samples, just females
#
# df_mnz <- df_corrected %>%
#   filter(site_description_name == "Minnesota Zoo", male_m_female_f_unk == "f") %>%
#   dplyr::select(shortname, 13:42)
#
# write_tsv(df_mnz, "./data/data_output/mnz_geno.txt", col_names = FALSE, na = c("0"))
#
#
# # all males, check for inbreeding
#
# df_males <- df_corrected %>%
#   filter(male_m_female_f_unk == "m") %>%
#   dplyr::select(shortname, 13:42)
#
# # map(c("B124", "BTERN01", "BT28", "BT10", "B96", "BT30", "BTMS0081", "BTMS0066", "BTMS0083", "B126", "BTMS0062", "BTERN02", "BTMS0086", "BL13", "BTMS0059"),
#
# #This works but I don't really understand how
#
# #From: https://stackoverflow.com/questions/49838191/r-iterate-over-pairs-of-columns-in-dataframe-by-regex-match
#
#
# df_matched <- map(c("b124", "btern01", "bt28", "bt10", "b96", "bt30", "btms0081", "btms0066", "btms0083", "b126", "btms0062", "btern02", "btms0086", "bl13", "btms0059"), ~ df_males %>%
#                     mutate(
#                       !!(paste0(.x, "_ans", collapse = "")) :=
#                       UQ(rlang::sym(paste0(.x, "_1", collapse = ""))) == UQ(rlang::sym(paste0(.x, "_2", collapse = ""))) )) %>%
#   reduce(., left_join) %>%
#   dplyr::select(shortname, contains("_ans")) %>%
#   rowwise() %>%
#   mutate(het_count = sum(c_across(all_of(contains("_ans"))) == "FALSE", na.rm = TRUE))
#
# df_matched %>% group_by(het_count) %>% tally()
#
# df_male_match <- inner_join(df_corrected, dplyr::select(df_matched, shortname, het_count))
#
#
# # map check (run 2021-03-19-rpbb-tarsi-collections-to-date.rmd first)
#
# v_historic_states <- df_affinis_historic %>% distinct(state) %>% pull(state)
#
# # v_extant_states <- df_affinis_extant %>% distinct(state) %>% pull(state)
#
# v_extant_counties <- df_affinis_extant %>% distinct(county, .keep_all = TRUE) %>% mutate(county = gsub(x = county, pattern = ".", replacement = "", fixed = TRUE), state_county = paste0(state, ",", county)) %>% pull(state_county)
#
# bg_map <- st_as_sf(map("state", regions = v_historic_states, plot = FALSE, fill = TRUE))
#
# bg_map_extant <- st_as_sf(map("county", regions = v_extant_counties, plot = FALSE, fill = TRUE))
#
# #df_rpbb_county_summary <- df_meta_tarsi_cleannames %>% group_by(location_state, location_county) %>% summarise(total_collected = sum(number_collected), lat = mean(location_utm_lat), long = mean(location_utm_long))
#
# p_map <- ggplot() +
#   geom_sf(data = bg_map, fill = "antiquewhite", alpha = 0.4, color = "grey80", size = 0.4) +
#   geom_sf(data = bg_map_extant, fill = "grey", alpha = 0.5, color = "grey80", size = 0.4) +
#   geom_point(data = df_broadsites, aes(x = longitude, y = latitude, fill = reg_site), alpha = 0.5, color = "black", shape = 21, size = 2) +
#   #annotation_scale(location = "bl", width_hint = 0.4) +
#   #annotation_north_arrow(location = "bl", which_north = "true",
#   #                       pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
#   #                       style = north_arrow_fancy_orienteering) +
#   theme_bw() +
#   theme(panel.grid.major = element_line(color = gray(0.9),
#                                         linetype = "dashed",
#                                         size = 0.2),
#         panel.background = element_rect(fill = "white")) +
#   labs(x = "", y = "")
#
# ggplotly(p_map) #%>%
# # style(hoverlabel = list(bgcolor = "white")) %>%
# style(hoveron = "fill")
