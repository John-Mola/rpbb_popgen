# also junk draft for the presentation on 2022-12-05...

#assumes df_rpbb_fulldata gets loaded in from 02a

df_rpbb_fulldata

# all males, check for inbreeding


df_males <- df_rpbb_fulldata %>%
  filter(sex == "male") %>%
  dplyr::select(internal_barcode, 17:42)

# map(c("B124", "BTERN01", "BT28", "BT10", "B96", "BT30", "BTMS0081", "BTMS0066", "BTMS0083", "B126", "BTMS0062", "BTERN02", "BTMS0086", "BL13", "BTMS0059"),

#This works but I don't really understand how

#From: https://stackoverflow.com/questions/49838191/r-iterate-over-pairs-of-columns-in-dataframe-by-regex-match


df_matched <- map(c("btern01", "bt28", "b96", "bt30", "btms0081", "btms0066", "btms0083", "b126", "btms0062", "btern02", "btms0086", "bl13", "btms0059"), ~ df_males %>%
                    mutate(
                      !!(paste0(.x, "_ans", collapse = "")) :=
                      UQ(rlang::sym(paste0(.x, "_1", collapse = ""))) == UQ(rlang::sym(paste0(.x, "_2", collapse = ""))) )) %>%
  reduce(., left_join) %>%
  dplyr::select(internal_barcode, contains("_ans")) %>%
  rowwise() %>%
  mutate(het_count = sum(c_across(all_of(contains("_ans"))) == "FALSE", na.rm = TRUE))

df_matched %>% group_by(het_count) %>% tally()

df_male_match <- inner_join(df_rpbb_fulldata, dplyr::select(df_matched, internal_barcode, het_count))


df_male_match %>% 
  filter(str_detect(site, pattern = "Zoo")) %>% 
  dplyr::select(internal_barcode, year, loci_w_data, het_count) %>% 
  mutate(diploid = if_else(het_count > 3, "yes", "no")) %>% 
  count(year, diploid)
         