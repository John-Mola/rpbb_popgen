##%######################################################%##
#                                                          #
####             PREPARING FILES FOR COLONY             ####
#                                                          #
##%######################################################%##

#THIS INTENDED AS A BRANCH OF 02A TO DO THIS A BIT "AUTOMATED" AND ALL IN R

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


# COMMON FILTERS FOR ALL GENO SETS ----------------------------------------


# filter out individuals with less than 10 markers, females, only individuals not from known nests

df_rpbb_female_geno <- df_rpbb_fulldata %>% 
  filter(loci_w_data >= 10,
         #NOTE that this is ONLY females in this dataset now!!
         sex == "female",
         #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
         is.na(which_nest)) 



# GENOTYPE FILES FROM EACH CROSS OF CLUSTER X YEAR ------------------------

# saving the grouped dataframe for use in group_split (it's a little silly but whatever)
df_grp_clst_yr <- df_rpbb_female_geno %>%
  group_by(named_cluster100, year)

# group splitting the dataframe by the groups defined above, creates a list of separate dfs for each group
list_clst_year <- group_split(df_grp_clst_yr)
#group_keys(df_grp_clst_yr)

# getting the "group keys" with count of group membership, creating a helper column for easy filtering and whatnot
df_groups_key <- df_grp_clst_yr %>% count() %>% mutate(clst_yr = str_replace(paste0(named_cluster100,"_",year), pattern = " ", "_"))

#saving this as a helper to match in later steps
df_keys_match <- df_groups_key %>% 
  rowid_to_column("rowid") %>% 
  mutate(rowid = as.character(rowid))

#save the helper!
saveRDS(df_keys_match, file = "./analyses/inputs_colony/02a01_rpbb_batch_group_keys.Rdata")

# taking the list of dfs and selecting only the internal barcode and loci data, replacing NAs with 0s for COLONY's rules

list_geno_clst_year <- map(list_clst_year, .f = list(. %>% dplyr::select(internal_barcode, all_of(v_loci_kept)) %>% mutate(across(v_loci_kept, ~replace_na(.x, 0)))))


# MAKE HEADERS ------------------------------------------------------------


# read in the header generic file

colony_header <- read_file("./analyses/inputs_colony/header_files/colony_header_partial.txt")


# takes the key, pastes all the shit into one column because that was the dumb way I figured out to be able to write it to a file easily...., split splits it back into dfs (which is the same as groupsplit but for some reason i used base here instead of tidy...; then back to tidy by using map to select only the one needed column...)

list_group_keys <- df_groups_key %>% 
  mutate(header_top = paste(clst_yr, clst_yr, n, colony_header, sep = "\n")) %>% 
  ungroup() %>% 
  # dplyr::select(header_top) %>% 
  split(.$clst_yr) %>% 
  map(., .f = list(. %>% dplyr::select(header_top)))

# save it all as header files to a folder...no touch folder!!
mapply(function (x,y) write_tsv(x, file = paste0('./analyses/inputs_colony/header_files/batch_rpbb/', y, '.txt'), col_names = FALSE), list_group_keys, paste0(names(list_group_keys),"_header"))  



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



# ADDING GENOTYPE BELOW HEADERS -------------------------------------------

# reading the error rate file back in even though it was just generated above...seems dumb....but I was afraid it wouldn't behave right in the cat below if not "fresh"...dunno, brain dumb

err_rates <- read_file("./analyses/inputs_colony/err_files/output_02a_rpbb_err.txt")

# This writes all of the genotype files generated in a list above to their own files/directory...we can then load that directory at the same time that we load the headers

mapply(function (x,y) write_tsv(x, file = paste0('./analyses/inputs_colony/geno_files/batch_rpbb/', y, '.txt'), col_names = FALSE), list_geno_clst_year, paste0(names(list_group_keys),"_geno")) 

# save a vector of the paths to the genofiles and a vector of paths to the header files for our big dumb system call

v_genofiles <- list.files("./analyses/inputs_colony/geno_files/batch_rpbb/", full.names = TRUE)
v_headerfiles <- list.files("./analyses/inputs_colony/header_files/batch_rpbb/", full.names = TRUE)

# this just pastes together the command that we're going to send to the terminal. I don't know if there's a way to do something similar in R, but I already knew how to do it in bash...so this basically just creates a command line command that we then loop over below. Takes the form of cat header errorfile genofile footer > DAT file name

df_files <- data_frame(cat = "cat", header = v_headerfiles, error = "./analyses/inputs_colony/err_files/output_02a_rpbb_err.txt", geno = v_genofiles, footer = "./analyses/inputs_colony/header_files/colony_footer.txt", pipe = "> ./analyses/inputs_colony/DAT_files/batch_rpbb/", cluster = df_groups_key$clst_yr) %>% 
  mutate(sysstring = paste(cat, header, error, geno, footer, paste0(pipe, cluster,".DAT")))

#loop over all of the commands in the df_files$sysstring column to run them all. Saves a bunch of files to the DAT_files/batch_rpbb directory

for (x in 1:nrow(df_files)) {
  system(df_files$sysstring[x])
}




# MAKE EXCEL FILE FOR EASY RUNNING OF COLONY ------------------------------

# this is just to help me run them all...I'm sure I could run this using the batch fuction available in colony...but I wanted them in separate folders and honestly it doesn't take that long to copy-paste each COLONY command so whatever...it's all automated until the very last step apparently! lol

v_datfiles <- list.files("./analyses/inputs_colony/DAT_files/batch_rpbb/", full.names = TRUE)

df_helper <- data_frame(batch = df_groups_key$clst_yr,
                        begin = "~/Colony2_Mac_01_02_2022/run_colony.out",
                        middle1 = "IFN:../../../.",
                        middle2 = v_datfiles,
                        end =  " ; afplay ~/Downloads/Dry\ Town.mp3") %>% 
  mutate(colcommand = paste(begin, paste0(middle1, middle2), end))


write_csv(df_helper, "./analyses/inputs_colony/rpbb_batch_jan_2023_helper_list.csv")

# NEXT STEPS ARE DONE OUTSIDE OF R ----------------------------------------

# In command line, copy-paste all of the commands from the last column of the helper file




