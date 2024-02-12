##%######################################################%##
#                                                          #
####           FIGURE MAKING COLLECTIONS MAP            ####
#                                                          #
##%######################################################

# PURPOSE - this script makes the main map of collections and putative "populations" (i.e. 100km clusters)


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(sf)
library(maps)
library(ggspatial)
# library(svglite)
# library(ggalt)
library(ggnewscale)
library(pals)


# DATA --------------------------------------------------------------------


# Merged genotype and metadata generated in 01d

df_rpbb_fulldata <- readRDS("./data/data_output/output_01d_merged_genotypes.Rdata")

#Metadata; historic affinis states for mapping

df_affinis_historic <- read_csv("./data/data_raw/meta_rpbb_external/rpbb_historic_counties.csv")  %>% clean_names()

# big aff database from Clint

df_all_aff <- read_csv("./data/data_raw/RPBB_gdb_updated_skinny_xy_20220516.csv") %>% clean_names()

# PLOTTING ----------------------------------------------------------------

df_modern_aff <- df_all_aff %>% 
  filter(year > 2015) 


v_historic_states <- df_affinis_historic %>% distinct(state) %>% pull(state)

# bg_map <- st_as_sf(map("state", regions = v_historic_states, plot = FALSE, fill = TRUE))

focused_bg_map <- st_as_sf(map("state", regions = c("Minnesota", "Wisconsin", "Iowa", "Illinois", "Indiana", "Michigan", "Ohio", "Kentucky", "West Virginia", "Virginia"), plot = FALSE, fill = TRUE)) %>% 
  mutate(cons_unit = case_when(
    ID %in% c("minnesota", "wisconsin") ~ "CU 1",
    ID %in% c("iowa", "illinois") ~ "CU 2",
    ID %in% c("indiana", "michigan", "ohio") ~ "CU 3", 
    ID %in% c("virginia", "west virginia", "kentucky") ~ "CU 4"
  ))

p_map <- ggplot() +
  geom_sf(data = focused_bg_map, alpha = 0.5, color = "grey60", size = 0.4, aes(fill = cons_unit)) +
  scale_fill_manual(values = hcl.colors(n = 4, palette = "Pastel 1")) + 
  new_scale_fill() +
  # geom_encircle(data = df_modern_aff, aes(x = longitude, y = latitude), alpha = 0.5, fill = "red", expand = 0.01) +
  # geom_point(data = df_modern_aff, aes(x = longitude, y = latitude), color = "black", alpha = 0.1, size = 2) +
  geom_bin2d(data = df_modern_aff, aes(x = longitude, y = latitude), bins = 50, alpha = 0.8) +   
  scale_fill_gradient(trans = "log", guide = "none", low = "grey90", high = "grey10") +
  new_scale_fill() +
  geom_point(data = df_rpbb_fulldata, aes(x = longitude, y = latitude, fill = as.factor(named_cluster100)), alpha = 0.7, size = 3, shape = 21) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_line(color = gray(0.9),
                                        linetype = "dashed",
                                        size = 0.2),
        panel.background = element_rect(fill = "white"),
        # legend.position = c(0.8, 0.7),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_rect(fill = FALSE)) +
  labs(x = "", y = "", color = "", fill = "") +
  scale_fill_manual(values = unname(alphabet(13)), breaks = c("Twin Cities", "SE Minnesota", "Iowa City", "Decorah",  "Quad Cities", "Madison", "Central Wisconsin",  "North Illinois","Milwaukee",  "Chicago", "North Milwaukee", "Green Bay",  "Appalachian"))
  
  

p_map


# SAVE OUTPUT -------------------------------------------------------------

#svg version -- need to manually move legends and re-color the UP of Michigan to get final figure. 
ggsave(p_map, filename = "./figures/figures_output/04a_collections_map.svg")

