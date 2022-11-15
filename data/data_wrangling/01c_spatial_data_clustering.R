##%######################################################%##
#                                                          #
####             SPATIAL DATA CLUSTERING OF             ####
####            POINTS TO ASSIGN POPULATIONS            ####
#                                                          #
##%######################################################%##

#PURPOSE - to take the metadata with lat-long and assign them to clusters based on some threshold rule of distance. This seems preferable to assigning to population based on some more arbitrary metric (like county or just...my decision)


# Things to do in this document

#TODO - ensure specimens have real lat-long and not obscured, if obscured get real and merge
#TODO - plot everything on a map and make sure they are as expected (go back and make corrections as necessary)
#TODO group nearby specimens into "populations" (not sure of the smartest way to do this...)
### Going to try to follow the StackOverflow page here: https://gis.stackexchange.com/questions/17638/clustering-spatial-data-in-r
#TODO - clean this thing up!

# Things to do in the next document

#TODO - merge genotype data with the output of this document
#TODO - add a loci count and other quality check columns as needed
#TODO - ensure output is ready to be converted to required files for COLONY


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(sp)
library(sf)
library(maps)
library(rgdal)
library(geosphere)
library(janitor)
library(plotly)


# DATA --------------------------------------------------------------------

df_meta <- readRDS("./data/data_output/output_01b_df_rpbb_metadata.Rdata")


# CONVERT METADATA FRAME INTO A SPATIAL DATA FRAME ------------------------

df_meta2=rename(df_meta, x=longitude,y=latitude) %>% 
  #TODO - !!!!!!!! NOTE!!!!!!!!! I HAD TO FILTER OUT THE SPECIMENS WITH MISSING LAT-LONG!!!!!!
  filter(!is.na(x))

# Transforming mw15 into UTM -----
geo_meta <- st_as_sf(df_meta2, coords = c("x", "y"), crs = 4326)
geo_meta2<-st_transform(x = geo_meta, crs = 32616) # ask Amy about doing these steps
geo_meta2$lon<-st_coordinates(geo_meta2)[,1] # get coordinates
geo_meta2$lat<-st_coordinates(geo_meta2)[,2] # get coordinates
# geo_meta3<-st_set_geometry(geo_meta2, NULL)


# # convert data to a SpatialPointsDataFrame object
# xy <- SpatialPointsDataFrame(
#   matrix(c(x,y), ncol=2), data.frame(ID=seq(1:length(geo_meta))),
#   proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# use the distm function to generate a geodesic distance matrix in meters
# mdist <- distm(geo_meta)

mtx_distance <- st_distance(geo_meta, geo_meta)

# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mtx_distance), method="complete")

# define the distance threshold, in this case 40 m
d=10000

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
geo_meta$cluster <- cutree(hc, h=d)


df_geo <-st_set_geometry(geo_meta, NULL) %>% 
  inner_join(., df_meta) %>% 
  group_by(cluster) %>% 
  add_tally()

df_geo_check <- df_geo %>% 
  mutate(lonely = if_else(n < 3, "lonely", "grouped")) %>% 
  filter(sex == "female")
  

# EXPORT RESULT -----------------------------------------------------------

saveRDS(df_geo, "./data/data_output/output_01c_df_rpbb_clustered.Rdata")



# CHECKING STUFF ----------------------------------------------------------



df_affinis_historic <- read_csv("./data/data_raw/meta_rpbb_external/rpbb_historic_counties.csv")  %>% clean_names()


v_historic_states <- df_affinis_historic %>% distinct(state) %>% pull(state)

bg_map <- st_as_sf(map("state", regions = v_historic_states, plot = FALSE, fill = TRUE))

p_map <- ggplot() +
  geom_sf(data = bg_map, fill = "antiquewhite", alpha = 0.4, color = "grey80", size = 0.4) +
  #geom_sf(data = bg_map_extant, fill = "grey", alpha = 0.5, color = "grey80", size = 0.4, aes(text = ID)) +
  geom_jitter(data = df_geo_check, aes(x = longitude, y = latitude, fill = cluster), alpha = 0.5, color = "black", shape = 21, size = 2) +
  # annotation_scale(location = "bl", width_hint = 0.4) +
  # annotation_north_arrow(location = "bl", which_north = "true",
  #                        pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
  #                        style = north_arrow_fancy_orienteering) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = gray(0.9),
                                        linetype = "dashed",
                                        size = 0.2),
        panel.background = element_rect(fill = "white"),
        legend.position = "right")
  # labs(x = "", y = "", fill = "Region")

ggplotly(p_map)

df_geo %>% count(cluster) %>% View
