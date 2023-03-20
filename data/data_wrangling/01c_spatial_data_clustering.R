##%######################################################%##
#                                                          #
####             SPATIAL DATA CLUSTERING OF             ####
####            POINTS TO ASSIGN POPULATIONS            ####
#                                                          #
##%######################################################%##

#PURPOSE - to take the metadata with lat-long and assign them to clusters based on some threshold rule of distance. This seems preferable to assigning to population based on some more arbitrary metric (like county or just...my decision)


### Going to try to follow the StackOverflow page here: https://gis.stackexchange.com/questions/17638/clustering-spatial-data-in-r


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
  #This next filtering step no longer does anything because there are no longer specimens lacking lat-long data...but I am afraid to change it. 
  filter(!is.na(x))

# Transforming into UTM -----
geo_meta <- st_as_sf(df_meta2, coords = c("x", "y"), crs = 4326)
geo_meta2<-st_transform(x = geo_meta, crs = 32616) # ask Amy about doing these steps
geo_meta2$lon<-st_coordinates(geo_meta2)[,1] # get coordinates
geo_meta2$lat<-st_coordinates(geo_meta2)[,2] # get coordinates


# use the st_distance function to generate a geodesic distance matrix in meters

mtx_distance <- st_distance(geo_meta2, geo_meta2)

# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mtx_distance), method="complete")

# define multiple distance thresholds
d05 = 5000
d10 = 10000
d25 = 25000
d50 = 50000
d100 = 100000

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
# geo_meta$cluster <- cutree(hc, h=d)
geo_meta2$cluster05 <- cutree(hc, h=d05)
geo_meta2$cluster10 <- cutree(hc, h=d10)
geo_meta2$cluster25 <- cutree(hc, h=d25)
geo_meta2$cluster50 <- cutree(hc, h=d50)
geo_meta2$cluster100 <- cutree(hc, h=d100)


#changed this to multiple thresholds on 2022-12-15...this is going to result in stuff downstream breaking as there is now no "cluster" column

df_geo <-st_set_geometry(geo_meta2, NULL) %>% 
  inner_join(., df_meta) %>% 
  add_count(cluster05, name = "n_cluster05") %>% 
  add_count(cluster10, name = "n_cluster10") %>% 
  add_count(cluster25, name = "n_cluster25") %>% 
  add_count(cluster50, name = "n_cluster50") %>% 
  add_count(cluster100, name = "n_cluster100")

df_geo_named <- df_geo %>% 
  mutate(named_cluster100 = case_when(
    cluster100 == 1 ~ "Appalachian",
    cluster100 == 2 ~ "Iowa City",
    cluster100 == 3 ~ "Quad Cities",
    cluster100 == 4 ~ "Decorah",
    cluster100 == 5 ~ "Twin Cities",
    cluster100 == 6 ~ "SE Minnesota",
    cluster100 == 7 ~ "Milwaukee",
    cluster100 == 8 ~ "Madison",
    cluster100 == 9 ~ "North Illinois",
    cluster100 == 10 ~ "Chicago",
    cluster100 == 11 ~ "North Milwaukee",
    cluster100 == 12 ~ "Green Bay",
    cluster100 == 13 ~ "Central Wisconsin"
  ))

# df_geo_check <- df_geo %>% 
#   mutate(lonely = if_else(n < 3, "lonely", "grouped")) %>% 
#   filter(sex == "female")
  

# EXPORT RESULT -----------------------------------------------------------

saveRDS(df_geo_named, "./data/data_output/output_01c_df_rpbb_clustered.Rdata")


# FINDING CENTROIDS OF CLUSTERS -------------------------------------------

df_cluster100_centroids <- df_geo_named %>% 
  group_by(named_cluster100) %>% 
  summarise(latitude_center = mean(latitude),
            longitude_center = mean(longitude),
            lat_m_center = mean(lat),
            long_m_center = mean(lon))

saveRDS(df_cluster100_centroids, "./data/data_output/output_01c_df_cluster100_centroids.Rdata")


# FINDING PAIRWISE DISTANCE OF ALL CENTROIDS ------------------------------

df_cluster100_centroids2 <- mutate(df_cluster100_centroids, FamilyID = "one")

df_joined_clusters <- left_join(df_cluster100_centroids2, df_cluster100_centroids2, by = "FamilyID", suffix=c(".A", ".B")) %>%
  filter(named_cluster100.A != named_cluster100.B) %>%
  rowwise %>%
  mutate(name = toString(sort(c(named_cluster100.A, named_cluster100.B)))) %>%
  distinct(name, .keep_all=T) %>%
  mutate(distance = sqrt((long_m_center.A - long_m_center.B)^2+(lat_m_center.A - lat_m_center.B)^2))


saveRDS(df_joined_clusters, "./data/data_output/output_01c_df_cluster100_pw_distances.Rdata")


# crossing(df_cluster100_centroids, df_cluster100_centroids)


# CHECKING STUFF ----------------------------------------------------------


#
# df_affinis_historic <- read_csv("./data/data_raw/meta_rpbb_external/rpbb_historic_counties.csv")  %>% clean_names()
# 
# 
# v_historic_states <- df_affinis_historic %>% distinct(state) %>% pull(state)
# 
# bg_map <- st_as_sf(map("state", regions = v_historic_states, plot = FALSE, fill = TRUE))
# 
# p_map <- ggplot() +
#   geom_sf(data = bg_map, fill = "antiquewhite", alpha = 0.4, color = "grey80", size = 0.4) +
#   #geom_sf(data = bg_map_extant, fill = "grey", alpha = 0.5, color = "grey80", size = 0.4, aes(text = ID)) +
#   geom_jitter(data = filter(df_geo_named, sex == "female"), aes(x = longitude, y = latitude, fill = as.factor(named_cluster100), shape = as.factor(year)), alpha = 0.5, color = "black", size = 2) +
#   # annotation_scale(location = "bl", width_hint = 0.4) +
#   # annotation_north_arrow(location = "bl", which_north = "true",
#   #                        pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
#   #                        style = north_arrow_fancy_orienteering) +
#   geom_point(data = df_cluster100_centroids, aes(x = longitude_center, y = latitude_center)) +
#   theme_bw() +
#   theme(panel.grid.major = element_line(color = gray(0.9),
#                                         linetype = "dashed",
#                                         size = 0.2),
#         panel.background = element_rect(fill = "white"),
#         legend.position = "right")
#   # labs(x = "", y = "", fill = "Region")
# 
# ggplotly(p_map)
# 
# # df_geo %>% count(cluster) %>% View
