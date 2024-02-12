##%######################################################%##
#                                                          #
####          MAKING STRUCTURE OUTPUT FIGURES           ####
#                                                          #
##%######################################################%##

# PURPOSE - this script is to make the evanno plots, typical structure barplot, and the map. I then merge the map and structure barplot in Inkscape (rather than trying to wrangle that here)


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(scatterpie)
library(janitor)
library(pophelper)
library(gridExtra)
library(sf)
library(maps)
library(ggnewscale)

#affinis colors
v_affinis_colors <- c("#060200", "#F4D75D", "#D17E23")

# DATA --------------------------------------------------------------------

# GENIND Object generated in 03a; this actually should contain filtered siblings now so it's a misnomer; check and correct
gen_rpbb <- readRDS("./analyses/analyses_output/03a_rpbb_femaleNOknown_NOSibs_genind.Rdata")

#Metadata; historic affinis states for mapping
df_affinis_historic <- read_csv("./data/data_raw/meta_rpbb_external/rpbb_historic_counties.csv")  %>% clean_names()

# centroids between 100km sites
df_centroids <- readRDS("./data/data_output/output_01c_df_cluster100_centroids.Rdata")

#STRUCTURE stuff
evannosummary <- readRDS("./analyses/analyses_output/03e_structure_evanno.Rdata")
mergedk3 <- readRDS("./analyses/analyses_output/03e_structure_mergedk3.Rdata")


# EVANNO PLOTS ------------------------------------------------------------

p_evanno <- evannoMethodStructure(data=evannosummary,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p_evanno)



# STRUCTURE BARPLOT -------------------------------------------------------

# pop list

df_pops <- tibble(pop = as.character(gen_rpbb@pop))

#making dataframe from the merged runs at k=3, adding ID columns for population and individual. Pivoting to make compatible with ggplot language. Factoring to rearrange for approximate left-to-right ordering by longitude, group by individual ID for creating a column of Clust3 values for arranging (to make the plot prettier, the calculation itself is not biologically meaningful, really)

df_admixed_labeled <- bind_rows(mergedk3) %>% 
  mutate(pop = as.character(gen_rpbb@pop),
         indid = row.names(gen_rpbb@tab)) %>% 
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "cluster",
               values_to = "proportion") %>% 
  mutate(pop = factor(pop, levels = c("Twin Cities", "SE Minnesota", "Iowa City", "Decorah",  "Quad Cities", "Madison", "Central Wisconsin",  "North Illinois","Milwaukee",  "Chicago", "North Milwaukee", "Green Bay",  "Appalachian"))) %>% 
  group_by(indid) %>% 
  mutate(clust3prop = if_else(cluster == "Cluster3", proportion, 0)) %>% 
  arrange(desc(clust3prop))


# plot it

(p_structure <- ggplot(df_admixed_labeled, aes(reorder(factor(indid), -clust3prop), proportion, fill = factor(cluster))) +
  geom_col(color = "grey", linewidth = 0.01, width = 1) +
  facet_grid(~pop, scales = "free", space = "free") +
  theme_minimal() + labs(x = "", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    strip.text.x = element_text(size = 19, angle = 90, hjust = 0)
  ) +
  scale_fill_manual(values = v_affinis_colors)
)


# PIE-MAP -----------------------------------------------------------------

#create dataframe summarizing values from merged3k
df_combined_admix <- bind_rows(mergedk3) %>% 
  mutate(pop = as.character(gen_rpbb@pop),
         indid = row.names(gen_rpbb@tab)) %>% 
  group_by(pop) %>% 
  summarise(cluster1 = mean(Cluster1),
            cluster2 = mean(Cluster2),
            cluster3 = mean(Cluster3),
            n = n())

#join with centroid location data
df_cent_admix <- inner_join(df_centroids, df_combined_admix, by = c("named_cluster100" = "pop"))

# focused map of known states
focused_bg_map <- st_as_sf(map("state", regions = c("Minnesota", "Wisconsin", "Iowa", "Illinois", "Indiana", "Michigan", "Ohio", "Kentucky", "West Virginia", "Virginia"), plot = FALSE, fill = TRUE)) %>% 
  mutate(cons_unit = case_when(
    ID %in% c("minnesota", "wisconsin") ~ "CU 1",
    ID %in% c("iowa", "illinois") ~ "CU 2",
    ID %in% c("indiana", "michigan", "ohio") ~ "CU 3", 
    ID %in% c("virginia", "west virginia", "kentucky") ~ "CU 4"
  ))

p_pie_map <- ggplot() +
  geom_sf(data = focused_bg_map, alpha = 0.5, color = "grey60", size = 0.4, aes(fill = cons_unit)) +
  scale_fill_manual(values = hcl.colors(n = 4, palette = "Pastel 1")) + 
  new_scale_fill() +
  geom_point(data = df_centroids, aes(x = longitude_center, y = latitude_center)) +
  geom_scatterpie(aes(x=longitude_center, y=latitude_center, group = named_cluster100, r = log(n)/4), data = df_cent_admix, cols = c("cluster1", "cluster2", "cluster3"), alpha = 0.8) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_line(color = gray(0.9),
                                        linetype = "dashed",
                                        size = 0.2),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", fill = "Cluster") +
  scale_fill_manual(values = v_affinis_colors)

# ggplotly(p_map)
p_pie_map



# SAVE OUTPUTS ------------------------------------------------------------

saveRDS(p_evanno, "./figures/figures_output/04d_figure_evanno.Rdata")

#Need to overlay the structure plot with the map; also need to manually jitter the pie locations; need to manually repair the color of UP with #FFC5D0; then this gets merged with DAPC figure made in 04e

ggsave(p_structure, filename =  "./figures/figures_output/04d_figure_structure.svg", width = 10.5, height = 4)

ggsave(p_pie_map, filename = "./figures/figures_output/04d_figure_pie_structure_map.svg", height = 10, width = 10)
