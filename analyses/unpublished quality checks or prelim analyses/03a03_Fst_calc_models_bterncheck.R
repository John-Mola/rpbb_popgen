##%######################################################%##
#                                                          #
####              CALCULATING F STATISTICS              ####
#                                                          #
##%######################################################%##

# Purpose - this script is used to calculate basic F statistics from the gen_ind object. Then saving some outputs for making figures/reporting results. 



# PACKAGES ----------------------------------------------------------------
# not all used, this is just copied from 03a02 for laziness...

library(adegenet)
library(pegas)
library(tidyverse)
library(lattice)
library(poppr)
library(hierfstat)
library(reshape2)
library(lme4)
library(MuMIn)
library(multcomp)
library(ecodist)



# DATA --------------------------------------------------------------------

gen_rpbb_flt <- readRDS("./analyses/analyses_output/03a01_gen_rpbb_flt.Rdata")

toRemove= c(1)
gen_rpbb_flt =gen_rpbb_flt[loc=-toRemove]

# for whatever reason, was having weird downstream issues with Green Bay being last in the list despite being properly ordered elsehwere...not usre...but this fixes it!
gen_rpbb_flt$pop <- factor(gen_rpbb_flt$pop, levels = c("Twin Cities", "SE Minnesota", "Iowa City", "Decorah",  "Quad Cities", "Madison", "Central Wisconsin",  "North Illinois","Milwaukee",  "Chicago", "North Milwaukee", "Green Bay",  "Appalachian"))

df_joined_clusters <- readRDS("./data/data_output/output_01c_df_cluster100_pw_distances.Rdata")

df_centroids <- readRDS("./data/data_output/output_01c_df_cluster100_centroids.Rdata") %>% 
  arrange(factor(named_cluster100, levels = c("Twin Cities", "SE Minnesota", "Iowa City", "Decorah",  "Quad Cities", "Madison", "Central Wisconsin",  "North Illinois","Milwaukee",  "Chicago", "North Milwaukee", "Green Bay",  "Appalachian")))

# FINDING PAIRWISE FST VALUES ---------------------------------------------

dist_fst <- genet.dist(gen_rpbb_flt, method = "WC84")

# Convert dist object to data.frame
fst.matrix <- as.matrix(dist_fst)
ind <- which( lower.tri(fst.matrix, diag = FALSE), arr.ind = TRUE)
fst.df <- data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))

# MATCHING TO CLUSTERS

df_fst <- fst.df %>% 
  rowwise() %>%
  mutate(pair_name = paste0(sort(c(Site1, Site2)), collapse = ", "))

df_flt_joined <- dplyr::select(df_joined_clusters, name, distance)

df_fst_join <- left_join(df_fst, df_flt_joined, by = c("pair_name" = "name"))
# 
# df_fst_join %>%
#   ggplot(., aes(x = distance/1000, y = Fst)) +
#   geom_jitter(size = 3.5, alpha = 0.5) +
#   geom_smooth(method = "lm") +
#   # geom_smooth(color = "red") +
#   theme_classic(base_size = 15) +
#   labs(x = "Distance Between Centroids (km)", y = "Fst")


# FITTING SIMPLE LINEAR MODELS --------------------------------------------

# lm(Fst ~ distance, data = df_fst_join) %>% summary()
# lm(Fst ~ log(distance), data = df_fst_join) %>% summary()


# MANTEL TEST -------------------------------------------------------------

  
dist_sites <- dist(cbind(df_centroids$lat_m_center, df_centroids$long_m_center))
  
#because of how data are distributed, it might make more sense to do the test on a log scale for distance...another approach is to use the spearman method instead
(mod_mantel_pearson <- vegan::mantel(dist_fst, dist_sites))
(mod_mantel_log <- vegan::mantel(dist_fst, log(dist_sites)))
(mod_mantel_spearman <- vegan::mantel(dist_fst, dist_sites, method = "spearman"))



# SAVE OUTPUTS ------------------------------------------------------------

# dataframe to be able to do scatterplot
saveRDS(df_fst_join, "./analyses/analyses_output/03a03_df_fst_joined.Rdata")

# save model outputs
saveRDS(mod_mantel_pearson, "./analyses/analyses_output/03a03_output_mod_mantel_pearson.Rdata")
saveRDS(mod_mantel_log, "./analyses/analyses_output/03a03_output_mod_mantel_log.Rdata")
saveRDS(mod_mantel_spearman, "./analyses/analyses_output/03a03_output_mod_mantel_spearman.Rdata")

# Fst matrix for making pairwise matrix plot

saveRDS(fst.matrix, "./analyses/analyses_output/03a03_fst_matrix.Rdata")

