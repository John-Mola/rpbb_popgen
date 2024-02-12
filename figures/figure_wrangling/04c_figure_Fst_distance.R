##%######################################################%##
#                                                          #
####   FIGURE OF PAIRWISE FST BY GEOGRAPHIC DISTANCE    ####
#                                                          #
##%######################################################%##

#PURPOSE - this script is used to make a plot of pairwise Fst values by geographic distance


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(patchwork)

v_affinis_colors <- c("#060200", "#F4D75D", "#D17E23")


# DATA --------------------------------------------------------------------

df_fst <- readRDS("./analyses/analyses_output/03a03_df_fst_joined.Rdata")

fst.matrix <- readRDS("./analyses/analyses_output/03a03_fst_matrix.Rdata")

mod_mantel_spearman <- readRDS("./analyses/analyses_output/03a03_output_mod_mantel_spearman.Rdata")

# PLOT FST BY GEOGRAPHIC DISTANCE -----------------------------------------

(p_fst <- df_fst %>%
  ggplot(., aes(x = distance/1000, y = Fst)) +
  geom_jitter(size = 3.5, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
   annotate(x = 900, y = 0, geom = "text", label = "Spearman Rank Correlation \n Mantel statistic r = 0.58 \n P = 0.004") +
  # geom_smooth(color = "red") +
  theme_classic(base_size = 15) +
  labs(x = "Distance Between Centroids (km)", y = expression(F[ST])) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250))
 
 )


(p_fst_log <- df_fst %>%
    ggplot(., aes(x = log(distance/1000), y = Fst)) +
    geom_jitter(size = 3.5, alpha = 0.5) +
    geom_smooth(method = "lm", color = "black") +
    annotate(x = 5, y = 0.10, geom = "text", label = "Spearman Rank Correlation \n Mantel statistic r = 0.58 \n P = 0.004") +
    # geom_smooth(color = "red") +
    theme_classic(base_size = 15) +
    labs(x = "Log of Distance Between Centroids (km)", y = expression(F[ST])))


# FST MATRIX --------------------------------------------------------------

# Convert dist object to data.frame
ind = which( lower.tri(fst.matrix, diag = FALSE), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


upper_tri <- get_upper_tri(fst.matrix)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Convert minus values to zero
melted_cormat$value[melted_cormat$value < 0] = 0

melted_cormat2 <- melted_cormat %>%
  mutate(value = ifelse(Var1 != Var2, value, NA))

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(melted_cormat2$value, na.rm = TRUE) / 2

# # Heatmap

# level_order <- c("Twin Cities", "SE Minnesota", "Iowa City", "Decorah",  "Quad Cities", "Madison", "Central Wisconsin",  "North Illinois","Milwaukee",  "Chicago", "North Milwaukee", "Green Bay",  "Appalachian")

(p_fst_mat <- ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = v_affinis_colors[1], mid = "#7D6D2F", high = v_affinis_colors[2], midpoint = mid, name = fst.label, limits = c(0, max(melted_cormat2$value)), na.value = NA) +
  theme_minimal(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed() +
  labs(x = "", y = "")
)


# COMBINE -----------------------------------------------------------------

(p_fst + p_fst_mat) + plot_annotation(tag_levels = 'A')

(p_fst_log + p_fst_mat) + plot_annotation(tag_levels = 'A')

# SAVE FIGURES ------------------------------------------------------------

ggsave(p_fst, filename = "./figures/figures_output/04c_figure_fst_by_distance.svg")

ggsave(p_fst_log, filename = "./figures/figures_output/04c_figure_Fst_by_log_distance.svg")

ggsave(p_fst_mat, filename = "./figures/figures_output/04c_figure_Fst_matrix.svg")





