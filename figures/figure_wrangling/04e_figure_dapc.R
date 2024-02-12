##%######################################################%##
#                                                          #
####                    DAPC FIGURE                     ####
#                                                          #
##%######################################################%##

#PURPOSE - this script plots the DAPC output


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(poppr)
library(pals)


# DATA --------------------------------------------------------------------

dapc1 <- readRDS("./analyses/analyses_output/03g_dapc_output.Rdata")


# PLOTTING ----------------------------------------------------------------

#colored by scale

# svg(filename = "./figures/figures_output/04e_figure_dapc_scatter.svg")
# p_dapc_rpbb <- scatter(dapc1, posi.da="bottomleft", bg="white", cstar=0, solid=0.8, cex=2, clabel = 0.8, col = viridis::cividis(13), cellipse = 1, legend = TRUE, posi.leg = "topleft")
# dev.off()


df_color <- tibble(names = c("Twin Cities", "SE Minnesota", "Iowa City", "Decorah",  "Quad Cities", "Madison", "Central Wisconsin",  "North Illinois","Milwaukee",  "Chicago", "North Milwaukee", "Green Bay",  "Appalachian"), 
                   colors = c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00" ,"#9DCC00", "#C20088")) %>% 
  arrange(names)

v_color_sort <- df_color %>% pull(colors)

#basically have to resort manually to get the correct order, DAPC scatter not very flexible...so this makes the colors match and then manually rearrange in inkscape for aesthetic consistency

svg(filename = "./figures/figures_output/04e_figure_dapc_scatter_rainbow.svg")
p_dapc_rpbb <- scatter(dapc1, posi.da="bottomleft", bg="white", cstar=0, solid=0.8, cex=2, clabel = 0.8, col = v_color_sort, cellipse = 1, legend = TRUE, posi.leg = "topleft")
dev.off()



