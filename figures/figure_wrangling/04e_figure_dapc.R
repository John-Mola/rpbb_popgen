##%######################################################%##
#                                                          #
####                    DAPC FIGURE                     ####
#                                                          #
##%######################################################%##

#PURPOSE - this script plots the DAPC output


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(poppr)


# DATA --------------------------------------------------------------------

dapc1 <- readRDS("./analyses/analyses_output/03g_dapc_output.Rdata")


# PLOTTING ----------------------------------------------------------------

svg(filename = "./figures/figures_output/04e_figure_dapc_scatter.svg")
p_dapc_rpbb <- scatter(dapc1, posi.da="bottomleft", bg="white", cstar=0, solid=0.8, cex=2, clabel = 0.8, col = viridis::cividis(13), cellipse = 1, legend = TRUE, posi.leg = "topleft")
dev.off()

