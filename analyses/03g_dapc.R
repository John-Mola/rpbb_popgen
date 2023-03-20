
##%######################################################%##
#                                                          #
####   DAPC ON POPULATIONS TO DETERMINE STRUCTURING     ####
#                                                          #
##%######################################################%##

# This script largely follows the advice from the tutorial here: https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-dapc.pdf


# PACKAGES ----------------------------------------------------------------

library(adegenet)


# DATA --------------------------------------------------------------------

genind_rpbb <- readRDS("./analyses/analyses_output/03a_rpbb_femaleNOknown_NOSibs_genind.Rdata")


# TRANSFORM TO PCA --------------------------------------------------------

#Not used here. This results in same/similar result to STRUCTURE; k=3 is the optimal solution with k=4 being similar
# grp <- find.clusters(genind_rpbb) #retained all PCs; 

# for whatever reason, was having weird downstream issues with Green Bay being last in the list despite being properly ordered elsehwere...not usre...but this fixes it!
genind_rpbb$pop <- factor(genind_rpbb$pop, levels = c("Appalachian", "Central Wisconsin", "Chicago", "Decorah", "Green Bay", "Iowa City", "Madison", "Milwaukee", "North Illinois", "North Milwaukee", "Quad Cities", "SE Minnesota", "Twin Cities"))

# run main function; the pca count is determined from xvalDapc (below); n.da = makes it so all eigenvalues are retained
dapc1 <- dapc(genind_rpbb, n.pca = 60, n.da = 20)


# PLOTTING ----------------------------------------------------------------

#my goal here is to demonstrate how the Twin Cities and Appalachian are the predominant clusters besides the smudge of all others in the middle. So, I chose two opposing colors and then a gradient within. The intention is not necessarily that these are super easy to distinguish within the smear, but that the 3 clusers are evident


# p_dapc_rpbb <- scatter(dapc1, posi.da="bottomleft", bg="white", cstar=0, solid=0.8, cex=2, clabel = 0.8, col = viridis::cividis(13), cellipse = 1, legend = TRUE, posi.leg = "topleft")

set.seed(999)
rpbbx <- xvalDapc(tab(genind_rpbb, NA.method = "mean"), pop(genind_rpbb))

#not typically run, it takes a crazy long time!!
# system.time(rpbbx <- xvalDapc(tab(genind_rpbb, NA.method = "mean"), pop(genind_rpbb),
#                               n.pca = 50:70, n.rep = 1000,
#                               parallel = "multicore", ncpus = 4L))

# SAVE RESULTS ------------------------------------------------------------

saveRDS(dapc1, "./analyses/analyses_output/03g_dapc_output.Rdata")


