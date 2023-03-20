##%######################################################%##
#                                                          #
####   RUNNING AMOVA TO DETERMINE SOURCE OF VARIATION   ####
#                                                          #
##%######################################################%##

# This script is intended to run AMOVA using the poppr package to help explain how much genetic variation is due to variation within individuals vs. among populations


# PACKAGES ----------------------------------------------------------------

library(tidyverse) #not sure if needed
library(poppr)


# DATA --------------------------------------------------------------------

genind_rpbb <- readRDS("./analyses/analyses_output/03a_rpbb_femaleNOknown_NOSibs_genind.Rdata")


# ANY WRANGLING -----------------------------------------------------------

# Set regional categories - I'd like to test populations within regions; only need these two cuz all others get lumped into the middle
v_mn <- c("SE Minnesota", "Twin Cities")
v_app <- c("Appalachian")

#set population to strata
strata(genind_rpbb) <- data.frame(pop = genind_rpbb$pop) %>% 
  mutate(region = if_else(pop %in% v_app, "appalachian", if_else(pop %in% v_mn, "north west", "central")))

#checking
table(strata(genind_rpbb, ~region/pop))  # Populations; a bit goofy as the regions are really imbalanced and there's only one population within the appalachian region...and only two within the northwest with one of them being kind of giant compared to the other

#in the end, this doesn't really change the results below and given the funkiness of this (imbalanced; and sort of biased selection of regions bc I know the results of STRUCTURE) I just proceed with only pop below


# FITTING AMOVA MODEL -----------------------------------------------------

amova_rpbb <- poppr.amova(genind_rpbb, ~pop)
amova_rpbb


# SIGNIFICANCE TESTING ----------------------------------------------------

set.seed(1999)
amova_signif_rpbb   <- randtest(amova_rpbb, nrepet = 999)

plot(amova_signif_rpbb)
amova_signif_rpbb


# RANDOMIZING ASSIGNMENTS FOR SANITY CHECK --------------------------------

set.seed(9001)
genind_shuf_rpbb <- genind_rpbb
strata(genind_shuf_rpbb) <- sample_frac(data.frame(pop = genind_shuf_rpbb$pop),1)

amova_shuf_rpbb         <- poppr.amova(genind_shuf_rpbb, ~pop)
amova_shuf_rpbb

amova_shuf_signif_rpbb    <- randtest(amova_shuf_rpbb, nrepet = 999)

# this results in variation within samples/between samples still being significant, but variations between populations no longer significant. I think that makes sense??
amova_shuf_signif_rpbb
plot(amova_shuf_signif_rpbb)



# # RE-RUNNING WITHOUT WITHIN SAMPLE VARIANCE -------------------------------
# 
# 
# # FITTING AMOVA MODEL -----------------------------------------------------
# 
# amova_rpbb_nowin <- poppr.amova(genind_rpbb, ~pop, within = FALSE)
# amova_rpbb_nowin
# 
# 
# # SIGNIFICANCE TESTING ----------------------------------------------------
# 
# set.seed(1997)
# amova_signif_rpbb_nowin   <- randtest(amova_rpbb_nowin, nrepet = 999)
# 
# plot(amova_signif_rpbb_nowin)
# amova_signif_rpbb_nowin


# SAVING SOME OF THE OUTPUT -----------------------------------------------

# primary AMOVA table

saveRDS(amova_rpbb, "./analyses/analyses_output/03f_amova_primary_output.Rdata")

# significance output table

saveRDS(amova_signif_rpbb, "./analyses/analyses_output/03f_amova_signif_output.RData")

