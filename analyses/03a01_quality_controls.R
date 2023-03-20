##%######################################################%##
#                                                          #
####               QUALITY CONTROL SCRIPT               ####
#                                                          #
##%######################################################%##

# PURPOSE: The purpose of this script is to do various "quality checks" on the genind object/loci. HWE; LD, missing loci, null alleles, etc. We then output a filtered genind object used for downstram analysis

# Various intermediate outputs are NOT saved at this time. Possibly what I'll do downstream is simply source this script into a Supplemental markdown file


# PACKAGES ----------------------------------------------------------------

library(adegenet)
library(pegas)
library(tidyverse)
library(lattice)
library(poppr)

# DATA --------------------------------------------------------------------

gen_rpbb <- readRDS("./analyses/analyses_output/03a_rpbb_femaleNOknown_NOSibs_genind.Rdata")


# HWE ---------------------------------------------------------------------

#See Shalene's paper on page 998 for inspiration on how to report this...Jha et al. 2015 Contemporary human-altered landscapes and oceanic barriers reduce bumble bee gene flow

#There is evidence that, globally, loci are out of HWE. But this is due to differences between subpopulations. Shalene handles this nicely in the paper noted above, so it's something to report...but not really anything to worry about downstream. 


#Standard test

hw.test(gen_rpbb, B=1000)

rpbbhwe.pop <- seppop(gen_rpbb) %>% lapply(hw.test, B = 0)
rpbbhwe.mat <- sapply(rpbbhwe.pop, "[", i = TRUE, j = 3) # Take the third column with all rows
alpha  <- 0.05
newmat <- rpbbhwe.mat
newmat[newmat > alpha] <- 1
levelplot(t(newmat), scales=list(x=list(rot=45)))

#correcting for multiple tests using this tutorial: https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html

# Chi-squared test: p-value
HWE.test <- data.frame(sapply(seppop(gen_rpbb), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))
HWE.test.chisq <- t(data.matrix(HWE.test))
{cat("Chi-squared test (p-values):", "\n")
  round(HWE.test.chisq,3)}

# Monte Carlo: p-value
HWE.test <- data.frame(sapply(seppop(gen_rpbb), 
                              function(ls) pegas::hw.test(ls, B=1000)[,4]))
HWE.test.MC <- t(data.matrix(HWE.test))
{cat("MC permuation test (p-values):", "\n")
  round(HWE.test.MC,3)}

#

alpha=0.05
Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                   MC=apply(HWE.test.MC<alpha, 2, mean))
#Prop.loci.out.of.HWE             # Type this line again to see results table

#

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean))
#Prop.pops.out.of.HWE             

#
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))
MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                 nrow=nrow(HWE.test.MC))

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                   MC=apply(HWE.test.MC<alpha, 1, mean),
                                   Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean),
                                   MC.fdr=apply(MC.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE             

# LD ----------------------------------------------------------------------

#LD between loci - not likely an issue here (p.rD and p.Ia > 0.05; no consistent linkage in the pairwise comparison either). 


poppr::ia(gen_rpbb, sample = 999)


# NULL ALLELES ------------------------------------------------------------

#Null allele frequences are within the range typical for pop gen analysis with microsatellites and unlikely to be an issue. Can cite https://www.nature.com/articles/6800545 and write it like they do in https://peerj.com/articles/13565/ where they say, "Nevertheless, the frequency of null alleles inferred with the methods of Brookfield and Chakraborty (0.09 and 0.13, respectively) is in line with values commonly reported in the literature and is unlikely to cause a major bias in downstream population structure analyses (Dakin & Avise, 2004)."

Null.alleles <- PopGenReport::null.all(gen_rpbb)
Null.alleles$null.allele.freq$summary2


# VARIOUS SMALL CHECKS ----------------------------------------------------

# MISSING LOCI; 
# There are no loci with <80% completeness

locimiss_rpbb = propTyped(gen_rpbb, by = "loc")
locimiss_rpbb[which(locimiss_rpbb < 0.80)] # print loci with < 80% complete genotypes


# INDIVIDUALS WITH POOR DATA 
# There are no individuals with poor data in the dataset...but that is because I also excluded them earlier in the process. So this is just a sanity check. 

indmiss_rpbb <- propTyped(gen_rpbb, by = "ind")
indmiss_rpbb[ which(indmiss_rpbb < 0.80) ] # print individuals with < 80% complete genotypes

# GENOTYPE COMPLETENESS
# remove individuals with less than 80% complete genotypes
gen_rpbb_flt <- missingno(gen_rpbb, type = "geno", cutoff = 0.20) #doesn't really do anything as they are pre-filtered already, but good to keep in for consistency

# CHECKING FOR DUPLICATES OR CLONES
#There are no duplicates or clones in the dataset. It's possible that there are prior to the COLONY filtering, so if re-running analysis keeping siblings in the dataset, double check. 

mlg(gen_rpbb_flt)

# CHECK LOCI ARE STILL POLYMORPHIC

isPoly(gen_rpbb_flt) %>% summary # they are all still polymorphic


# SAVING THE OUTPUT -------------------------------------------------------

saveRDS(gen_rpbb_flt, "./analyses/analyses_output/03a01_gen_rpbb_flt.Rdata")

