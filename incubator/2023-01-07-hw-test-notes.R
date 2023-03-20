library(pegas)
library(tidyverse)

genind <- read_rds("analyses/analyses_output/03a_rpbb_femaleNOknown_allSibs_genind.Rdata")

hw.test(genind, B=1000); beepr::beep()
