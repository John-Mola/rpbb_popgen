#prepping data for bottleneck


# PACKAGES ----------------------------------------------------------------

library(graph4lg)


# DATA --------------------------------------------------------------------

gen_rpbb_reg <- readRDS("./analyses/analyses_output/03a02_output_gen_rpbb_reg.Rdata")

genind_to_genepop(gen_rpbb_reg, output = "./analyses/analyses_output/output_03h_bottleneck.txt")
