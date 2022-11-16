##%######################################################%##
#                                                          #
####               CREATE GENEIND OBJECT                ####
#                                                          #
##%######################################################%##

#PURPOSE - this script is used to create a geneind object from our genotype file (really from the output of COLONY combined back with our genotype file!)


###CURRENT ISSUES WITH THIS FILE!!####
#TODO - I am not satisfied with the output of COLONY. I think because the population has low genetic diversity the results are funky. Try running this analysis both with exclusion of "siblings" and without any filtering. (So once on the "colonizer" output and once on the "raw" genotype file of the same batch)

# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(adegenet)
library(poppr)
library(hierfstat)
library(reshape2)

# DATA --------------------------------------------------------------------

#THIS CURRENT DATASET is females only and excluding specimens from known colonies; to match the 02a script where this filtering is done, I simply just copied below (as that script outputs only a txt file of a later step)

# pre-COLONY siblings data

df_rpbb_fulldata <- readRDS("./data/data_output/output_01d_merged_genotypes.Rdata") %>% 
  filter(loci_w_data >= 10,
         #NOTE that this is ONLY females in this dataset now!!
         sex == "female",
         #NOTE excluding individuals from known nests! I believe this helps COLONY's assumptions
         #!!!!!WHICHNEST IS NOT A RELIABLE FILTER!!!!!
         is.na(which_nest)) %>% mutate(cluster = paste0(state, " (",cluster,")"))


# post-COLONY siblings data

df_rpbb_colonizer <- read_csv("./analyses/outputs_colony/r_colonizer/rpbb_femalesNOknown_2022-11-14-colonizeR.csv") %>% mutate(cluster = paste0(state, " (",cluster,")"))

v_rpbb_keepers <- readRDS("./analyses/outputs_colony/r_colonizer/02c_v_rpbb_keepers.Rdata")

# Loci used in analysis

v_loci_kept <- readRDS("./data/data_output/output_01d_vector_good_loci.Rdata")


# PREPARE DATA FOR CREATING GENIND OBJECT ---------------------------------

# Need to have the genotypes in dataframe with no other information, then a vector of individual IDs and a vector of site IDs. We then provide that to df2genind() to create the genind object


# select only the individual shortnames, sites, and genotypes from the dataframe
df_rpbb_simple <- df_rpbb_fulldata %>% 
  dplyr::select(internal_barcode, cluster, all_of(v_loci_kept)) #%>%
  # !!!!!!!! THIS IS WHERE THE REMOVAL OF SIBLINGS HAPPENS!!!!!!!!
  #filter(internal_barcode %in% v_rpbb_keepers)

# create vector of individual IDs
v_rpbb_shortnames <- pull(df_rpbb_simple, internal_barcode)

# create vector of site names
v_rpbb_sites <- pull(df_rpbb_simple, cluster)

# remove shortname, site names, then combine adjacent columns to get each allele into a single cell, last a little ditty to remove the "_1" in the merged columns
df_rpbb_onlyGeno <- df_rpbb_simple %>% 
  dplyr::select(-internal_barcode, -cluster)

df_rpbb_mergeGenos <- mapply(function(x, y) {
  paste(x, y, sep = ",")},
  df_rpbb_onlyGeno[ ,seq(1, ncol(df_rpbb_onlyGeno), by = 2)],
  df_rpbb_onlyGeno[ ,seq(2, ncol(df_rpbb_onlyGeno), by = 2)])

colnames(df_rpbb_mergeGenos) <- gsub(x = colnames(df_rpbb_mergeGenos), pattern = "_1", replacement = "", fixed = TRUE)  

df_rpbb_alleles <- as.data.frame(df_rpbb_mergeGenos) %>% 
  mutate_all(funs(str_replace_all(., "NA,NA", NA_character_)))



# CREATE GENIND OBJECT ----------------------------------------------------

#NOTE - one weird thing here is that there's 97 individuals in the Minneapolis region...also, no sibs are removed (yet...)
gen_rpbb = df2genind(df_rpbb_alleles, ploidy = 2, ind.names = v_rpbb_shortnames, pop = v_rpbb_sites, sep = ",")



# SAVE GENEIND OBJECT -----------------------------------------------------

saveRDS(gen_rpbb, "./analyses/analyses_output/03a_rpbb_femaleNOknown_allSibs_genind.Rdata")



# GONNA GET REAL MESSY BELOW RIGHT NOW ------------------------------------


# FILTERING LOCI, GENOTYPES, INDIVIDUALS, ETC -----------------------------


# MISSING LOCI 


locimiss_rpbb = propTyped(gen_rpbb, by = "loc")
locimiss_rpbb[which(locimiss_rpbb < 0.80)] # print loci with < 80% complete genotypes

# there are no loci with <80% completeness

# # Barplot
barplot(locimiss_rpbb, ylim = c(0,1), ylab = "Complete genotypes (proportion)", xlab = "Locus", las = 2, cex.names = 0.7)



# INDIVIDUALS WITH POOR DATA

indmiss_rpbb <- propTyped(gen_rpbb, by = "ind")
indmiss_rpbb[ which(indmiss_rpbb < 0.80) ] # print individuals with < 80% complete genotypes

# remove individuals with less than 80% complete genotypes
gen_rpbb_flt <- missingno(gen_rpbb, type = "geno", cutoff = 0.20) # this does nothing now because I filtered them in the COLONY run....


# CHECKING FOR DUPLICATES OR CLONES

mlg(gen_rpbb_flt)

# there might be one duplicate?


# CHECK LOCI ARE STILL POLYMORPHIC

isPoly(gen_rpbb_flt) %>% summary # they are all still polymorphic


### Summary Stats






# SUMMARY STATISTICS ------------------------------------------------------

# NUMBER OF SAMPLES PER SITE

v_rpbb_popN <- summary(gen_rpbb_flt$pop)

df_rpbb_popN <- tibble::enframe(v_rpbb_popN) %>% rename(N = value) 


# PRIVATE ALLELES PER SITE

v_rpbb_pa <- private_alleles(gen_rpbb_flt) %>% apply(MARGIN = 1, FUN = sum)


# v_rpbb_popN <- summary(gen_rpbb_flt$pop)

df_rpbb_pa <- tibble::enframe(v_rpbb_pa) %>% rename(PA = value) 


# ALLELIC RICHNESS PER SITE

# allelic.richness(genind2hierfstat(gen_rpbb_flt))$Ar %>%
#   apply(MARGIN = 2, FUN = mean) %>% 
#   round(digits = 3)

v_rpbb_AR <- allelic.richness(genind2hierfstat(gen_rpbb_flt))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)


# v_rpbb_popN <- summary(gen_rpbb_flt$pop)

df_rpbb_AR <- tibble::enframe(v_rpbb_AR) %>% rename(AR = value) 


df_sum_stats <- full_join(df_rpbb_popN, df_rpbb_pa) %>% full_join(., df_rpbb_AR)

df_sum_stats %>% 
  arrange(desc(N))


# HETEROZYGOSITY ----------------------------------------------------------

# CALCULDATE BASIC STATS using hierfstat
basic_rpbb <- basic.stats(gen_rpbb_flt, diploid = TRUE)

# CALCULATE OBSERVED HETEROZYGOSITY

Ho_rpbb <- apply(basic_rpbb$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
#Ho_rpbb

# CALCULATE EXPECTED HETEROZYGOSITY

He_rpbb <- apply(basic_rpbb$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
#He_rpbb

# CREATE DATAFRAME OF OBSERVED AND EXPECTED HET

df_het_rpbb <- data.frame(Site = names(Ho_rpbb), Ho = Ho_rpbb, He = He_rpbb) %>%
  melt(id.vars = "Site") %>% 
  full_join(., df_rpbb_popN, by = c("Site" = "name"))

# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])

ggplot(data = filter(df_het_rpbb, N > 3), aes(x = reorder(Site, -N), y = value, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge(), colour = "black")+ 
  geom_text(aes(y = 0.9, label = N)) +
  # scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
  scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(hetlab.o, hetlab.e))+
  labs(x = "Region", y = "Heterozygosity") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# F STATISTICS ------------------------------------------------------------

#this is kinda meaningless and very hard to understand without real names for the clusters...

# Fis PER SITE

fis_rpbb <- apply(basic_rpbb$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

# Pairwise Fst

fst_rpbb = genet.dist(gen_rpbb_flt, method = "WC84")
#fst_rpbb %>% round(digits = 3)

# this plot kinda janky below


# Convert dist object to data.frame
fst.matrix = as.matrix(fst_rpbb)
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
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Convert minus values to zero
melted_cormat$value[melted_cormat$value < 0] = 0

# melted_cormat2 <- melted_cormat %>% 
#   filter(Var1 != Var2)

melted_cormat2 <- melted_cormat %>% 
  mutate(value = ifelse(Var1 != Var2, value, NA))

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(melted_cormat2$value, na.rm = TRUE) / 2

# Heatmap

ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(melted_cormat2$value)), breaks = c(0, 0.05, 0.10, 0.15, 0.20), na.value = NA) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  labs(x = "", y = "")

