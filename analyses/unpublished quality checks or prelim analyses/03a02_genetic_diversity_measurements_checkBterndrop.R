##%######################################################%##
#                                                          #
####           GENETIC DIVERSITY MEASUREMENTS           ####
#                                                          #
##%######################################################%##

# PURPOSE - this script is used to calculate PAR, AR, He, Ho, Fis (i.e. genetic diversity measurements). We also fit models to test differences among sites/regions. We then save various outputs including a summary df of the measurements and model results to be used in the output markdown file.

# PACKAGES ----------------------------------------------------------------

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
library(report)

# DATA --------------------------------------------------------------------

gen_rpbb_flt <- readRDS("./analyses/analyses_output/03a01_gen_rpbb_flt.Rdata")


df_centroids <- readRDS("./data/data_output/output_01c_df_cluster100_centroids.Rdata")

toRemove= c(1)
gen_rpbb_flt =gen_rpbb_flt[loc=-toRemove]


# SUMMARY STATISTICS BY 100KM CLUSTERS ------------------------------------


# NUMBER OF SAMPLES PER SITE

v_rpbb_popN <- summary(gen_rpbb_flt$pop)

df_rpbb_popN <- tibble::enframe(v_rpbb_popN) %>% rename(N = value)


# PRIVATE ALLELES PER SITE

v_rpbb_pa <- private_alleles(gen_rpbb_flt) %>% apply(MARGIN = 1, FUN = sum)

df_rpbb_pa <- tibble::enframe(v_rpbb_pa) %>% rename(PA = value)


# ALLELIC RICHNESS PER SITE

v_rpbb_AR <- allelic.richness(genind2hierfstat(gen_rpbb_flt))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>%
  round(digits = 3)

df_rpbb_AR <- tibble::enframe(v_rpbb_AR) %>% rename(AR = value)

# JOIN THEM ALL TOGETHER

df_sum_stats <- full_join(df_rpbb_popN, df_rpbb_pa) %>% full_join(., df_rpbb_AR)
                                                                                                                      


# FINDING OBSERVED VS. EXPECTED HETEROZYGOSITY ----------------------------

# CALCULDATE BASIC STATS using hierfstat

basic_rpbb <- basic.stats(gen_rpbb_flt, diploid = TRUE)

# EXTRACT OBSERVED HETEROZYGOSITY

Ho_rpbb <- apply(basic_rpbb$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)

# EXTRACT EXPECTED HETEROZYGOSITY

He_rpbb <- apply(basic_rpbb$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)

# CREATE DATAFRAME OF OBSERVED AND EXPECTED HET

df_het_rpbb <- data.frame(Site = names(Ho_rpbb), Ho = Ho_rpbb, He = He_rpbb) %>%
  melt(id.vars = "Site") %>%
  full_join(., df_rpbb_popN, by = c("Site" = "name"))
                                                                                                                          
# FITTING A MODEL OF EXPECTED HETEROZYGOSITY ACROSS SITES -----------------

# EXTRACTING HE FOR ALL SITES, JOINING WITH CENTROITS TO HAVE LAT-LONG PER SITE
df_hs_site <- basic_rpbb$Hs %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "site", values_to = "he") %>% 
  full_join(., df_centroids, by = c("site" = "named_cluster100"))

# MODEL BY LAT-LONG
# No clear relationship, not informative

mod_he_ll <- lmer(he ~ long_m_center + lat_m_center + (1 | locus), data = df_hs_site, REML = FALSE)
summary(mod_he_ll)
mod_null_he_ll <- lmer(he ~ 1 + (1 | locus), data = df_hs_site, REML = FALSE)
anova(mod_he_ll, mod_null_he_ll)
r.squaredGLMM(mod_he_ll)
# posthoc_he <- glht(mod_he, linfct = mcp(region = "Tukey"))
# summary(posthoc_he)

#MODEL BY SITE
# This results in some significance due to differences between some sites...but it's so many comparisons and I'm not sure it's really informative in any biologically meaningful way

mod_he_site <- lmer(he ~ site + (1 | locus), data = df_hs_site, REML = FALSE)
summary(mod_he_site)
mod_null_he_site <- lmer(he ~ 1 + (1 | locus), data = df_hs_site, REML = FALSE) # not really needed but clearer this way
anova(mod_he_site, mod_null_he_site)
r.squaredGLMM(mod_he_site)
posthoc_he_site <- glht(mod_he_site, linfct = mcp(site = "Tukey"))
summary(posthoc_he_site)



# FITTING A MODEL OF ALLELIC RICHNESS BY SITE -----------------------------

# EXTRACTING HE FOR ALL SITES, JOINING WITH CENTROITS TO HAVE LAT-LONG PER SITE
df_ar_site <- allelic.richness(genind2hierfstat(gen_rpbb_flt))$Ar %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "site", values_to = "richness") %>% 
    full_join(., df_centroids, by = c("site" = "named_cluster100"))

# MODEL BY LAT-LONG
# No clear relationship, not informative

mod_ar_ll <- lmer(richness ~ long_m_center + lat_m_center + (1 | locus), data = df_ar_site, REML = FALSE)
summary(mod_ar_ll)
mod_null_ar_ll <- lmer(richness ~ 1 + (1 | locus), data = df_ar_site, REML = FALSE)
anova(mod_ar_ll, mod_null_ar_ll)
r.squaredGLMM(mod_ar_ll)
#posthoc <- glht(mod_ar, linfct = mcp(region = "Tukey"))
# summary(posthoc)

# MODEL SITE
# Same as above; This results in some significance due to differences between some sites...but it's so many comparisons and I'm not sure it's really informative in any biologically meaningful way

mod_ar_site <- lmer(richness ~ site + (1 | locus), data = df_ar_site, REML = FALSE)
summary(mod_ar_site)
mod_null_ar_site <- lmer(richness ~ 1 + (1 | locus), data = df_ar_site, REML = FALSE)
anova(mod_ar_site, mod_null_ar_site)
r.squaredGLMM(mod_ar_site)
posthoc_ar_site <- glht(mod_ar_site, linfct = mcp(site = "Tukey"))
summary(posthoc_ar_site)


# FITTING A MODEL OF FIS BY SITE ------------------------------------------

# EXTRACTING FIS FOR ALL SITES, JOINING WITH CENTROITS TO HAVE LAT-LONG PER SITE

df_fis_site <- basic_rpbb$Fis %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "site", values_to = "fis") %>% 
    full_join(., df_centroids, by = c("site" = "named_cluster100"))


# MODEL LAT-LONG
# Not informative

mod_fis_ll <- lmer(fis ~ lat_m_center + long_m_center + (1 | locus), data = df_fis_site, REML = FALSE)
summary(mod_fis_ll)
mod_null_fis_ll <- lmer(fis ~ 1 + (1 | locus), data = df_fis_site, REML = FALSE)
anova(mod_fis_ll, mod_null_fis_ll)
r.squaredGLMM(mod_fis_ll)
# posthoc_fis_ll <- glht(mod_fis_ll, linfct = mcp(site = "Tukey"))
# summary(posthoc_fis_site)


# MODEL SITE
# Same as above; This results in some significance due to differences between some sites...but it's so many comparisons and I'm not sure it's really informative in any biologically meaningful way

mod_fis_site <- lmer(fis ~ site + (1 | locus), data = df_fis_site, REML = FALSE)
summary(mod_fis_site)
mod_null_fis_site <- lmer(fis ~ 1 + (1 | locus), data = df_fis_site, REML = FALSE)
anova(mod_fis_site, mod_null_fis_site)
r.squaredGLMM(mod_fis_site)
posthoc_fis_site <- glht(mod_fis_site, linfct = mcp(site = "Tukey"))
summary(posthoc_fis_site)






# AGAIN,  BUT AT THE LEVEL OF REGION --------------------------------------

# Might be somewhat fraught, as the "region" determined here is sort of rationalized from the results of STRUCTURE/DAPC...but that seems fair compared to say using STRUCTURE to then determine some sort of a priori sorting for DAPC or whatever

# SUMMARY STATISTICS ------------------------------------------------------

#saving a version of gen_rpbb_flt for region-specific analysis
gen_rpbb_reg <- gen_rpbb_flt

# Set regional categories - I'd like to test populations within regions; only need these two cuz all others get lumped into the middle; this kind of not needed given way i assign below but whatever
v_mn <- c("Twin Cities")
v_app <- c("Appalachian")
                                                                                                                    #set population to strata
strata(gen_rpbb_reg) <- data.frame(pop = gen_rpbb_reg$pop) %>% 
  mutate(region = if_else(pop %in% v_app, "appalachian", if_else(pop %in% v_mn, "north west", "central")))

# sets the region as the population (rather than subpops/"sites")
setPop(gen_rpbb_reg) <- ~region

# NUMBER OF SAMPLES PER SITE
v_rpbb_regN <- summary(gen_rpbb_reg$strata$region)
df_rpbb_regN <- tibble::enframe(v_rpbb_regN) %>% rename(N = value)

# PRIVATE ALLELES PER SITE
v_rpbb_pa_reg <- private_alleles(gen_rpbb_reg) %>% apply(MARGIN = 1, FUN = sum)
df_rpbb_pa_reg <- tibble::enframe(v_rpbb_pa_reg) %>% rename(PA = value)

# ALLELIC RICHNESS PER SITE
v_rpbb_AR_reg <- allelic.richness(genind2hierfstat(gen_rpbb_reg))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>%
  round(digits = 3)
df_rpbb_AR_reg <- tibble::enframe(v_rpbb_AR_reg) %>% rename(AR = value)

# JOIN THEM ALL TOGETHER

df_sum_stats_reg <- full_join(df_rpbb_regN, df_rpbb_pa_reg) %>% full_join(., df_rpbb_AR_reg)



# FITTING A MODEL OF HETEROZYGOSITY ACROSS REGIONS ------------------------

basic_rpbb_reg <- basic.stats(gen_rpbb_reg, diploid = TRUE)

# ASSIGNING HE (HS) TO REGIOn BY LOCUS
df_hs_reg <- basic_rpbb_reg$Hs %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "region", values_to = "he")

# FIT MODEL OF HE ~ REGION WITH LOCUS AS RANDOM
# No significant difference/relationship
mod_he_reg <- lmer(he ~ region + (1 | locus), data = df_hs_reg, REML = FALSE)
summary(mod_he_reg)
mod_null_he_reg <- lmer(he ~ 1 + (1 | locus), data = df_hs_reg, REML = FALSE)
anova(mod_he_reg, mod_null_he_reg)
r.squaredGLMM(mod_he_reg)
posthoc_he_reg <- glht(mod_he_reg, linfct = mcp(region = "Tukey"))
summary(posthoc_he_reg)
report(mod_he_reg)                                                                                                                          


# FITTING A MODEL OF ALLELIC RICHNESS ACROSS REGIONS ----------------------

# EXTRACT VALUES BY REGION
df_ar_reg <- allelic.richness(genind2hierfstat(gen_rpbb_reg))$Ar %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "region", values_to = "ar")

# FIT MODEL, RUN COMPARISONS
# Significant with Appalachian having lower allelic richness than others -- Appalachian collections though also represent a much smaller physical space!
mod_ar_reg <- lmer(ar ~ region + (1 | locus), data = df_ar_reg, REML = FALSE)
summary(mod_ar_reg)
mod_null_ar_reg <- lmer(ar ~ 1 + (1 | locus), data = df_ar_reg, REML = FALSE)
anova(mod_ar_reg, mod_null_ar_reg)
r.squaredGLMM(mod_ar_reg)
posthoc_ar_reg <- glht(mod_ar_reg, linfct = mcp(region = "Tukey"))
summary(posthoc_ar_reg)
# ggplot(df_ar_reg, aes(x = region, y = ar, fill = region)) + geom_boxplot()
report(mod_ar_reg)
                                                                                                                       


# MODELING FIS BY REGION --------------------------------------------------

# EXTRACT VALUES

df_fis_reg <- basic_rpbb_reg$Fis %>% as_tibble(rownames = "locus") %>% pivot_longer(!locus, names_to = "region", values_to = "fis")

# FIT A MODEL, COMPARE
# Appalachian substantially lower than Central or Northwest (twin cities, really); sort of a contrasting result in my mind with AR as lower Fis value would be considered "good"...observed # of heterozygotes approaching our expected value...

mod_fis_reg <- lmer(fis ~ region + (1 | locus), data = df_fis_reg, REML = FALSE)
summary(mod_fis_reg)
mod_null_fis_reg <- lmer(fis ~ 1 + (1 | locus), data = df_fis_reg, REML = FALSE)
anova(mod_fis_reg, mod_null_fis_reg)
r.squaredGLMM(mod_fis_reg)
posthoc_fis_reg <- glht(mod_fis_reg, linfct = mcp(region = "Tukey"))
summary(posthoc_fis_reg)
# ggplot(df_fis_reg, aes(x = region, y = fis, fill = region)) + geom_boxplot()
report(mod_fis_reg)


                                                                                              
# SAVED OUTPUTS -----------------------------------------------------------

# SAVE OTHER THINGS AS NEEDED; FOR NOW JUST SAVING THE SUMMARY DATAFRAMES

# Dataframe of cluster and region level summary statistics

saveRDS(df_sum_stats, "./analyses/analyses_output/03a02_output_df_sum_stats_by_site.Rdata")
saveRDS(df_sum_stats_reg, "./analyses/analyses_output/03a02_output_df_sum_stats_by_region.Rdata")

# Dataframe for making a He/Ho figure

saveRDS(df_het_rpbb, "./analyses/analyses_output/03a02_output_df_he_ho_by_site.Rdata")


# Basic Stats output for site and region level

saveRDS(basic_rpbb, "./analyses/analyses_output/03a02_output_basic_rpbb.Rdata")
saveRDS(basic_rpbb_reg, "./analyses/analyses_output/03a02_output_basic_rpbb_reg.Rdata")

# saving the genind region

saveRDS(gen_rpbb_reg, "./analyses/analyses_output/03a02_output_gen_rpbb_reg.Rdata")

# saving df ar, he, fis for plotting

saveRDS(df_ar_reg, "./analyses/analyses_output/03a02_output_df_ar_reg.Rdata")
saveRDS(df_hs_reg, "./analyses/analyses_output/03a02_output_df_hs_reg.Rdata")
saveRDS(df_fis_reg, "./analyses/analyses_output/03a02_output_df_fis_reg.Rdata")
