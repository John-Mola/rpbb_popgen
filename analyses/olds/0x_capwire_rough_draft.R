# this rough draft assumes I ran 02c01 before this and it's still in the env

# library(devtools) 
# install_github("mwpennell/capwire") #update capwire following fix by Matt

# install.packages("../../../../../../../../Downloads/capwire_1.3.tar.gz", repos = NULL, type="source")

## this should become functions to do all this shit rather than copy-pasting...

library(capwire)

v_ncol_mnz20 <- df_mnz20_colonizer %>% count(family_index) %>% arrange(desc(n)) %>% pull(n) %>% as.numeric()
v_ncol_mnz21 <- df_mnz21_colonizer %>% count(family_index) %>% arrange(desc(n)) %>% pull(n) %>% as.numeric()
v_ncol_turt21 <- df_turt21_colonizer %>% count(family_index) %>% arrange(desc(n)) %>% pull(n) %>% as.numeric()

captable_mnz20 <- buildClassTable(v_ncol_mnz20)
captable_mnz21 <- buildClassTable(v_ncol_mnz21)
captable_turt21 <- buildClassTable(v_ncol_turt21)

ft_mnz20 <- fitTirm(data=captable_mnz20, max.pop=10000)
ft_mnz21 <- fitTirm(data=captable_mnz21, max.pop=10000)
ft_turt21 <- fitTirm(data=captable_turt21, max.pop=10000)

rboot_mnz20 <- bootstrapCapwire(fit = ft_mnz20, bootstraps = 100)
rboot_mnz21 <- bootstrapCapwire(fit = ft_mnz21, bootstraps = 100)
rboot_turt21 <- bootstrapCapwire(fit = ft_turt21, bootstraps = 100)

(output_mnz20 <- tibble(ml.colony.num =  rboot_mnz20[[1]], CI.lower= rboot_mnz20[[2]][[1]], CI.upper = rboot_mnz20[[2]][[2]]))
(output_mnz21 <- tibble(ml.colony.num =  rboot_mnz21[[1]], CI.lower= rboot_mnz21[[2]][[1]], CI.upper = rboot_mnz21[[2]][[2]]))
(output_turt21 <- tibble(ml.colony.num =  rboot_turt21[[1]], CI.lower= rboot_turt21[[2]][[1]], CI.upper = rboot_turt21[[2]][[2]]))


# captable <- buildClassTable(mnzoo_cols)

# res.tirm <- fitTirm(data=captable, max.pop=1000)
# res.bootstrap <- bootstrapCapwire(res.tirm, bootstraps = 100)

output <- data_frame(ml.colony.num =  res.tirm[[3]], CI.lower= res.bootstrap[[2]][[1]], CI.upper = res.bootstrap[[2]][[2]])
output
beepr::beep()

save(output,file = "./data/data_output/mnzoo_capwire_2021-12-1.Rdata")
