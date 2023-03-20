#RESAMPLING VOSNESENSKII TO SEE FAMILY DISTRIBUTION

library(tidyverse)

# datasets

df_sra <- read_csv("data/data_raw/misc/sra_captures_post_COLONY.csv")

#MNZoo 2020: 9 families from 18 individuals
#MNZoo 2021: 13 families from 16 individuals
#Turtle Valley: 27 families from 42 individuals

# for mnz20 RUN THESE IN ORDER I DID NOT CHANGE NAMES THIS IS RISKY!

type_shuffle_f_m = function(df, v_sites) {

  df %>% 
    filter(site %in% v_sites) %>% 
    sample_n(size = 18) %>% 
    distinct(ClusterIndex) %>% 
    summarise(nrows = n())

}

sra_sites <- c("JM01", "JM02", "JM03", "JM04", "JM05", "JM06", "JM07", "JM10", "JMX", "JMY", "forest_north", "forest_south", "mdw_north", "mdw_south", "north_east")

df_sra15_vos <- df_sra %>% filter(species == "vosnesenskii", year == 2015)
df_sra15_bif <- df_sra %>% filter(species == "bifarius", year == 2015)
df_sra18_bif <- df_sra %>% filter(species == "bifarius", year == 2018)

df_shuf_vos15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites))) %>% 
  bind_rows() %>% mutate(set = "vos15")

df_shuf_bif15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites))) %>% 
  bind_rows() %>% mutate(set = "bif15")

df_shuf_bif18 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites))) %>% 
  bind_rows() %>% mutate(set = "bif18")


df_all_shuf <- bind_rows(df_shuf_bif15, df_shuf_bif18, df_shuf_vos15)


(p_mnz20 <- df_all_shuf %>% 
  group_by(nrows, set) %>% 
  tally() %>% 
  ggplot(., aes(x = nrows, y = n, fill = set)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_vline(xintercept = 9, color = "red", linetype = 2, size = 1.5) +
  scale_x_continuous(limits = c(8,20), breaks = c(8,10,12,14,16,18,20)) +
  theme_bw(base_size = 15) +
  labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset") +
  scale_fill_viridis_d())


# for mnz21 RUN THESE IN ORDER I DID NOT CHANGE NAMES THIS IS RISKY!

type_shuffle_f_m = function(df, v_sites) {
  
  df %>% 
    filter(site %in% v_sites) %>% 
    sample_n(size = 16) %>% 
    distinct(ClusterIndex) %>% 
    summarise(nrows = n())
  
}

sra_sites <- c("JM01", "JM02", "JM03", "JM04", "JM05", "JM06", "JM07", "JM10", "JMX", "JMY", "forest_north", "forest_south", "mdw_north", "mdw_south", "north_east")

df_sra15_vos <- df_sra %>% filter(species == "vosnesenskii", year == 2015)
df_sra15_bif <- df_sra %>% filter(species == "bifarius", year == 2015)
df_sra18_bif <- df_sra %>% filter(species == "bifarius", year == 2018)

df_shuf_vos15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites))) %>% 
  bind_rows() %>% mutate(set = "vos15")

df_shuf_bif15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites))) %>% 
  bind_rows() %>% mutate(set = "bif15")

df_shuf_bif18 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites))) %>% 
  bind_rows() %>% mutate(set = "bif18")


df_all_shuf <- bind_rows(df_shuf_bif15, df_shuf_bif18, df_shuf_vos15)


(p_mnz21 <- df_all_shuf %>% 
  group_by(nrows, set) %>% 
  tally() %>% 
  ggplot(., aes(x = nrows, y = n, fill = set)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_vline(xintercept = 13, color = "red", linetype = 2, size = 1.5) +
  scale_x_continuous(limits = c(10,20), breaks = c(10,12,14,16,18,20)) +
  theme_bw(base_size = 15) +
  labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset") +
  scale_fill_viridis_d())


# for mnz20 RUN THESE IN ORDER I DID NOT CHANGE NAMES THIS IS RISKY!

type_shuffle_f_m = function(df, v_sites) {
  
  df %>% 
    filter(site %in% v_sites) %>% 
    sample_n(size = 42) %>% 
    distinct(ClusterIndex) %>% 
    summarise(nrows = n())
  
}

sra_sites <- c("JM01", "JM02", "JM03", "JM04", "JM05", "JM06", "JM07", "JM10", "JMX", "JMY", "forest_north", "forest_south", "mdw_north", "mdw_south", "north_east")

df_sra15_vos <- df_sra %>% filter(species == "vosnesenskii", year == 2015)
df_sra15_bif <- df_sra %>% filter(species == "bifarius", year == 2015)
df_sra18_bif <- df_sra %>% filter(species == "bifarius", year == 2018)

df_shuf_vos15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_vos, sra_sites))) %>% 
  bind_rows() %>% mutate(set = "vos15")

df_shuf_bif15 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra15_bif,sra_sites))) %>% 
  bind_rows() %>% mutate(set = "bif15")

df_shuf_bif18 <- lapply(1:1000, function(i) data.frame(iteration = i, type_shuffle_f_m(df_sra18_bif, sra_sites))) %>% 
  bind_rows() %>% mutate(set = "bif18")


df_all_shuf <- bind_rows(df_shuf_bif15, df_shuf_bif18, df_shuf_vos15)


(p_turt21 <- df_all_shuf %>% 
  group_by(nrows, set) %>% 
  tally() %>% 
  ggplot(., aes(x = nrows, y = n, fill = set)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_vline(xintercept = 27, color = "red", linetype = 2, size = 1.5) +
  scale_x_continuous(limits = c(25,45), breaks = c(25,30,35,40,45)) +
  theme_bw(base_size = 15) +
  labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset") +
  scale_fill_viridis_d())


p_mnz20
p_mnz21
p_turt21

(4+2+33+3)/1000 # =0.042; observed lower than ~96% of all resamples
