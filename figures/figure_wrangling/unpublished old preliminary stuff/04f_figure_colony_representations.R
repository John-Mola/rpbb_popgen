##%######################################################%##
#                                                          #
####   PLOTTING EXAMPLES OF THE RESHUFFLING PROCEDURE   ####
#                                                          #
##%######################################################%##

#PURPOSE - this is to make a series of figures showing examples of the reshuffling procedure compared to representative datasets


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(patchwork)

# DATA --------------------------------------------------------------------

df_all_shuf_n9 <- readRDS("./analyses/analyses_output/03d_shuffle_results9.Rdata")
df_all_shuf_n11 <- readRDS("./analyses/analyses_output/03d_shuffle_results11.Rdata")
df_all_shuf_n19 <- readRDS("./analyses/analyses_output/03d_shuffle_results19.Rdata")
df_all_shuf_n42 <- readRDS("./analyses/analyses_output/03d_shuffle_results42.Rdata")



# PLOTTING FOR N=9 --------------------------------------------------------

(p_seedsave9 <- df_all_shuf_n9 %>% 
   group_by(nrows, set) %>% 
   tally() %>% 
   ggplot(., aes(x = nrows, y = n, fill = set)) + 
   geom_bar(stat = "identity", position = position_dodge()) +
   geom_vline(xintercept = 4, color = "red", linetype = 2, size = 1.5) +
   # scale_x_continuous(limits = c(25,45), breaks = c(25,30,35,40,45)) +
   theme_bw(base_size = 15) +
   labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset", title = "Seed Savers 2021", subtitle = "N=9") +
   scale_fill_viridis_d())


# PLOTTING FOR N=11 --------------------------------------------------------

(p_illin11 <- df_all_shuf_n11 %>% 
   group_by(nrows, set) %>% 
   tally() %>% 
   ggplot(., aes(x = nrows, y = n, fill = set)) + 
   geom_bar(stat = "identity", position = position_dodge()) +
   geom_vline(xintercept = 10, color = "red", linetype = 2, size = 1.5) +
   # scale_x_continuous(limits = c(25,45), breaks = c(25,30,35,40,45)) +
   theme_bw(base_size = 15) +
   labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset", title = "Illiniwek 2021", subtitle = "N=11") +
   scale_fill_viridis_d())


# PLOTTING FOR N=19 --------------------------------------------------------

(p_mnzoo19 <- df_all_shuf_n19 %>% 
   group_by(nrows, set) %>% 
   tally() %>% 
   ggplot(., aes(x = nrows, y = n, fill = set)) + 
   geom_bar(stat = "identity", position = position_dodge()) +
   geom_vline(xintercept = 13, color = "red", linetype = 2, size = 1.5) +
   # scale_x_continuous(limits = c(25,45), breaks = c(25,30,35,40,45)) +
   theme_bw(base_size = 15) +
   labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset", title = "Minnesota Zoo 2021", subtitle = "N=19") +
   scale_fill_viridis_d())


# PLOTTING FOR N=42 --------------------------------------------------------

(p_turtle42 <- df_all_shuf_n42 %>% 
   group_by(nrows, set) %>% 
   tally() %>% 
   ggplot(., aes(x = nrows, y = n, fill = set)) + 
   geom_bar(stat = "identity", position = position_dodge()) +
   geom_vline(xintercept = 27, color = "red", linetype = 2, size = 1.5) +
   # scale_x_continuous(limits = c(25,45), breaks = c(25,30,35,40,45)) +
   theme_bw(base_size = 15) +
   labs(x = "Unique Families in Reshuffle", y = "Count", fill = "Dataset", title = "Turtle Valley 2021", subtitle = "N=42") +
   scale_fill_viridis_d())


# COMBINE PLOTS -----------------------------------------------------------

(p_combined <- (p_seedsave9 + p_illin11) / (p_mnzoo19 + p_turtle42) + plot_layout(guides = 'collect'))



# SAVE OUTPUT -------------------------------------------------------------

ggsave(p_combined, filename = "./figures/figures_output/04f_figure_reshuffle_combined.svg")
