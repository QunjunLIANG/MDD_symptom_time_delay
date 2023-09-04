###############################
#
# This script is used to analyze
# the network-symptom similarity
#
# Liang Qunjun 2023/08/31

library(tidyverse)
library(ggeasy)
library(patchwork)
library(tidyverse)
library(yacca)
library(coin)
library(boot)
library(ggsci)
library(ggeasy)
library(patchwork)
library(cowplot)

# load the data
dat_td <- rio::import("inputs/REST_3_LPAgroup.xlsx") 

########################################################################
#
# network-level similarity analysis 
#
#######################################################################
set.seed(1234)

# define function to calculate spearman correlation fit the bootstrap
spearman_corr <- function(data, indices){
  sampled_data <- data[indices,]
  res_cor <- cor(sampled_data$x, sampled_data$y, method = "spearman")
  return(res_cor)
}

# obtain item-level HAMD score
mat_hamd <- dat_td %>% select(starts_with("HAMD_item")) %>%
  select(1:17) %>% t() %>% cor()
p_mat_hamd <- ggcorrplot::ggcorrplot(mat_hamd, type = "lower") + 
  ggtitle("HAMD score") +
  scale_fill_viridis_c(option = "G") + 
  theme_void() + easy_center_title() + easy_add_legend_title("HAMD similarity")
p_mat_hamd

Cairo::Cairo(file = "outputs/REST_4_HAMD_similarity_matrix.png")
print(p_mat_hamd)
dev.off()

mat_hamd_vec <- mat_hamd[upper.tri(mat_hamd)] # transform as vector
dat_boot <- data.frame(x = mat_hamd_vec, y = 1) # initial data frame for bootstrap

# network-level TD mean ------------------------------------------
mat_td <- dat_td %>% 
  select(visual_mean,somMot_mean,ATN_mean,
           salience_mean,FPN_mean,DMN_mean) %>%
  t() %>% cor()

p_mat_mean <- ggcorrplot::ggcorrplot(mat_td, type = "lower") + 
  ggtitle("Mean of TD") +
  scale_fill_viridis_c(option = "G") + 
  theme_void() + easy_center_title() + easy_add_legend_title("Mean similarity")
p_mat_mean

Cairo::Cairo(file = "outputs/REST_4_mean_similarity_matrix.png")
print(p_mat_mean)
dev.off()

mat_td_vec <- mat_td[upper.tri(mat_td)]

cor.test(mat_td_vec,mat_hamd_vec, method = "spearman")
res_perm <- spearman_test(mat_td_vec ~ mat_hamd_vec) # permutation test 

## bootstrap for confidential interval
# dat_boot$y <- mat_td_vec
# res_boot <- boot::boot(data = dat_boot, statistic = spearman_corr, R = 5000)
# res_boot_ci <- boot::boot.ci(res_boot, type = "perc")

## collect results in this part
dat_res_mean <- data.frame(
  level = "network",
  index = "mean",
  cor = res_boot$t0, perm_p = coin::pvalue(res_perm),
  lower_interval = res_boot_ci$percent[4], upper_interval = res_boot_ci$percent[5])

# network-level TD sd ------------------------------------------
mat_td <- dat_td %>% select(visual_sd,somMot_sd,ATN_sd,
                            salience_sd,FPN_sd,DMN_sd) %>%
  t() %>% cor()

p_mat_sd <- ggcorrplot::ggcorrplot(mat_td, type = "lower") + 
  ggtitle("SD of TD") +
  scale_fill_viridis_c(option = "G") + 
  theme_void() +easy_center_title() + easy_add_legend_title("SD similarity")
p_mat_sd

Cairo::Cairo(file = "outputs/REST_4_sd_similarity_matrix.png")
print(p_mat_sd)
dev.off()

mat_td_vec <- mat_td[upper.tri(mat_td)]

cor.test(mat_td_vec,mat_hamd_vec, method = "spearman")
res_perm <- spearman_test(mat_td_vec ~ mat_hamd_vec) # permutation test 
res_perm

## bootstrap for confidential interval
# dat_boot$y <- mat_td_vec
# res_boot <- boot::boot(data = dat_boot, statistic = spearman_corr, R = 5000)
# res_boot_ci <- boot::boot.ci(res_boot, type = "perc")

## collect results in this part
dat_res_sd <- data.frame(
  level = "network",
  index = "sd",
  cor = res_boot$t0, perm_p = coin::pvalue(res_perm),
  lower_interval = res_boot_ci$percent[4], upper_interval = res_boot_ci$percent[5])

# collect results for all whole-brain results ----------------------------------
dat_intersbj_simi <- rbind(dat_res_mean, dat_res_sd)

p_simi_bar <- ggplot(dat_intersbj_simi, aes(x = index, y = cor)) +
  geom_bar(stat = "identity", position = position_dodge2(), 
           width = .6, fill = "steelblue", alpha = .6) +
  geom_errorbar(aes(ymin = lower_interval, ymax = upper_interval),
                width = 0.2, size = 1, position = position_dodge(.7)) +
  ylab("Inter-subject similarity") + 
  theme_classic() + easy_text_size(15) + 
  easy_move_legend(to = "top") + easy_remove_x_axis(what = 'title')
p_simi_bar

Cairo::Cairo(file = 'outputs/REST_4_similarity_barplot.png')
print(p_simi_bar)
dev.off()
