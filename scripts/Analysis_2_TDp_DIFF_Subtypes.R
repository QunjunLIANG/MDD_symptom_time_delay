##########################
#
# Time delay projection  
# map for MDD-v.s.-HC contrast
#
# This script is used to explore 
# group difference of time delay
# projection map
#
# Liang Qunjun 2023/06/07

library(tidyverse)
library(bruceR)
library(ggstatsplot)
library(ggridges)
library(psych)
library(brainconn)
library(RColorBrewer)
library(ggeasy)
library(ggsci)
library(cowplot)
library(scales)
library(ggsignif)
library(patchwork)
source('scripts/function_ObtainNetID.R')
source('scripts/function_ObtainBrainData.R')
source('scripts/function_AddNetworkScore.R')

# indicate the path to the inputs and participant information
path_td <- "inputs/timeseries_TD/time_lag_estimation/"
sbj_info <- rio::import('inputs/Analysis1_subject_table.xlsx')
net_anna <- rio::import("inputs/Power264_Yeo7.xlsx")
outfile_name <- "inputs/Analysis2_TDp_collection.xlsx"

if (file.exists(outfile_name)) {
  dat_td <- rio::import(outfile_name)
}else{
  file_names <- list.files(path_td, pattern = 'sub-[0-9]*_projection_map_weighted', full.names = T)
  dat_td_raw <- ObtainBrainData(sbj_use = sbj_info$participant_id, type="TDp",
                                file_list = file_names)
  colnames(dat_td_raw)[2:ncol(dat_td_raw)] <- paste0("R",formatC(1:214, width = 3, flag = "0"))
  ## add network-level mean TDp
  dat_td_net <- AddNetworkScore_Power(dat_td_raw, net_anna)
  ## add group index
  dat_td <- dat_td_net %>%
    left_join(sbj_info %>% select(participant_id, gender, age, 
                                  educations, QC_bold_fd_mean, LPAgroup, serverity)) 
  
  VIM::matrixplot(dat_td)
}


## visualization
dat_td_matrix <- dat_td[,2:214] %>% as.matrix()
p_heatmap <- lattice::levelplot(dat_td_matrix %>% t(), 
                   xlab = 'Atlas ROIs', ylab = 'participants')
Cairo::Cairo(width = 800, height = 400, 
             file = 'outputs/Anay2_heatmap_weighed_projection_map.png')
print(p_heatmap)
dev.off()

#############################################################################
#
# Visualize using brainconn
#
############################################################################

net_anna["ROI.Name"] <- paste0("R",formatC(1:214, digits = 2, flag = "0"))
net_anna <- net_anna %>%
  rename(x.mni = "x", y.mni="y", z.mni="z")
net_anna$x.mni <- as.integer(net_anna$x.mni)
net_anna$y.mni <- as.integer(net_anna$y.mni)
net_anna$z.mni <- as.integer(net_anna$z.mni)

check_atlas(net_anna)
x <- example_unweighted_undirected
brainconn(atlas ="schaefer300_n7", conmat=x, node.size = 3, view="ortho")

net_anna <- net_anna %>% mutate(net_color = ifelse(network_new == "somMot", "cyan",
                                                     ifelse(network_new == "ATN","green",
                                                            ifelse(network_new == "DMN", "red",
                                                                   ifelse(network_new == "salience", "purple",
                                                                          ifelse(network_new == "visual", "blue",
                                                                                 ifelse(network_new == "FPN", "orange", "black"))))))) %>%
  rename(network = "network_new")

p_power264 <- brainconn(atlas = net_anna, 
          conmat = matrix(0, nrow = 214, ncol = 214),
          node.size =2.5,
          show.legend = T,
          all.nodes = T,
          view= "ortho") 
p_power264 + scale_color_aaas()

#############################################################################
#
# Testing TD difference between sever and moderated groups
#
############################################################################

# testing network-based difference between sever and moderated groups ---------------------
dat_td %>% 
  select(participant_id, serverity, ends_with("mean"), QC_bold_fd_mean, age, gender, educations) %>%
  pivot_longer(cols = 3:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "serverity", within = 'network',
                 covariate = c("age","gender","QC_bold_fd_mean"))

# focusing on B and D ------------------------------------------------------
dat_td %>% filter(serverity == "server") %>%
  select(participant_id, LPAgroup, ends_with("mean"), QC_bold_fd_mean, age, gender, educations) %>%
  pivot_longer(cols = 3:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "LPAgroup", within = 'network',
                 covariate = c("age","gender","QC_bold_fd_mean"))%>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network")
# ───────────────────────────────────────────────────────────────────────────────────────
# MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ───────────────────────────────────────────────────────────────────────────────────────
# LPAgroup                   0.004 0.003   1  65 1.464  .231       .022 [.000, .112] .003
# age                        0.009 0.003   1  65 3.457  .068 .     .050 [.000, .160] .007
# gender                     0.012 0.003   1  65 4.886  .031 *     .070 [.004, .188] .010
# QC_bold_fd_mean            0.002 0.003   1  65 0.875  .353       .013 [.000, .093] .002
# network                    0.012 0.003   6 390 4.411 <.001 ***   .064 [.020, .094] .055
# LPAgroup * network         0.006 0.003   6 390 2.297  .034 *     .034 [.001, .055] .030
# age * network              0.008 0.003   6 390 2.914  .009 **    .043 [.006, .067] .037
# gender * network           0.001 0.003   6 390 0.305  .934       .005 [.000, .002] .004
# QC_bold_fd_mean * network  0.007 0.003   6 390 2.695  .014 *     .040 [.004, .063] .035
# ───────────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
#   ────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E. df      t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────────
# D - B ATN_mean        -0.000 (0.014) 60 -0.032  .975     -0.007 [-0.430,  0.417]
# D - B DMN_mean        -0.003 (0.011) 60 -0.261  .795     -0.041 [-0.359,  0.276]
# D - B FPN_mean         0.012 (0.016) 60  0.732  .467      0.177 [-0.307,  0.662]
# D - B salience_mean    0.016 (0.019) 60  0.866  .390      0.240 [-0.315,  0.795]
# D - B somMot_mean     -0.028 (0.013) 60 -2.172  .034 *   -0.423 [-0.813, -0.034]
# D - B TD_mean         -0.008 (0.005) 60 -1.596  .116     -0.119 [-0.269,  0.030]
# D - B visual_mean     -0.035 (0.013) 60 -2.715  .009 **  -0.526 [-0.913, -0.138]
# ────────────────────────────────────────────────────────────────────────────────────
# Pooled SD for computing Cohen’s d: 0.067

# focusing on A and C ------------------------------------------------------
dat_td %>% filter(serverity == "moderated") %>%
  select(participant_id, LPAgroup, ends_with("mean"), QC_bold_fd_mean, age, gender, educations) %>%
  pivot_longer(cols = 4:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "LPAgroup", within = 'network',
                 covariate = c("age","gender","QC_bold_fd_mean"))%>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network")


#############################################################################
#
# plot the mean-difference testing results
#
############################################################################

## plot the contrast of somatomotor network
show_col(pal_lancet("lanonc")(4)) # obtain the color for B, D

f_log1 <- dat_td %>% filter(LPAgroup == "D" | LPAgroup == "B") %>% 
  .$LPAgroup
f_log1 <- c(.75, 1.75)[as.integer(factor(f_log1))]

p_log_dot_somMot <- dat_td %>% filter(LPAgroup == "D" | LPAgroup == "B") %>% 
  ggplot(aes(x = LPAgroup, y = somMot_mean)) +
  geom_boxplot(
    fill = c("#ED0000FF","#0099B4FF"),
    width = .3, alpha = .6, outlier.alpha = 0) +
  geom_point(aes(x = f_log1),  
             position = position_jitter(.03), 
             size = 3, alpha =.5) +
  geom_signif(annotations = c("*"), 
              textsize = 7, 
              vjust = .7,
              y_position = c(.12), 
              xmin = c(1), 
              xmax = c(2),
              tip_length = 0) +
  coord_flip() +
  ylim(c(-0.11, 0.13)) +
  ylab("Time delay") +  xlab("LPA group") +
  ggtitle("Somatomotor network") +
  theme_classic() +
  easy_text_size(13)
p_log_dot_somMot

## plot the visual network
f_log2 <- dat_td %>% filter(LPAgroup == "B" | LPAgroup == "D") %>% 
  .$LPAgroup
f_log2 <- c(.75, 1.75)[as.integer(factor(f_log2))]

p_log_dot_visual <- dat_td %>% filter(LPAgroup == "B" | LPAgroup == "D") %>% 
  ggplot(aes(x = LPAgroup, y = visual_mean)) +
  geom_boxplot(
    fill = c("#ED0000FF","#0099B4FF"),
    width = .3, alpha = .6, outlier.alpha = 0) +
  geom_point(aes(x = f_log2),  
             position = position_jitter(.03), 
             size = 3, alpha =.5) +
  geom_signif(annotations = c("**"), 
              textsize = 7, 
              vjust = .7,
              y_position = .11, 
              xmin = 1, xmax = 2,
              tip_length = 0) +
  coord_flip() +
  ylim(c(-.11, .12)) + xlab("LPA group") +
  ylab("Time delay") + ggtitle("Visual network") +
  theme_classic() +
  easy_text_size(13)
p_log_dot_visual

############################################################################
#
# The correlation between symptom and TD
#
############################################################################

dat_use <- dat_td %>% left_join(sbj_info)
rio::export(dat_use, file = outfile_name)

## HAMD wave1 total score
cor_mat <- dat_use %>% 
  select(HAMD_wave1_total, ends_with("mean")) %>%
  select(1:8) %>%
  bruceR::Corr(p.adjust = "fdr")
# Pearson's r and 95% confidence intervals:
# p values and 95% CIs are adjusted using the "fdr" method.
# ─────────────────────────────────────────────────────────────────
#                                     r      [95% CI]     p       N
# ─────────────────────────────────────────────────────────────────
# HAMD_wave1_total-TD_mean        -0.01 [-0.28, 0.26]  .929     128
# HAMD_wave1_total-somMot_mean    -0.17 [-0.42, 0.11]  .153     128
# HAMD_wave1_total-visual_mean    -0.15 [-0.40, 0.13]  .218     128
# HAMD_wave1_total-salience_mean   0.07 [-0.21, 0.33]  .668     128
# HAMD_wave1_total-ATN_mean        0.23 [-0.05, 0.47]  .044 *   128

## correlation between SMN ,VN and HAMD depression subfactor in group B and D
dat_hamd_cor_lpa_SMN <- data.frame()
for (lpa in c("B","D")) {
  cor_tmp <- dat_use %>% filter(LPAgroup == lpa) %>%
    cor.test(~ HAMD_wave1_depression + somMot_mean, data = .)
  dat_tmp <- data.frame(
    LPAgroup = lpa,
    network = "somMot",
    r_value = cor_tmp$estimate,
    conf_low = cor_tmp$conf.int[1],
    conf_up = cor_tmp$conf.int[2],
    p_value = cor_tmp$p.value
  )
  dat_hamd_cor_lpa_SMN <- rbind(dat_hamd_cor_lpa_SMN, dat_tmp)
}
dat_hamd_cor_lpa_SMN["p_fdr"] <- p.adjust(dat_hamd_cor_lpa_SMN$p_value, method = 'fdr')
# LPAgroup    r_value   conf_low    conf_up     p_value      p_fdr
# cor         B  0.1202648 -0.1831550  0.4027648 0.436799377 0.43679938
# cor1        D -0.5140404 -0.7517149 -0.1581855 0.007223762 0.01444752

dat_hamd_cor_lpa_VN <- data.frame()
for (lpa in c("B","D")) {
  cor_tmp <- dat_use %>% filter(LPAgroup == lpa) %>%
    cor.test(~ HAMD_wave1_depression + visual_mean, data = .)
  dat_tmp <- data.frame(
    LPAgroup = lpa,
    network = "visual",
    r_value = cor_tmp$estimate,
    conf_low = cor_tmp$conf.int[1],
    conf_up = cor_tmp$conf.int[2],
    p_value = cor_tmp$p.value
  )
  dat_hamd_cor_lpa_VN <- rbind(dat_hamd_cor_lpa_VN, dat_tmp)
}
dat_hamd_cor_lpa_VN["p_fdr"] <- p.adjust(dat_hamd_cor_lpa_VN$p_value, method = 'fdr')
# LPAgroup    r_value   conf_low     conf_up    p_value      p_fdr
# cor         B  0.0612300 -0.2400133  0.35171686 0.69296675 0.69296675
# cor1        D -0.4113547 -0.6889328 -0.02855228 0.03681651 0.07363302

dat_hamd_cor_collect <- rbind(
  dat_hamd_cor_lpa_SMN, dat_hamd_cor_lpa_VN
)
# LPAgroup network    r_value   conf_low     conf_up     p_value      p_fdr
# cor          B  somMot  0.1202648 -0.1831550  0.40276476 0.436799377 0.43679938
# cor1         D  somMot -0.5140404 -0.7517149 -0.15818545 0.007223762 0.01444752
# cor2         B  visual  0.0612300 -0.2400133  0.35171686 0.692966745 0.69296675
# cor11        D  visual -0.4113547 -0.6889328 -0.02855228 0.036816508 0.07363302

#############################################################################
#
# Plot the correlation results
#
############################################################################

## plot the ATN-HAMD total score correlation 
p_dot_ATN_hamd <- ggplot(dat_use, aes(y = ATN_mean, x = HAMD_wave1_total)) +
  geom_point(size = 3, alpha =.6) +
  geom_smooth(method = 'lm') + 
  xlab("HAMD score") + ylab("Attention network") +
  theme_classic() + easy_text_size(13)
p_dot_ATN_hamd

## plot the LPA group specific correlation
p_hamd_lpa_net <- ggplot(dat_hamd_cor_collect, aes(x = network, y = r_value, fill = LPAgroup)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  scale_fill_manual(values = c("#ED0000FF","#0099B4FF")) +
  geom_hline(yintercept = 0) +
  ylab("Correlation") +
  theme_classic() + easy_move_legend(to = "top") +
  easy_text_size(13) + easy_remove_x_axis(what = "title")

#############################################################################
#
# combine the plots
#
############################################################################

layout_use <- "
AAAADD
AAAAEE
BBBCCC
"

p_combine <- p_power264 + p_log_dot_somMot + 
  p_log_dot_visual + p_dot_ATN_hamd + p_hamd_lpa_net + 
  plot_layout(design = layout_use) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(width = 2700, height = 2300, dpi = 300,
             file = "outputs/Anay2_plot_combine.png")
print(p_combine)
dev.off()
