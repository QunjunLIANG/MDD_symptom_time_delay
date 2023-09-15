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
library(emmeans)
library(ggeasy)
library(ggsci)
library(cowplot)
library(scales)
library(ggsignif)
library(patchwork)
library(sjPlot)
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
    left_join(sbj_info) 
}

## factorize the LPAgroup
dat_use <- dat_td %>%
  mutate(LPAgroup = factor(LPAgroup, levels = c("I-MDD","NA-MDD","NI-MDD","A-MDD")))
rio::export(dat_use, file = outfile_name)

## check for missing value 
VIM::matrixplot(dat_td)

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
Cairo::Cairo(width = 1600, height = 1400, dpi = 300, file = "outputs/Power264.png")
print(p_power264)
dev.off()


#############################################################################
#
# Testing TD difference among networks in subgroups
#
############################################################################

## testing the difference between B and D, which is our interest
dat_use %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(participant_id, LPAgroup, ends_with("mean"), QC_bold_fd_mean, age, gender) %>%
  select(-TD_mean) %>%
  pivot_longer(cols = 3:8, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "LPAgroup", within = 'network',
                 covariate = c("age","gender","QC_bold_fd_mean"))%>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network", p.adjust = "fdr")
# ───────────────────────────────────────────────────────────────────────────────────────
#                                MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ───────────────────────────────────────────────────────────────────────────────────────
# LPAgroup                   0.003 0.002   1  65 1.255  .267       .019 [.000, .105] .002
# age                        0.008 0.002   1  65 3.459  .067 .     .051 [.000, .160] .007
# gender                     0.011 0.002   1  65 4.963  .029 *     .071 [.004, .189] .009
# QC_bold_fd_mean            0.002 0.002   1  65 0.930  .339       .014 [.000, .095] .002
# network                    0.014 0.003   5 325 4.436 <.001 ***   .064 [.018, .100] .056
# LPAgroup * network         0.007 0.003   5 325 2.306  .044 *     .034 [.000, .060] .030
# age * network              0.009 0.003   5 325 2.928  .013 *     .043 [.005, .072] .038
# gender * network           0.001 0.003   5 325 0.301  .912       .005 [.000, .004] .004
# QC_bold_fd_mean * network  0.009 0.003   5 325 2.707  .021 *     .040 [.003, .068] .035
# ───────────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
# ──────────────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E. df      t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────────────────────────
# (A-MDD) - (NA-MDD) ATN_mean         0.001 (0.014) 65  0.101  .920      0.018 [-0.339,  0.375]
# (A-MDD) - (NA-MDD) DMN_mean        -0.004 (0.010) 65 -0.369  .713     -0.048 [-0.306,  0.210]
# (A-MDD) - (NA-MDD) FPN_mean         0.010 (0.017) 65  0.577  .566      0.126 [-0.311,  0.564]
# (A-MDD) - (NA-MDD) salience_mean    0.022 (0.018) 65  1.192  .237      0.282 [-0.190,  0.754]
# (A-MDD) - (NA-MDD) somMot_mean     -0.026 (0.012) 65 -2.087  .041 *   -0.333 [-0.652, -0.014]
# (A-MDD) - (NA-MDD) visual_mean     -0.037 (0.012) 65 -3.038  .003 **  -0.483 [-0.800, -0.165]
# ──────────────────────────────────────────────────────────────────────────────────────────────

## controlling analysis: testing if the same effect could be found in A and C
dat_td %>% filter(LPAgroup == "I-MDD" | LPAgroup == "NI-MDD") %>%
  select(participant_id, LPAgroup, ends_with("mean"), 
         QC_bold_fd_mean, age, gender, educations) %>%
  select(-TD_mean) %>%
  pivot_longer(cols = 3:8, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "LPAgroup", within = 'network',
                 covariate = c("age","gender","QC_bold_fd_mean"),
                 sph.correction="GG")
# ─────────────────────────────────────────────────────────────────────────────────────────────
#                              MS   MSE   df1     df2     F     p     η²p [90% CI of η²p]  η²G
# ─────────────────────────────────────────────────────────────────────────────────────────────
# LPAgroup                   0.005 0.004 1.000  53.000 1.343  .252       .025 [.000, .130] .004
# age                        0.001 0.004 1.000  53.000 0.309  .581       .006 [.000, .082] .001
# gender                     0.001 0.004 1.000  53.000 0.231  .633       .004 [.000, .075] .001
# QC_bold_fd_mean            0.001 0.004 1.000  53.000 0.330  .568       .006 [.000, .083] .001
# network                    0.021 0.005 3.977 210.789 4.220  .003 **    .074 [.016, .124] .063
# LPAgroup * network         0.006 0.005 3.977 210.789 1.247  .292       .023 [.000, .050] .019
# age * network              0.006 0.005 3.977 210.789 1.243  .294       .023 [.000, .050] .019
# gender * network           0.007 0.005 3.977 210.789 1.385  .241       .025 [.000, .054] .021
# QC_bold_fd_mean * network  0.009 0.005 3.977 210.789 1.788  .133       .033 [.000, .066] .028
# ─────────────────────────────────────────────────────────────────────────────────────────────

#############################################################################
#
# plot the mean-difference testing results
#
############################################################################

## plot the contrast of somatomotor network
show_col(pal_lancet("lanonc")(4)) # obtain the color for B, D

f_log1 <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>% 
  .$LPAgroup
f_log1 <- c(.75, 1.75)[as.integer(factor(f_log1))]

p_log_dot_somMot <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>% 
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
  ylab("Time delay") +  xlab("Subtype") +
  ggtitle("Somatomotor") +
  theme_classic() +
  easy_text_size(13)
p_log_dot_somMot

## plot the visual network
f_log2 <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>% 
  .$LPAgroup
f_log2 <- c(.75, 1.75)[as.integer(factor(f_log2))]

p_log_dot_visual <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>% 
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
  ylim(c(-.11, .12)) + xlab("Subtype") +
  ylab("Time delay") + ggtitle("Visual") +
  theme_classic() +
  easy_text_size(13)
p_log_dot_visual

############################################################################
#
# The correlation between symptom and TD
#
############################################################################

## linear model for the predictive effect of TD in sensory networks on depression factor 
model_lm <- dat_use %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  mutate_at(vars(somMot_mean, visual_mean), ~ scale(., scale = F)) %>% # demean the time delay
  lm(HAMD_wave1_total ~ somMot_mean*LPAgroup + visual_mean*LPAgroup + 
       age + gender + QC_bold_fd_mean, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ──────────────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ──────────────────────────────────────────────────────────────────────────────────────────────
# somMot_mean                 0.023 (0.145)  0.162  .872     [-0.266,  0.313]      0.021   0.018
# LPAgroupA-MDD               0.060 (0.131)  0.457  .650     [-0.202,  0.321]      0.058   0.051
# visual_mean                 0.054 (0.196)  0.273  .786     [-0.338,  0.445]      0.035   0.031
# age                         0.232 (0.119)  1.953  .055 .   [-0.005,  0.469]      0.243   0.219
# gendermale                  0.107 (0.125)  0.857  .395     [-0.143,  0.358]      0.109   0.096
# QC_bold_fd_mean            -0.064 (0.125) -0.514  .609     [-0.314,  0.186]     -0.066  -0.058
# somMot_mean:LPAgroupA-MDD  -0.312 (0.143) -2.175  .033 *   [-0.598, -0.025]     -0.268  -0.244
# LPAgroupA-MDD:visual_mean  -0.154 (0.197) -0.779  .439     [-0.548,  0.241]     -0.099  -0.087
# ──────────────────────────────────────────────────────────────────────────────────────────────

### post-hoc analysis for the predictive effect of somMot
emtrends(model_lm, pairwise ~ LPAgroup, var="somMot_mean")
# $emtrends
# LPAgroup somMot_mean.trend   SE df lower.CL upper.CL
# NA-MDD                1.98 12.2 61    -22.5    26.47
# A-MDD               -45.86 18.9 61    -83.8    -7.97
# 
# Results are averaged over the levels of: gender 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast           estimate SE df t.ratio p.value
# (NA-MDD) - (A-MDD)     47.8 22 61   2.175  0.0335
# 
# Results are averaged over the levels of: gender 

interactions::sim_slopes(model_lm, johnson_neyman = F,
                         pred = "somMot_mean", modx = "LPAgroup")
# SIMPLE SLOPES ANALYSIS 
# 
# Slope of somMot_mean when LPAgroup = A-MDD: 
#   
#   Est.    S.E.   t val.      p
# -------- ------- -------- ------
#   -45.86   18.95    -2.42   0.02
# 
# Slope of somMot_mean when LPAgroup = NA-MDD: 
#   
#   Est.    S.E.   t val.      p
# ------ ------- -------- ------
#   1.98   12.25     0.16   0.87

## test if the same effect would be found in I-MDD and NI-MDD
model_lm <- dat_use %>% filter(LPAgroup == "I-MDD" | LPAgroup == "NI-MDD") %>%
  mutate_at(vars(somMot_mean, visual_mean), ~ scale(., scale = F)) %>% # demean the time delay
  lm(HAMD_wave1_total ~ somMot_mean*LPAgroup + visual_mean*LPAgroup + 
       age + gender + QC_bold_fd_mean, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ───────────────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ───────────────────────────────────────────────────────────────────────────────────────────────
# somMot_mean                 -0.297 (0.222) -1.334  .188     [-0.744,  0.150]     -0.187  -0.167
# LPAgroupNI-MDD               0.035 (0.130)  0.267  .790     [-0.227,  0.296]      0.038   0.033
# visual_mean                  0.033 (0.172)  0.192  .849     [-0.312,  0.378]      0.027   0.024
# age                         -0.147 (0.149) -0.989  .328     [-0.446,  0.152]     -0.140  -0.124
# gendermale                  -0.191 (0.141) -1.350  .183     [-0.475,  0.093]     -0.189  -0.169
# QC_bold_fd_mean             -0.288 (0.140) -2.060  .045 *   [-0.568, -0.007]     -0.282  -0.258
# somMot_mean:LPAgroupNI-MDD   0.060 (0.214)  0.282  .779     [-0.369,  0.490]      0.040   0.035
# LPAgroupNI-MDD:visual_mean  -0.183 (0.172) -1.064  .293     [-0.527,  0.162]     -0.150  -0.133
# ───────────────────────────────────────────────────────────────────────────────────────────────

############################################################################
#
# Plot the correlation results
#
############################################################################

## plot the LPA group specific correlation
p_hamd_lpa_cor <- dat_use %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
ggplot(., aes(x = somMot_mean, y = HAMD_wave1_total, color = LPAgroup)) +
  geom_point(size = 3, alpha = .4) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#ED0000FF","#0099B4FF")) +
  ylab("HAMD score") + xlab("Time delay") +
  theme_classic() + 
  easy_text_size(13) + easy_add_legend_title("Subtype")
p_hamd_lpa_cor

p_lm_emm <- plot_model(model_lm, type = "emm", terms=c("somMot_mean","LPAgroup")) +
  ylab("Predicted values") + xlab("Time delay") + ggtitle("Somatomotor") +
  theme_blank() + easy_remove_legend()
p_lm_emm

############################################################################
#
# treatment response prediction with logistic regression
#
############################################################################

## construct the data for modelling
dat_lm_logi <- dat_use %>%
  mutate(treatEff_recode = ifelse(treatGroup == "positive", 1, 0)) %>%
  mutate(treatEff_recode = factor(treatEff_recode, levels = c(0,1),
                                  labels = c("negative","positive")))

## build the model
model_logi <- dat_lm_logi %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  mutate_at(vars(somMot_mean, visual_mean, HAMD_wave1_total), ~ scale(., scale = F)) %>% 
  glm(treatEff_recode ~  somMot_mean*LPAgroup + visual_mean*LPAgroup + 
        age + gender + QC_bold_fd_mean + HAMD_wave1_total, 
      data = .,
      family = binomial())

## check the result
bruceR::GLM_summary(model_logi)
# Model Fit:
#   AIC = 86.331
# BIC = 108.816
# χ²(9) = 20.82, p = 0.013 *  
#   ─────── Pseudo-R² ───────
# McFadden’s R²   = 0.23887 (= 1 - logLik(model)/logLik(null.model))
# Nagelkerke’s R² = 0.36127 (= Cragg-Uhler’s R², adjusts Cox & Snell’s)
# 
# Unstandardized Coefficients:
#   Outcome Variable: treatEff_recode (family: binomial; link function: logit)
# N = 70
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#                                b     S.E.      z     p          [95% CI of b]                              OR   VIF
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# (Intercept)                  1.037 ( 1.404)  0.739  .460     [ -1.715,   3.788]                           2.820      
# somMot_mean                -13.473 ( 7.978) -1.689  .091 .   [-29.110,   2.164]                           0.000 1.815
# LPAgroupA-MDD                1.900 ( 1.352)  1.406  .160     [ -0.749,   4.549]                           6.685 3.863
# visual_mean                 10.095 ( 9.353)  1.079  .280     [ -8.236,  28.427]                       24225.727 1.756
# age                          0.047 ( 0.029)  1.614  .107     [ -0.010,   0.105]                           1.048 1.222
# gendermale                   0.750 ( 0.878)  0.854  .393     [ -0.970,   2.471]                           2.117 1.378
# QC_bold_fd_mean            -20.584 (10.913) -1.886  .059 .   [-41.973,   0.804]                           0.000 1.257
# HAMD_wave1_total            -0.014 ( 0.080) -0.181  .857     [ -0.170,   0.142]                           0.986 1.318
# somMot_mean:LPAgroupA-MDD   61.522 (29.122)  2.113  .035 *   [  4.444, 118.601] 523405233545091460204486600.000 5.326
# LPAgroupA-MDD:visual_mean  -39.845 (20.533) -1.940  .052 .   [-80.090,   0.400]                           0.000 2.611
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

## post-hoc anlaysis
emtrends(model_logi, pairwise ~ LPAgroup, var="somMot_mean")
# $emtrends
# LPAgroup somMot_mean.trend    SE  df asymp.LCL asymp.UCL
# NA-MDD               -13.5  7.98 Inf    -29.11      2.16
# A-MDD                 48.0 27.78 Inf     -6.39    102.49
# 
# Results are averaged over the levels of: gender 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast           estimate   SE  df z.ratio p.value
# (NA-MDD) - (A-MDD)    -61.5 29.1 Inf  -2.113  0.0346
# 
# Results are averaged over the levels of: gender 

interactions::sim_slopes(model_logi, johnson_neyman = F,
                         pred = "somMot_mean", modx = "LPAgroup")
# SIMPLE SLOPES ANALYSIS 
# 
# Slope of somMot_mean when LPAgroup = A-MDD: 
#   
#   Est.    S.E.   z val.      p
# ------- ------- -------- ------
#   48.05   27.78     1.73   0.08
# 
# Slope of somMot_mean when LPAgroup = NA-MDD: 
#   
#   Est.   S.E.   z val.      p
# -------- ------ -------- ------
#   -13.47   7.98    -1.69   0.09

############################################################################
#
# plot the result for prediction
#
############################################################################

p_logi_somMot <- dat_lm_logi %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(participant_id, LPAgroup, treatEff_recode, somMot_mean) %>%
  ggplot(aes(x = LPAgroup, y = somMot_mean, fill = treatEff_recode)) +
  geom_boxplot(alpha = .5) + 
  xlab('Subtypes') + ylab("Time delay") + ggtitle("Somatomotor") +
  scale_fill_manual(values = c("#2E2A2BFF", "#CF4E9CFF")) + ylim(-.15, .15) +
  theme_classic() + 
  easy_add_legend_title("Response") +
  easy_text_size(15)  
p_logi_somMot

p_logi_emm <- plot_model(model_logi, type = "emm", terms=c("somMot_mean","LPAgroup")) +
  xlab("Time delay") + ylab("Predicted probabilityes") +
  ggtitle("Somatomotor") +
  theme_blank() + easy_text_size(13) + easy_add_legend_title("Subtype")
p_logi_emm

#############################################################################
#
# combine the plots
#
############################################################################
layout_use <- "
AACC
BBCC
DDEE
DDEE
"

p_combine <- p_log_dot_somMot + p_log_dot_visual + p_lm_emm + 
  p_logi_emm + p_logi_somMot +
  plot_layout(guides = "collect", design = layout_use) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(width = 2700, height = 2300, dpi = 300,
             file = "outputs/Fig2.png")
print(p_combine)
dev.off()
