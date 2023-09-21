###############################
#
# This script is used to analyze
# the TD among LPA groups
#
# Liang Qunjun 2023/08/31

library(tidyverse)
library(ggeasy)
library(patchwork)
library(ggsci)
library(RColorBrewer)
library(scales)
library(cowplot)
library(sjPlot)
library(emmeans)
library(MatchIt)
library(neuroCombat)
library(ggsignif)

# load the data
dat_td <- rio::import("inputs/REST_3_LPAgroup.xlsx")
## summary the data
table1::table1(~ age + gender + center + HAMD + HAMA | LPAgroup, data = dat_td)

###############################################################################
#
# group difference in network
#
###############################################################################

## testing the difference in network
dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(ID, ends_with("mean"), LPAgroup, age, gender, fd_mean, HAMD) %>%
  pivot_longer(cols = 3:7, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "ID", dv = 'td', between = "LPAgroup",
                 within = "network", covariate = c("fd_mean","age","gender")) %>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network")
# ────────────────────────────────────────────────────────────────────────────────
# MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ────────────────────────────────────────────────────────────────────────────────
# LPAgroup            0.002 0.006   1 234 0.315  .575       .001 [.000, .020] .000
# fd_mean             0.010 0.006   1 234 1.579  .210       .007 [.000, .035] .001
# age                 0.047 0.006   1 234 7.370  .007 **    .031 [.005, .075] .007
# gender              0.002 0.006   1 234 0.280  .597       .001 [.000, .019] .000
# network             0.058 0.006   4 936 9.633 <.001 ***   .040 [.019, .059] .032
# LPAgroup * network  0.014 0.006   4 936 2.396  .049 *     .010 [.000, .020] .008
# fd_mean * network   0.011 0.006   4 936 1.795  .128       .008 [.000, .016] .006
# age * network       0.028 0.006   4 936 4.555  .001 **    .019 [.005, .033] .015
# gender * network    0.004 0.006   4 936 0.677  .608       .003 [.000, .007] .002
# ────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
# ───────────────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD) ATN_mean        -0.025 (0.010) 234 -2.600  .010 **  -0.231 [-0.407, -0.056]
# (NA-MDD) - (A-MDD) FPN_mean        -0.006 (0.011) 234 -0.522  .602     -0.052 [-0.247,  0.144]
# (NA-MDD) - (A-MDD) salience_mean    0.000 (0.011) 234  0.034  .973      0.004 [-0.204,  0.211]
# (NA-MDD) - (A-MDD) somMot_mean     -0.003 (0.011) 234 -0.294  .769     -0.030 [-0.234,  0.173]
# (NA-MDD) - (A-MDD) visual_mean      0.020 (0.009) 234  2.187  .030 *    0.186 [ 0.018,  0.354]
# ───────────────────────────────────────────────────────────────────────────────────────────────

###############################################################################
#
# TD & HAMD correlation
#
###############################################################################

## linear model for the predictive effect of TD in sensory networks on depression factor 
model_lm <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  mutate_at(vars(somMot_mean, visual_mean), ~ scale(., scale = F)) %>% 
  lm(HAMD ~ somMot_mean*LPAgroup + visual_mean*LPAgroup + 
       age + gender + fd_mean, data = .)        
bruceR::GLM_summary(model_lm) # check the results
# ───────────────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ───────────────────────────────────────────────────────────────────────────────────────────────
# somMot_mean                  0.018 (0.083)  0.223  .824     [-0.145,  0.181]      0.015   0.013
# LPAgroupNA-MDD              -0.371 (0.063) -5.876 <.001 *** [-0.495, -0.247]     -0.361  -0.354
# visual_mean                 -0.125 (0.079) -1.569  .118     [-0.281,  0.032]     -0.103  -0.095
# age                          0.060 (0.066)  0.905  .367     [-0.070,  0.189]      0.060   0.055
# gendermale                   0.001 (0.061)  0.015  .988     [-0.119,  0.121]      0.001   0.001
# fd_mean                     -0.059 (0.061) -0.977  .330     [-0.179,  0.060]     -0.064  -0.059
# somMot_mean:LPAgroupNA-MDD   0.043 (0.081)  0.533  .594     [-0.116,  0.202]      0.035   0.032
# LPAgroupNA-MDD:visual_mean   0.109 (0.079)  1.374  .171     [-0.047,  0.264]      0.090   0.083
# ───────────────────────────────────────────────────────────────────────────────────────────────

###############################################################################
#
# PSM for NA-MDD and A-MDD
#
###############################################################################
set.seed(123)
## slice the interest part
dat_focus <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  mutate(LPAgroup_recore = ifelse(LPAgroup == "A-MDD", 0, 1)) %>%
  mutate(LPAgroup_recore = factor(LPAgroup_recore))

## call PSM algorithm 
dat_psm <- matchit(LPAgroup_recore ~ center + age + gender + HAMD, 
                   method='optimal', distance = 'glm', 
                   ratio=1, replace=F, data=dat_focus)

## obtain the post-PSM data
dat_psm_use <- match.data(dat_psm)  

## summary the data
table1::table1(~ age + gender + center + HAMD + HAMA | LPAgroup, data = dat_psm_use)

dat_td <- dat_psm_use

###############################################################################
#
# group difference in network - PSM data
#
###############################################################################

## testing the difference in network
dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(ID, ends_with("mean"), LPAgroup, age, gender, fd_mean, HAMD) %>%
  pivot_longer(cols = 3:7, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "ID", dv = 'td', between = "LPAgroup",
                 within = "network", covariate = c("fd_mean","age","gender")) %>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network")
# ────────────────────────────────────────────────────────────────────────────────
# MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ────────────────────────────────────────────────────────────────────────────────
# LPAgroup            0.003 0.007   1 207 0.484  .488       .002 [.000, .026] .001
# fd_mean             0.007 0.007   1 207 1.119  .291       .005 [.000, .034] .001
# age                 0.034 0.007   1 207 5.221  .023 *     .025 [.002, .070] .006
# gender              0.000 0.007   1 207 0.042  .839       .000 [.000, .011] .000
# network             0.041 0.006   4 828 7.331 <.001 ***   .034 [.014, .053] .027
# LPAgroup * network  0.014 0.006   4 828 2.533  .039 *     .012 [.000, .023] .009
# fd_mean * network   0.007 0.006   4 828 1.177  .320       .006 [.000, .013] .004
# age * network       0.014 0.006   4 828 2.421  .047 *     .012 [.000, .022] .009
# gender * network    0.004 0.006   4 828 0.724  .576       .003 [.000, .008] .003
# ────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
# ───────────────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD) ATN_mean        -0.025 (0.010) 207 -2.532  .012 *   -0.239 [-0.425, -0.053]
# (NA-MDD) - (A-MDD) FPN_mean        -0.010 (0.011) 207 -0.910  .364     -0.094 [-0.297,  0.110]
# (NA-MDD) - (A-MDD) salience_mean    0.001 (0.012) 207  0.071  .944      0.008 [-0.210,  0.225]
# (NA-MDD) - (A-MDD) somMot_mean     -0.004 (0.012) 207 -0.354  .723     -0.040 [-0.261,  0.182]
# (NA-MDD) - (A-MDD) visual_mean      0.021 (0.009) 207  2.341  .020 *    0.197 [ 0.031,  0.363]
# ───────────────────────────────────────────────────────────────────────────────────────────────

###############################################################################
#
# TD & HAMD correlation
#
###############################################################################

## linear model for the predictive effect of TD in sensory networks on depression factor 
model_lm <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  mutate_at(vars(somMot_mean, visual_mean), ~ scale(., scale = F)) %>% 
  lm(HAMD ~ somMot_mean*LPAgroup + visual_mean*LPAgroup + 
       age + gender + fd_mean, data = .)        
bruceR::GLM_summary(model_lm) # check the results
# ───────────────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ───────────────────────────────────────────────────────────────────────────────────────────────
# somMot_mean                  0.095 (0.093)  1.019  .309     [-0.089,  0.278]      0.071   0.067
# LPAgroupNA-MDD              -0.281 (0.068) -4.168 <.001 *** [-0.415, -0.148]     -0.281  -0.273
# visual_mean                 -0.269 (0.097) -2.768  .006 **  [-0.461, -0.078]     -0.191  -0.181
# age                         -0.031 (0.069) -0.448  .654     [-0.166,  0.105]     -0.031  -0.029
# gendermale                   0.013 (0.066)  0.192  .848     [-0.118,  0.144]      0.013   0.013
# fd_mean                     -0.056 (0.066) -0.845  .399     [-0.186,  0.075]     -0.059  -0.055
# somMot_mean:LPAgroupNA-MDD  -0.016 (0.092) -0.170  .865     [-0.198,  0.166]     -0.012  -0.011
# LPAgroupNA-MDD:visual_mean   0.228 (0.097)  2.346  .020 *   [ 0.036,  0.419]      0.162   0.153
# ───────────────────────────────────────────────────────────────────────────────────────────────

### post-hoc analysis for the predictive effect of somMot
emtrends(model_lm, pairwise ~ LPAgroup, var="visual_mean")
# $emtrends
# LPAgroup visual_mean.trend   SE  df lower.CL upper.CL
# A-MDD               -13.39 4.84 203   -22.93    -3.85
# NA-MDD                2.27 4.57 203    -6.73    11.28
# 
# Results are averaged over the levels of: gender 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast           estimate   SE  df t.ratio p.value
# (A-MDD) - (NA-MDD)    -15.7 6.68 203  -2.346  0.0200
# 
# Results are averaged over the levels of: gender

interactions::sim_slopes(model_lm, johnson_neyman = F,
                         pred = "visual_mean", modx = "LPAgroup")
# SIMPLE SLOPES ANALYSIS 
# 
# Slope of visual_mean when LPAgroup = A-MDD: 
#   
#   Est.   S.E.   t val.      p
# -------- ------ -------- ------
#   -13.39   4.84    -2.77   0.01
# 
# Slope of visual_mean when LPAgroup = NA-MDD: 
#   
#   Est.   S.E.   t val.      p
# ------ ------ -------- ------
#   2.27   4.57     0.50   0.62

p_lm_hamd <- plot_model(model_lm, type = "emm", terms=c("visual_mean","LPAgroup")) +
  ylab("Predicted values") + xlab("Time delay") + ggtitle("HAMD") +
  theme_blank() 
p_lm_hamd

###############################################################################
#
# result visualization
#
###############################################################################

show_col(pal_lancet("lanonc")(4))

## boxplot for the visual network difference between LPA group B and D
dat_bd <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD")
f_2 <- dat_bd$LPAgroup
f_2 <- c(.7, 1.7)[as.integer(factor(f_2))]
p_visual_bd <- dat_bd %>% 
  ggplot(aes(x = LPAgroup, y = visual_mean)) +
  geom_boxplot(fill = c("#0099B4FF","#ED0000FF"),
               alpha = .6, width = .4) +
  geom_point(aes(x = f_2), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = "*", 
              textsize = 7, vjust = .7,
              y_position = .2, 
              xmin = 1, 
              xmax = 2,
              tip_length = 0) +
  ylim(-.25, .21) +
  ylab("Time delay") + xlab("Subtypes") +
  ggtitle("Visual") +
  theme_classic()
p_visual_bd

###############################################################################
#
# combine the plots
#
###############################################################################

p_combine <- p_visual_bd + p_lm_hamd +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(width = 2700, height = 1800, file = "outputs/Fig4.png", dpi = 300)
print(p_combine)
dev.off()
