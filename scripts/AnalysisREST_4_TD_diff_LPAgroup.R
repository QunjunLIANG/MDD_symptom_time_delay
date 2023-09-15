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

# load the data
dat_td <- rio::import("inputs/REST_3_LPAgroup.xlsx") %>%
  mutate(LPAgroup = factor(LPAgroup, levels = c("severe","NA-MDD","mild","A-MDD")))

###############################################################################
#
# group difference in network
#
###############################################################################

## testing the difference in network
dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(ID, ends_with("mean"), LPAgroup, age, gender, fd_mean) %>%
  pivot_longer(cols = 3:7, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "ID", dv = 'td', between = "LPAgroup",
                 within = "network",
                 covariate = c("age","gender","fd_mean")) %>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network")
# ────────────────────────────────────────────────────────────────────────────────
#                         MS MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ────────────────────────────────────────────────────────────────────────────────
# LPAgroup            0.002 0.006   1 234 0.315  .575       .001 [.000, .020] .000
# age                 0.047 0.006   1 234 7.370  .007 **    .031 [.005, .075] .007
# gender              0.002 0.006   1 234 0.280  .597       .001 [.000, .019] .000
# fd_mean             0.010 0.006   1 234 1.579  .210       .007 [.000, .035] .001
# network             0.058 0.006   4 936 9.633 <.001 ***   .040 [.019, .059] .032
# LPAgroup * network  0.014 0.006   4 936 2.396  .049 *     .010 [.000, .020] .008
# age * network       0.028 0.006   4 936 4.555  .001 **    .019 [.005, .033] .015
# gender * network    0.004 0.006   4 936 0.677  .608       .003 [.000, .007] .002
# fd_mean * network   0.011 0.006   4 936 1.795  .128       .008 [.000, .016] .006
# ────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
# ─────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────────
# D - B ATN_mean         0.025 (0.010) 234  2.600  .010 **   0.231 [ 0.056,  0.407]
# D - B FPN_mean         0.006 (0.011) 234  0.522  .602      0.052 [-0.144,  0.247]
# D - B salience_mean   -0.000 (0.011) 234 -0.034  .973     -0.004 [-0.211,  0.204]
# D - B somMot_mean      0.003 (0.011) 234  0.294  .769      0.030 [-0.173,  0.234]
# D - B visual_mean     -0.020 (0.009) 234 -2.187  .030 *   -0.186 [-0.354, -0.018]
# ─────────────────────────────────────────────────────────────────────────────────────

###############################################################################
#
# TD & HAMD correlation
#
###############################################################################

## HAMD wave1 total score
cor_mat <- dat_td %>% filter(LPAgroup == "B" | LPAgroup == "D") %>%
  select(HAMD, ends_with("mean")) %>%
  select(1:8) %>% select(-TD_mean) %>%
  bruceR::Corr(p.adjust = "fdr")
# ────────────────────────────────────────────────────────────
#                               r      [95% CI]     p       N
# ────────────────────────────────────────────────────────────
# HAMD-somMot_mean            0.02 [-0.18, 0.21]  .796     239
# HAMD-visual_mean           -0.09 [-0.28, 0.10]  .297     239
# HAMD-salience_mean         -0.05 [-0.24, 0.15]  .622     239
# HAMD-ATN_mean               0.18 [-0.01, 0.36]  .029 *   239
# HAMD-FPN_mean              -0.06 [-0.25, 0.13]  .542     239
# HAMD-DMN_mean               0.04 [-0.16, 0.23]  .670     239

## linear model for the predictive effect of TD in sensory networks on depression factor 
model_lm <- dat_td %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  mutate_at(vars(somMot_mean, visual_mean), ~ scale(., scale = F)) %>% 
  lm(HAMD ~ somMot_mean*LPAgroup + visual_mean*LPAgroup + 
       age + gender + fd_mean, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ─────────────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p       [95% CI of β] r(partial) r(part)
# ─────────────────────────────────────────────────────────────────────────────────────────────
# somMot_mean                 0.084 (0.093)  0.899  .370     [-0.100, 0.267]      0.059   0.054
# LPAgroupA-MDD               0.371 (0.063)  5.876 <.001 *** [ 0.247, 0.495]      0.361   0.354
# visual_mean                 0.049 (0.097)  0.501  .617     [-0.143, 0.241]      0.033   0.030
# age                         0.060 (0.066)  0.905  .367     [-0.070, 0.189]      0.060   0.055
# gendermale                  0.001 (0.061)  0.015  .988     [-0.119, 0.121]      0.001   0.001
# fd_mean                    -0.059 (0.061) -0.977  .330     [-0.179, 0.060]     -0.064  -0.059
# somMot_mean:LPAgroupA-MDD  -0.049 (0.092) -0.533  .594     [-0.230, 0.132]     -0.035  -0.032
# LPAgroupA-MDD:visual_mean  -0.134 (0.098) -1.374  .171     [-0.327, 0.058]     -0.090  -0.083
# ─────────────────────────────────────────────────────────────────────────────────────────────

###############################################################################
#
# result visualization
#
###############################################################################

show_col(pal_lancet("lanonc")(4))

## barplot for the ATN mean difference among LPA groups
f <- dat_td$LPAgroup
f <- c(.7, 1.7, 2.7, 3.7, 4.7, 5.7)[as.integer(factor(f))]
p_atn_bar <- dat_td %>% 
  ggplot(aes(x = LPAgroup, y = ATN_mean)) +
  geom_boxplot(fill = c("#00468BFF","#ED0000FF",
                        "#42B540FF","#0099B4FF"),
               alpha = .6, width = .4) +
  geom_point(aes(x = f), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = c("*","***","*","*","***"), 
              textsize = 7, vjust = .7,
              y_position = c(.12,.135,.15,.165,.18), 
              xmin = c(1, 1, 2, 2, 3), 
              xmax = c(2, 3, 3, 4, 4),
              tip_length = 0) +
  ylab("Attention network") +
  ylim(-.25, .19) +
  theme_classic()
p_atn_bar

p_td_bar <- dat_td %>% 
  ggplot(aes(x = LPAgroup, y = TD_mean)) +
  geom_boxplot(fill = c("#00468BFF","#ED0000FF",
                        "#42B540FF","#0099B4FF"),
               alpha = .6, width = .4) +
  geom_point(aes(x = f), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = c("*","**","*"), 
              textsize = 7, vjust = .7,
              y_position = c(.09, .1, .11), 
              xmin = c(1, 1, 3), 
              xmax = c(2, 3, 4),
              tip_length = 0) +
  ylab("Whole-brain time delay") +
  theme_classic()
p_td_bar

## boxplot for the visual network difference between LPA group B and D
dat_bd <- dat_td %>% filter(LPAgroup == "B" | LPAgroup == "D")
f_2 <- dat_bd$LPAgroup
f_2 <- c(.7, 1.7)[as.integer(factor(f_2))]
p_visual_bd <- dat_bd %>% 
  ggplot(aes(x = LPAgroup, y = visual_mean)) +
  geom_boxplot(fill = c("#ED0000FF","#0099B4FF"),
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
  ylab("Visual network") +
  theme_classic()
p_visual_bd

## plot the ATN-HAMD total score correlation 
p_dot_ATN_hamd <- ggplot(dat_td, aes(y = ATN_mean, x = HAMD)) +
  geom_point(size = 3, alpha =.6) +
  geom_smooth(method = 'lm') + 
  xlab("HAMD score") + ylab("Attention network") +
  theme_classic() + easy_text_size(13)
p_dot_ATN_hamd

## plot the LPA group specific correlation
p_hamd_lpa_net <- ggplot(dat_hamd_cor_collect, aes(x = network, y = r_value, fill = LPAgroup)) +
  geom_bar(stat = "identity", alpha = .6,
           position = position_dodge2()) +
  scale_fill_manual(values = c("#ED0000FF","#0099B4FF")) +
  geom_hline(yintercept = 0) +
  ylab("Correlation") +
  theme_classic() + easy_move_legend(to = "top") +
  easy_text_size(13) + easy_remove_x_axis(what = "title")

###############################################################################
#
# combine the plots
#
###############################################################################

p_combine <- (p_td_bar + p_atn_bar) / ( p_visual_bd + p_dot_ATN_hamd + p_hamd_lpa_net) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(width = 900, height = 700, file = "outputs/REST_4_plot_combine.png")
print(p_combine)
dev.off()
