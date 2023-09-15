##########################
#
# Subjects information description
#
# This script is used to summary 
# the selected participant obtained
# from REST 
#
# Liang Qunjun 2023/09/15

library(tidyverse)
library(NbClust)
library(ggiraphExtra)
library(ggsci)
library(ggeasy)
library(tidyLPA)
library(psych)
library(ggsignif)
library(patchwork)
library(RColorBrewer)
library(scales)
library(ggsignif)
library(patchwork)

# load the data
dat_mdd <- rio::import("inputs/REST_2_timeDelay.xlsx")
dat_cohor1 <- rio::import("inputs/Analysis2_TDp_collection.xlsx")

#############################################################
#
# calculate the factor and total score
#
#############################################################

# define the items
item_depression <- c(1,2,7,8,13)
item_anxiety <- c(9,10,11)
item_sleepness <- c(4,5,6)
fact_depression_total <- length(item_depression)*4
fact_anxiety_total <- length(item_anxiety)*4
fact_sleepness_total <- 3*2

# calculate the factor and total score
hamd_wave1 <- dat_mdd %>% select(ID, contains('HAMD_item')) 
hamd_wave1 <- hamd_wave1 %>%
  mutate(HAMD_depression = rowSums(hamd_wave1[, 1+item_depression])) %>%
  mutate(HAMD_anxiety = rowSums(hamd_wave1[, 1+item_anxiety])) %>%
  mutate(HAMD_sleepness = rowSums(hamd_wave1[, 1+item_sleepness]))

# merger the factor to the main data table
dat_mdd_factor <- dat_mdd %>% left_join(hamd_wave1) 

##############################################################
#
# Latent profile analysis - LPA
#
##############################################################

# explore the data-driven clusters
dat_mdd %>% select(starts_with("HAMD_item")) %>%
  select(1:17) %>% 
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1) %>%
  compare_solutions(statistics = c('AIC','BIC',"AWE", "CLC", "KIC"))
# Compare tidyLPA solutions:
#   
#   Model Classes AIC      BIC      AWE      CLC      KIC     
# 1     2       7712.478 7880.496 8306.970 7610.021 7767.478
# 1     3       7641.088 7867.266 8441.672 7502.861 7714.088
# 1     4       7599.961 7884.299 8607.031 7425.566 7690.961
# 
# Best model according to AIC is Model 1 with 4 classes.
# Best model according to BIC is Model 1 with 3 classes.
# Best model according to AWE is Model 1 with 2 classes.
# Best model according to CLC is Model 1 with 4 classes.
# Best model according to KIC is Model 1 with 4 classes.
# 
# An analytic hierarchy process, based on the fit indices AIC, AWE, BIC, CLC,
# and KIC (Akogul & Erisoglu, 2017), 
# suggests the best solution is Model 1 with 4 classes.

# call mixture model
model_lpa <- dat_mdd %>% select(starts_with("HAMD_item")) %>%
  select(1:17) %>% 
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1)
model_lpa


##############################################################
#
# HAMD difference among LPA groups
#
##############################################################

## merge the data
dat_mdd_factor_lpa <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "B", 
                           ifelse(Class == 2, "A",
                                  ifelse(Class == 3, "D", "C")))) %>% 
  mutate(LPAgroup = factor(LPAgroup, levels = c("A","B","C","D"),
                           labels = c("severe","NA-MDD","mild","A-MDD"))) %>%
  select(Class, LPAgroup) %>% 
  cbind(dat_mdd_factor)

## summary the demograpihc information among LPA groups 
table1::table1(~ age + gender + center + HAMD + HAMA | LPAgroup,
               dat_mdd_factor_lpa)

## test the difference in HAMD score
bruceR::MANOVA(dat_mdd_factor_lpa, dv = "HAMD", between = "LPAgroup", 
               covariate = c("age","gender","educations")) %>%
  bruceR::EMMEANS(effect = "LPAgroup", p.adjust = 'fdr')
# ──────────────────────────────────────────────────────────────────────────────
#                   MS    MSE df1 df2       F     p     η²p [90% CI of η²p]  η²G
# ──────────────────────────────────────────────────────────────────────────────
# LPAgroup    3026.612 10.677   3 295 283.469 <.001 ***   .742 [.705, .773] .742
# age            0.374 10.677   1 295   0.035  .852       .000 [.000, .007] .000
# gender         1.488 10.677   1 295   0.139  .709       .000 [.000, .001] .000
# educations     8.192 10.677   5 295   0.767  .574       .013 [.000, .028] .013
# ──────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
# ──────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df       t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - severe   -11.479 (0.631) 295 -18.192 <.001 *** -3.513 [-4.026, -3.000]
# mild - severe       -23.480 (0.819) 295 -28.655 <.001 *** -7.186 [-7.852, -6.520]
# mild - (NA-MDD)     -12.001 (0.694) 295 -17.281 <.001 *** -3.673 [-4.237, -3.108]
# (A-MDD) - severe     -8.736 (0.607) 295 -14.397 <.001 *** -2.674 [-3.167, -2.180]
# (A-MDD) - (NA-MDD)    2.743 (0.440) 295   6.233 <.001 ***  0.839 [ 0.482,  1.197]
# (A-MDD) - mild       14.744 (0.688) 295  21.436 <.001 ***  4.512 [ 3.953,  5.071]
# ──────────────────────────────────────────────────────────────────────────────────

## test the difference in HAMA score between B and D
dat_mdd_factor_lpa %>% filter(LPAgroup == "B" | LPAgroup == "D") %>%
  bruceR::TTEST(y = "HAMA", x = "LPAgroup", digits = 3)
# Results of t-test:
# ─────────────────────────────────────────────────────────────────────────────
# t df     p      Difference [95% CI]   Cohen’s d [95% CI]      BF10
# ─────────────────────────────────────────────────────────────────────────────
# 3.170 90  .002 **  4.800 [1.791, 7.808] 0.665 [0.248, 1.082] 1.600e+01
# ─────────────────────────────────────────────────────────────────────────────

##############################################################
#
# Plot the results form LPA 
#
##############################################################

# plot the profile of the subtypes
p_lpa_line <- model_lpa$model_1_class_4$estimates %>%
  filter(Category == "Means") %>% 
  mutate(item = rep(paste0("item", formatC(1:17,width = 2, flag = "0")), times = 4)) %>%
  mutate(LPAgroup = ifelse(Class == 1, "B", 
                           ifelse(Class == 2, "A",
                                  ifelse(Class == 3, "D", "C")))) %>% 
  mutate(LPAgroup = factor(LPAgroup, levels = c("A","B","C","D"),
                           labels = c("severe","NA-MDD","mild","A-MDD"))) %>%
  ggplot(aes(x = item, y = Estimate, group = LPAgroup, color = LPAgroup)) +
  geom_line() + geom_point() +
  scale_color_lancet() +
  xlab("HAMD-17 items") +
  theme_classic() + easy_remove_x_axis(what = "title") +
  easy_text_size(13) + easy_rotate_x_labels(angle = -45) +
  easy_add_legend_title("Subtype")
p_lpa_line

# bar plot for the number of patients
p_lpa_bar <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "B", 
                           ifelse(Class == 2, "A",
                                  ifelse(Class == 3, "D", "C")))) %>% 
  ggplot(aes(x = LPAgroup, fill = LPAgroup)) +
  geom_bar(width = .7, alpha = .5) +
  xlab("LPA subgroups") +
  scale_fill_lancet() + 
  coord_flip() +
  theme_classic() + easy_text_size(15)
p_lpa_bar

## visualize HAMD difference
f <- dat_mdd_factor_lpa$LPAgroup
f <- c(.7, 1.7, 2.7, 3.7)[as.integer(factor(f))]

p_lpa_boxplot <- ggplot(dat_mdd_factor_lpa, aes(x = LPAgroup, y = HAMD, fill = LPAgroup)) +
  geom_boxplot(alpha = .5, width = .3) +
  geom_point(aes(x = f), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = rep("***", times = 6), 
              textsize = 7, vjust = .7,
              y_position = c(40,42,44,46,48,50), 
              xmin = c(1, 1, 1, 2,2,3), 
              xmax = c(2, 3, 4, 3,4,4),
              tip_length = 0) +
  scale_fill_lancet() + ylab("HAMD score") + xlab("Subtype") +
  theme_classic() +
  easy_text_size(15)
p_lpa_boxplot

## visualize the HAMA difference
f1 <- dat_mdd_factor_lpa %>% filter(LPAgroup == "NA-MDD" | LPAgroup == "A-MDD") %>% .$LPAgroup
f1 <- c(.7, 1.7)[as.integer(factor(f1))]

p_bd_hama <- ggplot(dat_mdd_factor_lpa %>% filter(LPAgroup == "NA-MDD" | LPAgroup == "A-MDD"),
                    aes(x = LPAgroup, y = HAMA)) +
  geom_boxplot(alpha = .5, width = .3, fill = c("#ED0000FF","#0099B4FF")) +
  geom_point(aes(x = f1), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = "**", 
              textsize = 7, vjust = .7,
              y_position = 44, 
              xmin = 1, 
              xmax = 2,
              tip_length = 0) +
  scale_fill_lancet() + ylab("HAMA score") + xlab("Subtype") +
  theme_classic() +
  easy_text_size(15)
p_bd_hama

#############################################################
#
# comparing HAMD and HAMA socre between the main dataset and 
# the validation dataset
#
#############################################################

data.frame(dataset = c(rep("main", times = nrow(dat_cohor1)),rep("validation", times = nrow(dat_mdd_factor))),
           HAMD = c(dat_cohor1$HAMD_wave1_total, dat_mdd_factor$HAMD)) %>%
  bruceR::TTEST(y = "HAMD", x = "dataset", digits = 3)
# Results of t-test:
# ────────────────────────────────────────────────────────────────────────────────
# t  df     p         Difference [95% CI]      Cohen’s d [95% CI]      BF10
# ────────────────────────────────────────────────────────────────────────────────
# -7.254 432 <.001 *** -5.000 [-6.355, -3.645] -0.764 [-0.970, -0.557] 3.560e+09
# ────────────────────────────────────────────────────────────────────────────────

data.frame(dataset = c(rep("main", times = nrow(dat_cohor1)),rep("validation", times = nrow(dat_mdd_factor))),
           HAMD = c(dat_cohor1$HAMA_wave1_total, dat_mdd_factor$HAMA)) %>%
  bruceR::TTEST(y = "HAMD", x = "dataset", digits = 3)
# Results of t-test:
# ──────────────────────────────────────────────────────────────────────────────
# t  df     p         Difference [95% CI]      Cohen’s d [95% CI]      BF10
# ──────────────────────────────────────────────────────────────────────────────
# -3.401 263 <.001 *** -3.517 [-5.553, -1.481] -0.418 [-0.660, -0.176] 3.013e+01
# ──────────────────────────────────────────────────────────────────────────────

## plot the comparison between HAMD
plot_compa_hamd <- data.frame(HAMD = c(dat_cohor1$HAMD_wave1_total,
                                       validation = dat_mdd_factor$HAMD),
                              dataset = c(rep("main", times = length(dat_cohor1$HAMD_wave1_total)),
                                          rep("validation", times = length(dat_mdd_factor$HAMD)))) %>%
  ggplot(aes(x = HAMD, fill = dataset)) +
  geom_histogram(aes(y =..density..), bins = 40, alpha = .6, color = 'white') +
  geom_density(aes(y =-..density..), alpha = .6, color = "white") +
  xlab("HAMD total score") +
  theme_classic() + easy_text_size(13)
plot_compa_hamd

## plot the comparison between HAMA
plot_compa_hama <- data.frame(HAMD = c(dat_cohor1$HAMA_wave1_total,
                                       validation = dat_mdd_factor$HAMA),
                              dataset = c(rep("main", times = length(dat_cohor1$HAMD_wave1_total)),
                                          rep("validation", times = length(dat_mdd_factor$HAMD)))) %>%
  ggplot(aes(x = HAMD, fill = dataset)) +
  geom_histogram(aes(y =..density..), bins = 40, alpha = .6, color = 'white') +
  geom_density(aes(y =-..density..), alpha = .6, color = "white") +
  xlab("HAMA total score") +
  theme_classic() + easy_text_size(13) + easy_remove_legend()
plot_compa_hama

##############################################################
#
# HAMD & HAMA correlation
#
##############################################################

cor.test(dat_mdd$HAMA, dat_mdd$HAMD)
# data:  dat_mdd$HAMA and dat_mdd$HAMD
# t = 6.71, df = 64, p-value = 5.988e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.4743840 0.7655869
# sample estimates:
#   cor 
# 0.642632

p_hama_hamd <- ggplot(data = dat_mdd, aes(x = HAMA, y = HAMD)) +
  geom_point(size = 3, alpha = .6) +
  geom_smooth(method = 'lm') +
  theme_classic() + 
  easy_text_size(15)
p_hama_hamd

##############################################################
#
# combine the plots
#
##############################################################
layout_use <- "
AAAAA
AAAAA
BBBCC
BBBCC
"
p_combine <- p_lpa_line + p_lpa_boxplot + p_bd_hama +
  plot_layout(guides = "collect", design = layout_use) +
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(width = 2600, height = 2000, 
             file = "outputs/Fig3.png", dpi = 300)
print(p_combine)
dev.off()

rio::export(dat_mdd_factor_lpa, file = "inputs/REST_3_LPAgroup.xlsx")



