##########################
#
# Subjects information description
#
# This script is used to summary 
# the selected participant obtained
# from REST 
#
# Liang Qunjun 2023/08/30

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

#############################################################
#
# comparing HAMD and HAMA socre between the main dataset and 
# the validation dataset
#
#############################################################

t.test(dat_cohor1$HAMD_wave1_total, dat_mdd_factor$HAMD)
# data:  dat_cohor1$HAMD_wave1_total and dat_mdd_factor$HAMD
# t = 7.0474, df = 224.89, p-value = 2.207e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   3.588614 6.374435
# sample estimates:
#   mean of x mean of y 
# 26.35938  21.37785 

t.test(dat_cohor1$HAMA_wave1_total, dat_mdd_factor$HAMA)
# data:  dat_cohor1$HAMA_wave1_total and dat_mdd_factor$HAMA
# t = 3.4284, df = 255.79, p-value = 0.0007075
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1.496753 5.536892
# sample estimates:
#   mean of x mean of y 
# 20.86719  17.35036 

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

# plot the result
p_lpa_line <- model_lpa$model_1_class_4$estimates %>%
  filter(Category == "Means") %>% 
  mutate(item = rep(paste0("item", formatC(1:17,width = 2, flag = "0")), times = 4)) %>%
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                           ifelse(Class == 2, "B",
                                  ifelse(Class == 3, "D", "C")))) %>% 
  ggplot(aes(x = item, y = Estimate, group = LPAgroup, color = LPAgroup)) +
  geom_line() + geom_point() +
  scale_color_lancet() +
  xlab("HAMD-17 items") +
  theme_classic() + easy_remove_x_axis(what = "title") +
  easy_text_size(20) + easy_rotate_x_labels(angle = -45)
p_lpa_line

p_lpa_bar <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                           ifelse(Class == 2, "B",
                                  ifelse(Class == 3, "D", "C")))) %>% 
  ggplot(aes(x = LPAgroup, fill = LPAgroup)) +
  geom_bar(width = .7, alpha = .5) +
  xlab("LPA subgroups") +
  scale_fill_lancet() + 
  theme_classic() + easy_text_size(15)
p_lpa_bar


##############################################################
#
# HAMD difference among LPA groups
#
##############################################################

## merge the data
dat_mdd_factor_lpa <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                           ifelse(Class == 2, "B",
                                  ifelse(Class == 3, "D", "C")))) %>% 
  select(Class, LPAgroup) %>% 
  cbind(dat_mdd_factor)

## summary the demograpihc information among LPA groups 
table1::table1(~ age + gender + educations + center + HAMD + HAMA | LPAgroup,
               dat_mdd_factor_lpa)

## test the difference in HAMD score
bruceR::MANOVA(dat_mdd_factor_lpa, dv = "HAMD", between = "LPAgroup", 
               covariate = c("age","gender","educations")) %>%
  bruceR::EMMEANS(effect = "LPAgroup", p.adjust = 'fdr')

# ──────────────────────────────────────────────────────────────────────────────
# MS    MSE df1 df2       F     p     η²p [90% CI of η²p]  η²G
# ──────────────────────────────────────────────────────────────────────────────
# LPAgroup    3046.429 10.525   3 296 289.453 <.001 ***   .746 [.709, .776] .746
# ──────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
#   ────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df       t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────
# B - A  -11.305 (0.618) 296 -18.285 <.001 *** -3.485 [-3.991, -2.979]
# C - A  -23.183 (0.802) 296 -28.922 <.001 *** -7.146 [-7.802, -6.490]
# C - B  -11.878 (0.680) 296 -17.480 <.001 *** -3.661 [-4.218, -3.105]
# D - A   -8.536 (0.599) 296 -14.250 <.001 *** -2.631 [-3.122, -2.141]
# D - B    2.770 (0.437) 296   6.330 <.001 ***  0.854 [ 0.495,  1.212]
# D - C   14.647 (0.675) 296  21.708 <.001 ***  4.515 [ 3.962,  5.067]
# ────────────────────────────────────────────────────────────────────────
# Pooled SD for computing Cohen’s d: 3.244
# Results are averaged over the levels of: gender, educations
# P-value adjustment: FDR method for 6 tests.

## visualize the among-group difference in HAMD --------------------------------------
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
  scale_fill_lancet() + ylab("HAMD score") + xlab("LPA subgroups") +
  theme_classic() +
  easy_text_size(15)
p_lpa_boxplot

##############################################################
#
# combine the plots
#
##############################################################

p_combine <- (plot_compa_hama + plot_compa_hamd)/p_lpa_line/(p_hama_hamd + p_lpa_bar + p_lpa_boxplot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(width = 4300, height = 3000, 
             file = "outputs/REST_2_combine.png", dpi = 300)
print(p_combine)
dev.off()

rio::export(dat_mdd_factor_lpa, file = "inputs/REST_3_LPAgroup.xlsx")
