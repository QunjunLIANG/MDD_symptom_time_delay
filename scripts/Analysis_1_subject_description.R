##########################
#
# Subjects information description
#
# This script is used to summary 
# the demographic information about 
# subjects between MDD and HC
#
# Liang Qunjun 2023/09/15

library(tidyverse)
library(ggiraphExtra)
library(ggsci)
library(ggeasy)
library(tidyLPA)
library(psych)
library(mclust)
library(ggsignif)
library(patchwork)
library(scales)
library(forcats) 
source("scripts/function_PvalueForTable1.R")

# load the data
dat_sbj <- rio::import('inputs/subject_information.xlsx')

# fix the education 
dat_sbj <- dat_sbj %>% mutate(educations = ifelse(education == 1, 'Illiterate',
                                                  ifelse(education == 2, 'Primary education',
                                                         ifelse(education == 3, 'Junior high school',
                                                                ifelse(education == 4, 'Senior high school',
                                                                       ifelse(education == 5, 'Undergraduate', 'Graduate'))))))

# summary the clinical scales
dat_mdd <- dat_sbj %>% filter(!is.na(HAMD_wave3_Date))

#############################################################
#
# calculate the factor and total score
#
#############################################################

# define the items
item_depression <- c(1,2,7,8,13)
item_anxiety <- c(9,10,11)
item_sleepness <- c(4,5,6)

# calculate the factor and total score
hamd_wave1 <- dat_mdd %>% select(participant_id, contains('HAMD_wave1_item')) 
hamd_wave1 <- hamd_wave1 %>%
  mutate(HAMD_wave1_total = rowSums(hamd_wave1[,2:ncol(hamd_wave1)])) %>%
  mutate(HAMD_wave1_depression = rowSums(hamd_wave1[, 1+item_depression])) %>%
  mutate(HAMD_wave1_anxiety = rowSums(hamd_wave1[, 1+item_anxiety])) %>%
  mutate(HAMD_wave1_sleepness = rowSums(hamd_wave1[, 1+item_sleepness]))

hamd_wave2 <- dat_mdd %>% select(participant_id, contains('HAMD_wave2_item')) 
hamd_wave2 <- hamd_wave2 %>%
  mutate(HAMD_wave2_total = rowSums(hamd_wave2[,2:ncol(hamd_wave2)])) %>%
  mutate(HAMD_wave2_depression = rowSums(hamd_wave2[, 1+item_depression])) %>%
  mutate(HAMD_wave2_anxiety = rowSums(hamd_wave2[, 1+item_anxiety])) %>%
  mutate(HAMD_wave2_sleepness = rowSums(hamd_wave2[, 1+item_sleepness]))

hamd_wave3 <- dat_mdd %>% select(participant_id, contains('HAMD_wave3_item')) 
hamd_wave3 <- hamd_wave3 %>%
  mutate(HAMD_wave3_total = rowSums(hamd_wave3[,2:ncol(hamd_wave3)])) %>%
  mutate(HAMD_wave3_depression = rowSums(hamd_wave3[, 1+item_depression])) %>%
  mutate(HAMD_wave3_anxiety = rowSums(hamd_wave3[, 1+item_anxiety])) %>%
  mutate(HAMD_wave3_sleepness = rowSums(hamd_wave3[, 1+item_sleepness]))

# calculate HAMA wave1 total score
hama_wave1 <- dat_mdd %>% select(participant_id, contains('HAMA_wave1_item')) 
hama_wave1 <- hama_wave1 %>% mutate(HAMA_wave1_total = rowSums(hama_wave1[,2:ncol(hama_wave1)]) )

# merger the factor to the main data table
dat_mdd_factor <- dat_mdd %>% left_join(hamd_wave1) %>% left_join(hamd_wave2) %>%
  left_join(hamd_wave3) %>% left_join(hama_wave1)

##############################################################
#
# Latent profile analysis - LPA
#
##############################################################
set.seed(123)

dat_mdd_factor %>% select(starts_with("HAMD_wave1_item")) %>%
  select(1:17) %>% 
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1) %>%
  compare_solutions(statistics = c('AIC','BIC',"AWE", "CLC", "KIC"))
# Compare tidyLPA solutions:
#   
#   Model Classes AIC      BIC      AWE      CLC      KIC     
# 1     2       5148.582 5296.888 5703.214 5046.561 5203.582
# 1     3       5059.561 5259.203 5807.014 4921.392 5132.561
# 1     4       4987.296 5238.275 5927.348 4813.202 5078.296
# 
# Best model according to AIC is Model 1 with 4 classes.
# Best model according to BIC is Model 1 with 4 classes.
# Best model according to AWE is Model 1 with 2 classes.
# Best model according to CLC is Model 1 with 4 classes.
# Best model according to KIC is Model 1 with 4 classes.
# 
# An analytic hierarchy process, 
# based on the fit indices AIC, AWE, BIC, CLC, and KIC (Akogul & Erisoglu, 2017), 
# suggests the best solution is Model 1 with 4 classes.

# call mixture model
model_lpa <- dat_mdd_factor %>% select(starts_with("HAMD_wave1_item")) %>%
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1)

## merge the data
dat_mdd_factor_lpa <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                           ifelse(Class == 2, "B",
                                  ifelse(Class == 3, "C", "D")))) %>% 
  mutate(LPAgroup = factor(LPAgroup, levels = c("A","B","C","D"),
                              labels = c("I-MDD","NA-MDD","NI-MDD","A-MDD"))) %>%
  select(Class, LPAgroup) %>% 
  cbind(dat_mdd_factor)

##############################################################
#
# Difference of HAMD among LPA groups
#
##############################################################

## MANOVA with bruceR
bruceR::MANOVA(
  data = dat_mdd_factor_lpa,
  subID = "participant_id",
  between = "LPAgroup",
  dv = "HAMD_wave1_total",
  covariate = c("age","gender")
) %>%
  bruceR::EMMEANS(effect = "LPAgroup")
# ───────────────────────────────────────────────────────────────────────────
# MS    MSE df1 df2      F     p     η²p [90% CI of η²p]  η²G
# ───────────────────────────────────────────────────────────────────────────
# LPAgroup  1160.551 18.557   3 122 62.540 <.001 ***   .606 [.516, .672] .606
# age         24.506 18.557   1 122  1.321  .253       .011 [.000, .060] .011
# gender      14.368 18.557   1 122  0.774  .381       .006 [.000, .050] .006
# ───────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
#   ───────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df       t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────────────
# (I-MDD) - (A-MDD)    -12.092 (1.244) 122  -9.721 <.001 *** -2.807 [-3.581, -2.033]
# (NA-MDD) - (A-MDD)    -2.084 (1.094) 122  -1.905  .355     -0.484 [-1.165,  0.197]
# (NA-MDD) - (I-MDD)    10.007 (1.116) 122   8.965 <.001 ***  2.323 [ 1.628,  3.018]
# (NI-MDD) - (A-MDD)   -11.550 (1.140) 122 -10.132 <.001 *** -2.681 [-3.391, -1.971]
# (NI-MDD) - (I-MDD)     0.542 (1.159) 122   0.468 1.000      0.126 [-0.596,  0.847]
# (NI-MDD) - (NA-MDD)   -9.465 (0.980) 122  -9.655 <.001 *** -2.197 [-2.808, -1.587]
# ───────────────────────────────────────────────────────────────────────────────────

## retardation difference between A-MDD and NA-MDD
dat_mdd_factor_lpa %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(participant_id, LPAgroup, HAMA_wave1_item8) %>%
  mutate(LPAgroup_recode = ifelse(LPAgroup == "A-MDD", "A-MDD", "NA-MDD")) %>%
  bruceR::TTEST(y = "HAMA_wave1_item8", x = "LPAgroup_recode", digits = 3)
# Results of t-test: HAMA_wave1_item13: LPAgroup_recode (NA-MDD - A-MDD)  
# ──────────────────────────────────────────────────────────────────────────────
#      t df     p         Difference [95% CI]      Cohen’s d [95% CI]    BF10
# ──────────────────────────────────────────────────────────────────────────────
# -2.629 68  .011 *   -0.752 [-1.322, -0.181] -0.650 [-1.144, -0.157] 4.431e+00
# ──────────────────────────────────────────────────────────────────────────────

## somatic difference between A-MDD and NA-MDD
dat_mdd_factor_lpa %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(participant_id, LPAgroup, HAMA_wave1_item13) %>%
  mutate(LPAgroup_recode = ifelse(LPAgroup == "A-MDD", "A-MDD", "NA-MDD")) %>%
  bruceR::TTEST(y = "HAMA_wave1_item13", x = "LPAgroup_recode", digits = 3)
# Results of t-test: HAMA_wave1_item13: LPAgroup_recode (NA-MDD - A-MDD)  
# ──────────────────────────────────────────────────────────────────────────────
#       t df     p         Difference [95% CI]      Cohen’s d [95% CI]   BF10
# ──────────────────────────────────────────────────────────────────────────────
# -4.161 68 <.001 *** -0.703 [-1.040, -0.366] -1.029 [-1.523, -0.536] 2.399e+02
# ──────────────────────────────────────────────────────────────────────────────

## HAMA difference between A-MDD and NA-MDD
dat_mdd_factor_lpa %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD") %>%
  select(participant_id, LPAgroup, HAMA_wave1_total) %>%
  mutate(LPAgroup_recode = ifelse(LPAgroup == "A-MDD", "A-MDD", "NA-MDD")) %>%
  bruceR::TTEST(y = "HAMA_wave1_total", x = "LPAgroup_recode", digits = 3)
# Results of t-test:
# ──────────────────────────────────────────────────────────────────────────────
#       t df     p         Difference [95% CI]      Cohen’s d [95% CI]   BF10
# ──────────────────────────────────────────────────────────────────────────────
# -2.372 68  .021 *   -3.816 [-7.027, -0.606] -0.587 [-1.080, -0.093] 2.633e+00
# ──────────────────────────────────────────────────────────────────────────────

## HAMA difference between I-MDD and NI-MDD
dat_mdd_factor_lpa %>% filter(LPAgroup == "I-MDD" | LPAgroup == "NI-MDD") %>%
  select(participant_id, LPAgroup, HAMA_wave1_total) %>%
  mutate(LPAgroup_recode = ifelse(LPAgroup == "I-MDD", "I-MDD", "NI-MDD")) %>%
  bruceR::TTEST(y = "HAMA_wave1_total", x = "LPAgroup_recode", digits = 3)
# Results of t-test:
# ─────────────────────────────────────────────────────────────────────────────
#      t df     p       Difference [95% CI]    Cohen’s d [95% CI]      BF10
# ─────────────────────────────────────────────────────────────────────────────
# 0.102 56  .919  0.171 [-3.204, 3.547]   0.027 [-0.510, 0.565]   2.721e-01
# ─────────────────────────────────────────────────────────────────────────────

##############################################################
#
# Plot the results in LPA part
#
##############################################################

# lienplot for the estimation of item-level LPA
p_lpa_line <- model_lpa$model_1_class_4$estimates %>%
  filter(Category == "Means") %>% 
  mutate(item = rep(paste0("item", formatC(1:17,width = 2, flag = "0")), times = 4)) %>%
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                           ifelse(Class == 2, "B",
                                  ifelse(Class == 3, "C", "D")))) %>% 
  mutate(LPAgroup = factor(LPAgroup, levels = c("A","B","C","D"),
                           labels = c("I-MDD","NA-MDD","NI-MDD","A-MDD"))) %>%
  ggplot(aes(x = item, y = Estimate, group = LPAgroup, color = LPAgroup)) +
  geom_line() + geom_point() +
  scale_color_lancet() +
  xlab("HAMD-17 items") +
  theme_classic() + easy_remove_x_axis(what = "title") +
  easy_text_size(20) + easy_rotate_x_labels(angle = -45) +
  easy_add_legend_title("Subtypes")
p_lpa_line

## number of patient assign to each subtypes
p_lpa_bar_count <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                           ifelse(Class == 2, "B",
                                  ifelse(Class == 3, "C", "D")))) %>% 
  mutate(LPAgroup = factor(LPAgroup, levels = c("A","B","C","D"),
                           labels = c("I-MDD","NA-MDD","NI-MDD","A-MDD"))) %>%
  ggplot(aes(x = LPAgroup, fill = LPAgroup)) +
  geom_bar(width = .7, alpha = .5) +
  xlab("LPA subgroups") +
  scale_fill_lancet() + coord_flip() +
  theme_classic() + easy_text_size(20)
p_lpa_bar_count

## boxlplot to show the HAMD wave1 difference
f <- dat_mdd_factor_lpa$LPAgroup
f <- c(.7, 1.7, 2.7, 3.7)[as.integer(factor(f))]

p_lpa_hamd_diff <- ggplot(dat_mdd_factor_lpa, aes(x = LPAgroup, y = HAMD_wave1_total, fill = LPAgroup)) +
  geom_boxplot(alpha = .5, width = .4) +
  geom_point(aes(x = f), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = rep("***", times = 4), 
              textsize = 7, vjust = .7,
              y_position = c(46.5, 45, 43.5, 42), 
              xmin = c(3, 2, 1, 1), 
              xmax = c(4, 3, 4, 2),
              tip_length = 0) +
  scale_fill_lancet() + ylab("HAMD score") + xlab("Subtypes") +
  theme_classic() +
  easy_text_size(20) + easy_remove_legend()
p_lpa_hamd_diff

## HAMA difference between A-MDD and NA-MDD
dat_hama <- dat_mdd_factor_lpa %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD")
f_1 <- dat_hama$LPAgroup
f_1 <- c(.7, 1.7)[as.integer(factor(f_1))]

p_lpa_hama_diff <- ggplot(dat_hama, aes(x = LPAgroup, y = HAMD_wave1_total)) +
  geom_boxplot(alpha = .5, width = .3, fill = c("#ED0000FF","#0099B4FF")) +
  geom_point(aes(x = f_1), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = c("**"), 
              textsize = 7, vjust = .7,
              y_position = c(42), 
              xmin = c(1), 
              xmax = c(2),
              tip_length = 0) +
  coord_flip() +
  ylab("HAMA score") + xlab("Subtypes") +
  theme_classic() +
  easy_text_size(20)
p_lpa_hama_diff

#############################################################
#
# patients treat response - poor or well?
#
#############################################################

# summary the demographic info for treat group
dat_mdd_factor_lpa_treat <- dat_mdd_factor_lpa %>% 
  mutate(treatEff = HAMD_wave3_total/HAMD_wave1_total) %>%
  mutate(HAMD_decline = HAMD_wave1_total - HAMD_wave3_total) %>%
  mutate(treatGroup = ifelse((HAMD_wave1_total > 17 & treatEff <= .5) , 'positive', 
                             ifelse(HAMD_wave1_total < 17, "undefine", "negative")))

## and frequency
dat_mdd_factor_lpa_treat %>% filter(!is.na(HAMD_wave3_Date)) %>%
  .$treatGroup %>%
  table()
# negative positive undefine 
# 52       67        9 

#############################################################
#
# Plot the results for treatment response
#
#############################################################

## obtain colors
p_treteffect_bar <- dat_mdd_factor_lpa_treat %>%
  filter(treatGroup != "undefine") %>%
  ggplot(aes(x = treatGroup)) +
  geom_bar(fill = c("#2E2A2BFF", "#CF4E9CFF"), width = .7, alpha = .6) +
  xlab("Treat response") + 
  coord_flip() +
  #ggtitle("Amount of patients in different treatment effect") +
  theme_classic() + easy_text_size(20) + easy_plot_title_size(20)
p_treteffect_bar

### summary the LPA group 
table1::table1(data = dat_mdd_factor_lpa_treat, 
               ~ age + gender + HAMD_wave1_total + HAMD_wave3_total + 
                 treatGroup + HAMA_wave1_total | LPAgroup)

## proportion of treat response
p_lpa_treatGroup <- dat_mdd_factor_lpa_treat %>% 
  filter(treatGroup != "undefine" & !is.na(treatGroup)) %>%
  ggplot(aes(x = LPAgroup, fill = treatGroup)) +
  geom_bar(position = position_fill(), alpha = .6) +
  scale_fill_cosmic() + xlab("LPA subgroups") + ylab("Percentage") +
  coord_flip() + 
  theme_classic() +easy_add_legend_title("Response")
p_lpa_treatGroup

### summary the LPA group 
table1::table1(data = dat_mdd_factor_lpa_treat %>% filter(LPAgroup == "A-MDD" | LPAgroup == "NA-MDD"), 
               ~ age + gender + HAMD_wave1_total + HAMD_wave3_total + 
                 treatGroup + HAMA_wave1_total | LPAgroup,
               overall = F, extra.col=list("P-value"=pvalue)) 

#############################################################
#
# correlation - HAMD - HAMA - wave1
#
#############################################################

p_cor_hamd_hama <- ggplot(data = dat_mdd_factor, aes(x = HAMD_wave1_total, y = HAMA_wave1_total)) +
  geom_point(size = 3, alpha = .6) +
  geom_smooth(method = "lm") +
  xlab("HAMD score") + ylab("HAMA score") +
  theme_classic()
p_cor_hamd_hama
cor.test(dat_mdd_factor$HAMD_wave1_total, dat_mdd_factor$HAMA_wave1_total)
# data:  dat_mdd_factor$HAMD_wave1_total and dat_mdd_factor$HAMA_wave1_total
# t = 7.8057, df = 117, p-value = 2.769e-12
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.4528868 0.6922599
# sample estimates:
#   cor 
# 0.5851799 

##############################################################
#
# Combine all plots
#
##############################################################

# define layout
layout_use <- "
AAAA
BBCC
BBDD
"
# combine the plots
p_combine <- p_lpa_line + p_lpa_hamd_diff + p_lpa_hama_diff + p_lpa_treatGroup  +
  plot_layout(design = layout_use, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(file = 'outputs/Fig1.png', height = 2000, width = 2800, dpi = 300)
print(p_combine)
dev.off()

##############################################################
#
# Export the new subject information
#
##############################################################
rio::export(dat_mdd_factor_lpa_treat, file = "inputs/Analysis1_subject_table.xlsx")
