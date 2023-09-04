##########################
#
# Subjects information description
#
# This script is used to summary 
# the demographic information about 
# subjects between MDD and HC
#
# Liang Qunjun 2023/07/19

library(tidyverse)
library(NbClust)
library(ggiraphExtra)
library(ggsci)
library(ggeasy)
library(tidyLPA)
library(psych)
library(mclust)
library(ggsignif)
library(patchwork)
library(scales)

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

# use table1 to summary the gender, age and education between groups 
table1::table1(data = dat_mdd, ~ age + educations | gender)

#############################################################
#
# calculate the factor and total score
#
#############################################################

# define the items
item_depression <- c(1,2,7,8,10,13)
item_anxiety <- c(9,10,11)
item_sleepness <- c(4,5,6)
fact_depression_total <- length(item_depression)*4
fact_anxiety_total <- length(item_anxiety)*4
fact_sleepness_total <- 3*2

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

# weighted the factor score
hamd_wave1["HAMD_wave1_depression_weighted"] <- hamd_wave1$HAMD_wave1_depression/fact_depression_total
hamd_wave1["HAMD_wave1_anxiety_weighted"] <- hamd_wave1$HAMD_wave1_anxiety/fact_anxiety_total
hamd_wave1["HAMD_wave1_sleepness_weighted"] <- hamd_wave1$HAMD_wave1_sleepness/fact_sleepness_total

hamd_wave2["HAMD_wave2_depression_weighted"] <- hamd_wave2$HAMD_wave2_depression/fact_depression_total
hamd_wave2["HAMD_wave2_anxiety_weighted"] <- hamd_wave2$HAMD_wave2_anxiety/fact_anxiety_total
hamd_wave2["HAMD_wave2_sleepness_weighted"] <- hamd_wave2$HAMD_wave2_sleepness/fact_sleepness_total

hamd_wave3["HAMD_wave3_depression_weighted"] <- hamd_wave3$HAMD_wave3_depression/fact_depression_total
hamd_wave3["HAMD_wave3_anxiety_weighted"] <- hamd_wave3$HAMD_wave3_anxiety/fact_anxiety_total
hamd_wave3["HAMD_wave3_sleepness_weighted"] <- hamd_wave3$HAMD_wave3_sleepness/fact_sleepness_total

# calculate HAMA wave1 total score
hama_wave1 <- dat_mdd %>% select(participant_id, contains('HAMA_wave1_item')) 
hama_wave1 <- hama_wave1 %>% mutate(HAMA_wave1_total = rowSums(hama_wave1[,2:ncol(hama_wave1)]) )

# merger the factor to the main data table
dat_mdd_factor <- dat_mdd %>% left_join(hamd_wave1) %>% left_join(hamd_wave2) %>%
  left_join(hamd_wave3) %>% left_join(hama_wave1)

# use table1 to summary HAMD score between testing waves divided by gender
table1::table1(data = dat_mdd_factor,
               ~ age + educations + HAMD_wave1_total + HAMD_wave1_depression + HAMD_wave2_total + HAMD_wave3_total | gender)

#############################################################
#
# testing the HAMD change after 2-week period of treatment
#
#############################################################

t.test(dat_mdd_factor$HAMD_wave1_total, dat_mdd_factor$HAMD_wave3_total,)
# data:  dat_mdd_factor$HAMD_wave1_total and dat_mdd_factor$HAMD_wave3_total
# t = 19.281, df = 235.96, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   13.81438 16.95873
# sample estimates:
#   mean of x mean of y 
# 27.29412  11.90756 


# result visualization
ComparisonDotBarChart <- function(data, ylab_title, sig_y_position,
                                  sig_TextSize=7, sig_annotation = "p < .001",
                                  bar_width=0.5, dot_size = 4, dot_alpha = 0.7,
                                  xlab_title='survey wave', 
                                  main_title=NULL, textSize=20){
  library(ggeasy)
  library(ggsignif)
  
  f <- data$stage
  f <- c(1.2,1.8)[as.integer(factor(f))]
  
  p <- ggplot(data, aes(x = stage, y = value)) +
    stat_summary(fun = mean, geom = 'bar', aes(x = stage, y = value, fill = stage),
                 width = bar_width, position = position_identity(), alpha = .55) +
    geom_line(aes(x = f, group = participant_id), size = 1, color = 'grey') +
    geom_point(aes(x = f), size = dot_size, position = position_dodge(width = 1), alpha = dot_alpha) +
    geom_signif(annotations = c(sig_annotation), textsize = sig_TextSize, vjust = .7,
                y_position = sig_y_position, xmin = 1, xmax = 2, tip_length = 0) +
    xlab(xlab_title) + ylab(ylab_title) + ggtitle(main_title) +
    theme_classic() + easy_text_size(textSize) + easy_remove_legend() + 
    easy_plot_title_size(textSize)
  
  return(p)
}

# total score
p_hamd_change <- dat_mdd_factor %>% filter(!is.na(HAMD_wave3_Date)) %>%
  select(1, HAMD_wave1_total, HAMD_wave3_total) %>%
  pivot_longer(cols = 2:3) %>% 
  mutate(stage = ifelse(str_detect(name, pattern ='wave1'), 'wave1', 'wave3')) %>%
  ComparisonDotBarChart(data = ., ylab_title = 'HAMD score', sig_annotation="***", 
                        sig_y_position = 42, dot_alpha = 0.4
                        )

Cairo::Cairo(file = 'outputs/Anay1_HAMD_change.png', width = 580, height = 680)
print(p_hamd_change)
dev.off()

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

#############################################################
#
# patients treat response - poor or well?
#
#############################################################

# summary the demographic info for treat group
dat_mdd_factor_treat <- dat_mdd_factor %>% 
  mutate(treatEff = HAMD_wave3_total/HAMD_wave1_total) %>%
  mutate(HAMD_decline = HAMD_wave1_total - HAMD_wave3_total) %>%
  mutate(treatGroup = ifelse((HAMD_wave1_total > 17 & treatEff <= .5) , 'positive', 
                             ifelse(HAMD_wave1_total < 17, "undefine", "negative")))

## and frequency
dat_mdd_factor_treat %>% filter(!is.na(HAMD_wave3_Date)) %>%
  .$treatGroup %>%
  table()
# negative positive 
# 52       67 

## obtain colors
p_treteffect_bar <- dat_mdd_factor_treat %>%
  filter(treatGroup != "undefine") %>%
  ggplot(aes(x = treatGroup)) +
    geom_bar(fill = c("#2E2A2BFF", "#CF4E9CFF"), width = .7, alpha = .6) +
    xlab("Treat response") + 
    coord_flip() +
  #ggtitle("Amount of patients in different treatment effect") +
    theme_classic() + easy_text_size(20) + easy_plot_title_size(20)
p_treteffect_bar

# use table1 to summary HAMD score between testing waves divided by gender
table1::table1(data = dat_mdd_factor_treat, 
               ~ age + gender + educations + HAMD_wave1_total | treatGroup)

##############################################################
#
# Latent profile analysis - LPA
#
##############################################################
dat_mdd_factor_treat %>% select(starts_with("HAMD_wave1_item")) %>%
  select(1:17) %>% 
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1) %>%
  compare_solutions(statistics = c('AIC','BIC',"AWE", "CLC", "KIC"))
# Model Classes AIC      BIC      AWE      CLC      KIC     
# 1     2       4896.397 5040.911 5443.861 4793.961 4951.397
# 1     3       4680.790 4875.329 5418.038 4542.620 4753.790
# 1     4       4618.053 4862.616 5545.285 4443.947 4709.053
# 
# Best model according to AIC is Model 1 with 4 classes.
# Best model according to BIC is Model 1 with 4 classes.
# Best model according to AWE is Model 1 with 3 classes.
# Best model according to CLC is Model 1 with 4 classes.
# Best model according to KIC is Model 1 with 4 classes.
# 
# An analytic hierarchy process, 
# based on the fit indices AIC, AWE, BIC, CLC, and KIC (Akogul & Erisoglu, 2017), 
# suggests the best solution is Model 1 with 4 classes.

# call mixture model
model_lpa <- dat_mdd_factor_treat %>% select(starts_with("HAMD_wave1_item")) %>%
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1)

# plot the result
p_lpa_line <- model_lpa$model_1_class_4$estimates %>%
  filter(Category == "Means") %>% 
  mutate(item = rep(paste0("item", formatC(1:17,width = 2, flag = "0")), times = 4)) %>%
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                        ifelse(Class == 2, "B",
                               ifelse(Class == 3, "C", "D")))) %>%
  ggplot(aes(x = item, y = Estimate, group = LPAgroup, color = LPAgroup)) +
  geom_line() + geom_point() +
  scale_color_lancet() +
  xlab("HAMD-17 items") +
  theme_classic() + easy_remove_x_axis(what = "title") +
  easy_text_size(20) + easy_rotate_x_labels(angle = -45)
p_lpa_line

Cairo::Cairo(file = 'outputs/Anay1_LPA_lineplot.png')
print(p_lpa_line)
dev.off()

p_lpa_bar_count <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                        ifelse(Class == 2, "B",
                               ifelse(Class == 3, "C", "D")))) %>%
  ggplot(aes(x = LPAgroup, fill = LPAgroup)) +
  geom_bar(width = .7, alpha = .5) +
  xlab("LPA subgroups") +
  scale_fill_lancet() + coord_flip() +
  theme_classic() + easy_text_size(20)
p_lpa_bar_count

Cairo::Cairo(file = 'outputs/Anay1_LPA_barplot_count.png')
print(p_lpa_bar_count)
dev.off()

##########################################################
#
# Merge the LPA group and testing HAMD and HAMA difference
#
##########################################################
## merge the data
dat_mdd_factor_treat_lpa <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "A", 
                        ifelse(Class == 2, "B",
                               ifelse(Class == 3, "C", "D")))) %>% 
  mutate(serverity = ifelse(LPAgroup == "A" | LPAgroup == "C", "moderated", "server")) %>%
  select(Class, LPAgroup, serverity) %>% 
  cbind(dat_mdd_factor_treat)
### summary the LPA group 
table1::table1(data = dat_mdd_factor_treat_lpa, 
               ~ age + gender + educations + treatGroup + HAMD_wave1_total | LPAgroup) 

## MANOVA with bruceR
bruceR::MANOVA(
  data = dat_mdd_factor_treat_lpa,
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

## boxlplot to show the HAMD wave1 difference
f <- dat_mdd_factor_treat_lpa$LPAgroup
f <- c(.7, 1.7, 2.7, 3.7)[as.integer(factor(f))]

p_lpa_hamd_diff <- ggplot(dat_mdd_factor_treat_lpa, aes(x = LPAgroup, y = HAMD_wave1_total, fill = LPAgroup)) +
  geom_boxplot(alpha = .5, width = .3) +
  geom_point(aes(x = f), 
             position = position_jitter(.02), size = 3, alpha =.35) +
  geom_signif(annotations = rep("***", times = 4), 
              textsize = 7, vjust = .7,
              y_position = c(46.5, 45, 43.5, 42), 
              xmin = c(3, 2, 1, 1), 
              xmax = c(4, 3, 4, 2),
              tip_length = 0) +
  scale_fill_lancet() + ylab("HAMD score") + xlab("LPA subgroups") +
  theme_classic() +
  easy_text_size(20)
p_lpa_hamd_diff

## HAMA difference between B and D --------------------------------
dat_mdd_factor_treat_lpa %>% filter(serverity == "server") %>%
t.test(HAMA_wave1_total ~ LPAgroup , data = .)
# data:  HAMA_wave1_total by LPAgroup
# t = -2.6468, df = 48.926, p-value = 0.0109
# alternative hypothesis: true difference in means between group B and group D is not equal to 0
# 95 percent confidence interval:
#   -7.515368 -1.028268
# sample estimates:
#   mean in group B mean in group D 
# 22.56818        26.84000 

dat_hama <- dat_mdd_factor_treat_lpa %>% filter(serverity == "server")
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
  ylab("HAMA score") + xlab("LPA subgroups") +
  theme_classic() +
  easy_text_size(20)
p_lpa_hama_diff

## proportion of treat response
p_lpa_treatGroup <- dat_mdd_factor_treat_lpa %>% filter(treatGroup != "undefine" & !is.na(treatGroup)) %>%
  ggplot(aes(x = LPAgroup, fill = treatGroup)) +
  geom_bar(position = position_fill(), alpha = .6) +
  scale_fill_cosmic() + xlab("LPA subgroups") + ylab("Percentage") +
  coord_flip() +
  theme_classic()
p_lpa_treatGroup

##############################################################
#
# Combine all plots
#
##############################################################

# define layout
layout_use <- "
AABB
AACC
DDDD
EEFF
EEGG
"

# combine the plots
p_combine <- p_hamd_change + p_cor_hamd_hama + p_treteffect_bar + p_lpa_line + p_lpa_hamd_diff + p_lpa_bar_count + p_lpa_hama_diff  +
  plot_layout(design = layout_use, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))
p_combine

Cairo::Cairo(file = 'outputs/Anay1_plot_combine.png', height = 1000, width = 800)
print(p_combine)
dev.off()

##############################################################
#
# Export the new subject information
#
##############################################################
rio::export(dat_mdd_factor_treat_lpa, file = "inputs/Analysis1_subject_table.xlsx")
