#########################################
#
# This script is used to analyze the 
# relationship between baseline TDp
# and the treat effect.
#
#
# Liang Qunjun   2023/08/10

library(tidyverse)
library(ggeasy)
library(ggsci)
library(ggsignif)
library(patchwork)
library(caret)
library(ggsignif)

# load the data
dat_td   <- rio::import(file = "inputs/Analysis2_TDp_collection.xlsx") 

# combine HAMD and TD data
dat_use <- dat_td %>% filter(treatGroup != "undefine")

# the distribution of treatment response in LPA groups
p_lpa_treatGroup <- dat_use %>% filter(treatGroup != "undefine" & !is.na(treatGroup)) %>%
  ggplot(aes(x = LPAgroup, fill = treatGroup)) +
  geom_bar(position = position_fill(), alpha = .6) +
  scale_fill_manual(values = c("#2E2A2BFF", "#CF4E9CFF")) + xlab("LPA subgroups") + ylab("Percentage") +
  coord_flip() + 
  theme_classic() + easy_add_legend_title("Response") 
p_lpa_treatGroup

###############################################################################
#
# Profile of the treat-different groups
#
###############################################################################

# the baseline HAMD total score in different treat groups
bruceR::MANOVA(data = dat_use, subID = "participant_id", 
               dv = "HAMD_wave1_total", between = "treatGroup",
               covariate = c("age", "gender"))
# ───────────────────────────────
# "treatGroup"   Mean    S.D.  n
# ───────────────────────────────
# negative 25.529 (6.307) 51
# positive 28.672 (5.685) 67
# ───────────────────────────────
# Covariate(s):               age, gender, educations
# ───────────────────────────────────────────────────────────────────────────
# MS    MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ───────────────────────────────────────────────────────────────────────────
# treatGroup  234.760 35.154   1 115 6.678  .011 *     .055 [.007, .136] .055
# age          75.430 35.154   1 115 2.146  .146       .018 [.000, .078] .018
# gender        3.509 35.154   1 115 0.100  .753       .001 [.000, .028] .001
# ───────────────────────────────────────────────────────────────────────────

# the baseline TD in different treat groups
bruceR::MANOVA(data = dat_use, subID = "participant_id", 
               dv = "TD_mean", between = "treatGroup",
               covariate = c("age", "gender", "QC_bold_fd_mean"))
# ───────────────────────────────
# "treatGroup"   Mean    S.D.  n
# ───────────────────────────────
# negative -0.012 (0.014) 51
# positive -0.011 (0.014) 67
# ───────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ─────────────────────────────────────────────────────────────────────────────
# treatGroup       0.000 0.000   1 114 0.754  .387       .007 [.000, .053] .007
# ─────────────────────────────────────────────────────────────────────────────

# Chi-square test for treat group
ggstatsplot::ggbarstats(dat_use, x = treatGroup, y = LPAgroup)

table(dat_use %>% filter(LPAgroup == "D") %>% .$treatGroup) %>%
  chisq.test()
# Chi-squared test for given probabilities
# 
# data:  .
# X-squared = 9.8462, df = 1, p-value = 0.001702

f_log1 <- dat_use %>% .$treatGroup
f_log1 <- c(.75, 1.75)[as.integer(factor(f_log1))]

p_baseHAMD <- ggplot(dat_use, aes(x = treatGroup, y = HAMD_wave1_total)) +
  geom_boxplot(fill = c("#2E2A2BFF", "#CF4E9CFF"), alpha = .6, width = .3) +
  geom_point(aes(x = f_log1), position = position_jitter(.02), size = 3, alpha = .4) +
  geom_signif(annotations = c("**"), 
              textsize = 7, 
              vjust = .7,
              y_position = c(42), 
              xmin = c(1), 
              xmax = c(2),
              tip_length = 0) +
  ylab("baseline HAMD score") + xlab("Treatment response") +
  coord_flip() +
  theme_classic() + easy_text_size(15)
p_baseHAMD

Cairo::Cairo(file = "outputs/Anay4_baselineHAMD.png")
print(p_baseHAMD)
dev.off()

###############################################################################
#
# Logistic regression -- LPA subgroups
#
###############################################################################

dat_lm_logi <- dat_use %>%
  mutate(treatEff_recode = ifelse(treatGroup == "positive", 1, 0))
dat_lm_logi$treatEff_recode <- factor(dat_lm_logi$treatEff_recode, levels = c(0,1),
                                      labels = c("negative","positive"))

## model for all LPA group ---------------------------------------------------------
model_logi <- glm(treatEff_recode ~  visual_mean + somMot_mean + ATN_mean +
                       salience_mean + FPN_mean + DMN_mean +
                       age + gender + QC_bold_fd_mean + HAMD_wave1_total, 
                     data = dat_lm_logi,
                     family = binomial())
bruceR::GLM_summary(model_logi)
# ─────────────────────────────────────────────────────────────────────────────────
# b    S.E.      z     p         [95% CI of b]      OR   VIF
# ─────────────────────────────────────────────────────────────────────────────────
# salience_mean      4.618 (2.710)  1.704  .088 .   [ -0.694,  9.929] 101.251 1.078
# FPN_mean          -8.368 (3.434) -2.437  .015 *   [-15.099, -1.638]   0.000 1.169

## model for B and D LPA group - somatomotor ---------------------------------------------------------
dat_lm_logi %>% filter(LPAgroup == "B" | LPAgroup == "D") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "somMot_mean", 
                 between = c("LPAgroup","treatEff_recode"),
                 covariate = c("age","gender","QC_bold_fd_mean")) %>%
  bruceR::EMMEANS("treatEff_recode", by = "LPAgroup")
# ────────────────────────────────────────────────────────────────────────────────────────
# MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ────────────────────────────────────────────────────────────────────────────────────────
# LPAgroup                    0.020 0.002   1  63 9.446  .003 **    .130 [.028, .265] .130
# treatEff_recode             0.000 0.002   1  63 0.023  .879       .000 [.000, .028] .000
# age                         0.001 0.002   1  63 0.618  .435       .010 [.000, .085] .010
# gender                      0.001 0.002   1  63 0.385  .537       .006 [.000, .074] .006
# QC_bold_fd_mean             0.016 0.002   1  63 7.730  .007 **    .109 [.018, .240] .109
# LPAgroup * treatEff_recode  0.014 0.002   1  63 6.858  .011 *     .098 [.013, .227] .098
# ────────────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "treatEff_recode":
# ────────────────────────────────────────────────────────────────────────────────────────────
# Contrast "LPAgroup" Estimate    S.E. df      t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────────────────
# positive - negative          B   -0.033 (0.015) 63 -2.255  .028 *   -0.728 [-1.373, -0.083]
# positive - negative          D    0.037 (0.023) 63  1.642  .106      0.820 [-0.178,  1.817]
# ────────────────────────────────────────────────────────────────────────────────────────────

## model for B and D LPA group - visual ---------------------------------------------------------
dat_lm_logi %>% filter(LPAgroup == "B" | LPAgroup == "D") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "visual_mean", 
                 between = c("LPAgroup","treatEff_recode"),
                 covariate = c("age","gender","QC_bold_fd_mean"))

# visualize the logistic regression results
dat_lm_logi$LPAgroup <- factor(dat_lm_logi$LPAgroup, levels = c("A","C","B","D"))
label_use <- c(FPN_mean = "Fronto-parietal Net.",
               salience_mean = "Salience Net.",
               somMot_mean = "Somatosensory Net.")

p_logstic <- dat_lm_logi %>% 
  select(participant_id, LPAgroup, treatEff_recode, 
         somMot_mean, salience_mean, FPN_mean) %>%
  pivot_longer(cols = 4:6, names_to = "network", values_to = "td") %>%
  ggplot(aes(x = LPAgroup, y = td, fill = treatEff_recode)) +
  geom_boxplot(alpha = .5) + 
  facet_grid(.~network, labeller = labeller(network = label_use)) +
  xlab('LPA groups') + ylab("Mean time delay") +
  scale_fill_manual(values = c("#2E2A2BFF", "#CF4E9CFF")) + ylim(-.15, .15) +
  theme_bruce() + 
  bruceR::easy_add_legend_title("Response") +
  easy_text_size(15)
p_logstic


p_box_fpn <- dat_lm_logi %>% 
  select(participant_id, LPAgroup, treatEff_recode, FPN_mean) %>%
  ggplot(aes(x = LPAgroup, y = FPN_mean, fill = treatEff_recode)) +
  geom_boxplot(alpha = .5) + 
  xlab('LPA groups') + ylab("Time delay") + ggtitle("FPN") +
  scale_fill_manual(values = c("#2E2A2BFF", "#CF4E9CFF")) + ylim(-.15, .15) +
  theme_classic() + 
  easy_add_legend_title("Response") +
  easy_text_size(15)
p_box_fpn

p_box_somMot <- dat_lm_logi %>% filter(LPAgroup == "B" | LPAgroup == "D") %>%
  select(participant_id, LPAgroup, treatEff_recode, somMot_mean) %>%
  ggplot(aes(x = LPAgroup, y = somMot_mean, fill = treatEff_recode)) +
  geom_boxplot(alpha = .5) + 
  geom_signif(annotations = c("*"), 
              textsize = 7, 
              vjust = .3,
              y_position = .12, 
              xmin = 0.8, xmax = 1.2,
              tip_length = 0) +
  xlab('LPA groups') + ylab("Time delay") + ggtitle("somMot") +
  scale_fill_manual(values = c("#2E2A2BFF", "#CF4E9CFF")) + ylim(-.15, .15) +
  coord_flip() +
  theme_classic() + 
  easy_add_legend_title("Response") + easy_remove_legend() +
  easy_text_size(15)  
p_box_somMot
###############################################################################
#
# combine plots
#
###############################################################################

layout_use <- "
AACC
AACC
BBBB
BBBB
BBBB
"

p_combine <- p_lpa_treatGroup +  p_box_fpn + p_box_somMot +
  plot_layout(design = layout_use, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 15),
        axis.title.x = element_text(size =10),
        axis.text =  element_text(size =13),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_text(size =13),
        plot.tag = element_text(size = 18, face = "bold"))

p_combine

Cairo::Cairo(file = "outputs/Anay3_combine_plot.png", width = 2300, height = 1800, dpi = 300)
print(p_combine)
dev.off()
