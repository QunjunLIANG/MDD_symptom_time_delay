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
dat_td <- rio::import("inputs/REST_3_LPAgroup.xlsx")

# summary the left subject with TD 
table1::table1(~ gender + age + educations + center | LPAgroup, data=dat_td)

###############################################################################
#
# group difference in network
#
###############################################################################

dat_td %>%
  select(ID, LPAgroup, ends_with("mean"), age, gender, educations, fd_mean) %>%
  pivot_longer(cols = 3:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "ID", dv = 'td', between = "LPAgroup", within = "network",
                 covariate = c("age","gender","fd_mean"), sph.correction="GG") %>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network",p.adjust =  'fdr')
# ─────────────────────────────────────────────────────────────────────────────────────────
# MS   MSE    df1      df2      F     p     η²p [90% CI of η²p]  η²G
# ─────────────────────────────────────────────────────────────────────────────────────────
# LPAgroup * network  0.011 0.006 13.268 1326.807  1.738  .047 *     .017 [.000, .020] .014
# ─────────────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
#   ─────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────────
# B - A ATN_mean        -0.033 (0.013) 300 -2.493  .019 *   -0.350 [-0.722,  0.023]
# C - A ATN_mean        -0.068 (0.017) 300 -3.972 <.001 *** -0.719 [-1.199, -0.238]
# C - B ATN_mean        -0.035 (0.014) 300 -2.423  .019 *   -0.369 [-0.774,  0.036]
# D - B ATN_mean         0.024 (0.009) 300  2.529  .019 *    0.250 [-0.013,  0.513]
# D - C ATN_mean         0.058 (0.014) 300  4.082 <.001 ***  0.620 [ 0.216,  1.023]
# B - A TD_mean         -0.013 (0.005) 300 -2.438  .033 *   -0.133 [-0.278,  0.012]
# C - A TD_mean         -0.023 (0.007) 300 -3.452  .004 **  -0.243 [-0.431, -0.056]
# C - B TD_mean         -0.010 (0.006) 300 -1.856  .077 .   -0.110 [-0.268,  0.048]
# D - A TD_mean         -0.009 (0.005) 300 -1.903  .077 .   -0.101 [-0.242,  0.040]
# D - C TD_mean          0.013 (0.006) 300  2.408  .033 *    0.142 [-0.015,  0.300]
# ─────────────────────────────────────────────────────────────────────────────────────

dat_td %>% filter(LPAgroup == "B" | LPAgroup == "D") %>%
  select(ID, LPAgroup, ends_with("mean"), age, gender, educations, fd_mean) %>%
  pivot_longer(cols = 3:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "ID", dv = 'td', between = "LPAgroup", within = "network",
                 covariate = c("age","gender","fd_mean"))%>%
  bruceR::EMMEANS(effect = "LPAgroup", by = "network",p.adjust =  'fdr')
# ─────────────────────────────────────────────────────────────────────────────────
# MS   MSE df1  df2     F     p     η²p [90% CI of η²p]  η²G
# ─────────────────────────────────────────────────────────────────────────────────
# LPAgroup * network  0.009 0.004   6 1398 2.107  .050 *     .009 [.000, .015] .007
# ─────────────────────────────────────────────────────────────────────────────────
# Pairwise Comparisons of "LPAgroup":
#   ─────────────────────────────────────────────────────────────────────────────────────
# Contrast     "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────────
# D - B ATN_mean         0.024 (0.010) 233  2.441  .015 *    0.257 [ 0.050,  0.463]
# D - B visual_mean     -0.021 (0.009) 233 -2.273  .024 *   -0.223 [-0.417, -0.030]
# ─────────────────────────────────────────────────────────────────────────────────────
# Pooled SD for computing Cohen’s d: 0.093
# Results are averaged over the levels of: gender
# No need to adjust p values.

###############################################################################
#
# TD & HAMD correlation
#
###############################################################################

## HAMD wave1 total score
cor_mat <- dat_td %>% 
  select(HAMD, ends_with("mean")) %>%
  select(1:8) %>%
  bruceR::Corr(p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────
# r       [95% CI]     p       N
# ─────────────────────────────────────────────────────────────
# HAMD-TD_mean                0.17 [-0.01,  0.34]  .007 **  307
# HAMD-somMot_mean            0.09 [-0.08,  0.27]  .166     307
# HAMD-visual_mean           -0.02 [-0.20,  0.16]  .790     307
# HAMD-salience_mean         -0.03 [-0.20,  0.15]  .784     307
# HAMD-ATN_mean               0.24 [ 0.06,  0.40] <.001 *** 307
# HAMD-FPN_mean               0.02 [-0.15,  0.20]  .790     307
# HAMD-DMN_mean               0.06 [-0.11,  0.24]  .388     307


## correlation between SMN ,VN and HAMD depression subfactor in group B and D
dat_hamd_cor_lpa_SMN <- data.frame()
for (lpa in c("B","D")) {
  cor_tmp <- dat_td %>% filter(LPAgroup == lpa) %>%
    cor.test(~ HAMD_depression + somMot_mean, data = .)
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

dat_hamd_cor_lpa_VN <- data.frame()
for (lpa in c("B","D")) {
  cor_tmp <- dat_td %>% filter(LPAgroup == lpa) %>%
    cor.test(~ HAMD_depression + visual_mean, data = .)
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

dat_hamd_cor_collect <- rbind(
  dat_hamd_cor_lpa_SMN, dat_hamd_cor_lpa_VN
)
dat_hamd_cor_collect
# LPAgroup network     r_value    conf_low    conf_up    p_value     p_fdr
# cor          B  somMot  0.01855429 -0.17014124 0.20593727 0.84813729 0.8481373
# cor1         D  somMot -0.03626959 -0.20782107 0.13744649 0.68322368 0.8481373
# cor2         B  visual  0.16105056 -0.02789624 0.33888611 0.09433279 0.1045477
# cor11        D  visual -0.14357039 -0.30876228 0.03002922 0.10454769 0.1045477

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
