########################################
#
# This script is used to select participant
# from REST.
#
# Liang Qunjun 2023/08/30

library(tidyverse)
library(patchwork)
library(ggeasy)

######################################
#
# Load the meta information
#
######################################

# load the demographic information
dat_mdd_demog <- rio::import("inputs/REST-meta-MDD-PhenotypicData_WithHAMDSubItem_V4.xlsx",
                             sheet = 1)
# load the image QC information
dat_mdd_qc <- rio::import("inputs/REST_demog_imag_fd.xlsx")


#####################################
#
# Filtering the participants with 
# principals

# 1. remove age < 18
dat_mdd_age <- dat_mdd_demog %>% filter(Age >= 18)

# 2. remove missing item-score
dat_mdd_age_hamd <- dat_mdd_age[which(complete.cases(dat_mdd_age[,10:26])),]

dat_mdd_age_hamd_episode <- dat_mdd_age_hamd 

# add the qc to the filtered data
dat_mdd_age_hamd_episode_combine <- dat_mdd_age_hamd_episode %>% 
  left_join(dat_mdd_qc %>% select(ID, center, `QC Score`, fd_mean, timepoint, fd_precent))

# 4. remove QC score < 3
dat_mdd_age_hamd_episode_qc <- dat_mdd_age_hamd_episode_combine %>% 
  filter(`QC Score` > 2)

# 5. remove fd mean > 0.25
dat_mdd_age_hamd_episode_qc_fd <- dat_mdd_age_hamd_episode_qc %>%
  filter(fd_mean <= .25)

# 6. FD percent <= 20%
dat_mdd_age_hamd_episode_qc_fd_percent <- dat_mdd_age_hamd_episode_qc_fd %>%
  filter(fd_precent<=.2)

############################################
#
# summary the filtered subjects
#
###########################################

dat_filtered <- dat_mdd_age_hamd_episode_qc_fd_percent %>% 
  mutate(educations = ifelse(`Education (years)` == 0, "Illiterate",
                            ifelse(`Education (years)` > 0 & `Education (years)` <= 6, "Primary education",
                                   ifelse(`Education (years)` > 6 & `Education (years)` <= 9, "Junior high school",
                                          ifelse(`Education (years)` > 9 & `Education (years)` <= 12, "Senior high school",
                                                 ifelse(`Education (years)` > 12 & `Education (years)` <= 16, "Undergraduate", "Graduate")))))) %>%
  mutate(gender = ifelse(Sex == 1, "male", "female")) %>%
  rename(age = Age)
dat_filtered$HAMD <- as.numeric(dat_filtered$HAMD)
dat_filtered$HAMA <- as.numeric(dat_filtered$HAMA)
colnames(dat_filtered)[10:33] <- paste0("HAMD_item",formatC(1:24, digits = 1, flag = "0"))

table1::table1(data = dat_filtered, ~ age + educations + HAMD + HAMA + center  | gender)

############################################
#
# Export to the disk

rio::export(dat_filtered, file = "inputs/REST_1_subject_selection.xlsx")
