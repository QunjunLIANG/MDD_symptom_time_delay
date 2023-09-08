#########################################
#
# This script is used to collect the 
# TD of available subjects and merge 
# their demograpihc information together.
#
# The merged data will be exported. 
#
# Liang Qunjun 2023/08/31

library(tidyverse)

# load the demographic information
dat_info <- rio::import("inputs/REST_1_subject_selection.xlsx")

# load the TD data -----------------------------------------------------
path_td <- "REST/time_lag_estimation/"
file_list <- list.files(path_td, pattern = "^S.*_projection_map_weighted.csv")
dat_td <- data.frame()
for (i in file_list) {
  sbj_tmp <- str_extract(i, pattern = "^S.*[0-9]{4}")
  td_tmp <- rio::import(file.path(path_td,i))
  dat_tmp <- data.frame(ID = sbj_tmp, td_tmp)
  dat_td <- rbind(dat_td, dat_tmp)
}
colnames(dat_td)[2:215] <- paste0("R", formatC(1:214, digits = 2, flag = "0"))

# load the network identification
net_id <- rio::import("inputs/Power264_Yeo7.xlsx")

dat_td_net <- dat_td %>%
  mutate(TD_mean = rowMeans(across(all_of(2:215)), na.rm = T)) %>%
  mutate(somMot_mean = rowMeans(across(all_of(1+which(net_id$network_new == "somMot"))), na.rm = T)) %>%
  mutate(visual_mean = rowMeans(across(all_of(1+which(net_id$network_new == "visual"))), na.rm = T)) %>%
  mutate(salience_mean = rowMeans(across(all_of(1+which(net_id$network_new == "salience"))), na.rm = T)) %>%
  mutate(ATN_mean = rowMeans(across(all_of(1+which(net_id$network_new == "ATN"))), na.rm = T)) %>%
  mutate(FPN_mean = rowMeans(across(all_of(1+which(net_id$network_new == "FPN"))), na.rm = T)) %>%
  mutate(DMN_mean = rowMeans(across(all_of(1+which(net_id$network_new == "DMN"))), na.rm = T)) %>%
  rowwise() %>%
  mutate(visual_sd = sd(across(all_of(1+which(net_id$network_new == "visual"))), na.rm = T)) %>%
  mutate(somMot_sd = sd(across(all_of(1+which(net_id$network_new == "somMot"))), na.rm = T)) %>%
  mutate(salience_sd = sd(across(all_of(1+which(net_id$network_new == "salience"))), na.rm = T)) %>%
  mutate(ATN_sd = sd(across(all_of(1+which(net_id$network_new == "ATN"))), na.rm = T)) %>%
  mutate(FPN_sd = sd(across(all_of(1+which(net_id$network_new == "FPN"))), na.rm = T)) %>%
  mutate(DMN_sd = sd(across(all_of(1+which(net_id$network_new == "DMN"))), na.rm = T)) %>%
  mutate(TD_sd =  sd(across(all_of(2:215)), na.rm = T))

# merge the demographic and TD ---------------------------------------------
VIM::matrixplot(dat_td_net)
dat_all <- dat_td_net %>% left_join(dat_info)

rio::export(dat_all, file = "inputs/REST_2_timeDelay.xlsx")

##########################################################################
#
# Summary the main and validation dataset
#
##########################################################################

dat_validation <- rio::import("inputs/REST_2_timeDelay.xlsx") %>%
  select(ID, age, gender, educations, HAMD, HAMA) %>%
  rename(participant_id = "ID") %>%
  mutate(dataset = "Validation")
dat_main <- rio::import("inputs/Analysis1_subject_table.xlsx") %>%
  select(participant_id, age, gender, educations, HAMD_wave1_total, HAMA_wave1_total) %>%
  rename(HAMD = "HAMD_wave1_total", HAMA = "HAMA_wave1_total")%>%
  mutate(dataset = "Main")
dat_collect <- rbind(dat_validation, dat_main)

table1::table1(~ age + gender + educations + HAMD + HAMA | dataset, data = dat_collect)
