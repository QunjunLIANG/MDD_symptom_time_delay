################################
#
# This script is used to obtain
# TD from REST
#

library(tidyverse)
library(R.matlab)

# load the participant information 
dat_sbj <- rio::import('inputs/REST_1_subject_selection.xlsx')

##########################################################
#
# network identification
#
##########################################################

power_doc <- rio::import("inputs/Power2011ROI_Module-2hv15xc.xls")

# filter out undefine and subcortical network
power_doc_filter <- power_doc %>% filter(`Master Assignment`!=-1 & 
                                         `Master Assignment`!=6 & 
                                         `Master Assignment`!=13 &
                                         `Master Assignment`!=10)

use_region <- as.numeric(power_doc_filter$...1)

## tidy the network annotation
power_doc_filter_mod <- power_doc_filter %>% 
  select(7:9, 32, 35, 37)
colnames(power_doc_filter_mod) <- c("x","y","z",
                                    "clust_ori","color","network_ori")
power_doc_filter_mod$x <- as.numeric(power_doc_filter_mod$x) %>%
  round(digits = 0)
power_doc_filter_mod$y <- as.numeric(power_doc_filter_mod$y) %>%
  round(digits = 0)
power_doc_filter_mod$z <- as.numeric(power_doc_filter_mod$z) %>%
  round(digits = 0)

## re-assign the networks into 6
### DMN, FPN(fronto-parietal network), Salience, Attention(ATN), somatosensory (SMN), visual, auditory
### ATN contains dorsal and ventral attention networks, and cingulo-opercular network
### SMN contains all somatosensory-related and auditory networks
### reference (Watanabe & Rees, 2017, Nat. Comm.)
### 1 - visual | 2 - somMot | 3 - salience | 4 - ATN 
### 5 - FPN | 6 - DMN
power_doc_filter_mod_reassign <- power_doc_filter_mod %>%
  mutate(clust_new = ifelse(clust_ori == 1 | clust_ori == 2 , 2,
                            ifelse(clust_ori == 7, 1,
                                   ifelse(clust_ori == 9, 3,
                                          ifelse(clust_ori == 3| clust_ori == 11 | clust_ori == 12, 4,
                                                        ifelse(clust_ori == 8, 6, 7))))))

power_doc_filter_mod_reassign_rename <- power_doc_filter_mod_reassign %>%
  mutate(network_new = ifelse(clust_new == 1, "visual",
                              ifelse(clust_new == 2, "somMot",
                                     ifelse(clust_new == 3, "salience",
                                            ifelse(clust_new == 4, "ATN",
                                                   ifelse(clust_new == 5, "Aud",
                                                          ifelse(clust_new == 6, "FPN", "DMN")))))))

rio::export(power_doc_filter_mod_reassign_rename, file = "inputs/Power264_Yeo7.xlsx")


##########################################################
#
# Deal with the left td !!!!!!!!S19 did not contain Power264
#
###########################################################
# make target directory
if (!dir.exists("REST/time_lag_estimation")) {
  dir.create("REST/time_lag_estimation", recursive = T)
}
# collect the time delay 
dat_sbj_name <- dat_sbj$ID

# the original path
path <- "/media/lqj/research_data/REST/ROISignals_FunImgARCWF"

# find the left 64 subjects
file_list <- list.files(path)

# obtain the timeseries
subject_utility <- c()
for (i in dat_sbj_name) {
  print(i)
  test <- R.matlab::readMat(file.path(path, bruceR::Glue("ROISignals_{i}.mat")))
  dat <- test$ROISignals
  if (dim(dat)[2]<1833) {
    print(bruceR::Glue("{i} did not contain Power atlas, skip this subject"))
    next
  }
  dat <- data.frame(test$ROISignals[,c(1569 + use_region)])
  # test if any column value equals zero
  if (mean(apply(dat, 2, function(x) all(x==0)))) {
    print(bruceR::Glue("{i} has a zero column, skip this subject"))
    next
  }
  subject_utility <- c(subject_utility, i)
  write_csv(dat, file = bruceR::Glue("REST/{i}_power264.csv"), col_names = F)
}

# obtain the FD file
path_fd <- "/media/lqj/research_data/REST/FD"
file_list <- list.files(path_fd, pattern = "FD_Power_", full.names = T)
# match the subject data in use
file_use <- file_list[grep(paste(subject_utility, collapse = "|"), file_list)]
# move the td maps
file.copy(file_use, to = "REST")


write_tsv(data.frame(subject_utility), file = "inputs/REST_subject_list.tsv", col_names = F)

# *********** The left subjects' TD should be estimated in MATLAB ************