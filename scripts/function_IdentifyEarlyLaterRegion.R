#######################################
#
# This function used to extract the 
# early and later regions in MDD and HC.
#
# The standard for early and later regions are
# 15% (400x.15=60) according to previous study.
#
# Liang Qunjun 2023/08/04

IdentifyEarlyLaterRegion <- function(dat, percentage){
  
  library(tidyverse)
  # divide two groups
  dat_mdd <- dat %>% filter(group == 'MDD') %>% select(2:401)
  dat_hc <- dat %>% filter(group == 'HC') %>% select(2:401)
  
  # define the cut-off
  cut_off <- 400*percentage
  
  dat_early <- c()
  dat_later <- c()
  for (i in 1:nrow(dat_mdd)) {
    dat_tmp <- dat_mdd[i,] %>% t()%>% as.vector()
    early_ind <- sort(dat_tmp)[cut_off]
    later_ind <- sort(dat_tmp)[400 - cut_off]
    
    dat_tmp_early <- rep(0, times = length(dat_tmp))
    dat_tmp_early[which(dat_tmp < early_ind)] <- 1
    
    dat_tmp_later <- rep(0, times = length(dat_tmp))
    dat_tmp_later[which(dat_tmp > later_ind)] <- 1
    
    if (i == 1) {
      dat_early <- dat_tmp_early
      dat_later <- dat_tmp_later
    }else{
      dat_early <- (dat_early+dat_tmp_early)/2
      dat_later <- (dat_later+dat_tmp_later)/2
    }
  }
  
  # combine the early and later region in a whole
  dat_whole <- dat_early - dat_later
  dat_whole_bin <- dat_whole
  dat_whole_bin[which(dat_whole>0)] <- 1
  dat_whole_bin[which(dat_whole<0)] <- -1
  
  # collect the result for mdd
  dat_res_mdd <- data.frame(
    mdd_early = dat_early,
    mdd_later = dat_later,
    mdd_whole = dat_whole,
    mdd_wholeBin = dat_whole_bin
  )
  
  dat_early <- c()
  dat_later <- c()
  for (i in 1:nrow(dat_hc)) {
    dat_tmp <- dat_hc[i,] %>% t()%>% as.vector()
    early_ind <- sort(dat_tmp)[cut_off]
    later_ind <- sort(dat_tmp)[400 - cut_off]
    
    dat_tmp_early <- rep(0, times = length(dat_tmp))
    dat_tmp_early[which(dat_tmp < early_ind)] <- 1
    
    dat_tmp_later <- rep(0, times = length(dat_tmp))
    dat_tmp_later[which(dat_tmp > later_ind)] <- 1
    
    if (i == 1) {
      dat_early <- dat_tmp_early
      dat_later <- dat_tmp_later
    }else{
      dat_early <- (dat_early+dat_tmp_early)/2
      dat_later <- (dat_later+dat_tmp_later)/2
    }
  }
  # combine the early and later region in a whole
  dat_whole <- dat_early - dat_later
  dat_whole_bin <- dat_whole
  dat_whole_bin[which(dat_whole>0)] <- 1
  dat_whole_bin[which(dat_whole<0)] <- -1
  
  # collect the result for HC
  dat_res_hc <- data.frame(
    hc_early = dat_early,
    hc_later = dat_later,
    hc_whole = dat_whole,
    hc_wholeBin = dat_whole_bin
  )
  
  dat_res_out <- cbind(dat_res_mdd, dat_res_hc)
  
  return(dat_res_out)
}
