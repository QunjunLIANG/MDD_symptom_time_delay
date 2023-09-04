#######################################
#
# This function is used to obtain the 
# time delay projecton maps from inputs
# while matching the selected participants
#
# This function can load the data of TDp and 
# ETS. Indicates what kinds of data you want to
# load in [type] parameter

ObtainBrainData <- function(sbj_use, file_list, type="TDp"){
  
  # match the subject data in use
  file_use <- file_list[grep(paste(sbj_use, collapse = "|"), file_list)]
  
  dat_td <- data.frame()
  if (type == "TDp") {
    for (i in file_use) {
      sbj_tmp <- str_extract(i, pattern = 'sub-[0-9]*')
      td_tmp <- readr::read_csv(i, col_names = F, show_col_types = FALSE)
      dat_tmp <- data.frame(participant_id = sbj_tmp, td_tmp)
      dat_td <- rbind(dat_td, dat_tmp)
    }
  }else{
    for (i in file_use) {
      sbj_tmp <- str_extract(i, pattern = 'sub-[0-9]*')
      td_tmp <- readr::read_csv(i, col_names = F, show_col_types = FALSE) %>%
        t()
      dat_tmp <- data.frame(participant_id = sbj_tmp, td_tmp)
      dat_td <- rbind(dat_td, dat_tmp)
    }
  }
  
  return(dat_td)
}