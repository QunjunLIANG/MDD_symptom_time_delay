#############################################
#
# These functions are used to estimate the 
# duration between peak and the hight after 
# ETS rss estimation

ETSmetrices <- function(sbj, rss_raw){
  rss <- rss_raw %>% 
    pivot_longer(cols = 2:(t_point+1), names_to = "time") %>%
    filter(participant_id == sbj) %>%
    .$value
  
  # initial vectors
  t_peak <- c()
  t_trough <- c()
  # detect the peak and trough in loop
  for (i in 2:(length(rss)-1)) {
    if (rss[i] > rss[i-1] & rss[i] > rss[i+1]) {
      t_peak <- c(t_peak, i) # obtain the position
    }
    if (rss[i] < rss[i-1] & rss[i] < rss[i+1]) {
      t_trough <- c(t_trough, i)
    }
  }
  # calculate the difference between trough
  peak_duration <- diff(t_trough)
  peak_high <- rss[t_peak]
  trough_low <- rss[t_trough]

    res_table <- data.frame(
    participant_id = sbj,
    
    duration_n = length(peak_duration),
    duration_mean = mean(peak_duration), 
    duration_sd = sd(peak_duration),
    duration_cv = sd(peak_duration)/mean(peak_duration), 
    duration_range = range(peak_duration),
    
    peak_n = length(peak_high),
    peak_mean = mean(peak_high), 
    peak_sd = sd(peak_high),
    peak_cv = sd(peak_high)/mean(peak_high),
    peak_range = range(peak_high),
    
    trough_n = length(trough_low),
    trough_mean = mean(trough_low), 
    trough_sd = sd(trough_low),
    trough_cv = sd(trough_low)/mean(trough_low),
    trough_range = range(trough_low),
    
    peak_duration_ratio = length(peak_high)/length(peak_duration)
    )
    
  return(res_table)
}
