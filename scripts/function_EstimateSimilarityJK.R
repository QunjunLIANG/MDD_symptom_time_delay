##########################################
#
# This function is used to estimate the 
# similarity with jack-knife method
#
# Liang Qunjun 2023/08/10

EstimateSimilarityJK <- function(network_name, index, mat_hamd, mat_td, n_boot = 5000){
  
  library(boot)
  library(coin)
  # define function to calculate spearman correlation fit the bootstrap
  spearman_corr <- function(data, indices){
    sampled_data <- data[indices,]
    res_cor <- cor(sampled_data$x, sampled_data$y, method = "spearman")
    return(res_cor)
  }
  
  # obtain item-level HAMD score
  mat_hamd_vec <- mat_hamd[upper.tri(mat_hamd)] # transform as vector
  mat_td_vec <- mat_td[upper.tri(mat_td)]
  res_perm <- spearman_test(mat_td_vec ~ mat_hamd_vec) # permutation test 
  
  ## bootstrap for confidential interval
  dat_boot <- data.frame(x = mat_hamd_vec, y = mat_td_vec) # initial data frame for bootstrap
  res_boot <- boot::boot(data = dat_boot, statistic = spearman_corr, R = n_boot)
  res_boot_ci <- boot::boot.ci(res_boot, type = "perc")
  
  ## collect results in this part
  dat_res <- data.frame(
    off_net = network_name,
    index = index,
    cor = res_boot$t0, 
    perm_p = coin::pvalue(res_perm),
    lower_interval = res_boot_ci$percent[4],
    upper_interval = res_boot_ci$percent[5])
}
