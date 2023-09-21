############################
#
# Identify statistic function 
# to test the difference in 
# demographic information
# 

pvalue <- function(x, ...){
  # construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times = sapply(x, length)))
  if(is.numeric(y)){
    # for numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # for categorical varibables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-valuem using on HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "$lt;", format.pval(p, digits = 3, eps = 0.001)))
}