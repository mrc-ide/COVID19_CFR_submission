# Summarise a single column of mcmc output
summarise_mcmc_single <- function(x, q = c(0.025, 0.975)){
  
  # compute kernel density and normalise
  d <- density(x)
  d$y <- d$y/sum(d$y)
  
  # get mode
  x_mode <- d$x[which.max(d$y)]
  x_mean = mean(x)
  x_median = median(x)
  x_lower = quantile(x, q[1], names = FALSE)
  x_upper = quantile(x, q[2], names = FALSE)
  
  return(c(mode = x_mode, mean = x_mean, median = x_median, x_0.025 = x_lower, x_0.975 = x_upper))
  
}
# Summarise multiple columns of MCMC output
summarise_mcmc <- function(x, q = c(0.025, 0.975)){
  if(is.vector(x)){
    x <- matrix(x)
  }
  data.frame(t(apply(x, 2,  summarise_mcmc_single, q = q)))
}
