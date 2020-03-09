
# ------------------------------------------------------------------
# compute mode and quantiles from kernel density of raw draws
posterior_summary <- function(x, prob = c(0.025, 0.975), digits = 3) {
  
  # initialise return object
  np <- length(prob)
  ret <- rep(NA, np + 1)
  
  # special case for all x identical
  if (all(x == x[1])) {
    ret <- rep(x[1], np + 1)
    names(ret) <- c("mode", sprintf("Q%s", prob))
    ret <- round(ret, digits = digits)
    return(ret)
  }
  
  # get mode of kernel density
  d <- density(x)
  ret[1] <- d$x[which.max(d$y)]
  
  # get quantiles
  ret[-1] <- quantile(x, probs = prob)
  
  # tidy up output
  names(ret) <- c("mode", sprintf("Q%s", prob))
  ret <- round(ret, digits = digits)
  
  return(ret)
}

# ------------------------------------------------------------------
# produce ggplot histogram of posterior draws
posterior_hist <- function(df_plot, name, breaks = NULL) {
  
  # set default breaks
  if (is.null(breaks)) {
    breaks <- seq(min(df_plot[[name]]), max(df_plot[[name]]), l = 50)
  }
  
  # produce plot
  ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = .data[[name]]),
                            breaks = breaks, fill = "darkblue", col = NA) + 
    ggplot2::xlab(name) + ggplot2::ylab("Count")
}
