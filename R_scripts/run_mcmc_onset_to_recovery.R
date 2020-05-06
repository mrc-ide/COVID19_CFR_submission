
# run_mcmc_onset_to_recovery.R
#
# Author: Bob Verity
# Date: 2020-03-09
#
# Purpose:
# Run MCMC analysis of onset-to-recovery distribution with Bayesian imputation
# of onset-to-report times.
#
# ------------------------------------------------------------------
# Load packages and functions
# NB. important to use drjacoby vesion1.0 as input formats likely to change in
# future versions

library(devtools)
#devtools::install_github("mrc-ide/drjacoby", ref = "version1.0")  # uncomment and run this line the first time through
library(drjacoby)
library(parallel)
library(ggplot2)
library(cowplot)

# source functions
source("source/functions_mcmc_likelihoods.R")
source("source/functions_utils.R")

# ------------------------------------------------------------------
# Data preperation

# read in processed international data
data <- read.csv("output/data_international_processed_mcmc.csv", stringsAsFactors = FALSE)

# subset to ensure locality info present
data <- subset(data, !is.na(local_transmission_TRUE_FALSE))

# subset to recoveries only
data <- subset(data, outcome == "recovery")

# define growth rate based on local vs. non-local
data$growth_rate <- ifelse(data$local_transmission_TRUE_FALSE, 0.14, 0.05)

# DEBUG - REDUCE DATA SIZE FOR TESTING
#set.seed(1)
#data <- data[sample(nrow(data), 50), ]

# get data into drjacoby format
x <- c(data$rel_date_onset,
       data$rel_date_report,
       data$rel_date_outcome,
       data$date_onset_imputed,
       data$growth_rate)

# ------------------------------------------------------------------
# MCMC parameters

# sampling parameters
burnin <- 1e3
samples <- 1e5
chains <- 5
run_parallel <- TRUE
n_cores <- 5

# specify key parameters
# m_or: mean onset to recovery
# s_or: SD to mean ratio onset to recovery
# m_op: mean onset to report
# s_op: SD to mean ratio onset to report
df_params <- rbind.data.frame(list("m_or", 0, 100, 10),
                              list("s_or", 0, 1, 0.5),
                              list("m_op", 0, 100, 10),
                              list("s_op", 0, 1, 0.5))
names(df_params) <- c("name", "min", "max", "init")

# get key parameter names
param_names <- as.character(df_params$name)

# extend params dataframe with one parameter per imputation point
n_impute <- sum(data$date_onset_imputed)
df_params <- rbind(df_params,
                   data.frame(name = sprintf("delta%s", 1:n_impute), min = 0, max = 50, init = 3))


# ------------------------------------------------------------------
# Run MCMC

# set to run in serial/parallel
cl <- NULL
if (run_parallel) {
  cl <- parallel::makeCluster(n_cores)
}

# run MCMC
set.seed(1)
t0 <- Sys.time()
mcmc <- run_mcmc(data = x,
                 df_params = df_params,
                 loglike = cpp_loglike_otr,
                 logprior = cpp_flat_prior,
                 chains = chains,
                 burnin = burnin,
                 samples = samples,
                 cluster = cl)

Sys.time() - t0
if (run_parallel) {
  parallel::stopCluster(cl)
}

# ------------------------------------------------------------------
# Post-processing

# posterior diagnostic plots
#plot_par(mcmc, show = param_names)

# subset to sampling iterations and desired parameters only
mcmc_samples <- subset(mcmc$output, stage == "sampling", select = param_names)

# posterior histograms
hist_otr_m_or <- posterior_hist(mcmc_samples, "m_or", breaks = seq(10,40,l=101)) + ggtitle("m_or")
hist_otr_s_or <- posterior_hist(mcmc_samples, "s_or", breaks = seq(0,1,l=101)) + ggtitle("s_or")
hist_otr_m_op <- posterior_hist(mcmc_samples, "m_op", breaks = seq(0,20,l=101)) + ggtitle("m_op")
hist_otr_s_op <- posterior_hist(mcmc_samples, "s_op", breaks = seq(0,1,l=101)) + ggtitle("s_op")

# combined plot
cp_otr <- cowplot::plot_grid(hist_otr_m_or,
                             hist_otr_s_or,
                             hist_otr_m_op,
                             hist_otr_s_op)
cp_otr

# posterior summaries
summary_otr <- apply(mcmc_samples, 2, posterior_summary)
summary_otr <- rbind(summary_otr, ESS = round(mcmc$diagnostics$ess[param_names], digits = 0))

summary_otr


# ------------------------------------------------------------------
# Write to file

if (FALSE) {  # safety catch to avoid accidental overwriting
  
  # save summary table
  write.csv(summary_otr, "output/summary_otr.csv", row.names = TRUE)
  
}

