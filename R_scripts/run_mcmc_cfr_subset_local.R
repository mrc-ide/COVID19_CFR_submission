
# run_mcmc_cfr_subset_local.R
#
# Author: Bob Verity
# Date: 2020-03-09
#
# Purpose:
# Perform similar analysis as run_mcmc_cfr_all.R, except subsetting to local vs.
# non-local cases.
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

# define subsetting, and create appropriate suffix for all outputs
subset_local <- FALSE
suff <- ifelse(subset_local, "_local", "_nonlocal")

# read in processed international data
data <- read.csv("output/data_international_processed_mcmc.csv", stringsAsFactors = FALSE)

# subset to ensure locality info present
data <- subset(data, !is.na(local_transmission_TRUE_FALSE))

# define growth rate based on local vs. non-local
data$growth_rate <- ifelse(data$local_transmission_TRUE_FALSE, 0.14, 0.05)

# subset based on locality info. Recoveries that do not pass this filter are
# still retained and used to estimate onset-to-recovery distributions,
# accounting for epidemic growth
if (subset_local) {
  data <- subset(data, local_transmission_TRUE_FALSE == TRUE | outcome == "recovery")
  data$recovery_only <- (data$local_transmission_TRUE_FALSE == FALSE)
} else {
  data <- subset(data, local_transmission_TRUE_FALSE == FALSE | outcome == "recovery")
  data$recovery_only <- (data$local_transmission_TRUE_FALSE == TRUE)
}

# DEBUG - REDUCE DATA SIZE FOR TESTING
#set.seed(1)
#data <- data[sample(nrow(data), 50), ]

# explore data
table(data$local_transmission_TRUE_FALSE, data$outcome)

# get data into drjacoby format
x <- c(data$rel_date_onset,
       data$rel_date_report,
       data$rel_date_outcome,
       match(data$outcome, c("death", "recovery", "other")),
       data$date_onset_imputed,
       data$recovery_only,
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
# m_od: mean onset to death
# s_od: SD to mean ratio onset to death
# m_or: mean onset to recovery
# s_or: SD to mean ratio onset to recovery
# m_op: mean onset to report
# s_op: SD to mean ratio onset to report
# cfr: Case fatality rate
# p: probability of recovery being correctly recorded
df_params <- rbind.data.frame(list("m_od", 0, 100, 10),
                              list("s_od", 0, 1, 0.5),
                              list("m_or", 0, 100, 10),
                              list("s_or", 0, 1, 0.5),
                              list("m_op", 0, 100, 10),
                              list("s_op", 0, 1, 0.5),
                              list("cfr", 0, 1, 0.5),
                              list("p", 0, 1, 0.5))
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
                 loglike = cpp_loglike_cfr,
                 logprior = cpp_otd_prior,
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
hist_cfr_m_od <- posterior_hist(mcmc_samples, "m_od", breaks = seq(10,40,l=101)) + ggtitle("m_od")
hist_cfr_s_od <- posterior_hist(mcmc_samples, "s_od", breaks = seq(0,1,l=101)) + ggtitle("s_od")
hist_cfr_m_or <- posterior_hist(mcmc_samples, "m_or", breaks = seq(10,40,l=101)) + ggtitle("m_or")
hist_cfr_s_or <- posterior_hist(mcmc_samples, "s_or", breaks = seq(0,1,l=101)) + ggtitle("s_or")
hist_cfr_m_op <- posterior_hist(mcmc_samples, "m_op", breaks = seq(0,20,l=101)) + ggtitle("m_op")
hist_cfr_s_op <- posterior_hist(mcmc_samples, "s_op", breaks = seq(0,1,l=101)) + ggtitle("s_op")
hist_cfr_cfr <- posterior_hist(mcmc_samples, "cfr", breaks = seq(0,1,l=101)) + ggtitle("cfr")
hist_cfr_p <- posterior_hist(mcmc_samples, "p", breaks = seq(0,1,l=101)) + ggtitle("p")

# combined plot
cp_cfr <- cowplot::plot_grid(hist_cfr_m_od,
                             hist_cfr_s_od,
                             hist_cfr_m_or,
                             hist_cfr_s_or,
                             hist_cfr_m_op,
                             hist_cfr_s_op,
                             hist_cfr_cfr,
                             hist_cfr_p)
cp_cfr

# posterior summaries
summary_cfr <- apply(mcmc_samples, 2, posterior_summary)
summary_cfr <- rbind(summary_cfr, ESS = round(mcmc$diagnostics$ess[param_names], digits = 0))

summary_cfr

# ------------------------------------------------------------------
# Write to file

if (FALSE) {  # safety catch to avoid accidental overwriting
  
  # save summary table
  write.csv(summary_cfr, sprintf("output/summary_cfr%s.csv", suff), row.names = TRUE)
  
}

