####################################################################################
##                                                                                ##
##  2) Model Fitting and Parameter Inference                                      ##
##                                                                                ##
##    The age-disaggregated incidence data for Wuhan and outside Wuhan is then    ##
##    integrated with age-disaggregated deaths from Wuhan and outside Wuhan over  ##
##    the same time period. We then fit an age-specific CFR model to these data,  ##
##    which explicitly includes considerations of under-reporting related to      ##
##    age-specific patterns of disease severity as well as the capacity of local  ##
##    health systems.                                                             ##
##                                                                                ##
####################################################################################

# Loading Libraries
set.seed(101092)
source("source/mcmc_and_likelihood_functions.R")
library(dplyr); library(ggplot2); library(tidyverse); library(drjacoby) 

# Loading in Age and Location Disaggregated Case Data
age_disaggregated_case_onset <- readRDS("data/age_disaggregated_onset_incidence_data.rds")
age_cases <- age_disaggregated_case_onset %>%
  group_by(age_groups) %>%
  summarise(counts = sum(cases))
data <- age_disaggregated_case_onset %>%
  arrange(location, age_groups, date)

# Setting Static Variables
age_groups <- c("0_9", "10_19", "20_29", "30_39", "40_49", "50_59", "60_69", "70_79", "80+")
n_age_bands <- length(unique(data$age_groups))
max_date <- as.numeric(max(data$date) - min(data$date)) # most recent date in the date, relative to date start
observed_deaths <- c(0, 1, 7, 18, 38, 130, 309, 312, 208) # (by age, Table 1 of China CDC report)
prop_deaths_wuhan <- 820 / 1113  # scraped from official China Health Commission report, by date reporting end = 11-02-2020
                                 # note that the exact numbers differs slightly from China CDC paper. Assume proportion is correct
                                 # and apply this proportion to the number of deaths detailed in the China CDC paper.
deaths_wuhan <- round(prop_deaths_wuhan * sum(observed_deaths))
deaths_outside <- (1 - prop_deaths_wuhan) * sum(observed_deaths)
death_observation_censoring <- as.numeric(as.Date("2020-01-21") - min(data$date)) # Assume no deaths detecte before 21st Jan in China CDC numbers
deaths_x <- round(prop_deaths_wuhan * sum(observed_deaths))
deaths_n <- sum(observed_deaths)
# WHO sit rep #43 (03-03-2020):
WHO_observed_deaths <- 2946
WHO_observed_cases <- 80304

# Calculating the Prevalence on Repatriation Flights Spanning 30th January - 1st February Inclusive
# Persons tested on arrival
repat_flight_data <- read.csv("data/repat_flights.csv")
repat_flight_data <- repat_flight_data[(repat_flight_data$number_repatriated == repat_flight_data$number_tested) &
                                         !is.na(repat_flight_data$number_tested) & repat_flight_data$point_tested == "arrival", ]
sum(repat_flight_data$number_positive_initial_test)/sum(repat_flight_data$number_tested)
repat_x <- sum(repat_flight_data$number_positive_initial_test) # repatriation flight data number positive (on arrival)
repat_n <- sum(repat_flight_data$number_tested) # repatriation flight data number tested

# Input data
x <- c(0, n_age_bands,
       death_observation_censoring,
       min(data[data$location == "Outside", "nici"]),
       deaths_x, deaths_n,
       repat_x, repat_n,
       WHO_observed_deaths, WHO_observed_cases,
       observed_deaths,
       data$cases, as.numeric(data$age_groups),
       data$date - min(data$date), as.numeric(data$location == "Wuhan"),
       data$nici)

# Input data for output = prediction
x_output <- x
x_output[1] <- 1
saveRDS(x_output, "data/predicted_x.rds")

# Parameters
#   z scales wuhan relative to outside
#   r growth rate (fixed)
#   D detection window (fixed)
#   y and extra scaling used to match repatriation flight prevalence - not used in this simplified likelihood
df_paramsRR <- data.frame(name = c("m_od", "s_od", "maxday", "z", "r", "D", "y", "cfr_80+", "RR_0_9", "RR_10_19", "RR_20_29", "RR_30_39", "RR_40_49", "RR_50_59", "RR_60_69", "RR_70_79"),
                          min = c(10, 0, max_date, 0, 0, 7, 0, rep(0, 9)),
                          max = c(Inf, Inf, max_date, 1, 0.1, 14, 1, 1, rep(Inf, 8)),
                          init = c(15, 0.35, max_date, 0.005, 0.05, 10, 0.1, 0.1, rep(1, 8)))
params <- df_paramsRR$init
lL_plane(df_paramsRR$init, x)

### Run MCMC ###################################################################
burnin <- 50000
samples <- 200000
mcmc_output <- run_mcmc(data = x,
                        df_params = df_paramsRR,
                        loglike = lL_plane,
                        logprior = lP_plane,
                        burnin = burnin,
                        samples = samples,
                        chains = 1)
# (Run with chains > 1 for convergence checks)
saveRDS(mcmc_output, "data/complete_output.rds")
saveRDS(mcmc_output$output, "data/MCMC_fitting_output.rds") # save the MCMC output
saveRDS(list(inputs = x, 
             parameters = df_paramsRR,
             samples = samples,
             burnin = burnin), "data/MCMC_inputs_and_parameters.rds") # save the data list used to run the MCMC
