
# process_mcmc_data.R
#
# Author: Bob Verity
# Date: 2020-03-09
#
# Purpose:
# Read in raw international data, process into convenient form and save into
# outputs folder.
#
# ------------------------------------------------------------------
# Load packages and functions

# source functions
source("source/functions_data.R")

# ------------------------------------------------------------------
# Data processing

# read in raw data
data_raw <- read.csv("data/20200225-184311-23838f15.exported_cases_cleaned.csv",
                     stringsAsFactors = FALSE)

# process data
data_processed <- prepare_data_international(data_raw)

# save to output folder
write.csv(data_processed, "output/data_international_processed_mcmc.csv", row.names = FALSE)
