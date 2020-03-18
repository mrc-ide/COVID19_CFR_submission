####################################################################################
##                                                                                ##
##  1) Data Collation and Preprocessing - An Overview                             ##
##                                                                                ##
##   This section loads in data on the demographic structure of the Chinese       ##
##   population, as well as daily case onset incidence data for Wuhan and rest of ##
##   China from the the recent WHO Mission Report. This is integrated with        ##
##   information on the Age distribution of cases from the recent Chinese CDC     ##
##   paper to calculate the age distribution of cases for Wuhan and the rest of   ##
##   China specifically.                                                          ##
##                                                                                ##
##   Assuming an invariant age-distribution of cases over time, this is then      ##
##   used to convert the aggregate daily incidence of cases (by symptom onset)    ##
##   into age-specific daily case incidences.                                     ##
##                                                                                ##
####################################################################################

# Loading Libraries
library(tidyverse)

# Loading in Demographic and Case Data
population <- read.csv("data/china_population_demography.csv") 
WHO_mission_report_onset <- read.csv("data/who_mission_confirmed_and_clinical_diagnosis_onset_data.csv")

# Setting Static Variables
smoothing_centre <- as.Date("2020-02-01", format = "%Y-%m-%d")
days_to_smooth_either_side <- 3
age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

####################################################################################
##                                                                                ##
## Initial Preprocessing of Data Including 1st Feb Incidence Spike Smoothing      ##
##                                                                                ##
##    Ensuring data variables are in correct format, and calculating              ##
##    the incidence of cases outside China over time.                             ##
##    Also contains code to correct for incidence spike that occurred on          ##
##    February 1st. Initial linear interpolation of case numbers between          ##
##    January 31st and 2nd February, followed by distribution of cases            ##
##    equally around 3 days either side.                                          ##
##                                                                                ##
####################################################################################

# Loading in Raw Case Onset Incidence Data
cases_by_onset <- WHO_mission_report_onset %>%
  mutate(Date = as.Date(x = Date, format = "%d/%m/%Y")) %>%
  select(Date, Cases_Inside_Wuhan, Cases_China) %>%
  mutate(Cases_Outside_Wuhan = Cases_China - Cases_Inside_Wuhan)

# Smoothing of Incidence Spike Around 1st February 
## Creating vector of dates to smooth over
total_days_smoothing <- 2 * days_to_smooth_either_side + 1
dates_vector <- c()
start_date <- smoothing_centre - days_to_smooth_either_side
for (i in 1:total_days_smoothing) {
  dates_vector[i] <- as.character(start_date)
  start_date <- start_date + 1
}

## Linear interpolation of cases on 1st Feb from cases on 31st Jan & 2nd Feb
cases_31st_Jan <- cases_by_onset[cases_by_onset$Date == "2020-01-31", 2:4]
cases_1st_Feb <- cases_by_onset[cases_by_onset$Date == "2020-02-01", 2:4]
cases_2nd_Feb <- cases_by_onset[cases_by_onset$Date == "2020-02-02", 2:4]
adj_cases_1st_Feb <- (cases_31st_Jan + cases_2nd_Feb)/2 # linear intepolation

## Redstributing extra cases to surrounding 3 days either side & 1st Feb equally 
cases_by_onset[cases_by_onset$Date == "2020-02-01", 2:4] <- adj_cases_1st_Feb
cases_to_distribute <- (cases_1st_Feb - adj_cases_1st_Feb)/total_days_smoothing
rows_for_adjustment <- cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4]
for (i in 1:nrow(rows_for_adjustment)) {
  rows_for_adjustment[i, ] <- rows_for_adjustment[i, ] + cases_to_distribute
}
cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4] <- rows_for_adjustment

####################################################################################
## Calculating Age Distribution of Cases Inside and Outside Wuhan                 ##
##                                                                                ##
##    Using information from Figures 1A (Wuhan case age distribution) and         ##
##    Figure 1C (all China cases age distribution, including Wuhan), the age-     ##
##    distribution of cases outside Wuhan is calculated.                          ##
##                                                                                ##
####################################################################################

# Age Distribution of China's Population
China_all_age <- population$proportion[population$location == "China"]
China_all_age <- China_all_age/100

# Age Distributions of Cases for Wuhan and Nationally (Including Wuhan)
Wuhan_case_age <- c(0.4000000000001, 0.399999999999, 4.5, 13.1, 15.6, 22.0, 26.5, 12.6, 5) 
Wuhan_case_age <- Wuhan_case_age/100 # convert percentage to proportion
China_case_age <- c(0.9, 1.2, 8.1, 17.0, 19.2, 22.4, 19.2, 8.8, 3.2) 
China_case_age <- China_case_age/100 # convert percentage to proportion

# Calculating the Age Distribution of Cases Outside Wuhan
Wuhan_Cases_Onset <- cases_by_onset$Cases_Inside_Wuhan # Wuhan only
China_Cases_Onset <- cases_by_onset$Cases_China # All of China, including Wuhan
Outside_Wuhan_Cases <- sum(China_Cases_Onset) - sum(Wuhan_Cases_Onset)
Wuhan_cases_by_age <- sum(Wuhan_Cases_Onset) * Wuhan_case_age # cases by age for Wuhan
Total_cases_by_age <- sum(China_Cases_Onset) * China_case_age # cases by age for China (inc. Wuhan)
Outside_Wuhan_Cases_by_age <- Total_cases_by_age - Wuhan_cases_by_age # cases by age outside Wuhan
Outside_Wuhan_Case_Age_Dist <- Outside_Wuhan_Cases_by_age/Outside_Wuhan_Cases # proportion cases by age outside Wuhan


####################################################################################
## Calculating Age-Disaggregated Onset Incidence Over Time, for Inside &          ##
##  Outside Wuhan                                                                 ##
##                                                                                ##
##    Using the daily onset information for Wuhan and outside Wuhan in            ##
##    conjunction with the age distribution of cases for each of these settings,  ##
##    we are able to calculate onset incidence over time disaggregated by         ##
##    age-group, for each location.                                               ##
##    We also calculate the population/number of cases for each age-group and     ##
##    location as this forms the basis for a demographic adjustment we make       ##
##    later on in the model. China's population is taken from their 2018 National ##
##    Statistics, whilst a population of 11.081 million is assumed for Wuhan,     ##
##    based on recent reports.                                                    ##
##                                                                                ##
####################################################################################

# Calculating the Incidence of Cases In Different Age Groups Over Time for Inside Wuhan
##  Create dataframe of Wuhan case distribution by age
##  Create grid with all combinations of age group and daily case incidence
##  Add date to the dataframe
##  Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_Wuhan <- data.frame(age_groups, Wuhan_case_age)
Wuhan_Age_Cases_Time <- expand.grid(Wuhan_case_age, cases_by_onset$Cases_Inside_Wuhan)
colnames(Wuhan_Age_Cases_Time) <- c("Age_Prop", "Cases")
Wuhan_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(Wuhan_cases_by_age))
Wuhan_Age_Cases_Time <- Wuhan_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_Wuhan, by = c("Age_Prop" = "Wuhan_case_age")) %>%
  mutate(Location = "Wuhan")

# Calculating the Incidence of Cases In Different Age Groups Over Time for Outside Wuhan
##  Create dataframe of outside Wuhan case distribution by age
##  Create grid with all combinations of age group and daily case incidence
##  Add date to the dataframe
##  Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_Outside <- data.frame(age_groups, Outside_Wuhan_Case_Age_Dist)
Outside_Age_Cases_Time <- expand.grid(Outside_Wuhan_Case_Age_Dist, cases_by_onset$Cases_Outside_Wuhan)
colnames(Outside_Age_Cases_Time) <- c("Age_Prop", "Cases")
Outside_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(Outside_Wuhan_Cases_by_age))
Outside_Age_Cases_Time <- Outside_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_Outside, by = c("Age_Prop" = "Outside_Wuhan_Case_Age_Dist")) %>%
  mutate(Location = "Outside")

# Combining the Dataframes Together
age_disaggregated_counts_df <- rbind(Wuhan_Age_Cases_Time, Outside_Age_Cases_Time)
age_disaggregated_counts_df <- age_disaggregated_counts_df[, c("Date", "Location", "age_groups", "Cases")]
colnames(age_disaggregated_counts_df) <- c("date", "location", "age_groups", "cases")

# Integrating Age-Disaggregated Cases With Population to Calculate Adjustment Factor
#   Used in downstream analyses
raw_adjustment_factor_df <- age_disaggregated_counts_df %>%
  group_by(age_groups, location) %>%
  summarise(cases = sum(cases)) %>%
  left_join(population, by = c("age_groups" = "age_groups", "location" = "location")) %>%
  mutate(nici = population/cases) %>%
  select(age_groups, location, nici)

# Joining the Case Incidence and Adjustment Factor Datasets Together
age_disaggregated_counts_df <- age_disaggregated_counts_df %>%
  left_join(raw_adjustment_factor_df, by = c("age_groups" = "age_groups", "location" = "location"))

# Saving Created Dataset
saveRDS(age_disaggregated_counts_df, file = "data/age_disaggregated_onset_incidence_data.rds")



