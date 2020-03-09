
# ------------------------------------------------------------------
# process international data
prepare_data_international <- function(data,
                                       first_reliable_onset = as.Date("2020-01-01"),
                                       date_format = "%Y-%m-%d") {
  
  # format dates
  data$date_onset <- as.Date(data$date_onset, date_format)
  data$date_outcome <- as.Date(data$date_outcome, date_format)
  data$date_report <- as.Date(data$date_report_dd_mm_yyyy, date_format)
  data$date_admission <- as.Date(data$date_hospitalized_dd_mm_yyyy, date_format)
  data$date_confirmed <- as.Date(data$date_confirmed_dd_mm_yyyy, date_format)
  
  # make dates relative to an index day
  data$rel_date_report <- as.numeric(data$date_report - first_reliable_onset)
  data$rel_date_onset <- as.numeric(data$date_onset - first_reliable_onset)
  data$rel_date_outcome <- as.numeric(data$date_outcome - first_reliable_onset)
  
  # ensure case has outcome date, along with either onset or report date
  data <- subset(data, !is.na(date_outcome) & !(is.na(data$date_onset) & is.na(data$date_report)))
  
  # specify which dates require imputation
  data$date_onset_imputed <- is.na(data$date_onset)
  
  # bundle local transmission criteria into yes/no
  yes_levels <- c("y", "y - confirm", "y - confirmed", "y - implied", "y -implied")
  no_levels <- c("n", "n - implied")
  data$local_transmission_TRUE_FALSE <- NA
  data$local_transmission_TRUE_FALSE[data$local_transmission_y_n %in% yes_levels] <- TRUE
  data$local_transmission_TRUE_FALSE[data$local_transmission_y_n %in% no_levels] <- FALSE
  
  # tidy up
  data <- subset(data, select = c("record_number", "rel_date_report", "rel_date_onset", "rel_date_outcome",
                                  "outcome", "date_onset_imputed", "age_years", "local_transmission_TRUE_FALSE"))
  row.names(data) <- NULL
  
  return(data)
}
