######################################################
# Functions needed when preparing/running casefat code.
######################################################

PrepareExportedCasesDataset <- function(exported){
  
  date_format <- "%Y-%m-%d"
  names(exported)[grep("date_outcome",names(exported))] <- "date_outcome_factor"
  exported$date_onset <- as.Date(as.character(exported$date_onset_symptoms_dd_mm_yyyy), date_format)
  exported$date_admission <- as.Date(as.character(exported$date_hospitalized_dd_mm_yyyy), date_format)
  exported$date_outcome <- as.Date(as.character(exported$date_outcome_factor), date_format)
  exported$date_confirmed <- as.Date(as.character(exported$date_confirmed_dd_mm_yyyy), date_format)
  #exported$date_travelled_from_wuhan <- as.Date(as.character(exported$date_travelled_from_wuhan_dd_mm_yyyy), date_format)
  exported$date_report <- as.Date(as.character(exported$date_report_dd_mm_yyyy), date_format)
  
  analysis_date <- max(exported$date_outcome,na.rm=T) #as.Date("26/02/2020", date_format) 
  
  # set outcome date equal to analysis_date for those who have not died or recovered
  exported$date_outcome[which(exported$outcome=='other')] <-analysis_date
  
  # where is time_onset used? check this. 
  exported$time_onset<-0
  exported$time_onset[which(is.na(exported$date_onset))]<-NA
  
  exported$time_onset_admission <- as.numeric(exported$date_admission - exported$date_onset)
  exported$time_onset_outcome <- as.numeric(exported$date_outcome - exported$date_onset)
  exported$time_onset_outcome_incl_report <- exported$time_onset_outcome
  
  exported$time_onset_report<-NA
  inds<-which(!is.na(exported$date_onset) & !is.na(exported$date_report))
  exported$time_onset_report[inds]<-as.numeric(exported$date_report[inds]-exported$date_onset[inds])
  exported$time_onset_report[which(exported$time_onset_report<0)]<-NA  ## one person was asymptomatic when reported then onset afterwards
  
  # for people without date onset, set onset as earliest among date report, date confirmed and date admission
  inds <- which(is.na(exported$date_onset) & 
                  (!is.na(exported$date_report) | 
                     !is.na(exported$date_admission) | 
                     !is.na(exported$date_confirmed)))
  
  exported$latest_onset_date <- NA
  
  temp_data <- exported[inds,c("date_report","date_admission","date_confirmed", "local_transmission_y_n")]
  exported$latest_onset_date[inds] <- apply(temp_data,1,min,na.rm=TRUE)
  exported$time_onset_outcome_incl_report[inds] <- as.numeric(as.Date(exported$date_outcome[inds]) - 
                                                                as.Date(exported$latest_onset_date[inds]))
  exported$time_onset_outcome_incl_report[which(exported$time_onset_outcome_incl_report<0)] <- NA
  
  # check why transformed to as.character
  exported$latest_onset_date[inds] <- as.character(as.Date(exported$latest_onset_date[inds], date_format))
  
  exported$local_transmission_y_n<-as.character(exported$local_transmission_y_n)
  exported$local_bin<-NA
  exported$local_bin[which(exported$local_transmission_y_n=="n" |
                             exported$local_transmission_y_n=="n - implied")] <-0
  exported$local_bin[which(exported$local_transmission_y_n=="y" |
                             exported$local_transmission_y_n=="y - implied" |
                             exported$local_transmission_y_n=="y -implied" |
                             exported$local_transmission_y_n=="y - confirmed")] <-1
  
  
  
  # define outcome as factor
  levels(exported$outcome)[which(levels(exported$outcome)=="death")]<-"Death"
  levels(exported$outcome)[which(levels(exported$outcome)=="other")]<-"Other"
  levels(exported$outcome)[which(levels(exported$outcome)=="recovery")]<-"Discharge"
  #names(exported)[which(names(exported)=="outcome_censored")]<-"outcome_interval_censored"
  
  exported$age_over60<-NA
  exported$age_over60[which(exported$age_years<60)] <-0
  exported$age_over60[which(exported$age_years>=60)] <-1
  
  
  # keep columns used in analysis
  exported <- exported[, c("country", "record_number", "date_report", "date_onset", "date_confirmed", "date_admission", 
                           "date_outcome", "latest_onset_date", "outcome", "time_onset_outcome", "time_onset_outcome_incl_report", 
                           "time_onset_report",  "local_transmission_y_n", "local_bin", "recovered_y_n",
                           "age_years", "age_over60")]
  
  return(exported)
}



##########################################################################
# IMPUTE ONSET DATES
##########################################################################

imputeOnsets<-function(exported, times_onset_report) {
  exported$date_report<-as.Date(exported$date_report)
  exported$date_outcome<-as.Date(exported$date_outcome)
  exported$date_onset<-as.Date(exported$date_onset)
  exported$time_onset_outcome_impute_from_report<-exported$time_onset_outcome_incl_report
  #times_onset_report<-exported$time_onset_report[which(!is.na(exported$time_onset_report))]
  inds<-which(is.na(exported$date_onset) & !is.na(exported$date_report))
  exported$time_onset_outcome_impute_from_report[inds]<-exported$time_onset_outcome_incl_report[inds] +
    sample(times_onset_report,size=length(inds),replace=T)
  #rgamma(length(inds),shape=shape_O2Rep_noR,scale=scale_O2Rep_noR)
  ## remove next line when data cleaned.
  exported$time_onset_outcome_impute_from_report[which(exported$time_onset_outcome_impute_from_report<0)]<-NA
  return(exported)
}

##########################################################################
# IMPUTE WHO HAS RECOVERED
##########################################################################
imputeRecoveries<-function(exported_incl_missing_onset) {
  exported_incl_missing_onset$outcome_impute<-exported_incl_missing_onset$outcome
  exported_incl_missing_onset$time_onset_outcome_impute_report_and_rec<-exported_incl_missing_onset$time_onset_outcome_impute_from_report
  countries<-as.character(exported_country_summ$country[which(exported_country_summ$unassigned_recoveries>0)])
  exported_incl_missing_onset$recovered_y_n<-as.character(exported_incl_missing_onset$recovered_y_n)
  countries<-as.character(exported_country_summ$country[which(exported_country_summ$unassigned_recoveries>0)])
  for(i in 1:length(countries)) {
    curr_country<-countries[i]  ## work by country
    n_recs<-exported_country_summ$unassigned_recoveries[which(exported_country_summ$country==curr_country)]
    inds<-which(exported_incl_missing_onset$country==curr_country & exported_incl_missing_onset$outcome=='Other'
                & (is.na(exported_incl_missing_onset$recovered_y_n) | 
                     exported_incl_missing_onset$recovered_y_n=="n - implied" | 
                     exported_incl_missing_onset$recovered_y_n=="y - implied"))
    if(length(inds)<n_recs) {
      #print(paste0("In ", curr_country," not enough eligible people to recover"))
      #print(n_recs)
      #print(length(inds))
      new_recs<-inds
    } else {
      weights<- pgamma(exported_incl_missing_onset$time_onset_outcome_impute_from_report[inds],shape =shape_O2R ,scale =scale_O2R ,lower.tail =T)
      # safer version of sample that can deal with vectors of length 1.
      resample <- function(x,size0,replace0=F,prob0=NULL) x[sample.int(length(x),size=size0,replace=replace0,prob=prob0)]
      new_recs<-resample(x=inds, size0=n_recs,replace0=F,prob0=weights)
    }
    exported_incl_missing_onset$outcome_impute[new_recs]<-'Discharge'
    exported_incl_missing_onset$time_onset_outcome_impute_report_and_rec[new_recs]<-
      rgamma(length(new_recs),shape =shape_O2R ,scale =scale_O2R)
  }
  return(exported_incl_missing_onset)
}
