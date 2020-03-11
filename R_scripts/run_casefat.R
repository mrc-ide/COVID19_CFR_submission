######################################################################
## NON PARAMETRIC
######################################################################

source("R_scripts/casefat.R")
source("R_scripts/utils_casefat.R")

################################################################################
# read in and reformat exported/international cases dataset
###############################################################################

exported_cases <- read.csv("data/20200225-184311-23838f15.exported_cases_cleaned.csv")
exported <- PrepareExportedCasesDataset(exported_cases)
exported_country_summ<-read.csv("data/summary_by_country_cleaned.csv")

exported_incl_missing_onset <- subset(exported, !is.na(exported$time_onset_outcome_incl_report))


#O2R for imputing recoveries
mean<-22.59
s<-0.325
shape_O2R <- 1/(s^2) 
scale_O2R <- 1/(1/(mean*s^2))

#####################################################################
# NON PARAMETRIC ANALYSIS
#####################################################################

########### Everyone.
data<-exported_incl_missing_onset
data$f<-0
data$f[which(data$outcome=="Death")]<-1
data$f[which(data$outcome=="Discharge")]<-2
dim(data)

times_onset_report_obs<-exported_incl_missing_onset$time_onset_report[which(!is.na(exported_incl_missing_onset$time_onset_report))]

nIter<-500
cf_res<-data.frame(logit_cfr=rep(NA,nIter),se_logit_cfr=NA)
for(i in 1:nIter) {
  ## must impute in full dataset or recovery numbers won't match country summaries.
  exported_incl_missing_onset<-imputeOnsets(exported_incl_missing_onset, times_onset_report_obs)
  exported_incl_missing_onset<-imputeRecoveries(exported_incl_missing_onset)
  data<-exported_incl_missing_onset
  data$f<-0
  data$f[which(data$outcome_impute=="Death")]<-1
  data$f[which(data$outcome_impute=="Discharge")]<-2
  cf<-casefat(data$time_onset_outcome_impute_report_and_rec,
              data$f)
  cf_res$logit_cfr[i]<-cf$logit_cfr
  cf_res$se_logit_cfr[i]<-cf$se_logit_cfr
}

mean_logit_cfr<-mean(cf_res$logit_cfr)
invlogit<-function(x) exp(x)/(1+exp(x))
invlogit(mean_logit_cfr)
#Within: V_W = 1/m sum_m (S.E_i)^2 where S.E_i is the standard error for imputed set i (calculated using formula 7) 
within_var<-1/nIter * sum(cf_res$se_logit_cfr^2)
#Between: V_B = 1/m-1 * (sum_m [cfr_i – cfr_hat]^2)
between_var<- 1/(nIter -1) * (sum((cf_res$logit_cfr-mean_logit_cfr)^2))
tot_var<-within_var+between_var
## lower bound
invlogit(mean_logit_cfr-1.96*sqrt(tot_var))
invlogit(mean_logit_cfr+1.96*sqrt(tot_var))
print(paste0(100*round(invlogit(mean_logit_cfr),3)," (95% CI ",100*round(invlogit(mean_logit_cfr-1.96*sqrt(tot_var)),3),"-"
             ,100*round(invlogit(mean_logit_cfr+1.96*sqrt(tot_var)),3),")"))




########### Travellers only
data<-exported_incl_missing_onset[which(exported_incl_missing_onset$local_bin==0),]
data$f<-0
data$f[which(data$outcome=="Death")]<-1
data$f[which(data$outcome=="Discharge")]<-2

times_onset_report_obs<-exported_incl_missing_onset$time_onset_report[which(!is.na(exported_incl_missing_onset$time_onset_report))]

nIter<-500
cf_res<-data.frame(logit_cfr=rep(NA,nIter),se_logit_cfr=NA)
for(i in 1:nIter) {
  ## must impute in full dataset or recovery numbers won't match country summaries.
  exported_incl_missing_onset<-imputeOnsets(exported_incl_missing_onset, times_onset_report_obs)
  exported_incl_missing_onset<-imputeRecoveries(exported_incl_missing_onset)
  data<-exported_incl_missing_onset[which(exported_incl_missing_onset$local_bin==0),]
  data$f<-0
  data$f[which(data$outcome_impute=="Death")]<-1
  data$f[which(data$outcome_impute=="Discharge")]<-2
  cf<-casefat(data$time_onset_outcome_impute_report_and_rec,
              data$f)
  cf_res$logit_cfr[i]<-cf$logit_cfr
  cf_res$se_logit_cfr[i]<-cf$se_logit_cfr
}

mean_logit_cfr<-mean(cf_res$logit_cfr)
invlogit<-function(x) exp(x)/(1+exp(x))
invlogit(mean_logit_cfr)
#Within: V_W = 1/m sum_m (S.E_i)^2 where S.E_i is the standard error for imputed set i (calculated using formula 7) 
within_var<-1/nIter * sum(cf_res$se_logit_cfr^2)
#Between: V_B = 1/m-1 * (sum_m [cfr_i – cfr_hat]^2)
between_var<- 1/(nIter -1) * (sum((cf_res$logit_cfr-mean_logit_cfr)^2))
tot_var<-within_var+between_var
## lower bound
invlogit(mean_logit_cfr-1.96*sqrt(tot_var))
invlogit(mean_logit_cfr+1.96*sqrt(tot_var))
print(paste0(100*round(invlogit(mean_logit_cfr),3)," (95% CI ",100*round(invlogit(mean_logit_cfr-1.96*sqrt(tot_var)),3),"-"
             ,100*round(invlogit(mean_logit_cfr+1.96*sqrt(tot_var)),3),")"))



### CFR Casefat local only 
data<-exported_incl_missing_onset[which(exported_incl_missing_onset$local_bin==1),]
data$f<-0
data$f[which(data$outcome=="Death")]<-1
data$f[which(data$outcome=="Discharge")]<-2

times_onset_report_obs<-exported_incl_missing_onset$time_onset_report[which(!is.na(exported_incl_missing_onset$time_onset_report))]

nIter<-500
cf_res<-data.frame(logit_cfr=rep(NA,nIter),se_logit_cfr=NA)
for(i in 1:nIter) {
  ## must impute in full dataset or recovery numbers won't match country summaries.
  exported_incl_missing_onset<-imputeOnsets(exported_incl_missing_onset, times_onset_report_obs)
  exported_incl_missing_onset<-imputeRecoveries(exported_incl_missing_onset)
  data<-exported_incl_missing_onset[which(exported_incl_missing_onset$local_bin==1),]
  data$f<-0
  data$f[which(data$outcome_impute=="Death")]<-1
  data$f[which(data$outcome_impute=="Discharge")]<-2
  cf<-casefat(data$time_onset_outcome_impute_report_and_rec,
              data$f)
  cf_res$logit_cfr[i]<-cf$logit_cfr
  cf_res$se_logit_cfr[i]<-cf$se_logit_cfr
}

mean_logit_cfr<-mean(cf_res$logit_cfr)
invlogit<-function(x) exp(x)/(1+exp(x))
invlogit(mean_logit_cfr)
#Within: V_W = 1/m sum_m (S.E_i)^2 where S.E_i is the standard error for imputed set i (calculated using formula 7) 
within_var<-1/nIter * sum(cf_res$se_logit_cfr^2)
#Between: V_B = 1/m-1 * (sum_m [cfr_i – cfr_hat]^2)
between_var<- 1/(nIter -1) * (sum((cf_res$logit_cfr-mean_logit_cfr)^2))
tot_var<-within_var+between_var
## lower bound
invlogit(mean_logit_cfr-1.96*sqrt(tot_var))
invlogit(mean_logit_cfr+1.96*sqrt(tot_var))
print(paste0(100*round(invlogit(mean_logit_cfr),3)," (95% CI ",100*round(invlogit(mean_logit_cfr-1.96*sqrt(tot_var)),3),"-"
             ,100*round(invlogit(mean_logit_cfr+1.96*sqrt(tot_var)),3),")"))


### CFR Casefat over 60 
data<-exported_incl_missing_onset[which(exported_incl_missing_onset$age_over60==1),]
data$f<-0
data$f[which(data$outcome=="Death")]<-1
data$f[which(data$outcome=="Discharge")]<-2
dim(data)

times_onset_report_obs<-exported_incl_missing_onset$time_onset_report[which(!is.na(exported_incl_missing_onset$time_onset_report))]

nIter<-500
cf_res<-data.frame(logit_cfr=rep(NA,nIter),se_logit_cfr=NA)
for(i in 1:nIter) {
  ## must impute in full dataset or recovery numbers won't match country summaries.
  exported_incl_missing_onset<-imputeOnsets(exported_incl_missing_onset, times_onset_report_obs)
  exported_incl_missing_onset<-imputeRecoveries(exported_incl_missing_onset)
  data<-exported_incl_missing_onset[which(exported_incl_missing_onset$age_over60==1),]
  data$f<-0
  data$f[which(data$outcome_impute=="Death")]<-1
  data$f[which(data$outcome_impute=="Discharge")]<-2
  cf<-casefat(data$time_onset_outcome_impute_report_and_rec,
              data$f)
  cf_res$logit_cfr[i]<-cf$logit_cfr
  cf_res$se_logit_cfr[i]<-cf$se_logit_cfr
}

mean_logit_cfr<-mean(cf_res$logit_cfr)
invlogit<-function(x) exp(x)/(1+exp(x))
invlogit(mean_logit_cfr)
#Within: V_W = 1/m sum_m (S.E_i)^2 where S.E_i is the standard error for imputed set i (calculated using formula 7) 
within_var<-1/nIter * sum(cf_res$se_logit_cfr^2)
#Between: V_B = 1/m-1 * (sum_m [cfr_i – cfr_hat]^2)
between_var<- 1/(nIter -1) * (sum((cf_res$logit_cfr-mean_logit_cfr)^2))
tot_var<-within_var+between_var
## lower bound
invlogit(mean_logit_cfr-1.96*sqrt(tot_var))
invlogit(mean_logit_cfr+1.96*sqrt(tot_var))
print(paste0(100*round(invlogit(mean_logit_cfr),3)," (95% CI ",100*round(invlogit(mean_logit_cfr-1.96*sqrt(tot_var)),3),"-"
             ,100*round(invlogit(mean_logit_cfr+1.96*sqrt(tot_var)),3),")"))


### CFR Casefat under 60 
data<-exported_incl_missing_onset[which(exported_incl_missing_onset$age_over60==0),]
data$f<-0
data$f[which(data$outcome=="Death")]<-1
data$f[which(data$outcome=="Discharge")]<-2

times_onset_report_obs<-exported_incl_missing_onset$time_onset_report[which(!is.na(exported_incl_missing_onset$time_onset_report))]

nIter<-500
cf_res<-data.frame(logit_cfr=rep(NA,nIter),se_logit_cfr=NA)
for(i in 1:nIter) {
  ## must impute in full dataset or recovery numbers won't match country summaries.
  exported_incl_missing_onset<-imputeOnsets(exported_incl_missing_onset, times_onset_report_obs)
  exported_incl_missing_onset<-imputeRecoveries(exported_incl_missing_onset)
  data<-exported_incl_missing_onset[which(exported_incl_missing_onset$age_over60==0),]
  data$f<-0
  data$f[which(data$outcome_impute=="Death")]<-1
  data$f[which(data$outcome_impute=="Discharge")]<-2
  cf<-casefat(data$time_onset_outcome_impute_report_and_rec,
              data$f)
  cf_res$logit_cfr[i]<-cf$logit_cfr
  cf_res$se_logit_cfr[i]<-cf$se_logit_cfr
}

mean_logit_cfr<-mean(cf_res$logit_cfr)
invlogit<-function(x) exp(x)/(1+exp(x))
invlogit(mean_logit_cfr)
#Within: V_W = 1/m sum_m (S.E_i)^2 where S.E_i is the standard error for imputed set i (calculated using formula 7) 
within_var<-1/nIter * sum(cf_res$se_logit_cfr^2)
#Between: V_B = 1/m-1 * (sum_m [cfr_i – cfr_hat]^2)
between_var<- 1/(nIter -1) * (sum((cf_res$logit_cfr-mean_logit_cfr)^2))
tot_var<-within_var+between_var
## lower bound
invlogit(mean_logit_cfr-1.96*sqrt(tot_var))
invlogit(mean_logit_cfr+1.96*sqrt(tot_var))
print(paste0(100*round(invlogit(mean_logit_cfr),3)," (95% CI ",100*round(invlogit(mean_logit_cfr-1.96*sqrt(tot_var)),3),"-"
             ,100*round(invlogit(mean_logit_cfr+1.96*sqrt(tot_var)),3),")"))
