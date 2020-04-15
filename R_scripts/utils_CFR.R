#####################################################
# define functions 
#####################################################

prepareDeathsDataset<-function(data) {
  date_format <- "%d/%m/%Y"
  data <- data[, c("death_id", "age", "date_onset", "date_death", "onset.to.death")]
  
  data$age_over60<-NA
  data$age_over60[which(data$age<60)] <-0
  data$age_over60[which(data$age>=60)] <-1
  
  #####################################################
  # format dates
  #first_reliable_death <- min(as.Date(data$date_death, date_format), na.rm  = TRUE)
  
  
  data$date_onset <- as.Date(data$date_onset, date_format)
  data$date_death <- as.Date(data$date_death, date_format)
  data$onset.to.death <- as.numeric(as.Date(data$date_death, date_format) - 
                                      as.Date(data$date_onset, date_format))
  
  # remove subjects with unknown onset.to.death
  data <- data[which(!is.na(data$onset.to.death)), ]
  
  # only keep subjects with reliable dates 
  data <- data[which(data$date_onset >= first_reliable_onset & data$date_death >= first_reliable_death
                     & data$date_death <= latest_reliable_death),]
  
  # remove subjects with onset.to.death > 40
  #data <- data[which(data$onset.to.death <= 40), ]
  
   
  # number of subjects used in the analysis
  message("number of subjects")
  message(dim(data)[1])
  return(data)
}

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

#####################################################
# calculate onset to outcome distributions
# as function of mean and s (cv = sd/mean)
# corrected for epidemic growth
#####################################################

calcOnsetToOutcome <- function(data,r=0,var_name_outcome = "onset.to.death",grid=grid){
  
  data<-data[!is.na(data[[var_name_outcome]]),]
  alpha <- 1/(grid$s*grid$s) # k in Wikipedia param
  beta <- 1/(1/(grid$mean*grid$s*grid$s)+r)
  prob <- alpha
  lprob <- 0
  
  for(i in 1:nrow(data))
  {
    lprob <- lprob+log(1e-100+(pgamma(data[[var_name_outcome]][i]+1, shape = alpha, scale = beta, lower.tail = TRUE)
                               -pgamma(data[[var_name_outcome]][i], shape = alpha, scale = beta,lower.tail = TRUE))
                       /pgamma(cens[i],shape = alpha, scale = beta, lower.tail = TRUE))
  }
  mlprob <- max(lprob)
  prob <- exp(lprob-mlprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
  
}


#####################################################
# calculate onset to outcome distributions
# as function of mean and s (cv = sd/mean)
# corrected for epidemic growth
#####################################################

calcOnsetToOutcome_2r <- function(data,r_trav=0,r_local=0.14,var_name_outcome = "onset.to.death",grid=grid){
  
  data<-data[!is.na(data[[var_name_outcome]]),]
  alpha <- 1/(grid$s*grid$s) # k in Wikipedia param
  beta_local <- 1/(1/(grid$mean*grid$s*grid$s)+r_local)
  beta_trav <- 1/(1/(grid$mean*grid$s*grid$s)+r_trav)
  prob <- alpha
  lprob <- 0
  
  for(i in 1:nrow(data))
  {
    if(data$local_bin[i]==1) {
      lprob <- lprob+log(1e-100+(pgamma(data[[var_name_outcome]][i]+1, shape = alpha, scale = beta_local, lower.tail = TRUE)
                                 -pgamma(data[[var_name_outcome]][i], shape = alpha, scale = beta_local,lower.tail = TRUE))
                         /pgamma(cens[i],shape = alpha, scale = beta_local, lower.tail = TRUE))
    } else if(data$local_bin[i]==0) {
      lprob <- lprob+log(1e-100+(pgamma(data[[var_name_outcome]][i]+1, shape = alpha, scale = beta_trav, lower.tail = TRUE)
                                 -pgamma(data[[var_name_outcome]][i], shape = alpha, scale = beta_trav,lower.tail = TRUE))
                         /pgamma(cens[i],shape = alpha, scale = beta_trav, lower.tail = TRUE))
    } else {
      lprob<-lprob
    }
  }
  mlprob <- max(lprob)
  prob <- exp(lprob-mlprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
  
}


#####################################################
# define function
# para[1] = mean, param[2]= s, param[3]= cfr
#####################################################
calcCFR_DeathsOnly<-function(data,postO2D,vcfr,var_name_time_outcome) {
  
  alpha <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  beta <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  lprob_onset_to_death <- log(postO2D[,3]+1e-100)
  prob <- vcfr
  
  for(j in 1:length(vcfr))
  {
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_death
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)+cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )
      
    }
    mtprob <- max(tprob)
    prob[j] <- log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
}

calcCFR_DeathsOnlyParallel<-function(data,postO2D,vcfr,var_name_time_outcome) {
  
  alpha <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  beta <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  lprob_onset_to_death <- log(postO2D[,3]+1e-100)
  prob <- vcfr
  
  prob=foreach(j = 1:length(vcfr), .combine = rbind) %dopar% {
    
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_death
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)+cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )
      
    }
    mtprob <- max(tprob)
    log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
}

calc_full_liks_DeathsOnly<-function(data,postO2D,vcfr,var_name_time_outcome) {
  
  alpha <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  beta <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  lprob_onset_to_death <- log(postO2D[,3]+1e-100)
  prob <- vcfr
  
  res<-data.frame(mean=rep(postO2D$mean,length(vcfr)),s=rep(postO2D$s,length(vcfr)))
  res$cfr<-rep(vcfr,each=nrow(postO2D))
  res$prior_OD<-rep(lprob_onset_to_death,length(vcfr))
  res$loglik<-NA
  res$loglik_no_prior<-NA
  
  for(j in 1:length(vcfr))
  {
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_death
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)+cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )
      
    }
    inds<-which((res$cfr+1e-10)==cfr)
    res$loglik[inds]<-tprob
    res$loglik_no_prior[inds]<-tprob-lprob_onset_to_death
    mtprob <- max(tprob)
    prob[j] <- log(sum(exp(tprob-mtprob)))+mtprob
    
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(res)
}

calcCFR_DeathsAndDischarges<-function(data,postO2D,preO2R,vcfr,var_name_time_outcome) {
  
  
  row_r <- sample(nrow(preO2R))
  preO2R_rand <- preO2R[row_r,]
  
  alpha <- 1/(preO2R_rand[,2]*preO2R_rand[,2]) # k in Wikipedia param
  beta <- preO2R_rand[,1]*preO2R_rand[,2]*preO2R_rand[,2]
  
  shape_death <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  scale_death <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  
  lprob_onset_to_discharge <- log(postO2D[,3]+preO2R_rand[,3]+1e-100)
  prob <- vcfr
  
  for(j in 1:length(vcfr)) {
    
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_discharge
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_death,scale=scale_death,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob <- tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                             pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=F)
                           +cfr*pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail = F) )
      
    }
    mtprob <- max(tprob)
    prob[j] <- log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
}


calcCFR_DeathsAndDischargesParallel<-function(data,postO2D,preO2R,vcfr,var_name_time_outcome) {
  
  row_r <- sample(nrow(preO2R))
  preO2R_rand <- preO2R[row_r,]
  
  alpha <- 1/(preO2R_rand[,2]*preO2R_rand[,2]) # k in Wikipedia param
  beta <- preO2R_rand[,1]*preO2R_rand[,2]*preO2R_rand[,2]
  
  shape_death <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  scale_death <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  
  lprob_onset_to_discharge <- log(postO2D[,3]+preO2R_rand[,3]+1e-100)
  prob <- vcfr
  
  prob <- foreach(j = 1:length(vcfr), .combine = rbind) %dopar% {
    
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_discharge
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_death,scale=scale_death,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob <- tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                             pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=F)
                           +cfr*pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail = F) )
      
    }
    mtprob <- max(tprob)
    log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
}

calcCFR_DeathsAndDischargesParallel_testBob<-function(data,postO2D,preO2R,vcfr,var_name_time_outcome) {
  
  row_r <- sample(nrow(preO2R))
  preO2R_rand <- preO2R[row_r,]
  
  alpha <- 1/(preO2R_rand[,2]*preO2R_rand[,2]) # k in Wikipedia param
  beta <- preO2R_rand[,1]*preO2R_rand[,2]*preO2R_rand[,2]
  
  shape_death <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  scale_death <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  
  lprob_onset_to_discharge <- rep(0,nrow(postO2D)) #log(postO2D[,3]*preO2R_rand[,3]+1e-100) # PRODUCT NOT SUM
  prob <- vcfr
  
  prob <- foreach(j = 1:length(vcfr), .combine = rbind) %dopar% {
    
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_discharge
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_death,scale=scale_death,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob <- tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                             pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=F)
                           +cfr*pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail = F) )
      
    }
    mtprob <- max(tprob)
    log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
}

calc_mOD_DeathsAndDischargesParallel<-function(data,postO2D,preO2R,vcfr,var_name_time_outcome) {
  
  row_r <- sample(nrow(preO2R))
  preO2R_rand <- preO2R[row_r,]
  
  alpha <- 1/(preO2R_rand[,2]*preO2R_rand[,2]) # k in Wikipedia param
  beta <- preO2R_rand[,1]*preO2R_rand[,2]*preO2R_rand[,2]
  
  shape_death <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  scale_death <- postO2D[,1]*postO2D[,2]*postO2D[,2]
  
  lprob_onset_to_discharge <- log(postO2D[,3]+preO2R_rand[,3]+1e-100)
  prob <- vcfr
  posterior_mOD<-rep(0,length(row_r))
  posterior_mOD_with_prior<-lprob_onset_to_discharge
  
  #prob <- foreach(j = 1:length(vcfr), .combine = rbind) %dopar% {
  for(j in 1:length(vcfr)) {
      
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_discharge
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death') 
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_death,scale=scale_death,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob <- tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                             pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=F)
                           +cfr*pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail = F) )
      
      posterior_mOD<-posterior_mOD+tprob
      posterior_mOD_with_prior<-posterior_mOD_with_prior+tprob
    }
    mtprob <- max(tprob)
    log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(cbind(posterior_mOD,posterior_mOD_with_prior))
}


calcCFR_DeathsAndDischargesOnce<-function(data,r=0) {
  
  var_name_time_outcome  <-  "time_onset_outcome_incl_report"
  
  m  <- 22.3
  s  <- 0.42
  alpha <- 1/(s*s) # k in Wikipedia param
  beta  <- 1/(1/(m*s*s)+r)
  
  shape_discharge <- alpha
  scale_discharge <- beta
  mean_d <- 22.2 #22.2  
  s_d <- 0.45  #0.45  
  scale_discharge <- 1/(1/(s_d*s_d*mean_d)+r)
  shape_discharge <- 1/(s_d*s_d)
  
  
  lprob_onset_to_death  <- log(prob_onset_to_death+1e-100)
  prob  <- vcfr
  
  for(j in 1:length(vcfr))
  {
    cfr <- vcfr[j]+1e-10
    tprob <- lprob_onset_to_death
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob <- tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                         pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob <- tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_discharge,scale=scale_discharge,lower.tail=T) - 
                                             pgamma(data[[var_name_time_outcome]][i],shape=shape_discharge,scale=scale_discharge,lower.tail=T)))
      else
        tprob <- tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=shape_discharge,scale=scale_discharge,lower.tail=F)
                           +cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )
      
    }
    mtprob <- max(tprob)
    prob[j] <- log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob <- max(prob)
  prob <- exp(prob-mprob)
  sprob <- sum(prob)
  prob <- prob/sprob
  return(prob)
}

#####################################################
# calculate distribution stats
#####################################################

calc_dist_stats <- function(prob,vs,vmean,plotm=F){
  
  # ncov marginal posterior distributions 
  matrix_prob <- matrix(prob, nrow = length(vmean),byrow = FALSE)
  marginal_mean_prob <- rowSums(matrix_prob)
  marginal_s_prob <- colSums(matrix_prob)
  
  if(plotm) {
    layout(matrix(c(1:2),nrow=1,ncol=2,byrow = TRUE))
    plot(marginal_mean_prob ~ vmean, type = "l")
    plot(marginal_s_prob ~ vs, type = "l")
  }
  layout(matrix(c(1),nrow=1,ncol=1,byrow = TRUE))
  
  prob_grid <- matrix(prob, ncol=length(vs), byrow = TRUE)
  
  # stats mean ncov
  mean_mean <- sum(vmean*marginal_mean_prob)
  mode_mean <- vmean[which(marginal_mean_prob == max(marginal_mean_prob))]
  df <- abs(cumsum(marginal_mean_prob) - 0.5)
  median_mean <- vmean[which(df == min(df))]
  df <- abs(cumsum(marginal_mean_prob) - 0.025)
  lower_bound_mean <- vmean[which(df == min(df))]
  df <- abs(cumsum(marginal_mean_prob)-0.975)
  upper_bound_mean <- vmean[which(df == min(df))]
  
  # stats s ncov
  mean_s <- sum(vs*marginal_s_prob)
  mode_s <- vs[which(marginal_s_prob == max(marginal_s_prob))]
  df <- abs(cumsum(marginal_s_prob) - 0.5)
  median_s <- vs[which(df == min(df))]
  df <- abs(cumsum(marginal_s_prob) - 0.025)
  lower_bound_s <- vs[which(df == min(df))]
  df <- abs(cumsum(marginal_s_prob)-0.975)
  upper_bound_s <- vs[which(df == min(df))]
  
  
  list(distribution = cbind(grid, 
                            prob),
       prob_grid = prob_grid, 
       marginal_mean = cbind(vmean, 
                             marginal_mean_prob),
       marginal_s = cbind(vs, 
                          marginal_s_prob),
       mean_stats = cbind(mean_mean, 
                          mode_mean, 
                          median_mean, 
                          lower_bound_mean, 
                          upper_bound_mean), 
       s_stats = cbind(mean_s, 
                       mode_s, 
                       median_s, 
                       lower_bound_s, 
                       upper_bound_s))
  
}

calc_dist_stats2 <- function(prob,grid,plotm=F){
  
  # ncov marginal posterior distributions 
  vmean<-unique(grid$mean)
  vs<-unique(grid$s)
  #matrix_prob <- matrix(prob, nrow = length(vmean),byrow = FALSE)
  #marginal_mean_prob <- rowSums(matrix_prob)
  marginal_mean_prob <- tapply(prob,grid$mean,sum)
  marginal_s_prob <- tapply(prob,grid$s,sum)
  
  if(plotm) {
    layout(matrix(c(1:2),nrow=1,ncol=2,byrow = TRUE))
    plot(marginal_mean_prob ~ vmean, type = "l")
    plot(marginal_s_prob ~ vs, type = "l")
  }
  
  #prob_grid <- matrix(prob, ncol=length(vs), byrow = TRUE)
  
  # stats mean ncov
  mean_mean <- sum(vmean*marginal_mean_prob)
  mode_mean <- vmean[which(marginal_mean_prob == max(marginal_mean_prob))]
  df <- abs(cumsum(marginal_mean_prob) - 0.5)
  median_mean <- vmean[which(df == min(df))]
  df <- abs(cumsum(marginal_mean_prob) - 0.025)
  lower_bound_mean <- vmean[which(df == min(df))]
  df <- abs(cumsum(marginal_mean_prob)-0.975)
  upper_bound_mean <- vmean[which(df == min(df))]
  
  # stats s ncov
  mean_s <- sum(vs*marginal_s_prob)
  mode_s <- vs[which(marginal_s_prob == max(marginal_s_prob))]
  df <- abs(cumsum(marginal_s_prob) - 0.5)
  median_s <- vs[which(df == min(df))]
  df <- abs(cumsum(marginal_s_prob) - 0.025)
  lower_bound_s <- vs[which(df == min(df))]
  df <- abs(cumsum(marginal_s_prob)-0.975)
  upper_bound_s <- vs[which(df == min(df))]
  
  
  list(distribution = cbind(grid, 
                            prob),
       #prob_grid = prob_grid, 
       marginal_mean = cbind(vmean, 
                             marginal_mean_prob),
       marginal_s = cbind(vs, 
                          marginal_s_prob),
       mean_stats = cbind(mean_mean, 
                          mode_mean, 
                          median_mean, 
                          lower_bound_mean, 
                          upper_bound_mean), 
       s_stats = cbind(mean_s, 
                       mode_s, 
                       median_s, 
                       lower_bound_s, 
                       upper_bound_s))
  
}


calc_cfr_stats <- function(prob,vcfr,plotm=F){
  
  mean_cfr <- sum(vcfr*prob)
  
  mode_cfr <- vcfr[which(prob == max(prob))]
  
  df <- abs(cumsum(prob) - 0.5)
  median_cfr <- vcfr[which(df == min(df))]
  
  df <- abs(cumsum(prob) - 0.025)
  lower_bound_cfr <- vcfr[which(df == min(df))]
  
  df <- abs(cumsum(prob)-0.975)
  upper_bound_cfr <- vcfr[which(df == min(df))]
  
  plot(prob ~ vcfr, type = "l")
  abline(v = mean_cfr, lty = 1)
  abline(v = lower_bound_cfr, lty = 2)
  abline(v = upper_bound_cfr, lty = 2)
  
  plot(cumsum(prob)~vcfr, type ="l")
  abline(h = 0.025, col="red")
  abline(h= 0.975, col="red")
  
  # cfr results 
  message("CFR mean and 95% CrI:")
  message(c(round(mean_cfr, 6), round(lower_bound_cfr,6), round(upper_bound_cfr, 6))) 
  
  message("CFR mean")
  message(mean_cfr) 
  
  message("CFR mode")
  message(mode_cfr) 
  
  message("CFR median")
  message(median_cfr)
  
  message("CFR lower bound")
  message(lower_bound_cfr) 
  
  message("CFR upper bound")
  message(upper_bound_cfr) 
  
  list(distribution = cbind(vcfr, prob), 
       stats = rbind(mean_cfr, mode_cfr, median_cfr, lower_bound_cfr, upper_bound_cfr))
  
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

##########################################################################
# IMPUTE WHO HAS RECOVERED - BY AGE>60, U60
##########################################################################
imputeRecoveriesAge<-function(exported_incl_missing_onset) {
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
      print(paste0("In ", curr_country," not enough eligible people to recover"))
      print(n_recs)
      print(length(inds))
      new_recs<-inds
    } else {
      inds_u60<-which(exported_incl_missing_onset$country==curr_country & exported_incl_missing_onset$outcome=='Other'
                  & exported_incl_missing_onset$age_over60==0 &
                    (is.na(exported_incl_missing_onset$recovered_y_n) | 
                       exported_incl_missing_onset$recovered_y_n=="n - implied" | 
                       exported_incl_missing_onset$recovered_y_n=="y - implied"))
      inds_o60<-which(exported_incl_missing_onset$country==curr_country & exported_incl_missing_onset$outcome=='Other'
                      & exported_incl_missing_onset$age_over60==1 &
                        (is.na(exported_incl_missing_onset$recovered_y_n) | 
                           exported_incl_missing_onset$recovered_y_n=="n - implied" | 
                           exported_incl_missing_onset$recovered_y_n=="y - implied"))
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

##########################################################################
# NON PARAMETRIC CFR ESTIMATATION - CASEFAT BY JAMIE GRIFFIN (ORIGINALLY A STATA MODULE)
##########################################################################

casefat = function(t, f){
  ###############################################################
  # Survivor function and variance for combined endpoint 
  
  c0 = survfit(Surv(t, factor(as.numeric(f>0), 0:1, labels=c("0", "1")))~1)
  si = which(c0$states!="1")
  S0 = c0$pstate[,si]
  V0 = (c0$std.err[,si])^2
  n = length(V0)
  if(V0[n]<1E-12){
    V0[n]=V0[n-1]
  }
  nrisk = c0$n.risk[,si]
  
  ###############################################################
  # CFR calculation
  
  c = survfit(Surv(t, factor(f, 0:2, labels=c("0", "1", "2")))~1)
  di = which(c$states=="1")
  ri = which(c$states=="2")
  
  # hazard contributions for each endpoint
  h1 = c$n.event[,di]/nrisk
  h2 = c$n.event[,ri]/nrisk
  
  # c$pstate is cumulative incidence function for each endpoint
  theta1 = max(c$pstate[,di])
  theta2 = max(c$pstate[,ri])
  
  cfr = theta1/(theta1+theta2)
  
  ###############################################################
  # Variances
  
  # Greenwood-like method
  M = diag(V0)
  for(j in 2:nrow(M)){
    for(k in 1:(j-1)){
      M[j, k] = V0[k]*S0[j]/S0[k]
      M[k, j] = M[j, k]	
    }
  }
  
  v1 = as.numeric(sum((S0)^2*h1/pmax(nrisk, 1)) + (h1 %*% M %*% h1))
  v2 = as.numeric(sum((S0)^2*h2/pmax(nrisk, 1)) + (h2 %*% M %*% h2))
  cov12 = as.numeric((h1 %*% M %*% h2))
  
  secfr = sqrt((theta2^2*v1 + theta1^2*v2 - 2*theta1*theta2*cov12))/(theta1+theta2)^2
  
  ###############################################################
  # logit scale for CI
  
  lc = log(cfr/(1-cfr))
  sel = sqrt(v1/theta1^2 + v2/theta2^2 - 2*cov12/(theta1*theta2))
  llc = lc-1.96*sel
  ulc = lc+1.96*sel
  lcfr = exp(llc)/(1+exp(llc))
  ucfr = exp(ulc)/(1+exp(ulc))
  
  ###############################################################
  # Two simple methods
  Nt = c$n
  Nd = sum(c$n.event[,di])
  Nr = sum(c$n.event[,ri])
  
  e1 = Nd/Nt
  see1 = sqrt(e1*(1-e1)/Nt)
  a1 = Nd+0.5
  b1 = Nt-Nd+0.5
  le1 = qbeta(0.025, shape1=a1, shape2=b1)
  ue1 = qbeta(0.975, shape1=a1, shape2=b1)
  
  e2 = Nd/(Nd+Nr)
  see2 = sqrt(e2*(1-e2)/(Nd+Nr))
  a2 = Nd+0.5
  b2 = Nr+0.5
  le2 = qbeta(0.025, shape1=a2, shape2=b2)
  ue2 = qbeta(0.975, shape1=a2, shape2=b2)
  
  ###############################################################
  
  # return(list(cfr=cfr, secfr=secfr, lcfr=lcfr, ucfr=ucfr, 
  # 			e1=e1, see1=see1, le1=le1, ue1=ue1, 
  # 			e2=e2, see2=see2, le2=le2, ue2=ue2))
  return(list(cfr=cfr, secfr=secfr, logit_cfr=lc, se_logit_cfr=sel, lcfr=lcfr, ucfr=ucfr, 
              e1=e1, see1=see1, le1=le1, ue1=ue1, 
              e2=e2, see2=see2, le2=le2, ue2=ue2))
  
}
