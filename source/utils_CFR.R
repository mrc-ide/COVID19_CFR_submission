#####################################################
# define functions 
# loglikelihood onset-to death as function of mean and s (cv = sd/mean)
# corrected for epidemic growth
#####################################################
calcOnsetToOutcome <- function(data,r=0,var_name_outcome = "onset.to.death"){
  
  alpha <- 1/(grid[,2]*grid[,2]) # k in Wikipedia param
  beta=1/(1/(grid[,1]*grid[,2]*grid[,2])+r)
  prob=alpha
  lprob=0

  for(i in 1:nrow(data))
  {
    lprob <- lprob+log(1e-100+(pgamma(data[[var_name_outcome]][i]+1, shape = alpha, scale = beta, lower.tail = TRUE)
                             -pgamma(data[[var_name_outcome]][i], shape = alpha, scale = beta,lower.tail = TRUE))
                     /pgamma(cens,shape = alpha, scale = beta, lower.tail = TRUE))
  }
  mlprob=max(lprob)
  prob=exp(lprob-mlprob)
  sprob=sum(prob)
  prob=prob/sprob
  return(prob)
  
}

#####################################################
# define function
# para[1] = mean, param[2]= s, param[3]= cfr
#####################################################
calcCFR_DeathsOnly<-function(data,postO2D,vcfr,var_name_time_outcome) {
  
  alpha <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  beta=postO2D[,1]*postO2D[,2]*postO2D[,2]
  lprob_onset_to_death=log(postO2D[,3]+1e-100)
  prob=vcfr
  
  for(j in 1:length(vcfr))
  {
    cfr=vcfr[j]+1e-10
    tprob=lprob_onset_to_death
    for(i in 1:nrow(data))
      {
      if(data$outcome[i]=='Death')
          tprob=tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
          pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
          tprob=tprob+log(1e-100+(1-cfr)+cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )

    }
    mtprob=max(tprob)
    prob[j]=log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob=max(prob)
  prob=exp(prob-mprob)
  sprob=sum(prob)
  prob=prob/sprob
  return(prob)
}

calcCFR_DeathsOnlyParallel<-function(data,postO2D,vcfr,var_name_time_outcome) {
  
  alpha <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  beta=postO2D[,1]*postO2D[,2]*postO2D[,2]
  lprob_onset_to_death=log(postO2D[,3]+1e-100)
  prob=vcfr
  
  prob=foreach(j = 1:length(vcfr), .combine = rbind) %dopar% {

    cfr=vcfr[j]+1e-10
    tprob=lprob_onset_to_death
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob=tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                      pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob=tprob+log(1e-100+(1-cfr)+cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )
      
    }
    mtprob=max(tprob)
    log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob=max(prob)
  prob=exp(prob-mprob)
  sprob=sum(prob)
  prob=prob/sprob
  return(prob)
}

calcCFR_DeathsAndDischarges<-function(data,postO2D,preO2R,vcfr,var_name_time_outcome) {
  
  
  row_r=sample(nrow(preO2R))
  preO2R_rand=preO2R[row_r,]
  
  alpha <- 1/(preO2R_rand[,2]*preO2R_rand[,2]) # k in Wikipedia param
  beta=preO2R_rand[,1]*preO2R_rand[,2]*preO2R_rand[,2]
  
  shape_death <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  scale_death=postO2D[,1]*postO2D[,2]*postO2D[,2]
  
  lprob_onset_to_discharge=log(postO2D[,3]+preO2R_rand[,3]+1e-100)
  prob=vcfr
  
  for(j in 1:length(vcfr)) {
    
    cfr=vcfr[j]+1e-10
    tprob=lprob_onset_to_discharge
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob=tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_death,scale=scale_death,lower.tail=T) - 
                                      pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob=tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                          pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob=tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=F)
                        +cfr*pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail = F) )
      
    }
    mtprob=max(tprob)
    prob[j]=log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob=max(prob)
  prob=exp(prob-mprob)
  sprob=sum(prob)
  prob=prob/sprob
  return(prob)
}


calcCFR_DeathsAndDischargesParallel<-function(data,postO2D,preO2R,vcfr,var_name_time_outcome) {
  

  row_r=sample(nrow(preO2R))
  preO2R_rand=preO2R[row_r,]
  
  alpha <- 1/(preO2R_rand[,2]*preO2R_rand[,2]) # k in Wikipedia param
  beta=preO2R_rand[,1]*preO2R_rand[,2]*preO2R_rand[,2]
  
  shape_death <- 1/(postO2D[,2]*postO2D[,2]) # k in Wikipedia param
  scale_death=postO2D[,1]*postO2D[,2]*postO2D[,2]
  
  lprob_onset_to_discharge=log(postO2D[,3]+preO2R_rand[,3]+1e-100)
  prob=vcfr
  
  prob=foreach(j = 1:length(vcfr), .combine = rbind) %dopar% {

    cfr=vcfr[j]+1e-10
    tprob=lprob_onset_to_discharge
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob=tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_death,scale=scale_death,lower.tail=T) - 
                                      pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob=tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                          pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else
        tprob=tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=F)
                        +cfr*pgamma(data[[var_name_time_outcome]][i],shape=shape_death,scale=scale_death,lower.tail = F) )
      
    }
    mtprob=max(tprob)
    log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob=max(prob)
  prob=exp(prob-mprob)
  sprob=sum(prob)
  prob=prob/sprob
  return(prob)
}


calcCFR_DeathsAndDischargesOnce<-function(data,r=0) {
  
  var_name_time_outcome = "time_onset_outcome_incl_report"
  
  m=22.3
  s=0.42
  alpha <- 1/(s*s) # k in Wikipedia param
  beta=1/(1/(m*s*s)+r)
  
  shape_discharge<-alpha
  scale_discharge<-beta
  mean_d<- 22.2 #22.2  
  s_d<- 0.45  #0.45  
  scale_discharge<-1/(1/(s_d*s_d*mean_d)+r)
  shape_discharge<-1/(s_d*s_d)
  
  
  lprob_onset_to_death=log(prob_onset_to_death+1e-100)
  prob=vcfr
  
  for(j in 1:length(vcfr))
  {
    cfr=vcfr[j]+1e-10
    tprob=lprob_onset_to_death
    for(i in 1:nrow(data))
    {
      if(data$outcome[i]=='Death')
        tprob=tprob+log(1e-100+cfr*(pgamma(data[[var_name_time_outcome]][i]+1,shape=alpha,scale=beta,lower.tail=T) - 
                                      pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail=T)))
      else if(data$outcome[i]=='Discharge')
        tprob=tprob+log(1e-100+(1-cfr)*(pgamma(data[[var_name_time_outcome]][i]+1,shape=shape_discharge,scale=scale_discharge,lower.tail=T) - 
                                          pgamma(data[[var_name_time_outcome]][i],shape=shape_discharge,scale=scale_discharge,lower.tail=T)))
      else
        tprob=tprob+log(1e-100+(1-cfr)*pgamma(data[[var_name_time_outcome]][i],shape=shape_discharge,scale=scale_discharge,lower.tail=F)
                        +cfr*pgamma(data[[var_name_time_outcome]][i],shape=alpha,scale=beta,lower.tail = F) )
      
    }
    mtprob=max(tprob)
    prob[j]=log(sum(exp(tprob-mtprob)))+mtprob
  }
  mprob=max(prob)
  prob=exp(prob-mprob)
  sprob=sum(prob)
  prob=prob/sprob
  return(prob)
}





#################################
#### save distribution stats ####
#################################

save_dist_stats <- function(nfile,prob,plotm=F){

  # ncov marginal posterior distributions 
  matrix_prob <- matrix(prob, nrow = length(vmean),byrow = FALSE)
  marginal_mean_prob <- rowSums(matrix_prob)
  marginal_s_prob <- colSums(matrix_prob)
  
  if(plotm) {
  plot(marginal_mean_prob ~ vmean, type = "l")
  }
  
  # output statistics ncov marginal posterior distributions
  write.csv(cbind(grid, prob), file=paste0(nfile,"_distributions.csv"))
  prob_grid <- matrix(prob, ncol=length(vs), byrow = TRUE)
  write.csv(prob_grid, file=paste0(nfile,"_grid.csv"))
  write.csv(cbind(vmean, marginal_mean_prob), file=paste0(nfile,"_marginal_distribution_mean.csv"))
  write.csv(cbind(vs, marginal_s_prob), file=paste0(nfile,"_marginal_distribution_s.csv"))
  
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
  
  write.csv(rbind(mean_mean,
                  mode_mean,
                  median_mean,
                  lower_bound_mean,
                  upper_bound_mean,
                  mean_s,
                  mode_s,
                  median_s,
                  lower_bound_s,
                  upper_bound_s),
              file=paste0(nfile,"_stats.csv"))
}

save_cfr_stats <- function(nfile,prob,vcfr,plotm=F){
  

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
  
  ################################################################################
  # cfr results for deaths only
  ################################################################################
  
  print("CFR mean and 95% CrI:")
  print(c(round(mean_cfr, 6), round(lower_bound_cfr,6), round(upper_bound_cfr, 6))) 
  
  print("CFR mean")
  print(mean_cfr) 
  
  print("CFR mode")
  print(mode_cfr) 
  
  print("CFR median")
  print(median_cfr)
  
  print("CFR lower bound")
  print(lower_bound_cfr) 
  
  print("CFR upper bound")
  print(upper_bound_cfr) 
  
  write.csv(cbind(vcfr, prob), file = paste0(nfile,"_distributions.csv"))
  write.csv(rbind(mean_cfr, mode_cfr, median_cfr, lower_bound_cfr, upper_bound_cfr),file=paste0(nfile,"_stats.csv"))
}
