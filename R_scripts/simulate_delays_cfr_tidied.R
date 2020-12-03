########################################################
# SIMULATE ONSET - DEATH AND ONSET - RECOVERY WITH EXPONENTIAL GROWTH
# FIT TO OUTCOMES OBSERVED MID-EPIDEMIC WITH AND WITHOUT CORRECTED FOR GROWTH
########################################################

# parameters
cfr<-0.05
## O2D
mean_death<-18.8  # mean days onset to death.
s_death<-0.45
shape_death <- 1/(s_death^2)
scale_death=mean_death*s_death^2
## O2R
mean_discharge<-25.8  # mean days onset to rec
s_discharge<-0.36
shape_discharge <- 1/(s_discharge^2) # k in Wikipedia param
scale_discharge=mean_discharge*s_discharge^2

r<-0.14  # growth rate in incidence
init_inc<-2
max_time<-55
inc<-init_inc*exp(r*1:max_time)
inc<-round(inc)

# doubling time
print(log(2)/r)

cfr<-0.05
N0<-sum(inc)
sim<-data.frame(id=1:N0)
sim$time_onset<-rep(1:max_time,inc)
sim$died<-runif(N0)<cfr
## draw death times
sim$time_outcome<-NA
sim$time_outcome[sim$died]<-sim$time_onset[sim$died]+rgamma(length(which(sim$died)),shape=shape_death,scale=scale_death)
## draw discharge times
sim$time_outcome[!sim$died]<-sim$time_onset[!sim$died]+rgamma(length(which(!sim$died)),shape=shape_discharge,scale=scale_discharge)
## true onset-outcome
sim$true_O2outcome<-sim$time_outcome-sim$time_onset

### now fit the simulated onset-outcome (for recoveries, but can also be done for deaths).
source("source/utils_CFR.R")
# define  grid
vmean <- seq(10, 100, 0.1)
vs <- seq(0.2, 0.8, 0.01)
grid <- as.data.frame(expand.grid(mean = vmean, s = vs))

## OBSERVED ONSET - OUTCOME DISTRIBUTION MID EPIDEMIC ######
### Date of observation - days into the epidemic
day_obs<-55 # day at which we observe epidemic
sim_trunc<-sim[which(sim$time_onset<day_obs),]
### whose outcome would be unobserved?
ind<-which(sim_trunc$time_outcome>day_obs)
sim_trunc$time_outcome[ind]<-NA  ## set to missing
sim_trunc$observed_O2outcome<-sim_trunc$true_O2outcome
sim_trunc$observed_O2outcome[ind]<-day_obs  ## set to day of censoring
sim_trunc$censored<-0  ## identify which have not been censored
sim_trunc$censored[ind]<-1

# refit to the observed values not allowing for growth
inds<-which(sim_trunc$censored==0 & !sim_trunc$died)
cens<-rep(1000,length(inds))  ## to switch off cens, set it to be a big number.
prob_O2R_m_s_obs_noR<-calcOnsetToOutcome(sim_trunc[inds,],r=0,var_name_outcome="observed_O2outcome")
save_dist_stats("output/prob_O2R_m_s_obs_noR",prob_O2R_m_s_obs_noR,plot=T)

res<-read.csv("output/prob_O2R_m_s_obs_noR_stats.csv")
mean_mean_fit_obs_noR<-res$V1[res$X=="mean_mean"]
mean_s_fit_obs_noR<-res$V1[res$X=="mean_s"]
shape_fit_obs_noR <- 1/(mean_s_fit_obs_noR^2)
scale_fit_obs_noR=mean_mean_fit_obs_noR*mean_s_fit_obs_noR^2

# refit to the observed values WITH GROWTH
inds<-which(sim_trunc$censored==0 & !sim_trunc$died)
cens<-rep(1000,length(inds))  ## to switch off cens, set it to be a big number.
prob_O2R_m_s_obs_withR<-calcOnsetToOutcome(sim_trunc[inds,],r=0.14,var_name_outcome="observed_O2outcome")
save_dist_stats("output/prob_O2R_m_s_obs_withR",prob_O2R_m_s_obs_withR,plot=T)

res<-read.csv("prob_O2R_m_s_obs_withR_stats.csv")
mean_mean_fit_obs_withR<-res$V1[res$X=="mean_mean"]
mean_s_fit_obs_withR<-res$V1[res$X=="mean_s"]
shape_fit_obs_withR <- 1/(mean_s_fit_obs_withR^2)
scale_fit_obs_withR=mean_mean_fit_obs_withR*mean_s_fit_obs_withR^2


#### COMPARE OUTPUTS.
layout(matrix(c(1),nrow=1,ncol=1, byrow=T))
plot(seq(0,60,1),length(which(!sim_trunc$died))*dgamma(seq(0,60,1),shape=shape_fit_obs_noR,scale=scale_fit_obs_noR),
     type="l",lwd=2, col="red",  main="", xlab="days since onset",ylab="frequency")
hist(sim_trunc$true_O2outcome[which(!sim_trunc$died)],breaks=seq(0,ceiling(max(sim_trunc$true_O2outcome[which(!sim_trunc$died)])),1),add=T)
inds<-which(sim_trunc$censored==0 & !sim_trunc$died)
hist(sim_trunc$observed_O2outcome[inds],breaks=seq(0,ceiling(max(sim_trunc$observed_O2outcome[inds])),1),add=T,col="red")
lines(seq(0,60,1),length(which(!sim_trunc$died))*dgamma(seq(0,60,1),shape=shape_fit_obs_withR,scale=scale_fit_obs_withR),
      col="blue",lwd=2)

