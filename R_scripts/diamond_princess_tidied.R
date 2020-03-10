##########################################################
# ESTIMATE IFR ON THE DIAMOND PRINCESS
##########################################################

library(binom)
logit<-function(x) log(x/(1-x))
invlogit<-function(x) exp(x)/(1+exp(x))

tests<-read.csv("data/test_dates_dp.csv")

date_format <- "%d/%m/%Y"
tests$date_report<-as.numeric(as.Date(as.character(tests$date_report),date_format))
tests<-tests[which(!is.na(tests$date_report)),]
tot_cases<-max(tests$cum_positive,na.rm=T)
tests$cum_pc_pos<-tests$cum_positive/tot_cases
tests$logit_cum_pos<-logit(tests$cum_pc_pos-0.5/tot_cases)
tests$time_elapsed<-tests$date_report-min(tests$date_report)

### FITTED ONSET TO DEATH PARAMETERS FROM MAIN ANALYSIS
mean_OD<-17.8324  #18.8
s_OD<-0.4226  #0.45
alpha <- 1/(s_OD^2) 
beta <- 1/(1/(mean_OD*s_OD^2))


### ANALYSIS
today<-as.Date("2020-03-05")
## generate weight data points from each day.
tests$p<-(tests$num_positive+1e-10)/tests$num_tested  ##   ## sums to 657, not the total 706 cases but the most we have.
tests$wgt<- 1/(tests$p*(1-tests$p)/tests$num_tested)  ## 1/variance of proportion.

cx<-lm(tests$logit_cum_pos ~tests$time_elapsed, weights=tests$wgt)
a<-cx$coefficients[1]
b<-cx$coefficients[2]
tests$predict_logit_pos<-a +tests$time_elapsed*b
tests$predict_pc_pos<-invlogit(tests$predict_logit_pos)
tests$predict_pc_pos<-tests$predict_pc_pos/max(tests$predict_pc_pos)

x<-seq(0,30,0.1)
predict_pos<-invlogit(a +x*b)
plot(tests$time_elapsed,tests$cum_pc_pos,xlab="days since 5th Feb",ylab="proportion positive",pch=19)
lines(x,predict_pos)

tests$predict_new<-tests$predict_pc_pos
tests$predict_new[2:nrow(tests)]<-tests$predict_pc_pos[2:nrow(tests)]-tests$predict_pc_pos[1:(nrow(tests)-1)]

### estimate proportion of deaths we expect to have occurred by now allowing for onset times
tests$elapsed<-as.numeric(as.Date("2020-03-05"))-tests$date_report
tests$prop_deaths_by_now<-pgamma(tests$elapsed,shape=alpha,scale=beta,lower.tail = T)
tests$prop_deaths_by_now_weighted<-tests$predict_new*tests$prop_deaths_by_now

prop_deaths_by_now<-sum(tests$prop_deaths_by_now_weighted)
print(prop_deaths_by_now)

### 20.6 total deaths are predicted, obtained by applying age-specific IFR estimates from China to the age distribution  
# of passengers on the diamond princess
prop_deaths_by_now*20.6  ## number of deaths expected by 5th March using the age specific IFR.

## 95% CI around deaths currently observed among cases on the ship.
binom.confint(x=7,n=706 , method='exact')
