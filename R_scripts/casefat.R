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
