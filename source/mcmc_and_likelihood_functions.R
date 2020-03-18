lL_plane <- function(params, x){
  ### Unpack parameters ########################################################
  m_od <- params[1]
  s_od <- params[2]
  
  alpha <- 1 / (s_od * s_od)
  beta <- 1 / (m_od * s_od * s_od)
  
  # Maximum day of data
  maxday <- params[3]
  # Proportion missed Wuhan
  z <- params[4]
  # Growth rate 
  r <- params[5]
  # Detection window
  D <- params[6]
  # Scaling
  y <- params[7]
  # CFR for each age band
  cfr <- params[8:length(params)]
  ##############################################################################
  cfr_adjustment <- cfr[1] * c(1, cfr[2:9])
  cfr <- c(cfr_adjustment[2:9], cfr_adjustment[1])
  
  ### Unpack data ##############################################################
  # Output type
  output <- x[1]
  # Number of age bands
  age_bands <- x[2]
  # Death observation cut off
  death_obv_cut <- x[3]
  # "Gold standard" ascertainment
  bestz <- x[4]
  # Total observed deaths Wuhan
  deaths_x <- x[5]
  deaths_n <- x[6]
  # Positive on repatriation
  repat_x <- x[7]
  # Tested on repatriation
  repat_n <- x[8]
  # WHO Observed Cases & Age
  WHO_observed_deaths <- x[9]
  WHO_observed_cases <- x[10]
  # Observed deaths by age
  observed_deaths <- x[11:(11 + age_bands - 1)]
  dm <- matrix(x[(11 + age_bands):length(x)], ncol = 5)
  # Cases
  cases_subset <- dm[,1]
  # Age band
  ageband <- dm[,2]
  # Time
  t <- dm[,3]
  # Indicator if wuhan or not
  wuhan <- dm[,4]
  # NiCi
  NiCi <- dm[,5]
  ##############################################################################
  
  ### Model prediction #########################################################
  # Adjust cases for ascertainment
  # This adjustment = 1 for the "gold standard" age group outside of Wuhan
  zNiCi <- ifelse(wuhan, z, 1 / bestz) * NiCi
  cases_subset_adjusted <- cases_subset * zNiCi
  
  # Death gamma density matrix (distributed our death over future days given the O2D)
  dmat <- estimate_deaths_mat(0:71, 72, alpha, beta)
  
  # Multiply cases by dmat for (i = wuhan and outside) and (j = agebands)
  ed <- list()
  k <- 1
  for(i in 0:1){
    for(j in 1:age_bands){
      index <- ageband == j & wuhan == i
      d2 <- Map("*", cases_subset_adjusted[index] * cfr[j], dmat) 
      ed[[k]] <- Reduce("+", d2)
      k <- k + 1
    }
  }
  pred_deaths <- unlist(ed) # pred deaths is all the deaths that should've occurred in the time window
  pred_deaths[t < death_obv_cut] <- 0 # Adjust for not observing deaths before cutoff
  
  # Crude CFR 
  predicted_deaths_all <- sum(cases_subset_adjusted * cfr[ageband]) # all deaths that should occur based on the cases we've estimated
  observed_cases <- sum(cases_subset)
  crude_cfr <- min(1, predicted_deaths_all/observed_cases)
  
  # Prevalence
  cases_subset_adjusted <- cases_subset_adjusted / y
  cases_wuhan_31st <- sum(cases_subset_adjusted[t == 60 & wuhan == 1])
  wuhan_pop <- 11081000
  prevalence <- min(1, ((cases_wuhan_31st / wuhan_pop) * (1 - exp(-r * D))) / r)
  
  # Sum cases and deaths by age group
  deaths_by_age <- tapply(pred_deaths, ageband, sum)
  prop_wuhan <- sum(pred_deaths * wuhan) / sum(pred_deaths)
  
  ##############################################################################
  
  # Output our predictions if specified
  if(output == 1){
    return(list(
      age_deaths = deaths_by_age,
      prop_wuhan = prop_wuhan,
      prevalence = prevalence,
      crude_cfr = crude_cfr
    ))
  }
  
  ### Likelihood calculation ###################################################
  ll <- sum(dpois(observed_deaths, deaths_by_age, log = TRUE)) +  
    dbinom(deaths_x, deaths_n, prop_wuhan, log = TRUE) +
    dbinom(repat_x, repat_n, prevalence, log = TRUE) + 
    dbinom(WHO_observed_deaths, WHO_observed_cases, crude_cfr, log = TRUE)
  return(ll)
  ##############################################################################
}

# Predict expected distribution of deaths in the future from a given timeperiod and O2D distribution
dp <- function(t, maxday, alpha, beta){
  c(rep(0, t), dgamma(1:(maxday - t), shape = alpha, rate = beta))
}
# Estimate the deaths matrix given the expected distribution of deaths from an O2D distribution 
estimate_deaths_mat <- function(times, maxday, alpha, beta){
  purrr::map(times, dp, maxday = maxday, alpha= alpha, beta = beta)
}


# Prior params from fit to 26 deaths dataset
lP_plane <- function(params){
  m_od <- params[1]
  s_od <- params[2]
  cfr <- params[8]
  
  alpha <- 1
  beta <- 17.6
  k <- 10
  shape1 <- 23.05
  shape2 <- 31.57
  
  if(is.nan(log(beta^alpha/gamma(alpha) * 1/(m_od - k)^(alpha + 1) * exp(-beta/(m_od - k))) )){
    message("m_od", m_od)
    message("s_od", s_od)
  }
  
  log(beta^alpha/gamma(alpha) * 1/(m_od - k)^(alpha + 1) * exp(-beta/(m_od - k))) +
    dbeta(s_od, shape1, shape2, log = TRUE) +
    8*log(cfr)
  
}
