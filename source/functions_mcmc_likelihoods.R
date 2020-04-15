
# ------------------------------------------------------------------
# flat prior
cpp_flat_prior <- "SEXP logprior(std::vector<double> params) {
  return Rcpp::wrap(0.0);
}"

# ------------------------------------------------------------------
# prior on onset-to-death from age-stratified CFR analysis
cpp_otd_prior <- "SEXP logprior(std::vector<double> params) {
  
  // extract parameters
  double m_od = params[0];
  double s_od = params[1];
  
  // define composite/fixed parameters
  double mu = 2.721402;
  double sigma = 0.0306896;
  double grad = 6.414599;
  double shape1 = 23.42268;
  double shape2 = 31.10047;
  
  // evaluate log-probability
  double z = m_od - grad * s_od;
  double l1 = R::dlnorm(z, mu, sigma, true);
  double l2 = R::dbeta(s_od, shape1, shape2, true);
  
  // return as SEXP
  double ret = l1 + l2;
  return Rcpp::wrap(ret);
}"

# ------------------------------------------------------------------
# loglikelihood onset-to-recovery analysis
# note that the C++ function Rcpp::stats::pgamma_1(x*b, a, true, true) is
# equivalent to the R function pgamma(x, shape = a, rate = b, lower.tail = TRUE,
# log.p = TRUE)
cpp_loglike_otr <- "SEXP loglike(std::vector<double> params, std::vector<double> x) {
  
  // define fixed quantities
  const double OVERFLO_DOUBLE = DBL_MAX/100;
  
  // extract data
  int di = 0;
  int n = x.size()/5;
  std::vector<int> t_onset(n);
  std::vector<int> t_report(n);
  std::vector<int> t_outcome(n);
  std::vector<int> imputed(n);  // 0 = not imputed, 1 = imputed
  std::vector<double> growth_rate(n);
  for (int i = 0; i < n; ++i) {
    t_onset[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    t_report[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    t_outcome[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    imputed[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    growth_rate[i] = x[di++];
  }
  
  // extract parameters
  int pi = 0;
  double m_or = params[pi++];
  double s_or = params[pi++];
  double m_op = params[pi++];
  double s_op = params[pi++];
  
  // define fixed or derived parameters
  double alpha_or = 1/(s_or*s_or);
  double beta_or = 1/(m_or*s_or*s_or);
  double alpha_op = 1/(s_op*s_op);
  double beta_op = 1/(m_op*s_op*s_op);
  
  // zero likelihood if outside truncation range
  double gamma_tail = R::qgamma(0.95, alpha_or, 1/beta_or, true, false);
  if (gamma_tail > 100) {
    return Rcpp::wrap(-OVERFLO_DOUBLE);
  }
  
  // impute onset dates
  for (int i = 0; i < n; ++i) {
    if (imputed[i] == 1) {
      t_onset[i] -= params[pi++];
    }
  }
  
  // sum log-likelihood over data
  double ret = 0;
  for (int i = 0; i < n; ++i) {
    
    // get growth rate
    double r = growth_rate[i];
    
    // get intervals
    double otr_interval = t_outcome[i] - t_onset[i];
    double otp_interval = t_report[i] - t_onset[i];
    
    // onset-to-recovery
    {
      // get separate parts that make up likelihood
      double l1 = Rcpp::stats::pgamma_1((otr_interval + 1)*(beta_or + r), alpha_or, true, true);
      double l2 = Rcpp::stats::pgamma_1(otr_interval*(beta_or + r), alpha_or, true, true);
      double l3 = Rcpp::stats::pgamma_1((t_outcome[i] + 1)*(beta_or + r), alpha_or, true, true);
      
      // add to running log-likelihood
      ret += log(exp(l1) - exp(l2)) - l3;
    }
    
    // onset-to-report
    // check for -1, indicating missing report date
    if (t_report[i] != -1) {
      
      // get separate parts that make up likelihood
      double l1 = Rcpp::stats::pgamma_1((otp_interval + 1)*(beta_op + r), alpha_op, true, true);
      double l2 = Rcpp::stats::pgamma_1(otp_interval*(beta_op + r), alpha_op, true, true);
      double l3 = Rcpp::stats::pgamma_1((t_report[i] + 1)*(beta_op + r), alpha_op, true, true);
      
      // add to running log-likelihood
      ret += log(exp(l1) - exp(l2)) - l3;
    }
    
  }
  
  // catch underflow in likelihood
  if (!std::isfinite(ret)) {
    return Rcpp::wrap(-OVERFLO_DOUBLE);
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"

# ------------------------------------------------------------------
# loglikelihood international CFR analysis
# note that the C++ function Rcpp::stats::pgamma_1(x*b, a, true, true) is
# equivalent to the R function pgamma(x, shape = a, rate = b, lower.tail = TRUE,
# log.p = TRUE)
cpp_loglike_cfr <- "SEXP loglike(std::vector<double> params, std::vector<double> x) {
  
  // define fixed quantities
  const double OVERFLO_DOUBLE = DBL_MAX/100;
  
  // extract data
  int di = 0;
  int n = x.size()/7;
  std::vector<int> t_onset(n);
  std::vector<int> t_report(n);
  std::vector<int> t_outcome(n);
  std::vector<int> outcome_type(n);   // 1 = death, 2 = recovery, 3 = other
  std::vector<int> imputed(n);        // 0 = not imputed, 1 = imputed
  std::vector<int> recovery_only(n);  // 0 = all outcomes, 1 = recovery only
  std::vector<double> growth_rate(n);
  for (int i = 0; i < n; ++i) {
    t_onset[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    t_report[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    t_outcome[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    outcome_type[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    imputed[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    recovery_only[i] = x[di++];
  }
  for (int i = 0; i < n; ++i) {
    growth_rate[i] = x[di++];
  }
  
  // extract parameters
  int pi = 0;
  double m_od = params[pi++];
  double s_od = params[pi++];
  double m_or = params[pi++];
  double s_or = params[pi++];
  double m_op = params[pi++];
  double s_op = params[pi++];
  double cfr = params[pi++];
  double p = params[pi++];
  
  // define fixed or derived parameters
  double alpha_od = 1/(s_od*s_od);
  double beta_od = 1/(m_od*s_od*s_od);
  double alpha_or = 1/(s_or*s_or);
  double beta_or = 1/(m_or*s_or*s_or);
  double alpha_op = 1/(s_op*s_op);
  double beta_op = 1/(m_op*s_op*s_op);
  
  // zero likelihood if outside truncation range
  double gamma_tail_od = R::qgamma(0.95, alpha_od, 1/beta_od, true, false);
  double gamma_tail_or = R::qgamma(0.95, alpha_or, 1/beta_or, true, false);
  double gamma_tail_op = R::qgamma(0.95, alpha_op, 1/beta_op, true, false);
  if ((gamma_tail_od > 100) || (gamma_tail_or > 100) || (gamma_tail_op > 100)) {
    return Rcpp::wrap(-OVERFLO_DOUBLE);
  }
  
  // impute onset dates
  for (int i = 0; i < n; ++i) {
    if (imputed[i] == 1) {
      t_onset[i] -= params[pi++];
    }
  }
  
  // sum log-likelihood over data
  double ret = 0.0;
  for (int i = 0; i < n; ++i) {
    
    // get growth rate
    double r = growth_rate[i];
    
    // get intervals
    double oto_interval = t_outcome[i] - t_onset[i];
    double otp_interval = t_report[i] - t_onset[i];
    
    // if recovery-only data then account for growth
    if (recovery_only[i] == 1) {
      
      // get separate parts that make up likelihood
      double l1 = Rcpp::stats::pgamma_1((oto_interval + 1)*(beta_or + r), alpha_or, true, true);
      double l2 = Rcpp::stats::pgamma_1(oto_interval*(beta_or + r), alpha_or, true, true);
      double l3 = Rcpp::stats::pgamma_1((t_outcome[i] + 1)*(beta_or + r), alpha_or, true, true);
      
      // add to running log-likelihood
      ret += log(exp(l1) - exp(l2)) - l3;
      
    } else {  // if all-outcomes data
      
      // onset-to-outcome
      if (outcome_type[i] == 1) {  // death
        
        // get separate parts that make up likelihood
        double l1 = Rcpp::stats::pgamma_1((oto_interval + 1)*beta_od, alpha_od, true, true);
        double l2 = Rcpp::stats::pgamma_1(oto_interval*beta_od, alpha_od, true, true);
        
        // add to running log-likelihood
        ret += log(cfr) + log(exp(l1) - exp(l2));
        
      } else if (outcome_type[i] == 2) {  // recovery
        
        // get separate parts that make up likelihood
        double l1 = Rcpp::stats::pgamma_1((oto_interval + 1)*beta_or, alpha_or, true, true);
        double l2 = Rcpp::stats::pgamma_1(oto_interval*beta_or, alpha_or, true, true);
        
        // add to running log-likelihood
        ret += log(p) + log(1.0 - cfr) + log(exp(l1) - exp(l2));
        
      } else if (outcome_type[i] == 3) {  // other
        
        // get separate parts that make up likelihood
        double l1 = Rcpp::stats::pgamma_1(oto_interval*beta_od, alpha_od, false, true);
        double l2 = Rcpp::stats::pgamma_1(oto_interval*beta_or, alpha_or, false, true);
        
        // add to running log-likelihood
        ret += log(cfr*exp(l1) + p*(1.0 - cfr)*exp(l2) + (1.0 - p)*(1.0 - cfr));
      }
      
      // onset-to-report
      // check for -1, indicating missing report date
      if (t_report[i] != -1) {
        
        // get separate parts that make up likelihood
        double l1 = Rcpp::stats::pgamma_1((otp_interval + 1)*(beta_op + r), alpha_op, true, true);
        double l2 = Rcpp::stats::pgamma_1(otp_interval*(beta_op + r), alpha_op, true, true);
        double l3 = Rcpp::stats::pgamma_1((t_report[i] + 1)*(beta_op + r), alpha_op, true, true);
        
        // add to running log-likelihood
        ret += log(exp(l1) - exp(l2)) - l3;
      }
    }
    
  }  // end loop over international data
  
  // catch underflow in likelihood
  if (!std::isfinite(ret)) {
    return Rcpp::wrap(-OVERFLO_DOUBLE);
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}"