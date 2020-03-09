
# ------------------------------------------------------------------
# flat prior
cpp_flat_prior <- "SEXP logprior(std::vector<double> params) {
  return Rcpp::wrap(0.0);
}"

# ------------------------------------------------------------------
# log likelihood
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
  
  // sum log-likelihood over onset-to-outcome data
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
    {
      // get separate parts that make up likelihood
      double l1 = Rcpp::stats::pgamma_1((otp_interval + 1)*(beta_op + r), alpha_op, true, true);
      double l2 = Rcpp::stats::pgamma_1(otp_interval*(beta_op + r), alpha_op, true, true);
      double l3 = Rcpp::stats::pgamma_1((t_report[i] + 1)*(beta_op + r), alpha_op, true, true);
      
      // add to running log-likelihood
      ret += log(exp(l1) - exp(l2)) - l3;
    }
    
  }
  
  // catch underflow in likelihood
  if (!isfinite(ret)) {
    return Rcpp::wrap(-OVERFLO_DOUBLE);
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"