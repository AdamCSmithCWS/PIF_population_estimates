
// function compute_E allows for efficient vectorization of the likelihood
// without saving n_counts*4 parameters (E, ste_)

functions {
  vector compute_E(
      vector log_beta,
      array[] int route,
      array[] int week,
      vector yeareffect_flat,
      array[] int strat_year_idx,
      matrix log_mean_rel_abund,
      real RHO,
      real sd_obs, vector obs_raw, array[] int observer,
      int use_pois, real sdnoise, vector noise_raw) {
    int n = size(route);
    vector[n] noise_effect;

    if (use_pois) {
      noise_effect = sdnoise * noise_raw;
    } else {
      noise_effect = rep_vector(0, n);
    }

    return log_beta[route]
           + yeareffect_flat[strat_year_idx]
           + sd_obs * obs_raw[observer]
           + RHO*log_mean_rel_abund[route,week]
           + noise_effect;
  }
}




// data here are the raw BBS counts
// adjustmen factors from PIF population estimates
// log mean eBird relative abundance values from within BBS route-path buffers
data {
  int<lower=1> n_routes; // number of BBS routes in species seasonal range
  int<lower=1> n_obs; // number of BBS observers in species seasonal range
  int<lower=1> n_weeks; // number of weeks of observations and relative abundance values
  int<lower=1> n_counts; // number of counts from those routes over last 10-years
  int<lower=2> n_years;
  int<lower=1> ebird_year; // index of year for the eBird relative abundance surface
  array[ebird_year-1] int<lower=1> yrev; // reverse year vector (ebird_year-1):1 (necessary for Stan indexing)
  int<lower=1> n_strata; // number of unique values of doy (max(doy))

  array[n_counts] int<lower=0> count; // Raw BBS counts
  array[n_counts] int<lower=1> route; // route indicators for counts
  array[n_counts] int<lower=1> observer; // observer indicators for counts
  array[n_counts] int<lower=1> year; // year indicators for counts
  array[n_counts] int<lower=1,upper=n_weeks> week; // week indicator for counts
  array[n_counts] int<lower=1> strata; // strata indicators for route (to support varying doy effect and trend)
  // predictor
  matrix[n_routes,n_weeks] log_mean_rel_abund; // log-transformed relative abund within BBS route buffers by week



// spatial varying time series components
  array[n_years] int<lower=0> y_2020; //indicators for 2020 = 0 if 2020 and missing if fixed_year
  vector[n_strata] zero_gammas;   // a vector of zeros to fill fixed log_beta values for fixed_year
  //data for spatial iCAR among strata
  int<lower=1> n_edges;
  array [n_edges] int<lower=1, upper=n_strata> node1;  // node1[i] adjacent to node2[i]
  array [n_edges] int<lower=1, upper=n_strata> node2;  // and node1[i] < node2[i]
  real<lower=0> scaling_factor;


// spline components for doy
  int<lower=1> n_knots_doy;  // number of knots in the basis function for year
  matrix[n_doy, n_knots_doy] doy_basis; // basis function matrix


// PIF adjustment factors
  real c_p; // mean of pair correction values
  real sd_c_p; // sd of pair correction values

  real c_t; // mean time of day correction value or availability correction
  real sd_c_t; //sd of time of day correction value or sd of availability correction
  int<lower=0,upper=1> use_availability; // indicator for calculating availability
          // if use_availability == 1, then use availability with 1/rnorm(c_d,sd_c_d)
          // if use_availability == 0, then use availability with rnorm(c_d,sd_c_d)


  real c_d; // mean distance adjustment value
  real sd_c_d; // sd of distance adjustment value
  real c_d_lower; // lower bound of distance adjustment value
  real c_d_upper; // upper bound of distance adjustment value
  int<lower=0,upper=1> use_edr; // indicator for calculating distance with edr
          // if use_edr == 1, then use edr with c_d and sd_c_d
          // if use_edr == 0, then use uniform dist between c_d_lower and c_d_upper

  int<lower=0,upper=1> use_pois; //indicator if count variation should be based on over-dispersed Poisson (if ==1) or Negative binomial (if == 0)
  int<lower=0,upper=1> use_pair; //indicator if pair-adjustment should be used
  int<lower=0,upper=1> use_t; //indicator if route variation should be based on heavy tailed t-distribution instead of a normal
  int<lower=0,upper=1> est_rho; //indicator if model should estimate slope of log-log relationship or leave it fixed at 1.0
  int<lower=0,upper=1> use_ppc; //indicator if model should generate posterior predictions
  int<lower=0,upper=1> use_likelihood; //indicator if model should use data if == 0, then only prior predictions

}

transformed data {


      int<lower=1> n_years_m1 = n_years - 1 ;

      // supporting the to_vector(yeareffect[s,t])
      // column-major order so strat_year_idx = c(1:nstrata year1, (nstrata+1):nstrata*year2,)
     array[n_counts] int<lower=1, upper=n_strata*n_years> strat_year_idx;

      for (i in 1:n_counts) {
        strat_year_idx[i] = (year[i] - 1) * n_strata + strata[i];
        }



}
//
parameters {

  real log_BETA; // mean route-level callibration log scale
  real<lower=0> sd_log_beta; // sd of calibration among routes
  real RHO_raw; // slope of the log-log relationship - used as a check, not for inference
  sum_to_zero_vector[n_routes] log_beta_raw; // route-specific calibrations
  sum_to_zero_vector[n_obs] obs_raw; // route-specific calibrations
  real<lower=0> sd_obs; // sd of calibration among routes
  real<lower=3> nu;
  real<lower=3> nu_obs;


  vector[n_counts*use_pois] noise_raw; // over-dispersion if use_pois == 1
  real<lower=0> sdnoise;    // sd of over-dispersion, if use_pois == 1


  real<lower=0> sd_gamma_raw;    // sd of annual changes among strata
  real<lower=0> sd_GAMMA;    // sd of overall annual changes

  vector[n_years_m1] GAMMA_raw;//_hyperparameter of overall annual change values - "differences" between years
  array[n_years_m1] sum_to_zero_vector[n_strata] gamma_raw;         // strata level parameters

  real GAMMA_raw_2020;//_hyperparameter of overall annual change value for the 2-year gap spanning 2020
  sum_to_zero_vector[n_strata] gamma_raw_2020;         // strata level parameters for hte two-year gap spanning 2020


}

transformed parameters{

  vector[n_routes] log_beta; // route-level intercepts (on linear scale these are multiplicative scaling factors for the predictor)
  real<lower=0> phi; //transformed sdnoise if use_pois == 0 (and therefore Negative Binomial)
  vector[n_years_m1] GAMMA;  // hyperparameter of annual differences
  matrix[n_strata,n_years_m1] gamma;         // strata-level mean differences (0-centered deviation from continental mean GAMMA)
  real GAMMA_2020;  // hyperparameter of annual differences
  vector[n_strata] gamma_2020;  // strata lelve annual differences for the 2020 span
  matrix[n_strata,n_years] yeareffect;  // matrix of estimated annual values of trajectory
  vector[n_years] YearEffect; // hyperparameter for population-wide annual trajectory
  real RHO; // slope of the log-log relationship - used as a check, not for inference
  real<lower=0> sd_gamma = sd_gamma_raw/sqrt(scaling_factor);


if(est_rho){
  RHO = RHO_raw;
}else{
  RHO = 1;
}


// overdispersion
   if(use_pois){
    phi = 1;
  }else{
    phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  }

// route level intercepts
  log_beta = sd_log_beta*log_beta_raw + log_BETA;//


// Time series
 GAMMA = sd_GAMMA * GAMMA_raw;
 // change values for the two-year span between 2019 and 2021
 // BBS was cancelled in 2020, so there are no data
 // assuming the same log-scale change value in the two years, and therefore
 // multiplying by 0.5
 GAMMA_2020 = (sd_GAMMA * 0.5) * GAMMA_raw_2020;
  gamma_2020 = (sd_gamma * 0.5) * gamma_raw_2020;

 yeareffect[,ebird_year] = zero_gammas; //fixed at zero
 YearEffect[ebird_year] = 0; //fixed at zero


// first half of time-series - runs backwards from ebird_year
  for(t in yrev){
  if(y_2020[t]){ // all years not equal to 2020 or 2021
    gamma[,t] = (sd_gamma * gamma_raw[t,]) + GAMMA[t];
        YearEffect[t] = YearEffect[t+1] - GAMMA[t]; // hyperparameter trajectory interesting to monitor but no direct inference

  }else{
    gamma[,t] = gamma_2020 + GAMMA_2020; // in 2020 and 2021, difference values represent 2-year change
      YearEffect[t] = YearEffect[t+1] - GAMMA_2020; // hyperparameter trajectory interesting to monitor but no direct inference
}
    yeareffect[,t] = yeareffect[,t+1] - gamma[,t]; //change from next year
  }
// second half of time-series - runs forwards from fixed_year
if(n_years - ebird_year){ // only used if ebird_year is not equal to n_years, i.e., ebird_year != final year of the BBS data
for(t in (ebird_year+1):n_years){
  if(y_2020[t]){ // all years not equal to 2019 or 2020
    gamma[,t-1] = (sd_gamma * gamma_raw[t-1,]) + GAMMA[t-1];//t-1 indicators to match dimensionality
        YearEffect[t] = YearEffect[t-1] + GAMMA[t-1];
      }else{
    gamma[,t-1] = gamma_2020 + GAMMA_2020; // in 2020 and 2019, difference values represent 2-year change
      YearEffect[t] = YearEffect[t-1] + GAMMA_2020;
      }
    yeareffect[,t] = yeareffect[,t-1] + gamma[,t-1]; // change from last year
  }
}



}
//
model {
  nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  nu_obs ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  log_BETA ~ student_t(3,0,2);
  RHO_raw ~ normal(1,0.5);
  if(use_t){
  log_beta_raw ~ student_t(nu,0,1);
  obs_raw ~ student_t(nu_obs,0,1);
  }else{
  log_beta_raw ~ normal(0,1);
  obs_raw ~ normal(0,1);
  }

  sd_log_beta ~ student_t(3,0,1);
  sd_obs ~ student_t(3,0,1);
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial

  sd_gamma_raw ~ gamma(2,10);
  //sd_gamma_raw ~ student_t(3,0,0.2); // prior on sd of annual differences among strata, see Smith et al. 2023 supplements for prior pred checks
  sd_GAMMA ~ student_t(3,0,0.1); // prior on sd of mean hyperparameter annual differences
  GAMMA_raw ~ std_normal();
  GAMMA_raw_2020 ~ std_normal();

//gamma_raw_2020 ~ icar_normal(n_strata, node1, node2);
target += -0.5 * dot_self(gamma_raw_2020[node1] - gamma_raw_2020[node2]);

for(t in 1:n_years_m1){

  if(y_2020[t]){ // y_2020 == 1 for all years not equal to 2020
     if(ebird_year-t){ // zero/false for the ebird year
   //gamma_raw[,t] ~ icar_normal(n_strata, node1, node2);
    target += -0.5 * dot_self(gamma_raw[t,node1] - gamma_raw[t,node2]); // ICAR prior
    }else{
     gamma_raw[t,] ~ std_normal(); // arbitrarily small because this is the ebird year and these are not included in the likelihood

     }
  }else{ //if year is on either side of 2020 - no spatial variance because there are no BBS observations that year
//These values are currently meaningless because they are not included in the likelihood
// for years with y_2020[t] == 0, these values are ignored in favour of the values
// from gamma_2020 + GAMMA_2020
    gamma_raw[t,] ~ std_normal();


  }
}


if(use_likelihood){

    vector[n_counts] E_m = compute_E(
      log_beta,
      route,
      week,
      to_vector(yeareffect),
      strat_year_idx,
      log_mean_rel_abund,
      RHO,
      sd_obs, obs_raw, observer,
      use_pois, sdnoise, noise_raw);


if(use_pois){
  count ~ poisson_log(E_m); //vectorized count likelihood with log-transformation
}else{
   count ~ neg_binomial_2_log(E_m,phi); //vectorized count likelihood with log-transformation

}
}


}

generated quantities {
  real cp;
  vector[2] cp_sel;
  real cp_draw;
  real ct;
  real p_avail;
  vector[2] p_avail_max;
  vector[2] p_avail_min;
  real cd;
  real c_area;
  real pred_count;
  real pred_count_alt1;
  real pred_count_alt2;
  real pred_count_mean;
  real pred_count_median;
  real calibration;
  real calibration_alt1;
  real calibration_alt2;
  real calibration_median;
  real calibration_mean;
  vector[n_routes] calibration_r;
  vector[n_routes] pred_count_r;
  real adj;
  array[n_counts*use_ppc] int<lower=0> y_rep;
  vector[n_counts*use_ppc] E; // this adjusts the size of E (length = 0 if use_ppc == 0)

  array[n_years,2] real raw_prediction;


if(use_t){
  adj = (1.422*(nu^0.906))/(1+(1.422*(nu^0.906)));
}else{
  adj = 1;
}

// posterior predictive check
if(use_ppc){

      E = compute_E(
      log_beta,
      route,
      week,
      to_vector(yeareffect),
      strat_year_idx,
      log_mean_rel_abund,
      RHO,
      sd_obs, obs_raw, observer,
      use_pois, sdnoise, noise_raw);


y_rep = neg_binomial_2_log_rng(E,phi);

}
// ridiculous work around to allow for the limits on pair correction

  if(use_pair){
  cp_draw = normal_rng(c_p,sd_c_p);
  cp_sel[1] = 1.0;
  cp_sel[2] = cp_draw;
  cp = max(cp_sel);
  cp_sel[1] = 2.0;
  cp_sel[2] = cp;
  cp = min(cp_sel);
  }else{
    cp = 1;
  }

  if(use_availability){
    p_avail = normal_rng(c_t,sd_c_t); // constrained between 0.001 and 0.999
  p_avail_min[2] = 1.0;
  p_avail_max[2] = 0.0001;
  p_avail_min[1] = p_avail;
  p_avail = min(p_avail_min);
  p_avail_max[1] = p_avail;
  p_avail = max(p_avail_max);

  ct = 1/p_avail;
  }else{
  ct = lognormal_rng(c_t,sd_c_t);
  p_avail = 1/ct;
  p_avail_min[2] = 1.0;
  p_avail_max[2] = 0.0001;
  p_avail_min[1] = p_avail;
  p_avail = min(p_avail_min);
  p_avail_max[1] = p_avail;
  p_avail = max(p_avail_max);
  }


  if(use_edr){
    cd = normal_rng(c_d,sd_c_d); // cd values in meters
  }else{
    cd = uniform_rng(c_d_lower,c_d_upper);
  }

  c_area = (50*3.14159*(cd^2))/1000000; // area in square km


// this calibration assumes the distribution of route-level log_betas is symetrical
  calibration = (exp(log_BETA + 0.5*(sd_log_beta/adj)^2) * cp * ct) / c_area;
  // this calibration assumes the distribution of route-level log_betas is log-normal
  calibration_alt1 = (exp(log_BETA + 0.5*(sd_log_beta)^2) * cp * ct) / c_area;
  // this calibration uses the median of the retransformed distribution - will underestimated
  calibration_alt2 = (exp(log_BETA) * cp * ct) / c_area;

  pred_count = exp(RHO*mean(log_mean_rel_abund) + log_BETA + 0.5*(sd_log_beta/adj)^2);
  pred_count_alt1 = exp(RHO*mean(log_mean_rel_abund) + log_BETA + 0.5*(sd_log_beta)^2);
  pred_count_alt2 = exp(RHO*mean(log_mean_rel_abund) + log_BETA);

// predictions (year = 2023 - year of the predicted surface) and calibrations for each route
  for(j in 1:n_routes){
    calibration_r[j] = (exp(log_beta[j]) * cp * ct) / c_area;
    pred_count_r[j] = (exp(log_beta[j] + RHO*log_mean_rel_abund[j]));
 }
  // this calibration does not assume anything about the distribution of log_beta[j]
  // but the following means and medians are estimated across the realised set of routes
  // it is also highly sensitive to long tails of the distribution of route-level variation
  calibration_mean = mean(calibration_r);
  calibration_median = quantile(calibration_r,0.5);
  pred_count_mean = mean(pred_count_r);
  pred_count_median = quantile(pred_count_r,0.5);


// useful mean predictions for each year to approximately visualise how the population changes with time average across all strata
for(y in 1:n_years){
  raw_prediction[y,1] = exp(RHO*min(log_mean_rel_abund) + log_BETA + 0.5*(sd_log_beta)^2 + YearEffect[y] );
  raw_prediction[y,2] = exp(RHO*max(log_mean_rel_abund) + log_BETA + 0.5*(sd_log_beta)^2 + YearEffect[y] );
}


}
