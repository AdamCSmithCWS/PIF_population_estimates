
// iCAR function, from Morris et al. 2019
// Morris, M., K. Wheeler-Martin, D. Simpson, S. J. Mooney, A. Gelman, and C. DiMaggio (2019).
// Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan.
// Spatial and Spatio-temporal Epidemiology 31:100301.

 functions {
   real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
     return -0.5 * dot_self(bb[n1] - bb[n2])
       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
  }
 }


// data here are the raw BBS counts
// adjustmen factors from PIF population estimates
// log mean eBird relative abundance values from within BBS route-path buffers
data {
  int<lower=1> n_routes; // number of BBS routes in species seasonal range
  int<lower=1> n_obs; // number of BBS observers in species seasonal range

  int<lower=1> n_counts; // number of counts from those routes over last 10-years
  int<lower=1> n_years;
  int<lower=1> ebird_year; // index of year for the eBird relative abundance surface
  array[ebird_year-1] int<lower=1> yrev; // reverse year vector (ebird_year-1):1 (necessary for Stan indexing)
  int<lower=1> n_doy; // number of unique values of doy (max(doy))
  int<lower=1> n_strata; // number of unique values of doy (max(doy))
  int<lower=1> mean_doy; // mean day of year across all BBS surveys

  array[n_counts] int<lower=0> count; // Raw BBS counts
  array[n_counts] int<lower=1> route; // route indicators for counts
  array[n_counts] int<lower=1> observer; // observer indicators for counts
  array[n_counts] int<lower=1> year; // year indicators for counts
  array[n_counts] int<lower=1> doy; // doy indicators for counts
  array[n_counts] int<lower=1> strata; // strata indicators for route (to support varying doy effect and trend)
  // predictor
  vector[n_routes] log_mean_rel_abund; // log-transformed mean relative abund within BBS route buffers



// spatial varying time series components
  array[n_years] int<lower=0> y_2020; //indicators for 2020 = 0 if 2020 and missing if fixed_year
  vector[n_strata] zero_gammas;   // a vector of zeros to fill fixed log_beta values for fixed_year
  //data for spatial iCAR among strata
  int<lower=1> n_edges;
  array [n_edges] int<lower=1, upper=n_strata> node1;  // node1[i] adjacent to node2[i]
  array [n_edges] int<lower=1, upper=n_strata> node2;  // and node1[i] < node2[i]


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
  int<lower=0,upper=1> use_season; //indicator if seasonal-adjustment should be used

}

//
parameters {

  real log_BETA; // mean route-level callibration log scale
  real<lower=0> sd_log_beta; // sd of calibration among routes
  real RHO_raw; // slope of the log-log relationship - used as a check, not for inference
  vector[n_routes] log_beta_raw; // route-specific calibrations
  vector[n_obs] obs_raw; // route-specific calibrations
  real<lower=0> sd_obs; // sd of calibration among routes
  real<lower=3> nu;
  real<lower=3> nu_obs;


  vector[n_counts*use_pois] noise_raw; // over-dispersion if use_pois == 1
  real<lower=0> sdnoise;    // sd of over-dispersion, if use_pois == 1


  real<lower=0> sd_gamma;    // sd of annual changes among strata
  real<lower=0> sd_GAMMA;    // sd of overall annual changes

  vector[n_years] GAMMA_raw;//_hyperparameter of overall annual change values - "differences" between years
  matrix[n_strata,n_years] gamma_raw;         // strata level parameters

  real GAMMA_raw_2020;//_hyperparameter of overall annual change value for the 2-year gap spanning 2020
  vector[n_strata] gamma_raw_2020;         // strata level parameters for hte two-year gap spanning 2020

  vector[n_knots_doy] DOY_raw;         // GAM hyper coefficients for doy
  matrix[n_strata,n_knots_doy] doy_raw;         // GAM strata level parameters
  vector<lower=0>[n_strata] sd_doy;    // sd of doy effects among strata
  real<lower=0> sd_DOY;    // sd of doy effect hyperparameter

}

transformed parameters{

  vector[n_routes] log_beta; // route-level intercepts (on linear scale these are multiplicative scaling factors for the predictor)
  vector[n_obs] obs; // observer effects
  vector[n_counts] E; // log-scale expected counts
  real<lower=0> phi; //transformed sdnoise if use_pois == 0 (and therefore Negative Binomial)
  matrix[n_doy,n_strata] doy_pred; // predictions for day of year effect in each stratum
  vector[n_doy] DOY_pred; // global mean for day of year effect
  vector[n_knots_doy] DOY_b;         // GAM coefficients hyperparameter on day of year effect
  matrix[n_strata,n_knots_doy] doy_b;         // GAM strata level parameters
  vector[n_years] GAMMA;  // hyperparameter of annual differences
  matrix[n_strata,n_years] gamma;         // strata-level mean differences (0-centered deviation from continental mean GAMMA)
  real GAMMA_2020;  // hyperparameter of annual differences
  vector[n_strata] gamma_2020;  // strata lelve annual differences for the 2020 span
  matrix[n_strata,n_years] yeareffect;  // matrix of estimated annual values of trajectory
  vector[n_years] YearEffect; // hyperparameter for population-wide annual trajectory
  real RHO; // slope of the log-log relationship - used as a check, not for inference


if(est_rho){
  RHO = RHO_raw;
}else{
  RHO = 1;
}

  // doy effects
  DOY_b = sd_DOY*DOY_raw;
  DOY_pred = doy_basis*DOY_b;

  for(k in 1:n_strata){
    doy_b[k,] = (sd_doy[k] * doy_raw[k,]) + transpose(DOY_b);
    doy_pred[,k] = doy_basis * transpose(doy_b[k,]);
  }

// overdispersion
   if(use_pois){
    phi = 1;
  }else{
    phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  }

// route level intercepts
  log_beta = sd_log_beta*log_beta_raw + log_BETA;//
// route level intercepts
  obs = sd_obs*obs_raw;//


// Time series
 GAMMA = sd_GAMMA * GAMMA_raw;
 // change values for the two-year span between 2019 and 2021
 // BBS was cancelled in 2020, so there are no data
 // assuming the same log-scale change value in the two years, and therefore
 // multiplying by 0.5
 GAMMA_2020 = sd_GAMMA * GAMMA_raw_2020 * 0.5;
  gamma_2020 = (sd_gamma * 0.5) * gamma_raw_2020;

  gamma[,ebird_year] = zero_gammas; //fixed at zero
  yeareffect[,ebird_year] = zero_gammas; //fixed at zero
  YearEffect[ebird_year] = 0; //fixed at zero


// first half of time-series - runs backwards from fixed_year
  for(t in yrev){
  if(y_2020[t]){ // all years not equal to 2020 or 2019
    gamma[,t] = (sd_gamma * gamma_raw[,t]) + GAMMA[t];
        YearEffect[t] = YearEffect[t+1] - GAMMA[t]; // hyperparameter trajectory interesting to monitor but no direct inference

  }else{
    gamma[,t] = gamma_2020 + GAMMA_2020; // in 2020 and 2019, difference values represent 2-year change
      YearEffect[t] = YearEffect[t+1] - GAMMA_2020; // hyperparameter trajectory interesting to monitor but no direct inference
}
    yeareffect[,t] = yeareffect[,t+1] - gamma[,t]; //change from next year
  }
// second half of time-series - runs forwards from fixed_year
if(n_years - ebird_year){ // only used if ebird_year is not equal to n_years, i.e., ebird_year != final year of the BBS data
for(t in (ebird_year+1):n_years){
  if(y_2020[t]){ // all years not equal to 2020 or 2021
    gamma[,t] = (sd_gamma * gamma_raw[,t]) + GAMMA[t-1];//t-1 indicators to match dimensionality
        YearEffect[t] = YearEffect[t-1] + GAMMA[t];
      }else{
    gamma[,t] = gamma_2020 + GAMMA_2020; // in 2020 and 2019, difference values represent 2-year change
      YearEffect[t] = YearEffect[t-1] + GAMMA_2020;
      }
    yeareffect[,t] = yeareffect[,t-1] + gamma[,t]; // change from last year
  }
}


    for(i in 1:n_counts){
       real noise;
       real doy_temp;

      if(use_pois){ // noise term only meaningful for observation-level extra-Poisson variance
        noise = sdnoise*noise_raw[i];
      }else{
        noise = 0; //otherwise the phi parameter in the Negative Binomial models overdispersion
      }

      if(use_season){ // using season in the likelihood if desired
        doy_temp = doy_pred[doy[i],strata[i]];
      }else{
        doy_temp = 0; //otherwise doy terms are only based on the priors
      }

    E[i] = log_beta[route[i]] + obs[observer[i]] + doy_temp + yeareffect[strata[i],year[i]] + RHO*log_mean_rel_abund[route[i]] + noise;
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
  sum(log_beta_raw) ~ normal(0,0.001*n_routes); // soft sum to zero constraint
  sum(obs_raw) ~ normal(0,0.001*n_obs); // soft sum to zero constraint

  sd_log_beta ~ student_t(3,0,1);
  sd_obs ~ student_t(3,0,1);
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
 // sd_gamma ~ gamma(2,10); //time-series annual variance in gamma
  sd_doy ~ gamma(2,10); // weak shrinkage prior with zero-avoidance
  sd_DOY ~ normal(0,1);
  DOY_raw ~ normal(0,1);         // GAM hyper coefficients for doy
// alternate non-spatial treatment of doy effects
//   for(s in 1:n_strata){
//  doy_raw[s,1:n_knots_doy] ~ normal(0,1);         // GAM strata level parameters for doy
// }

// spatial treatment of doy effects
  for(k in 1:n_knots_doy){
 doy_raw[1:n_strata,k] ~ icar_normal(n_strata, node1, node2);         // GAM strata level parameters for doy
}

  sd_gamma ~ student_t(3,0,0.2); // prior on sd of annual differences among strata, see Smith et al. 2023 supplements for prior pred checks
  sd_GAMMA ~ student_t(3,0,0.1); // prior on sd of mean hyperparameter annual differences
  GAMMA_raw ~ std_normal();
  GAMMA_raw_2020 ~ std_normal();

gamma_raw_2020 ~ icar_normal(n_strata, node1, node2);

for(t in 1:n_years){

  if(y_2020[t]){ // y_2020 == 1 for all years not equal to 2020
     if(ebird_year-t){ // zero/false for the ebird year
   gamma_raw[,t] ~ icar_normal(n_strata, node1, node2);

    }else{
     gamma_raw[,t] ~ normal(0,0.001); // arbitrarily small because this is the ebird year and these are not included in the likelihood

     }
  }else{ //if year is on either side of 2020 - no spatial variance because there are no BBS observations that year
//These values are currently meaningless because they are not included in the likelihood
// for years with y_2020[t] == 1, these values are ignored in favour of the values
// from gamma_2020 + GAMMA_2020
    gamma_raw[,t] ~ std_normal();
    sum(gamma_raw[,t]) ~ normal(0,0.001*n_strata);

  }
}




if(use_pois){
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
}else{
   count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
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
  array[n_counts] int<lower=0> y_rep;
  array[n_years,2] real raw_prediction;


if(use_t){
  adj = (1.422*(nu^0.906))/(1+(1.422*(nu^0.906)));
}else{
  adj = 1;
}

// posterior predictive check
if(use_ppc){
for(i in 1:n_counts){
  real doy_temp;
        if(use_season){ // using season in the likelihood if desired
        doy_temp = doy_pred[doy[i],strata[i]];
      }else{
        doy_temp = 0; //otherwise doy terms are only based on the priors
      }

y_rep[i] = neg_binomial_2_log_rng(log_beta[route[i]] + obs[observer[i]] + doy_temp + yeareffect[strata[i],year[i]] + RHO*log_mean_rel_abund[route[i]],phi);
}
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
