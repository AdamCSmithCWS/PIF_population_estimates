

// data here are the raw BBS counts
// density factors from PIF population estimates
// mean eBird relative abundance values from within BBS route-path buffers
data {
  int<lower=1> n_routes; // number of BBS routes in species seasonal range
  int<lower=1> n_counts; // number of counts from those routes over last 10-years
  int<lower=1> n_years;
  int<lower=1> ebird_year; // index of year for the eBird relative abundance surface
  array[ebird_year-1] int<lower=1> yrev; // reverse year vector (ebird_year-1):1

  array[n_counts] int<lower=0> count; // Raw BBS counts
  array[n_counts] int<lower=1> route; // route indicators for counts
  array[n_counts] int<lower=1> year; // year indicators for counts

  vector[n_routes] log_mean_rel_abund; // log-transformed mean relative abund within BBS route buffers

  real c_p; // mean of pair correction values
  real sd_c_p; // sd of pair correction values

  real c_t; // mean time of day correction value or availability correction
  real sd_c_t; //sd of time of day correction value or sd of availability correction
  int<lower=0,upper=1> use_availability; // indicator for calculating availability
          // if use_availability == 1, then use availability with 1/rnorm(c_d,sd_c_d)
          // if use_availability == 0, then use availability with rnorm(c_d,sd_c_d)


  real c_d; // mean distance correction value
  real sd_c_d; // sd of distance correction value
  real c_d_lower; // lower bound of distance correction value
  real c_d_upper; // upper bound of distance correction value
  int<lower=0,upper=1> use_edr; // indicator for calculating distance with edr
          // if use_edr == 1, then use edr with c_d and sd_c_d
          // if use_edr == 0, then use uniform dist between c_d_lower and c_d_upper

  int<lower=0,upper=1> use_pois; //indicator if count variation should be based on over-dispersed Poisson (if ==1) or Negative binomial (if == 0)
  int<lower=0,upper=1> use_pair; //indicator if pair-correction should be used
  int<lower=0,upper=1> use_t; //indicator if route variation should be based on heavy tailed t-distribution instead of a normal

}

//
parameters {

  real BETA; // mean route-level callibration log scale
  real<lower=0> sd_beta; // sd of calibration among routes
  vector[n_routes] beta_raw; // route-specific calibrations
  real<lower=3> nu;

  vector[n_counts*use_pois] noise_raw; // over-dispersion if use_pois == 1
  real<lower=0> sdnoise;    // sd of over-dispersion, if use_pois == 1
  vector[n_years] gamma_raw;
  real<lower=0> sd_gamma;    // sd of year effect

}

transformed parameters{

  vector[n_routes] beta; //
  vector[n_counts] E;
  real<lower=0> phi; //transformed sdnoise if use_pois == 0 (and therefore Negative Binomial)
  vector[n_years] gamma;

   if(use_pois){
    phi = 1;
  }else{
    phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  }

  beta = sd_beta*beta_raw + BETA;//

// replacing the missing fixed-year in gamma_raw
  for(y in 1:(ebird_year-1)){
    gamma[y] = gamma_raw[y];
  }
  gamma[ebird_year] = 0; // centers the year-effects on the year of the eBird abundance surface
  for(y in (ebird_year+1):n_years){
    gamma[y] = gamma_raw[y];
  }

    for(i in 1:n_counts){
       real noise;

      if(use_pois){
        noise = sdnoise*noise_raw[i];
      }else{
        noise = 0;
      }

//could consider an annual mean beta to adjust for yearly variation in the calibration
    E[i] = beta[route[i]] + gamma[year[i]] + log_mean_rel_abund[route[i]] + noise;
  }


}
//
model {
  nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  BETA ~ student_t(3,0,2);
  if(use_t){
  beta_raw ~ student_t(nu,0,1);
  }else{
  beta_raw ~ normal(0,1);
  }
  sum(beta_raw) ~ normal(0,0.001*n_routes); // soft sum to zero constraint

  sd_beta ~ student_t(3,0,1);
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  sd_gamma ~ gamma(2,10); //time-series annual variance in gamma


for(y in yrev){
  gamma_raw[y] ~ normal(gamma_raw[y+1],sd_gamma);
}
gamma_raw[ebird_year] ~ normal(0,0.001); // strong prior forcing this to 0, but it is exluded from likelihood
for(y in (ebird_year+1):n_years){
  gamma_raw[y] ~ normal(gamma_raw[y-1],sd_gamma);
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
  real<lower=0.001,upper=0.999> p_avail;
  real<lower=0.001,upper=0.999> ctp; //ensuring that the rng_normal below doesn't generate impossible availability rate values
  real cd;
  real c_area;
  real calibration;
  real calibration_alt;
  real calibration_median;
  real calibration_log;
  vector[n_routes] calibration_r;
  real adj;
  array[n_counts] int<lower=0> y_rep;
  array[n_years,2] real raw_prediction;

for(y in 1:n_years){
  raw_prediction[y,1] = exp(min(log_mean_rel_abund + BETA+gamma[y] + 0.5*sd_beta^2));
  raw_prediction[y,2] = exp(max(log_mean_rel_abund + BETA+gamma[y] + 0.5*sd_beta^2));
}
// posterior predictive check
for(i in 1:n_counts){
y_rep[i] = neg_binomial_2_log_rng(beta[route[i]] + gamma[year[i]] + log_mean_rel_abund[route[i]],phi);
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
    ctp = normal_rng(c_t,sd_c_t); // constrained between 0.001 and 0.999
  p_avail = ctp;// probability of at least one success
  // with binomial probability = ctp
  ct = 1/p_avail;
  }else{
  ct = normal_rng(c_t,sd_c_t);
  p_avail = 1/ct;
  }


  if(use_edr){
    cd = normal_rng(c_d,sd_c_d);
  }else{
    cd = uniform_rng(c_d_lower,c_d_upper);
  }

  c_area = (50*3.14159*(cd^2))/1000000; // area in square km

if(use_t){
  adj = (1.422*(nu^0.906))/(1+(1.422*(nu^0.906)));
}else{
  adj = 1;
}

// this calibration assumes the distribution of route-level betas is symetrical
  calibration = ((exp(BETA + 0.5*(sd_beta/adj)^2)) * cp * ct) / c_area;
  calibration_log = (BETA + 0.5*(sd_beta/adj)^2 + log(cp) + log(ct)) -log(c_area);

  for(j in 1:n_routes){
    calibration_r[j] = (exp(beta[j]) * cp * ct) / c_area;
  }
  // this calibration does not assume a normal distribution of beta[j]
  // but also assumes estimates a mean across the realised set of routes
  // it is also highly sensitive to long tails of the distribution of route-level variation
  calibration_alt = mean(calibration_r);
  calibration_median = quantile(calibration_r,0.5);
  // The difference between these two calibration metrics may be a good criterion
  // by which to identify strange variation among routes

}
