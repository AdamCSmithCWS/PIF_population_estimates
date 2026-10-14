
library(tidyverse)
library(bbsBayes2)
library(cmdstanr)
library(doParallel)
library(foreach)
#setwd("C:/Users/SmithAC/Documents/GitHub/PIF_population_estimates")

output_dir <- "output"
#output_dir <- "F:/PIF_pop_estimates/output"
yr_ebird <- 2023 # prediction year for eBird relative abundance

#selected_species <- readRDS("data/selected_species.rds")

re_fit <- TRUE

model_components <- data.frame(model = c("napops",
                                         "napops_edr",
                                         "napops_avail",
                                         "no_napops"),
                               use_avail = c(TRUE,FALSE,TRUE,FALSE),
                               use_edr = c(TRUE,TRUE,FALSE,FALSE),
                               use_tod = c(FALSE,TRUE,FALSE,TRUE),
                               use_mdd = c(FALSE,FALSE,TRUE,TRUE))


sp_sel <- "Chestnut-collared Longspur"
sp_aou <- 5380
sp_ebird <- "chclon"


#Parallel species loop ---------------------------------------------------
## run on machine with 64 cores, 128GB RAM
# n_cores <- 3 # requires total of 32 (8*4) cores 4 parallel chains for each of 8 species
# cluster <- makeCluster(n_cores, type = "PSOCK")
# registerDoParallel(cluster)
#
#
# test <- foreach(sp_sel = selected_species,
#                 .packages = c("bbsBayes2",
#                               "tidyverse",
#                               "cmdstanr",
#                               "ebirdst"),
#                 .errorhandling = "pass") %dopar%
#   {

#for(sp_sel in selected_species){

for(vers in model_components$model){

  if(!file.exists(paste0("stan_data/stan_data_",vers,"_",sp_aou,"_",sp_ebird,".rds"))){

    next}

  stan_data <- readRDS(paste0("stan_data/stan_data_",vers,"_",sp_aou,"_",sp_ebird,".rds"))
  # params_to_summarise <- c("nu",
  #                          "log_BETA",
  #                          "log_beta_raw",
  #                          "RHO",
  #                          "log_beta",
  #                          "sd_log_beta",
  #                          "sd_obs",
  #                          "obs",
  #                          "nu_obs",
  #                          "sdnoise",
  #                          "sd_gamma",
  #                          "GAMMA",
  #                          "sd_GAMMA",
  #                          "gamma",
  #                          "gamma_2020",
  #                          "GAMMA_2020",
  #                          "yeareffect",
  #                          "YearEffect",
  #                          "phi",
  #                          "cp",
  #                          "cp_sel",
  #                          "ct",
  #                          "cd",
  #                          "p_avail",
  #                          "c_area",
  #                          "calibration",
  #                          "calibration_alt2",
  #                          "calibration_median",
  #                          "calibration_alt1",
  #                          "calibration_mean",
  #                          "calibration_r",
  #                          "adj",
  #                          "raw_prediction",
  #                          "sd_doy",
  #                          "sd_DOY",
  #                          "doy_raw",
  #                          "DOY_raw",
  #                          "DOY_b",
  #                          "doy_b",
  #                          "DOY_pred",
  #                          "doy_pred",
  #                          "pred_count",
  #                          "pred_count_alt2",
  #                          "pred_count_median",
  #                          "pred_count_alt1",
  #                          "pred_count_mean",
  #                          "pred_count_r")

  #vers <- "seas_opt_"

  # if(vers == "seas_opt_"){ # alternate version of model that removes the seasonal component
  #   stop("Wrong model called check seas_opt_")
  #   stan_data[["use_season"]] <- 0
  #   model <- cmdstanr::cmdstan_model("models/ebird_rel_abund_calibration_spatial_year_doy_sep_obs_rte_rhoopt_seas_opt.stan")
  # }

    #model <- cmdstanr::cmdstan_model("models/ebird_rel_abund_calibration_spatial_year_doy_sep_obs_rte_rhoopt.stan")
    model <- cmdstanr::cmdstan_model("models/bbs_rel_abund_calibration_simple_season.stan",
                                     stanc_options = list("O1"))


stan_data$use_ppc <- 0

  fit <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500,
                    iter_warmup = 2000,
                    iter_sampling = 1000,
                    adapt_delta = 0.95,
                    max_treedepth = 11,
                    show_exceptions = TRUE)


summ <- fit$summary()


fit$save_object(paste0(output_dir,"/calibration_fit_alt_",vers,"_",sp_aou,"_",sp_ebird,".rds"))


saveRDS(summ,paste0("convergence/parameter_summary_alt_",vers,"_",sp_aou,"_",sp_ebird,".rds"))

}

}

parallel::stopCluster(cluster)



