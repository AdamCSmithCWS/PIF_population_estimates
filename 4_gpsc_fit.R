## script to run species analysis on GPSC


library(bbsBayes2)
library(cmdstanr)





sp_aou <- Sys.getenv("aou") # aou index for each species supplied by bash script
# location for saved output
out_loc <- "/gpfs/fs7/eccc/esrp/cws/acs001/PIF/output"

conv_loc <- "/gpfs/fs7/eccc/esrp/cws/acs001/PIF/convergence"


model_components <- data.frame(model = c("napops",
                                         "napops_edr",
                                         "napops_avail",
                                         "no_napops"),
                               use_avail = c(TRUE,FALSE,TRUE,FALSE),
                               use_edr = c(TRUE,TRUE,FALSE,FALSE),
                               use_tod = c(FALSE,TRUE,FALSE,TRUE),
                               use_mdd = c(FALSE,FALSE,TRUE,TRUE))



for(vers in model_components$model){


  # load prepared species data


if(file.exists(paste0("/home/acs001/PIF/stan_data/stan_data_",vers,"_",sp_aou,".rds"))){

  stan_data <- readRDS(paste0("/home/acs001/PIF/stan_data/stan_data_",vers,"_",sp_aou,".rds"))


  model <- cmdstanr::cmdstan_model("/home/acs001/PIF/models/bbs_rel_abund_calibration_simple_season.stan",
                                   stanc_options = list("O1"))


  stan_data$use_ppc <- 0

  fit <- model$sample(data = stan_data,
                      parallel_chains = 4,
                      refresh = 0,
                      iter_warmup = 1000,
                      iter_sampling = 1000,
                      adapt_delta = 0.95,
                      max_treedepth = 11,
                      show_exceptions = TRUE)


  summ <- fit$summary()


  fit$save_object(paste0(out_loc,"/calibration_fit_alt_",vers,"_",sp_aou,".rds"))


  saveRDS(summ,paste0(conv_loc,"/parameter_summary_alt_",vers,"_",sp_aou,".rds"))

}else{
  warning(paste0("data not found for ",sp_aou,"_",vers))
}
}

