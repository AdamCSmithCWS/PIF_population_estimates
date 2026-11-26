
library(sf)
library(tidyverse)
library(ebirdst)
library(patchwork)
library(terra)
library(bbsBayes2)
library(cmdstanr)
library(ggrepel)
library(tidyterra)
library(tidybayes)
library(HDInterval)
library(ggtext)
library(doParallel)
library(foreach)

# Script to summarise the fitted models -----------------------------------

## generates various graphs, summaries, and comparisons for each species


source("functions/hpdi.R")
source("functions/posterior_abundance.R")


yr_ebird <- 2025 # prediction year for eBird relative abundance


# mean counts to id routes on which species has been observed -------------
mean_counts_all <- readRDS("data/mean_counts_by_route.rds")

raw_counts_all <- readRDS("data/all_counts_by_route.rds")

#output_dir <- "output"

output_dir <- "E:/PIF_pop_estimates/output"


# adjustments <- readRDS("pop_estimates_control_file.rds")


routes_buf_all <- readRDS("data/all_routes_buffered_400m.rds")

strata <- bbsBayes2::load_map("bbs") |>
  mutate(bcr = as.integer(bcr))

strata_names <- strata %>%
  sf::st_drop_geometry()







# simple quantile wrapper that removes names and supplies 3 signif figures
quant <- function(x,q){
  y <- quantile(x,q,names = FALSE,digits = 3)
  return(y)
}
## skew function
skew <- function(x, na.rm = FALSE){
  if(na.rm) x <- x[!is.na(x)]
  n <- length(x)
  sum((x - mean(x))^3)/(n - 2)/var(x)^(3/2)
}

# base map for countries and continents -----------------------------------
country_codes <- readxl::read_xlsx("data/iso_codes.xlsx") %>%
  rename(country_name = Name) #

americas_codes <- country_codes %>%
  filter(CC %in% c("SA","NA"))

states <- rnaturalearthhires::states10 %>%
  left_join(country_codes,by = c("adm0_a3" = "a-3")) %>%
  filter(country_name %in% c("Canada","United States of America","Mexico")) %>%
  group_by(name_en,country_name,iso_3166_2,postal,name_fr,name_es) |>  #%>%
  summarise(do_union = FALSE) |>
  select(name_en,country_name,iso_3166_2,postal,name_fr,name_es)

countries <- rnaturalearthhires::countries10 %>%
  rename_with(.fn = str_to_lower) %>%
  left_join(country_codes,by = c("adm0_a3" = "a-3")) %>%
  select(CC,adm0_a3,country_name,name_en) %>%
  mutate(CC = ifelse(name_en %in% c("Western Sahara",
                                    "Somaliland",
                                    "South Sudan"),
                     "AF",CC),
         CC = ifelse(name_en %in% c("Palestine",
                                    "Turkish Republic of Northern Cyprus"),
                     "AS",CC),
         CC = ifelse(name_en %in% c("Kosovo"),
                     "EU",CC)) %>%
  filter(!is.na(country_name))


USACAN <- countries %>%
  filter(name_en %in% c("Canada","United States of America")) %>%
  ungroup() %>%
  summarise(do_union = FALSE)



continents <- countries %>%
  group_by(CC) %>%
  summarise(do_union = FALSE)
#plot(continents)

#Bird Studies Canada and NABCI.  2014.  Bird Conservation Regions.  Published by
#Bird Studies Canada on behalf of the North American Bird Conservation Initiative.
# https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions  Accessed:  Nov 8 2024
bcrs <- readRDS("data/bcr2026_statprov_simplified.rds") |>
  group_by(bcr_label, bcr_label_) |>
  summarise(do_union = FALSE) |>
  rename(bcr = bcr_label,
         bcr_name = bcr_label_)



estimate_rho <- FALSE



model_components <- data.frame(model = c("napops",
                                         "napops_edr",
                                         "napops_avail",
                                         "no_napops"),
                               use_avail = c(TRUE,FALSE,TRUE,FALSE),
                               use_edr = c(TRUE,TRUE,FALSE,FALSE),
                               use_tod = c(FALSE,TRUE,FALSE,TRUE),
                               use_mdd = c(FALSE,FALSE,TRUE,TRUE))



output_dir <- "E:/PIF_pop_estimates/output"
#output_dir <- "output"
 #

sps_list <- readRDS("pop_estimates_control_file.rds")

# species loop convergence ------------------------------------------------------------
summ_out <- NULL

 for(i in 1:nrow(sps_list)){


  sp_sel <- unname(unlist(sps_list[i,"common_name"]))

 sp_aou <- unname(unlist(sps_list[i,"sp_aou"]))
  sp_ebird <- unname(unlist(sps_list[i,"sp_ebird"]))





if(!file.exists(paste0("convergence/parameter_summary_alt_",sp_aou,".rds"))){
  next
}# skip if no convergence

summ <- readRDS(paste0("convergence/parameter_summary_alt_",sp_aou,".rds")) |>
  mutate(species = sp_sel)
# convergence summary -----------------------------------------------------


summ_out <-  bind_rows(summ_out,summ)

}




saveRDS(summ_out,"convergence_summary.rds")

problems <- summ_out |>
  filter(rhat > 1.02, ess_bulk < 200) |>
  mutate(variable_type = sub("\\[.*","",variable)) |>
  group_by(variable_type,species) |>
  summarise(n = n(),
            max_rhat = max(rhat),
            min_ess = min(ess_bulk))



sp_re_run <- problems |>
  ungroup() |>
  select(species) |>
  distinct()



aou_list2 <- sps_list |>
  filter(common_name %in% sp_re_run$species) |>
  select(sp_aou) |>
  distinct()

if(nrow(aou_list2) != 0){
write_tsv(aou_list2, "aou_list_rerun.txt",
          col_names = FALSE)
}else{
  warning("Congratulations! all species models have converged")
}


# end convergence check ---------------------------------------------------




# species loop summary ------------------------------------------------------------

## run on machine with
n_cores <- 4 #
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)
# library(sf)
# library(tidyverse)
# library(ebirdst)
# library(patchwork)
# library(terra)
# library(bbsBayes2)
# library(cmdstanr)
# library(ggrepel)
# library(tidyterra)
# library(tidybayes)
# library(HDInterval)
# library(ggtext)

test <- foreach(i = 1:nrow(sps_list),
                .packages = c("bbsBayes2",
                              "tidyverse",
                              "cmdstanr",
                              "ebirdst",
                              "sf",
                              "exactextractr",
                              "tidyterra",
                              "terra",
                              "tidybayes",
                              "ggtext",
                              "patchwork",
                              "HDInterval",
                              "ggrepel"),
                .errorhandling = "pass") %dopar%
  {
    summ_out <- NULL
    pop_ests_out_save <- NULL
    trajectories_out_all <- NULL

    sp_sel <- unname(unlist(sps_list[i,"common_name"]))

    sp_aou <- unname(unlist(sps_list[i,"sp_aou"]))
    sp_ebird <- unname(unlist(sps_list[i,"sp_ebird"]))



    if(!file.exists(paste0("data/species_relative_abundance_",yr_ebird,"/",sp_ebird,"_relative_abundance.rds"))){
      next} # skip if no saved data


if(!file.exists(paste0(output_dir,"/calibration_fit_alt_",
                       sp_aou,".rds"))){
  next} #skip of no saved model

# load saved fitted model, summary, and data ------------------------------
    if(!file.exists(paste0("convergence/parameter_summary_alt_",sp_aou,".rds"))){
      next
    }# skip if no convergence

    summ <- readRDS(paste0("convergence/parameter_summary_alt_",sp_aou,".rds")) |>
      mutate(species = sp_sel)
    # convergence summary -----------------------------------------------------


    summ_out <-  bind_rows(summ_out,summ)


  fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",sp_aou,".rds"))


  stan_data <- readRDS(paste0("stan_data/stan_data_",sp_aou,".rds"))

  combined <- readRDS(paste0("stan_data/dataframe_",sp_aou,".rds"))

  routes_buf <- readRDS(paste0("stan_data/routes_buffered_",sp_aou,"_",sp_ebird,".rds"))

  # tr_summ <- combined |>
  #   group_by(route_name,route) |>
  #   summarise(mean_obs = mean(count),
  #             .groups = "drop")
  #
  # tr_keep <- tr_summ |> filter(mean_obs > 0)
  # rt_k <- routes_buf |> filtercount# rt_k <- routes_buf |> filter(mean_count > 0)

# posterior of calibration ------------------------------------------------

# load the hyperparameter estimates,
# but only used to provide structure for the trimmed mean
# calculation below, not used in the final calibration
cali_alt_post <- fit$draws(variables = "log_BETA",
                           format = "df")


# trimmed mean function to exclude the outer 5% of routes with extreme
# calibration values. These extreme values have a strong influence on the
# mean of the retransformed values. particularly routes with very high calibration
# values which greatly increase the overall mean calibration and the total
# population estimate.
# Ideally, we'll figure out what is missing from the model and why
# we're getting these heavy tails, but until then this trimmed-mean generates
# much more robust population estimates.

  # ID routes with extreme betas -----------------------------------------------

  betas <- summ %>% filter(grepl("log_beta[",variable,fixed = TRUE))

  # tmp <- summ %>% filter(grepl("ct",variable,fixed = TRUE))
  # tmp <- summ %>% filter(variable == "cd")
  #
  w_excl <- c(which(betas$mean > quantile(betas$mean,0.975)), # identifying the routes in the outer 2.5% tails
              which(betas$mean < quantile(betas$mean,0.025)))

  #w_excl2 <- tr_summ |> filter(mean_obs == 0) |> select(route) |> unlist() |> unname()

  betas_keep <- c(1:stan_data$n_routes)[-w_excl] # dropping the routes in the tails

#posterior of the route-level calibration values for centre 95% of routes (trimmed mean)
cali_route_post <- fit$draws(variables = paste0("log_beta[",betas_keep,"]"),
                           format = "df")

# calculate the trimmed mean
cali_trim_post <- cali_route_post |>
  mutate(log_BETA_trim = rowMeans(pick(starts_with("log_beta")), na.rm = TRUE)) |>
  left_join(cali_alt_post,
            by = c(".chain", ".iteration", ".draw"))

# cali_route_all <- fit$draws(variables = paste0("log_beta[",1:stan_data$n_routes,"]"),
#                             format = "df") |>
#   mutate(log_BETA = rowMeans(pick(paste0("log_beta[",betas_keep,"]")), na.rm = TRUE))
#


use_approx_edr <- TRUE

if(unname(unlist(sps_list[i,"EDR_model"])) %in% c("Intercept",
                                  "Roadside") &
   !(unname(unlist(sps_list[i,"highly_mobile"])) == TRUE)){
edr <-  exp(rnorm(nrow(cali_trim_post),log(unname(unlist(sps_list[i,"edr"]))),unname(unlist(sps_list[i,"edr_sd"]))/unname(unlist(sps_list[i,"edr"]))))

c_area = (50*3.14159*(edr^2))/1000000; # area in square km

vers <- "NA-POPS edr"
}

# if(unname(unlist(sps_list[i,"highly_mobile"])) &
#    use_approx_edr){
#   edr <-  exp(rnorm(nrow(cali_trim_post),unname(unlist(sps_list[i,"approx_log_edr_median"])),unname(unlist(sps_list[i,"approx_log_edr_sd"]))))
#
#   c_area = (50*3.14159*(edr^2))/1000000; # area in square km
#
#   vers <- "approximate edr"
# }

if(unname(unlist(sps_list[i,"highly_mobile"])) &
   !use_approx_edr){
  edr <-  runif(nrow(cali_trim_post),unname(unlist(sps_list[i,"Dist.Lower"])),unname(unlist(sps_list[i,"Dist.Upper"])))

  c_area = (50*3.14159*(edr^2))/1000000; # area in square km

  vers <- "PIF MDD"
}

if((unname(unlist(sps_list[i,"highly_mobile"])) &
   use_approx_edr) |
   (!unname(unlist(sps_list[i,"EDR_model"])) %in% c("Intercept",
                                                   "Roadside") &
    use_approx_edr) ){
  edr <-  exp(rnorm(nrow(cali_trim_post),unname(unlist(sps_list[i,"approx_log_edr_median"])),unname(unlist(sps_list[i,"approx_log_edr_sd"]))))

  c_area = (50*3.14159*(edr^2))/1000000; # area in square km

  vers <- "approximate edr"
}




# availability ------------------------------------------------------------


if(unname(unlist(sps_list[i,"use_availability"])) &
   !unname(unlist(sps_list[i,"peak_detectability_outside_morning"]))){

  avail <- rnorm(nrow(cali_trim_post),unname(unlist(sps_list[i,"availability"])),unname(unlist(sps_list[i,"availability_sd"])))
  if(any(avail > 1)){
    avail[which(avail > 1)] <- 1
  }
  if(any(avail <= 0)){
    avail[which(avail <= 0)] <- 0.01
  }

  vers <- paste(vers,"NA-POPS Availability",collapse = " ; ")

}else{

  avail <- 1/rlnorm(nrow(cali_trim_post),unname(unlist(sps_list[i,"TimeAdj.meanlog"])),unname(unlist(sps_list[i,"TimeAdj.sdlog"])))
  if(any(avail > 1)){
    avail[which(avail > 1)] <- 0.99
  }

  vers <- paste(vers,"PIF ToD",collapse = " ; ")

}


# pair --------------------------------------------------------------------


pair <- rnorm(nrow(cali_trim_post),unname(unlist(sps_list[i,"Pair2"])),
              0.13) # original was 0.13

if(any(pair > 2)){
  pair[which(pair > 2)] <- 2
}
if(any(pair < 1)){
  pair[which(pair < 1)] <- 1
}


#calibration_r[j] = (exp(log_beta[j]) * cp * ct) / c_area;

cali_trim_post <- cali_trim_post |>
  mutate(edr_ = edr,
         c_area_ = c_area,
         avail_ = avail,
         pair_ = pair,
         calibration = (exp(log_BETA) * pair_ *(1/avail_) * (1/c_area_)))

saveRDS(cali_trim_post,paste0("estimates/calibration_posterior_",sp_aou,"_",sp_ebird,".rds"))



# calculating population --------------------------------------------------


# load relative abundance surface -----------------------------------------


# this is the range-wide (global) relative abundance surface for the
# selected set of weeks that overlap the BBS field season
breed_abundance <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",
                                  sp_ebird,
                                  "_derived_breeding_relative_abundance.rds"))

names(breed_abundance) <- "breeding_abundance"

breed_abundance_full <- breed_abundance






  cali_use <- cali_trim_post[,c(".draw",".iteration",".chain","calibration")]


### 3^2 scaling is to account for the area of each grid-cell 9-km^2
###


# BBS strata level estimates (to match PIF stratification) ----------------


## process is to summarise the relative abundance within given regions
## Then multiply those sums by the full posterior of the calibration
## to generate estimates with uncertainty
##
## One could apply the full posterior of the calibration to the
## pixel values (generating a full posterior of the mapped surface)
## one could also combine the uncertainty of the calibration with the
## uncertainty of the pixel=level relative abundances, but both of those
## approaches would have extreme computational costs.
##
## The first approach would be mathematically equivalent to what's done here
## because the uncertainty at a pixel level would be the same across each posterior
## draw i.e., sum(calibration * pixel-values) = calibration*sum(pixel-values)
##
## The second approach would account for the uncertainty in the relative abundance
## surface and is worth exploring, with full consideration of the spatial autocorrelation
## in the eBird uncertainty estimates (e.g., each pixel's uncertainty is not independent,
##  but how dependent is it?)
##

strata_proj <- st_transform(strata,
                            crs = st_crs(breed_abundance))


## extract the sum of the relative abundance surface in each region
abundance_in_strata <- terra::extract(breed_abundance_full,
                                      strata_proj,
                                       fun = sum,
                                       na.rm = TRUE,
                                       ID = TRUE,
                                       exact = FALSE)



strata_abund <- data.frame(region = strata_proj$strata_name,
                           region_type = "strata",
                           sum_abund = 3^2 * unlist(abundance_in_strata[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                                 fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird,
         strata_name = region) #%>%
  #left_join(.,strata_names)
  #filter(pop_mean != 0)
strata_sampled <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",sp_ebird,"_bbs_strata_relative_abundance.rds")) %>%
  select(strata_name,mean_ebird_abundance:log_ratio)

strata_abund <- strata_abund %>%
  left_join(strata_sampled,
            by = "strata_name")




# Annual strata-level abundances ------------------------------------------
# population size trajectories for each stratum with BBS data
# Excludes strata with population but no BBS time-series estimates
strats <- combined %>%
  select(strata_name,strata) %>%
  distinct() %>%
  arrange(strata)

n_yr <- max(combined$yr)
n_strat <- max(combined$strata)
fyr <- min(combined$year)


yrs <- data.frame(year = fyr:(fyr+n_yr-1),
                  yr = 1:n_yr)

time_series <- expand_grid(yrs,strats)


year_draws <- gather_draws(fit, yeareffect[strata,yr]) %>%
  inner_join(strats)

cali_use_y <- expand_grid(cali_use,yrs) %>%
    inner_join(year_draws,
               by = c("yr",".draw",".iteration",".chain"))






  strata_population_posterior <- data.frame(region = strata_proj$strata_name,
                             region_type = "strata",
                             sum_abund = 3^2 * unlist(abundance_in_strata[[2]])) %>%
    filter(!is.na(sum_abund)) %>%
    expand_grid(yrs) %>%
    inner_join(cali_use_y,
               by = c("region" = "strata_name",
                      "yr","year")) %>%
    mutate(population = sum_abund*calibration*exp(.value),
           species = sp_sel,
           species_ebird = sp_ebird,
           strata_name = region) %>%
    left_join(.,strata_names,
              by = "strata_name") |>
    mutate(country_code = str_sub(strata_name,1,2),
           prov_state = str_sub(strata_name,4,5),
           bcr = str_sub(strata_name,7))




 # saveRDS(strata_population_posterior,paste0("output/strata_population_posterior_",sp_aou,"_",sp_ebird,".rds"))




  strata_trajectories <- strata_population_posterior %>%
    group_by(region, region_type, year) %>%
  summarise(pop_mean = mean(population),
            pop_median = median(population),
            pop_lci_80 = interval_function_hpdi(population,probs = 0.1),
            pop_uci_80 = interval_function_hpdi(population,probs = 0.9),
            pop_lci_95 = interval_function_hpdi(population,probs = 0.025),
            pop_uci_95 = interval_function_hpdi(population,probs = 0.975))



  # traj_plot <- ggplot(data = strata_trajectories,
  #                aes(x = year,y = pop_median))+
  #   geom_ribbon(aes(ymin = pop_lci_80,ymax = pop_uci_80),
  #               alpha = 0.3)+
  #   geom_line()+
  #   facet_wrap(vars(region),
  #              scales = "free")

  #traj_plot

  #saveRDS(traj_plot,paste0("figures/saved_ggplots/abund_trajectory_strata_",sp_aou,"_",sp_ebird,".rds"))






  ## state province trajectories
  state_trajectories <- strata_population_posterior %>%
    group_by(prov_state,year, .draw) %>%
    summarise(population = sum(population)) %>%
    group_by(prov_state,year) %>%
    summarise(pop_mean = mean(population),
              pop_median = median(population),
              pop_lci_80 = interval_function_hpdi(population,probs = 0.1),
              pop_uci_80 = interval_function_hpdi(population,probs = 0.9),
              pop_lci_95 = interval_function_hpdi(population,probs = 0.025),
              pop_uci_95 = interval_function_hpdi(population,probs = 0.975)) |>
    mutate(region_type = "prov_state",
           region = prov_state)





  ## USA_CAN
  USA_CAN_trajectories <- strata_population_posterior %>%
    group_by(year, .draw) %>%
    summarise(population = sum(population)) %>%
    group_by(year) %>%
    summarise(pop_mean = mean(population),
              pop_median = median(population),
              pop_lci_80 = interval_function_hpdi(population,probs = 0.1),
              pop_uci_80 = interval_function_hpdi(population,probs = 0.9),
              pop_lci_95 = interval_function_hpdi(population,probs = 0.025),
              pop_uci_95 = interval_function_hpdi(population,probs = 0.975)) |>
    mutate(region_type = "USA_CAN",
           regino = "USA_CAN")

  #
  #
  # traj_plot <- ggplot(data = USA_CAN_trajectories,
  #                     aes(x = year,y = pop_median))+
  #   geom_ribbon(aes(ymin = pop_lci_80,ymax = pop_uci_80),
  #               alpha = 0.3)+
  #   geom_line()+
  #   scale_y_continuous(limits = c(0,NA),
  #                      labels = scales::unit_format(unit = "M", scale = 1e-6))+
  #   scale_x_continuous(breaks = seq(2015,2025,2))+
  #   xlab("")+
  #   ylab("Annual population size\nUnited States of America and Canada")+
  #   labs(title = sp_sel)+
  #   theme_bw()
  #
  #
  # traj_plot

  trajectories_out <- bind_rows(USA_CAN_trajectories,
                                state_trajectories,
                                strata_trajectories) |>
    mutate(species = sp_sel,
           version = vers)

  trajectories_out_all <- bind_rows(trajectories_out_all,
                                    trajectories_out)

#saveRDS(traj_plot,paste0("figures/saved_ggplots/abund_trajectory_",sp_aou,"_",sp_ebird,".rds"))
  # compare to traditional estimates ----------------------------------------
#
# strata_trad_j <- strata_trad %>%
#   mutate(med.PopEst = pop_median,
#          LCI80.PopEst = pop_lci_80,
#          UCI80.PopEst = pop_uci_80) %>%
#   select(st_abrev,
#          BCR,
#          med.PopEst,
#          LCI80.PopEst,
#          UCI80.PopEst)
#
# pal_labs <- c("< -1",
#               "-1, -0.5",
#               "-0.5, -0.15",
#               "-0.15, 0.15",
#               "0.15, 0.5",
#               "0.5, 1",
#               "> 1")
#
# strata_compare <- strata_abund %>%
#   full_join(.,strata_trad_j,
#             by = c("prov_state" = "st_abrev",
#                    "bcr" = "BCR")) %>%
#   mutate(ratio_cat = cut(log_ratio,breaks = c(-Inf,-1,-0.5,-0.15,0.15,0.5,1,Inf),
#                          levels = c("(-Inf,-1]","(-1,-0.5]",
#                          "(-0.5,-0.15]","(-0.15,0.15]",
#                          "(0.15,0.5]","(0.5,1]","(1, Inf]"),
#                          labels = pal_labs),
#          log_pop_ratio = log(pop_median/med.PopEst))
#
# div_pal <- viridisLite::viridis(n = 7,option = "turbo", direction = -1)
#
# names(div_pal) <- pal_labs#levels(strata_compare$ratio_cat)
#
#
# comp_trad_new_plot <- ggplot(data = strata_compare,
#                              aes(x = med.PopEst,
#                                  y = pop_median))+
#   # geom_abline(slope = 1, intercept = 0)+
#   # geom_abline(slope = 2, intercept = 0,linetype = 2)+
#   # geom_abline(slope = 5, intercept = 0,linetype = 3)+
#   # geom_abline(slope = 0.5, intercept = 0,linetype = 2)+
#   # geom_abline(slope = 0.2, intercept = 0,linetype = 3)+
#   geom_errorbar(aes(ymin = pop_lci_80,ymax = pop_uci_80),
#                 alpha = 0.3, width = 0)+
#   geom_errorbarh(aes(xmin = LCI80.PopEst,xmax = UCI80.PopEst),
#                 alpha = 0.3)+
#   geom_point(aes(colour = ratio_cat))+
#   geom_smooth(method = "lm", se = FALSE,
#               colour = grey(0.5))+
#   geom_text_repel(aes(label = strata_name),
#                   size = 2)+
#   scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
#   scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
#   scale_colour_manual(values = div_pal,
#                       name = "BBS sampling bias \n log(sampled/avail)",
#                       na.translate = FALSE)+
#   # scale_colour_viridis_c(option = "turbo", direction = -1,
#   #                        name = "BBS sampling bias \n log(sampled/avail)")+
#   xlab("Existing PIF population estimate")+
#   ylab("Revised PIF population estimate")+
#   labs(title = paste(sp_sel))+
#   theme_bw()
#
# png(filename = paste0("Figures/comp_trad_new_alt_",sp_aou,"_",sp_ebird,".png"),
#     res = 300,
#     height = 6,
#     width = 6,
#     units = "in")
# print(comp_trad_new_plot)
# dev.off()
#
# pdf(paste0("Final_figures/comp_trad_new_alt_",sp_aou,"_",sp_ebird,".pdf"),
#     height = 4,
#     width = 7)
# print(comp_trad_new_plot)
# dev.off()
#





  ## state province territory
  ##

  states_proj <- st_transform(states,
                              crs = st_crs(breed_abundance))

  ## extract the sum of the relative abundance surface in each region
  abundance_in_states <- terra::extract(breed_abundance_full,
                                        states_proj,
                                        fun = sum,
                                        na.rm = TRUE,
                                        ID = TRUE,
                                        exact = FALSE)


  state_abund <- data.frame(region = states_proj$name_en[abundance_in_states$ID],
                            region_type = "State_Prov",
                            sum_abund = 3^2 * unlist(abundance_in_states[[2]])) %>%
    mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                                 draws = cali_use$calibration,
                                 fun = "mean"),
           pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                   fun = "median",p = 0.5),
           pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                   fun = "hpdi",p = 0.1),
           pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                   fun = "hpdi",p = 0.9),
           pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                   fun = "hpdi",p = 0.025),
           pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                   fun = "hpdi",p = 0.975),
           species = sp_sel,
           species_ebird = sp_ebird)




# National summaries ------------------------------------------------------





countries_proj <- st_transform(countries,
                            crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_countries <- terra::extract(breed_abundance,
                                         countries_proj,
                                      fun = sum,
                                      na.rm = TRUE,
                                      ID = TRUE,
                                      exact = FALSE)

abundance_in_countries <- abundance_in_countries %>%
  filter(!is.na(breeding_abundance),
         breeding_abundance > 0)



country_abund <- data.frame(region = countries_proj$country_name[abundance_in_countries$ID],
                            region_type = "country",
                           sum_abund = 3^2 * unlist(abundance_in_countries[[2]])) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#country_abund # will be combined with other regional estimates below



# BCR summaries -----------------------------------------------------------






bcr_proj <- st_transform(bcrs,
                               crs = st_crs(breed_abundance))


abundance_in_bcr <- terra::extract(breed_abundance,
                                   bcr_proj,
                                         fun = sum,
                                         na.rm = TRUE,
                                         ID = TRUE,
                                         exact = FALSE)

abundance_in_bcr <- abundance_in_bcr %>%
  filter(!is.na(breeding_abundance),
         breeding_abundance > 0)



bcr_abund <- data.frame(region = as.character(bcr_proj$bcr[abundance_in_bcr$ID]),
                            region_type = "bcr",
                            sum_abund = 3^2 * unlist(abundance_in_bcr[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#bcr_abund
# will be combined with other regional estimates below


# mapping with BCR extent -------------------------------------------------



re_do_map <- TRUE #option to skip because they take a long time.
if(re_do_map){
  # abund_mapable <- breed_abundance %>%
  #   terra::project(.,"EPSG:9822")

  bcrs_w_pop <- bcr_abund %>%
    filter(pop_median > 0)

  bcrs_plot <- bcrs %>%
    st_transform(crs = st_crs(strata)) %>%
    filter(bcr %in% bcrs_w_pop$region)

countries_plot <- countries %>%
  st_transform(crs = st_crs(strata)) %>%
  filter(CC == "NA",
         country_name != "Greenland")

strata_w_pop <- strata_abund %>%
  filter(sum_abund > 0)

strata_bound <- strata %>%
  filter(strata_name %in% strata_w_pop$region)

bb <- sf::st_bbox(strata_bound)

  breed_abundance_plot <- breed_abundance
  names(breed_abundance_plot) <- "breeding"

  NA_clip <- countries_plot |>
    st_transform(st_crs(breed_abundance_plot))

  breed_abundance_plot <- breed_abundance_plot %>%
    terra::crop(NA_clip) |>
    mutate(breeding = ifelse(breeding == 0,NA,breeding),
           breeding = (breeding*median(cali_use$calibration))/9) # per km2 mean density

hist_df <- data.frame(breding = as.numeric(breed_abundance_plot$breeding)) |>
  drop_na() |>
  rename(birds_km2 = breeding)

  # hist_inset <- ggplot(data = hist_df)+
  #   geom_histogram(aes(x = birds_km2,
  #                      colour = birds_km2),
  #                  binwidth = 1)+
  #   scale_colour_viridis_b(direction = -1,
  #                          option = "B",
  #                          end = 0.9)
  # hist_inset


qs <- quantile(hist_df$birds_km2,
                c(0,0.1,0.3,0.5,0.75,0.9,0.95,0.99,1))
capt <- paste0("Colour scale binned to show breaks at the following percentiles of\nthe distribution of density values in the map: \n",
              paste((c(0,0.1,0.3,0.5,0.75,0.9,0.95,0.99)*100),collapse = "%,"),"%", "\n The lightest colour shows the grid cells with density in the lowest 10%",
              "\n and the darkest colour shows the highest 1%.")
names(qs) <- signif(qs,2)

min(hist_df$birds_km2)
max(hist_df$birds_km2)

# Mapping -----------------------------------------------------------------


  abund_map <- ggplot()+
    geom_sf(data = countries_plot)+
    geom_sf(data = bcrs_plot, fill = "white")+
    geom_spatraster(data = breed_abundance_plot,
                    maxcell = 16000000)+
    geom_sf(data = bcrs_plot, fill = NA)+
    # geom_sf_text(data = bcrs_plot,aes(label = bcr),
    #              size = 3)+
    #geom_sf(data = routes_buf,aes(colour = mean_count),fill = NA)+
    coord_sf(xlim = c(bb[c("xmin","xmax")]),
             ylim = c(bb[c("ymin","ymax")]))+
    #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
    scale_fill_viridis_b(direction = -1,
                         option = "B",
                         na.value = NA,
                         end = 0.98,
                         begin = 0.1,
                         breaks = qs,
                         name = paste0(sp_sel,"\n","Birds/km<sup>2</sup>"))+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(legend.title = element_markdown())+
    labs(title = paste(sp_sel),
         subtitle = paste0(bquote("Density within 3km^2 grid cells")),
         caption = capt)


  #

  png(filename = paste0("Figures/density_maps/density_map_",sp_aou,"_",sp_ebird,".png"),
      res = 400,
      height = 7,
      width = 6.5,
      units = "in")
  print(abund_map)
  dev.off()
  # save(list = c("bb",
  #               "breed_abundance_plot",
  #               "countries_plot",
  #               "bcrs_plot"),
  #      file = paste0("figures/saved_ggplots/abund_map_alt_",vers,"_",sp_aou,"_",sp_ebird,".rdata"))

}

#saveRDS(comp_trad_new_plot,paste0("figures/saved_ggplots/trad_vs_new_alt_",vers,"_",sp_aou,"_",sp_ebird,".rds"))
#saveRDS(strata_compare,paste0("estimates/strata_comparison_",vers,"_",sp_aou,"_",sp_ebird,".rds"))




# USACAN estimates -----------------------------------------



USACAN_proj <- st_transform(USACAN,
                               crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_USACAN <- terra::extract(breed_abundance,
                                     USACAN_proj,
                                         fun = sum,
                                         na.rm = TRUE,
                                         ID = TRUE,
                                         exact = FALSE)



USACAN_abund <- data.frame(region = "USACAN",
                          region_type = "USACAN",
                            sum_abund = 3^2 * unlist(abundance_in_USACAN[[2]])) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#USACAN_abund
# will be combined with other regional estimates below

# Continent estimates -----------------------------------------



continents_proj <- st_transform(continents,
                            crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_continents <- terra::extract(breed_abundance,
                                      continents_proj,
                                      fun = sum,
                                      na.rm = TRUE,
                                      ID = TRUE,
                                      exact = FALSE)

abundance_in_continents <- abundance_in_continents %>%
  filter(!is.na(breeding_abundance),
         breeding_abundance > 0)




continents_abund <- data.frame(region = as.character(continents_proj$CC[abundance_in_continents$ID]),
                           region_type = "continent",
                           sum_abund = 3^2 * unlist(abundance_in_continents[[2]])) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#continents_abund
# will be combined with other regional estimates below



abundance_in_global <- sum(values(breed_abundance_full),na.rm = TRUE)



global_abund <- data.frame(region = "global",
                           region_type = "global",
                           sum_abund = 3^2 * abundance_in_global) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#global_abund
# will be combined with other regional estimates in next line


pop_ests_out <- bind_rows(USACAN_abund,
                          continents_abund,
                          global_abund,
                          country_abund,
                          state_abund,
                          bcr_abund,
                          strata_abund) |>
  mutate(version = vers)

# Sum of population estimates

# write_excel_csv(pop_ests_out,
#                 paste0("estimates/pop_ests_alt_",vers,"_",sp_aou,sp_ebird,".csv"))

pop_ests_out_save <- bind_rows(pop_ests_out,
                               pop_ests_out_save)
#
# pop_ests_out_trad_sel <- pop_ests_out_trad %>%
#   filter(species == sp_sel)
#
#
# pop_compare_stack <- pop_ests_out %>%
#   select(region,region_type,model,
#          species,species_ebird,
#          pop_median, pop_lci_80, pop_uci_80,
#          pop_lci_95, pop_uci_95) %>%
#   bind_rows(pop_ests_out_trad_sel)
#
#
# pop_compare_stack_sel <- pop_compare_stack %>%
#   filter(region_type %in% c("USACAN","continent","global","bcr"))
#
#
# side_plot <- ggplot(data = pop_compare_stack_sel,
#               aes(y = region,
#                   x = pop_median,
#                   colour = model))+
#   geom_point(position = position_dodge(width = 0.5))+
#   geom_errorbarh(aes(xmin = pop_lci_80,
#                      xmax = pop_uci_80),
#                  height = 0,
#                  position = position_dodge(width = 0.5))+
#   scale_colour_viridis_d(end = 0.8)+
#   ylab("BCRs and broadscale estimates")+
#   xlab("Population estimate with 80% CI (Millions)")+
#   scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
#   theme_bw()
#
#
# saveRDS(side_plot,paste0("figures/saved_ggplots/side_plot_alt_",vers,"_",sp_aou,"_",sp_ebird,".rds"))
#

seas_plot <- FALSE
if(seas_plot){
# Exploring the estimated seasonal patterns by strata ------------------------------------------
strat_names <- combined %>%
  select(strata,strata_name,day_of_year) %>%
  distinct()

seasonal_strat <- summ %>% filter(grepl("doy_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = rep(c(1:stan_data$n_doy),times = stan_data$n_strata) + (min(combined$day_of_year)-1),
         strata = rep(c(1:stan_data$n_strata),each = stan_data$n_doy)) %>%
  inner_join(strat_names,by = c("strata",
                                "doy" = "day_of_year"))

seasonal <- summ %>% filter(grepl("DOY_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = c(1:stan_data$n_doy)+ (min(combined$day_of_year)-1),
         strata_name = "USA_CAN") %>%
  filter(doy %in% strat_names$day_of_year)

strat_labs1 <- seasonal_strat %>%
  group_by(strata_name) %>%
  summarise(doy = max(doy))

strat_labs <- strat_labs1 %>%
  inner_join(seasonal_strat,
             by = c("strata_name",
                    "doy"))

xlbl <- data.frame(dtes = c("2025/05/01",
                            "2025/05/15",
                             "2025/06/01",
                             "2025/06/15",
                             "2025/07/01"),
                   days = c("May-01",
                            "May-15",
                            "June-01",
                            "June-15",
                            "July-01")) %>%
  mutate(doy = yday(dtes))


vis_season <- ggplot(data = seasonal_strat,
                     aes(x = doy,y = median,
                         group = strata_name,
                         colour = strata_name))+
  geom_ribbon(data = seasonal,aes(x = doy,y = median,
                                  ymin = q5, ymax = q95),
              alpha = 0.1,
              inherit.aes = FALSE)+
  geom_line(data = seasonal,aes(x = doy,y = median),
            inherit.aes = FALSE)+
  geom_line(alpha = 0.5)+
  coord_cartesian(xlim = c(100,240))+
  ggrepel::geom_text_repel(data = strat_labs, aes(label = strata_name),
                           size = 1.5,
                           min.segment.length = 0,
                           xlim = c(200,240),
                           nudge_x = 20)+
  scale_colour_viridis_d()+
  scale_x_continuous(breaks = xlbl$doy, labels = xlbl$days)+
  xlab("")+
  ylab("Seasonal effect (log-scale)")+
  scale_y_continuous(transform = "exp")+
  theme_bw()+
  theme(legend.position = "none")


#saveRDS(seasonal_strat,paste0("output/seasonal_strat_",vers,"_",sp_aou,"_",sp_ebird,".rds"))

#saveRDS(vis_season,paste0("figures/saved_ggplots/season_",vers,"_",sp_aou,"_",sp_ebird,".rds"))
png(filename = paste0("Figures/seasonal_effect_",vers,"_",sp_aou,"_",sp_ebird,".png"),
    res = 300,
    height = 6,
    width = 6,
    units = "in")
print(vis_season)
dev.off()

}# end season plot


  saveRDS(summ_out,paste0("estimates/convergence_population_estimates_4_models_",sp_aou,".rds"))

  saveRDS(trajectories_out_all,paste0("estimates/trajectories_population_estimates_4_models_",sp_aou,".rds"))
  saveRDS(pop_ests_out_save,paste0("estimates/combined_population_estimates_4_models_",sp_aou,".rds"))
} # end species loop

stopCluster(cluster)




# End species loop --------------------------------------------------------
















#summ_out <- NULL
pop_ests_out_save <- NULL
trajectories_out_save <- NULL

# sp_sel <- sps_list[i,"common_name"]
#
# sp_aou <- sps_list[i,"sp_aou"]
# sp_ebird <- sps_list[i,"sp_ebird"]

for(i in 1:nrow(sps_list)){
  sp_aou <- sps_list[i,"sp_aou"]

  if(file.exists(paste0("estimates/combined_population_estimates_4_models_",sp_aou,".rds"))){


    pop_ests_out <- readRDS(paste0("estimates/combined_population_estimates_4_models_",sp_aou,".rds"))

  pop_ests_out_save <- bind_rows(pop_ests_out_save,
                                 pop_ests_out)

  trajectories_out <- readRDS(paste0("estimates/trajectories_population_estimates_4_models_",sp_aou,".rds"))

  trajectories_out_save <- bind_rows(trajectories_out_save,
                                 trajectories_out)


  }
}

pop_ests_out_save <- pop_ests_out_save |>
  mutate(version = factor(version, levels = c("NA-POPS edr NA-POPS Availability",
                                              "NA-POPS edr PIF ToD",
                                              "approximate edr NA-POPS Availability",
                                              "approximate edr PIF ToD" ),
                          ordered = TRUE))

tst_plot_data <- pop_ests_out_save |>
    filter(region %in% c("USACAN")) |>
  mutate(species = fct_reorder(species,
                                     .x = pop_median)) |>
  arrange(species)

  tst_plot1 <- ggplot()+
    geom_pointinterval(data = filter(tst_plot_data,pop_median > 15e6),
                       aes(x = species, y = pop_median,
                           ymin = pop_lci_95, ymax = pop_uci_95,
                           colour = version),
                       position = position_dodge(width = 0.5),
                       linewidth = 0.1)+
    scale_y_continuous(transform = "log10",
                       breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                       labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
    ylab("Population estimate Canada and USA")+
    theme(axis.text = element_text(size = 8),
          line = element_line(lineend = "round"))+
    theme_bw()+
    scale_colour_viridis_d(option = "viridis")+
    coord_flip()


  tst_plot2 <- ggplot()+
    geom_pointinterval(data = filter(tst_plot_data,pop_median <= 15e6 & pop_median > 1.2e6),
                       aes(x = species, y = pop_median,
                           ymin = pop_lci_95, ymax = pop_uci_95,
                           colour = version),
                       position = position_dodge(width = 0.5),
                       linewidth = 0.1)+
    scale_y_continuous(transform = "log10",
                       breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                       labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
    ylab("Population estimate Canada and USA")+
    theme(axis.text = element_text(size = 8),
          line = element_line(lineend = "round"))+
    theme_bw()+
    scale_colour_viridis_d(option = "viridis")+
    coord_flip()


  tst_plot3 <- ggplot()+
    geom_pointinterval(data = filter(tst_plot_data,pop_median <= 1.2e6),
                       aes(x = species, y = pop_median,
                           ymin = pop_lci_95, ymax = pop_uci_95,
                           colour = version),
                       position = position_dodge(width = 0.5),
                       linewidth = 0.1)+
    scale_y_continuous(transform = "log10",
                       breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                       labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
    ylab("Population estimate Canada and USA")+
    theme(axis.text = element_text(size = 8),
          line = element_line(lineend = "round"))+
    theme_bw()+
    scale_colour_viridis_d(option = "viridis")+
    coord_flip()




pdf(paste0("Population_estimates_",nrow(tst_plot_data),"_species.pdf"),
    width = 11,
    height = 17)
print(tst_plot1)
print(tst_plot2)
print(tst_plot3)
dev.off()






tst_plot_data <- pop_ests_out_save |>
  filter(region %in% c("Canada")) |>
  mutate(species = fct_reorder(species,
                               .x = pop_median)) |>
  arrange(species)

tst_plot1 <- ggplot()+
  geom_pointinterval(data = filter(tst_plot_data,pop_median > 15e6),
                     aes(x = species, y = pop_median,
                         ymin = pop_lci_95, ymax = pop_uci_95,
                         colour = version),
                     position = position_dodge(width = 0.5),
                     linewidth = 0.1)+
  scale_y_continuous(transform = "log10",
                     breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                     labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
  ylab("Population estimate Canada")+
  theme(axis.text = element_text(size = 8),
        line = element_line(lineend = "round"))+
  theme_bw()+
  scale_colour_viridis_d(option = "viridis")+
  coord_flip()


tst_plot2 <- ggplot()+
  geom_pointinterval(data = filter(tst_plot_data,pop_median <= 15e6 & pop_median > 1.2e6),
                     aes(x = species, y = pop_median,
                         ymin = pop_lci_95, ymax = pop_uci_95,
                         colour = version),
                     position = position_dodge(width = 0.5),
                     linewidth = 0.1)+
  scale_y_continuous(transform = "log10",
                     breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                     labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
  ylab("Population estimate Canada")+
  theme(axis.text = element_text(size = 8),
        line = element_line(lineend = "round"))+
  theme_bw()+
  scale_colour_viridis_d(option = "viridis")+
  coord_flip()


tst_plot3 <- ggplot()+
  geom_pointinterval(data = filter(tst_plot_data,pop_median <= 1.2e6),
                     aes(x = species, y = pop_median,
                         ymin = pop_lci_95, ymax = pop_uci_95,
                         colour = version),
                     position = position_dodge(width = 0.5),
                     linewidth = 0.1)+
  scale_y_continuous(transform = "log10",
                     breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                     labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
  ylab("Population estimate Canada")+
  theme(axis.text = element_text(size = 8),
        line = element_line(lineend = "round"))+
  theme_bw()+
  scale_colour_viridis_d(option = "viridis")+
  coord_flip()




pdf(paste0("Population_estimates_Canada_",nrow(tst_plot_data),"_species.pdf"),
    width = 11,
    height = 17)
print(tst_plot1)
print(tst_plot2)
print(tst_plot3)
dev.off()






tst_plot_data <- pop_ests_out_save |>
  filter(region %in% c("United States of America")) |>
  mutate(species = fct_reorder(species,
                               .x = pop_median)) |>
  arrange(species)

tst_plot1 <- ggplot()+
  geom_pointinterval(data = filter(tst_plot_data,pop_median > 15e6),
                     aes(x = species, y = pop_median,
                         ymin = pop_lci_95, ymax = pop_uci_95,
                         colour = version),
                     position = position_dodge(width = 0.5),
                     linewidth = 0.1)+
  scale_y_continuous(transform = "log10",
                     breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                     labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
  ylab("Population estimate USA")+
  theme(axis.text = element_text(size = 8),
        line = element_line(lineend = "round"))+
  theme_bw()+
  scale_colour_viridis_d(option = "viridis")+
  coord_flip()


tst_plot2 <- ggplot()+
  geom_pointinterval(data = filter(tst_plot_data,pop_median <= 15e6 & pop_median > 3e6),
                     aes(x = species, y = pop_median,
                         ymin = pop_lci_95, ymax = pop_uci_95,
                         colour = version),
                     position = position_dodge(width = 0.5),
                     linewidth = 0.1)+
  scale_y_continuous(transform = "log10",
                     breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                     labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
  ylab("Population estimate USA")+
  theme(axis.text = element_text(size = 8),
        line = element_line(lineend = "round"))+
  theme_bw()+
  scale_colour_viridis_d(option = "viridis")+
  coord_flip()


tst_plot3 <- ggplot()+
  geom_pointinterval(data = filter(tst_plot_data,pop_median <= 3e6),
                     aes(x = species, y = pop_median,
                         ymin = pop_lci_95, ymax = pop_uci_95,
                         colour = version),
                     position = position_dodge(width = 0.5),
                     linewidth = 0.1)+
  scale_y_continuous(transform = "log10",
                     breaks = c(1e3,1e4,1e5,1e6,1e7,1e8,1e9),
                     labels = c(paste(c(1e3,1e4,1e5,1e6,1e7,1e8)/1e6,"Million"),"1 Billion"))+
  ylab("Population estimate USA")+
  theme(axis.text = element_text(size = 8),
        line = element_line(lineend = "round"))+
  theme_bw()+
  scale_colour_viridis_d(option = "viridis")+
  coord_flip()




pdf(paste0("Population_estimates_United_States_",nrow(tst_plot_data),"_species.pdf"),
    width = 11,
    height = 17)
print(tst_plot1)
print(tst_plot2)
print(tst_plot3)
dev.off()


# Format population estimates and trajectories ----------------------------

signif_M <- function(x,sig = 3){
  signif(x,sig)/1e6
}

pop_ests_out_save <- pop_ests_out_save |>
  select(-c(sum_abund,strata_name,mean_ebird_abundance,total_ebird_abundance,sampled_abundance,
            sampled_avail_ratio,log_ratio,pop_mean)) |>
  mutate(across(starts_with("pop"),
                .fns = ~signif_M(.x))) |>
  rename_with(.fn = ~str_replace(.x,
                                 "pop_",
                                 "population_estimate_millions_"))



write_csv(pop_ests_out_save,
          "all_population_estimates.csv")
write_csv(trajectories_out_save,
          "all_population_estimate_trajectories.csv")


tmp2 <- pop_ests_out_save |>
  filter(species == "Gray-cheeked Thrush")




# NOTE: Friday something is wrong that are multiple (4) estimates -------------
# for each region in ped....

ped <- pop_ests_out_save  |>
  mutate(model = ifelse(model == "no_napops",
                        "PIF TOD and MDD",
                        model)) |>
  select(species,
         species_ebird,
         region,
         region_type,
         model,
         pop_median,
         pop_lci_80,
         pop_uci_80,
         pop_lci_95,
         pop_uci_95) |>
  mutate(across(starts_with("pop"),
                .fns = ~signif_M(.x)))


trajs <- trajectories_out_all |>
  select(species,
         #species_ebird,
         region,
         region_type,
         model,
         year,
         pop_median,
         pop_lci_80,
         pop_uci_80,
         pop_lci_95,
         pop_uci_95) |>
  mutate(across(starts_with("pop"),
                .fns = ~signif_M(.x)))

write_csv(ped,"ped_new.csv")


ped <- tst_plot_data |>
  select(species,
         species_ebird,
         region,
         region_type,
         model,
         pop_median,
         pop_lci_80,
         pop_uci_80,
         pop_lci_95,
         pop_uci_95) |>
  mutate(across(starts_with("pop"),
                .fns = ~signif_M(.x)))


write_csv(ped,"some_population_estimates.csv")


write_csv(trajs,"annual_peds_new.csv")



