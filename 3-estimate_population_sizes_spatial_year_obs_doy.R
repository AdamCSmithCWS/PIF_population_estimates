
library(sf)
library(tidyverse)
library(ebirdst)
library(patchwork)
library(terra)
library(bbsBayes2)
library(cmdstanr)
library(ggrepel)
library(tidyterra)
library(napops)
library(tidybayes)

source("functions/GAM_basis_function_mgcv.R")
source("functions/neighbours_define.R")
yr_ebird <- 2022 # prediction year for eBird relative abundance


## Example species

#	American Robin (widespread, generalist, with potential to for estimate to decrease based on eBird model calibration)
#	Western Meadowlark or Clay-colored Sparrow (or other … intent is to look at a widespread grassland species that has also been addressed in Barry’s recent work)
#	Sharp-tailed Grouse (or other species with a patch distribution)
#	Sora (or other waterbird with good BBS coverage)
#	Killdeer or Long-billed Curlew (shorebirds with good BBS coverage but lower previous population estimates based on shorebird analysis)
#	Yellow-rumped Warbler or Olive-sided Flycatcher (or other widespread boreal breeder with only moderate BBS coverage)
#	Purple Martin (or other species with substantial regional trend differences)
#	Barn Swallow (or other largely visually-detected species)
#	Verdin or Ash-throated Flycatcher or Black-throated Sparrow (to represent species with substantial breeding populations in both US and Mexico, to explore potential to extend estimates into Mexico)

sp_example <- c("American Robin",
                "Western Meadowlark","Clay-colored Sparrow",
                "Sharp-tailed Grouse",
                "Sora",
                "Killdeer","Long-billed Curlew",
                "Yellow-rumped Warbler","Olive-sided Flycatcher",
                "Purple Martin",
                "Barn Swallow",
                "Verdin","Ash-throated Flycatcher","Black-throated Sparrow",
                "Blue Jay","Varied Thrush","Veery","Wood Thrush","Chestnut-collared Longspur",
                "Bobolink","Savannah Sparrow","Grasshopper Sparrow",
                "Horned Lark","Baird's Sparrow",
                "Eastern Whip-poor-will", "Common Nighthawk")

# Load BBS data -----------------------------------------------------------


# mean counts to id routes on which species has been observed -------------
mean_counts_all <- readRDS("data/mean_counts_by_route.rds")

raw_counts_all <- readRDS("data/all_counts_by_route.rds")

# general summaries of time and season
surveys <- load_bbs_data(release = 2024)$routes %>%
  filter(year > 2012)  %>%  # 10 years from 2013 - 2023, missing 2020.
  mutate(doy = lubridate::yday(paste(year,month,day,sep = "-")),
         start_hr = floor(start_time/100),
         end_hr = floor(end_time/100),
         start_min = start_hr*60 + ((start_time/100)-start_hr)*100,
         end_min = end_hr*60 + ((end_time/100)-end_hr)*100,
         mid_time = ((end_min - start_min)/120)-0.5) %>%
  summarise(mean_time = mean(mid_time),
            mean_doy = mean(doy))

survey_dates <- load_bbs_data(release = 2024)$routes %>%
  filter(year > 2012)  %>%  # 10 years from 2013 - 2023, missing 2020.
  mutate(doy = lubridate::yday(paste(year,month,day,sep = "-")),
         start_hr = floor(start_time/100),
         end_hr = floor(end_time/100),
         start_min = start_hr*60 + ((start_time/100)-start_hr)*100,
         end_min = end_hr*60 + ((end_time/100)-end_hr)*100,
         mid_time = ((end_min - start_min)/120)-0.5)

mid_time = round(surveys$mean_time,1)
mean_doy = round(surveys$mean_doy)




# Add na-pops EDRs --------------------------------------------------------


Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")
napops_species <- napops::list_species() |>
  dplyr::select(Species,Common_Name)

re_napops <- FALSE # set to True to re-run the napops extraction

if(re_napops){
# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE)



ExpAdjs <- ExpAdjs |>
  dplyr::left_join(napops_species,
                   by = c("cn" = "Common_Name"))

for(i in 1:nrow(ExpAdjs)){
  if(is.na(ExpAdjs[i,"Species"])){
    print(paste("no napops edr for",ExpAdjs[i,"cn"]))
          next}

  sp <- ExpAdjs[i,"Species"]

  # choose better supported model between models 1 (intercept) and 2 (roadside)
  cefs <- napops::coef_distance(species = sp) |>
    dplyr::filter(Model %in% c(1,2)) |>
    dplyr::arrange(AIC)

  if(nrow(cefs) == 0){next} # skip if no napops output

  w_mod <- cefs[1,"Model"]

  edrt <- try(napops::edr(species = sp,
                  model = w_mod,
                  forest = 0.5,
                  road = TRUE,
                  quantiles = c(0.1625), # 1 sd below the mean
                  samples = 1000),
              silent = TRUE)

  if(edrt$EDR_est > 500){ # two species where the model 2 estimates of edr are silly (> 5 km!)
    edrt <- try(napops::edr(species = sp,
                    model = 1,
                    forest = 0.5,
                    road = TRUE,
                    quantiles = c(0.1625),# 1 sd below the mean
                    samples = 1000),
                silent = TRUE)
    w_mod <- 1
  }

  if(class(edrt) == "try-error"){
    warning(paste("edr extraction did not work for",sp))
    next
  }


  ExpAdjs[i,"edr"] <- as.numeric(edrt$EDR_est)

  if("EDR_16.25" %in% names(edrt)){
    ExpAdjs[i,"edr_sd"] <- as.numeric(edrt$EDR_est)-as.numeric(edrt$EDR_16.25)
  }

  ExpAdjs[i,"EDR_model"] <- ifelse(w_mod == 1,"Intercept","Roadside")

  rm(edrt, w_mod)
}

ExpAdjs <- ExpAdjs |>
  dplyr::mutate(edr = ifelse(edr == 1,NA,edr),
                use_edr = ifelse(is.na(edr),FALSE,TRUE))




# Availability from napops ------------------------------------------------



for(i in 1:nrow(ExpAdjs)){
  if(is.na(ExpAdjs[i,"Species"])){
    print(paste("no napops avail for",ExpAdjs[i,"cn"]))
    next}

  sp <- ExpAdjs[i,"Species"]

  # choose best supported model
  cefs <- napops::coef_removal(species = sp) |>
    dplyr::arrange(AIC)

  if(nrow(cefs) == 0){next} # skip if no napops output

  w_mod <- cefs[1,"Model"]

  availability <- try(napops::avail(species = sp,
                           time = 3,
                          model = w_mod,
                          od = mean_doy, # mean of the doy for all BBS surveys
                          tssr = mid_time, #mean of time of day for all BBS surveys
                          quantiles = c(0.1625), # 1 sd below the mean
                          samples = 1000),
              silent = TRUE)


  if(class(availability) == "try-error"){
    warning(paste("availability extraction did not work for",sp))
    next
  }


  ExpAdjs[i,"availability"] <- as.numeric(availability$p_est)

  if("p_16.25" %in% names(availability)){
    ExpAdjs[i,"availability_sd"] <- as.numeric(availability$p_est)-as.numeric(availability$p_16.25)
  }

  ExpAdjs[i,"availability_model"] <- w_mod

  rm(availability, w_mod)
}

ExpAdjs <- ExpAdjs |>
  dplyr::mutate(use_availability = ifelse(is.na(availability),FALSE,TRUE))


write_csv(ExpAdjs,"Species_correction_factors_w_edr_availability.csv")
}else{
  ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")
}

# ID species missing edrs or availability-rates
corrections_example_sp <- ExpAdjs %>%
  filter(cn %in% sp_example)

if(any(!corrections_example_sp$use_availability) |
   any(!corrections_example_sp$use_edr)){
  warning("edr or availability missing for some species")
}

routes_buf_all <- readRDS("data/all_routes_buffered.rds")

strata <- bbsBayes2::load_map("bbs_usgs")

strata_names <- strata %>%
  sf::st_drop_geometry()




# potentially important detail from Golden Eagle paper --------------------------------------
##  That paper USED MAX WEEKLY RELATIVE ABUNDANCE DURING KEY SEASON
# Using the 2019 relative abundance data product as a starting point, we calculated
# the maximum eBird relative abundance value for each 3x3-km pixel over the weeks
# of the aerial survey (August 17th–September 14th) because these units were
# expected to best match those of the golden eagle survey.



# load traditional estimates for later comparison -------------------------
strata_join <- strata_names %>%
  select(strata_name,country_code,
         prov_state, bcr) %>%
  st_drop_geometry() %>%
  distinct()

strata_trad_all <- read_csv("Stanton_2019_code/output/PSest_Prov_by_BCR__100iter.csv") %>%
  left_join(.,strata_join,
            by = c("st_abrev" = "prov_state",
                   "BCR" = "bcr")) %>%
  mutate(region = strata_name,
         region_type = "strata",
         version = "traditional",
         species = cn,
         aou = Aou,
         pop_median = med.PopEst,
         pop_lci_80 = LCI80.PopEst,
         pop_uci_80 = UCI80.PopEst,
         pop_lci_95 = LCI95.PopEst,
         pop_uci_95 = UCI95.PopEst)

USACAN_trad_all <- read_csv("Stanton_2019_code/output/PSest_Global_WH_NA__100iter.csv") %>%
  mutate(region = "USACAN",
         region_type = "USACAN",
         version = "traditional",
         species = cn,
         aou = Aou,
         pop_median = USCAN.med.PopEst,
         pop_lci_80 = USCAN.80LCI.PopEst,
         pop_uci_80 = USCAN.80UCI.PopEst,
         pop_lci_95 = USCAN.95LCI.PopEst,
         pop_uci_95 = USCAN.95UCI.PopEst)



global_trad_all <- read_csv("Stanton_2019_code/output/PSest_Global_WH_NA__100iter.csv") %>%
  mutate(region = "global",
         region_type = "global",
         version = "traditional",
         species = cn,
         aou = Aou,
         pop_median = GL.med.PopEst,
         pop_lci_80 = GL.80LCI.PopEst,
         pop_uci_80 = GL.80UCI.PopEst,
         pop_lci_95 = GL.95LCI.PopEst,
         pop_uci_95 = GL.95UCI.PopEst)



prov_trad_all <- read_csv("Stanton_2019_code/output/PSest_by_Prov__100iter.csv") %>%
  mutate(region = st_abrev,
         region_type = "prov_state",
         version = "traditional",
         species = cn,
         aou = Aou,
         pop_median = med.PopEst,
         pop_lci_80 = LCI80.PopEst,
         pop_uci_80 = UCI80.PopEst,
         pop_lci_95 = LCI95.PopEst,
         pop_uci_95 = UCI95.PopEst)

bcr_trad_all <- read_csv("Stanton_2019_code/output/PSest_by_BCR__100iter.csv")%>%
  mutate(region = as.character(BCR),
         region_type = "bcr",
         version = "traditional",
         species = cn,
         aou = Aou,
         pop_median = med.PopEst,
         pop_lci_80 = LCI80.PopEst,
         pop_uci_80 = UCI80.PopEst,
         pop_lci_95 = LCI95.PopEst,
         pop_uci_95 = UCI95.PopEst)




pop_ests_out_trad <- bind_rows(USACAN_trad_all,
                               global_trad_all,
                          #country_abund,
                          bcr_trad_all,
                          strata_trad_all) %>%
  select(region,region_type,version,
         species,aou,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95)


saveRDS(pop_ests_out_trad,"all_traditional_pop_estimates.rds")







# Running full model using eBird relative abundance ------------------------------------------------

#
sps_list <- readRDS("data/all_bbs_data_species_list.rds") #

sps_sel <- c("Tennessee Warbler","Bay-breasted Warbler",
             "Swainson's Thrush","Blue Jay",
             "Blue-headed Vireo","Bicknell's Thrush",
             "Wilson's Snipe","American Kestrel",
             "Steller's Jay")
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
log_calibration <- FALSE # conduct calibration on log transformed abundance

# base map for countries and continents -----------------------------------
country_codes <- readxl::read_xlsx("data/iso_codes.xlsx") %>%
  rename(country_name = Name) #

americas_codes <- country_codes %>%
  filter(CC %in% c("SA","NA"))

states <- rnaturalearthhires::states10 %>%
  left_join(country_codes,by = c("adm0_a3" = "a-3")) %>%
  select(CC,iso_a2,adm0_a3,country_name,name_en,postal,
         woe_name,gn_name)

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
bcrs <- read_sf("data/BCR_terrestrial/BCR_Terrestrial_master_International.shp") %>%
  group_by(BCR,Label) %>%
  summarise(do_union = FALSE)

adjs_out <- NULL
today <- as_date(Sys.Date())





# species loop -----------------------------------------------------
use_traditional <- FALSE
if(use_traditional){
  vers <- "trad_"
}else{
  vers <- ""
}

if(file.exists(paste0("adjs_out",vers,today,".csv"))){
  adjs_out <- read_csv(paste0("adjs_out",vers,today,".csv"))
  wh_drop_already <- which(sp_example %in% adjs_out$cn)
  wh_drop <- c(wh_drop_already,8)
}else{
  wh_drop <- 8
}


re_fit <- TRUE # set to true to fit model and re-run plotting and summaries

output_dir <- "D:/PIF_pop_estimates/output"

for(sp_sel in rev(sps_list$english[1:348])){#sp_example[-wh_drop]){#list$english){
#sp_sel = "Connecticut Warbler"
 sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)


  if(file.exists(paste0(output_dir,"/calibration_fit_alt_",
                        vers,sp_aou,"_",sp_ebird,".rds")) &
     !re_fit){next}


  mean_counts <- mean_counts_all %>%
    filter(english == sp_sel,
           mean_count > 0)

  raw_counts <- raw_counts_all %>%
    filter(grepl(paste0("^",sp_sel),english))

  nsp_names <- length(unique(raw_counts$english))

  if(nsp_names > 1){
    warning(paste(sp_sel,"has too many name matches"))
    next}
  if(nrow(raw_counts) == 0){
    warning(paste(sp_sel,"has no name matches"))
    next
  }

strata_trad <- strata_trad_all %>%
  filter(cn == sp_sel)

global_trad <- global_trad_all %>%
  filter(cn == sp_sel)

bcr_trad <- bcr_trad_all %>%
  filter(cn == sp_sel)





# filtering routes on which species wasn't observed.
routes_buf <- routes_buf_all %>%
  filter(route_name %in% mean_counts$route_name)

mn_c_join <- mean_counts %>%
  select(route_name,mean_count)

routes_buf <- routes_buf %>%
  left_join(mn_c_join, by = "route_name")



if(!file.exists(paste0("data/species_relative_abundance/",
                       sp_ebird,"_relative_abundance.rds"))){next}

rel_abund <- readRDS(paste0("data/species_relative_abundance/",
                            sp_ebird,"_relative_abundance.rds"))




combined <- rel_abund %>%
inner_join(.,raw_counts,
           by = "route_name") %>%
  mutate(logm = log(ebird_abund),
         day_of_year = doy,
         route = as.integer(factor(route_name)),
         prov_state = str_extract(pattern = ".+(?=[[:punct:]])",string = route_name),
         yr = year-(min(year)-1),
         observer = as.integer(factor(obs_n)),
         strata = as.integer(factor(strata_name)),
         doy = 1+(day_of_year-(min(day_of_year))),
         route_obs = as.integer(factor(paste(route_name,obs_n,sep = "-")))) # doy centered on earliest survey day for species


# spatial setup -----------------------------------------------------------

strata_incl <- combined %>%
  select(strata_name,strata) %>%
  distinct()

strata_map <- bbsBayes2::load_map("bbs_usgs") %>%
  inner_join(strata_incl) %>%
  arrange(strata)

neighbours <- neighbours_define(strata_map,
                                strat_indicator = "strata",
                                save_plot_data = FALSE)


all_years <- data.frame(year = c(min(combined$year):max(combined$year)))

y_2020 <- combined %>%
  select(year,yr) %>%
  distinct() %>%
  full_join(all_years, by = "year") %>%
  arrange(year) %>%
  mutate(y_2020 = ifelse(is.na(yr),0,1))



mean_abund <- combined %>%
  select(route_obs,logm) %>%
  distinct() %>%
  arrange(.,route_obs)


#adjustment scores
adjs <- ExpAdjs %>%
  filter(cn == sp_sel)
if(nrow(adjs) == 0){next}
if(!use_traditional & !(adjs$use_edr | adjs$use_availability)){next}

#library(cmdstanr)
#
#
# Fit stan model ----------------------------------------------------------

# setup for GAM of doy ----------------------------------------------------


doy_gam <- gam_basis(combined$doy, nknots = 7,
                     predpoints = c(1:max(combined$doy)),
                     sm_name = "doy")

model <- cmdstanr::cmdstan_model("models/ebird_rel_abund_calibration_spatial_year_doy_obs.stan")


stan_data <- list(n_route_obs = max(combined$route_obs),
                  n_counts = nrow(combined),
                  n_years = max(combined$yr),
                  n_knots_doy = 7,
                  n_doy = max(combined$doy),
                  mean_doy = as.integer(round(mean(combined$doy))),
                  n_strata = max(combined$strata),

                  count = combined$count,
                  route = combined$route_obs,
                  year = combined$yr,
                  strata = combined$strata,
                  doy = combined$doy,
                  doy_basis = doy_gam$doy_basispred,
                  ebird_year = yr_ebird-(min(combined$year)-1),
                  yrev = seq(from = (yr_ebird-(min(combined$year))),to = 1, by = -1),
                  log_mean_rel_abund = mean_abund$logm,
                  y_2020 = as.integer(unname(unlist(y_2020$y_2020))),
                  zero_gammas = rep(0,max(combined$strata)),
                  n_edges = neighbours$N_edges,
                  node1 = neighbours$node1,
                  node2 = neighbours$node2,

                  # PIF adjustments
                  c_p = adjs$Pair2,
                  sd_c_p = 0.13,
                  c_d = ifelse(!use_traditional & adjs$use_edr,adjs$edr,adjs$Dist2),
                  sd_c_d = ifelse(!use_traditional & adjs$use_edr,adjs$edr_sd,1),
                  c_t = ifelse(!use_traditional & adjs$use_availability,
                               adjs$availability,
                               adjs$TimeAdj.meanlog),
                  sd_c_t = ifelse(!use_traditional & adjs$use_availability,
                                  adjs$availability_sd,
                                  adjs$TimeAdj.sdlog),
                  c_d_lower = adjs$Dist.Lower,
                  c_d_upper = adjs$Dist.Upper,

                  # Conditional statments
                  use_edr = ifelse(!use_traditional & adjs$use_edr,1,0),
                  use_availability = ifelse(!use_traditional & adjs$use_availability,1,0),
                  use_pois = 0,
                  use_pair = 1,
                  use_t = 1
                  )


if(nrow(adjs) == 0){
  # stan_data$c_p = mean(ExpAdjs$Pair2,na.rm = T)
  # stan_data$c_t = exp(mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + 0.5*(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2))
  # stan_data$sd_c_t = exp(2*mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)*(exp(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)-1)
  # stan_data$c_d_lower = mean(ExpAdjs$Dist.Lower,na.rm = T)
  # stan_data$c_d_upper = mean(ExpAdjs$Dist.Upper,na.rm = T)

}
adjs$n_routes <- stan_data$n_routes


fit <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500,
                    iter_warmup = 2000,
                    iter_sampling = 2000,
                    adapt_delta = 0.8,
                    max_treedepth = 11,
                    show_exceptions = FALSE)



params_to_summarise <- c("nu",
                         "BETA",
                         "beta_raw",
                         "beta",
                         "sd_beta",
                         "sdnoise",
                         "sd_gamma",
                         "GAMMA",
                         "sd_GAMMA",
                         "gamma",
                         "yeareffect",
                         "YearEffect",
                         "phi",
                         "cp",
                         "cp_sel",
                         "ct",
                         "cd",
                         "p_avail",
                         "c_area",
                         "calibration",
                         "calibration_alt",
                         "calibration_median",
                         "calibration_log",
                         "calibration_r",
                         "adj",
                         "raw_prediction",
                         "sd_doy",
                         "sd_DOY",
                         "doy_raw",
                         "DOY_raw",
                         "DOY_b",
                         "doy_b",
                         "DOY_pred",
                         "doy_pred")

summ <- fit$summary(variables = params_to_summarise)
#shinystan::launch_shinystan(fit)


fit$save_object(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

saveRDS(summ,paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))





avail_correction_realised <- fit$summary(variable = "ct",
                                "mean",
                                "median",
                                "sd",
                                "rhat",
                                "ess_bulk",
                                q2_5 = ~quant(.x,q = 0.025),
                                q10 = ~quant(.x,q = 0.10),
                                q90 = ~quant(.x,q = 0.90),
                                q97_5 = ~quant(.x,q = 0.975)) %>%
  mutate(inference = "realised correction for time or availability")

p_avail_realised <- fit$summary(variable = "p_avail",
                                "mean",
                                "median",
                                "sd",
                                "rhat",
                                "ess_bulk",
                                q2_5 = ~quant(.x,q = 0.025),
                                q10 = ~quant(.x,q = 0.10),
                                q90 = ~quant(.x,q = 0.90),
                                q97_5 = ~quant(.x,q = 0.975)) %>%
  mutate(inference = "realised p availability 3-minutes")



edr_realised <- fit$summary(variable = "cd",
                            "mean",
                            "median",
                            "sd",
                            "rhat",
                            "ess_bulk",
                            q2_5 = ~quant(.x,q = 0.025),
                            q10 = ~quant(.x,q = 0.10),
                            q90 = ~quant(.x,q = 0.90),
                            q97_5 = ~quant(.x,q = 0.975)) %>%
  mutate(inference = "realised EDR")



surveyed_area_realised <- fit$summary(variable = "c_area",
                                      "mean",
                                      "median",
                                      "sd",
                                      "rhat",
                                      "ess_bulk",
                                      q2_5 = ~quant(.x,q = 0.025),
                                      q10 = ~quant(.x,q = 0.10),
                                      q90 = ~quant(.x,q = 0.90),
                                      q97_5 = ~quant(.x,q = 0.975)) %>%
  mutate(inference = "realised surveyed area per BBS route")


cali <- fit$summary(variable = "calibration",
                    "mean",
                    "median",
                    "sd",
                    "rhat",
                    "ess_bulk",
                    q2_5 = ~quant(.x,q = 0.025),
                    q10 = ~quant(.x,q = 0.10),
                    q90 = ~quant(.x,q = 0.90),
                    q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration assuming log normal route distribution")


cali_alt <- fit$summary(variable = "calibration_alt",
                    "mean",
                    "median",
                    "sd",
                    "rhat",
                    "ess_bulk",
                    q2_5 = ~quant(.x,q = 0.025),
                    q10 = ~quant(.x,q = 0.10),
                    q90 = ~quant(.x,q = 0.90),
                    q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration realised mean across groups")


cali_log <- fit$summary(variable = "calibration_log",
                        "mean",
                        "median",
                        "sd",
                        "rhat",
                        "ess_bulk",
                        q2_5 = ~quant(.x,q = 0.025),
                        q10 = ~quant(.x,q = 0.10),
                        q90 = ~quant(.x,q = 0.90),
                        q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration log-scale")


cali_median <- fit$summary(variable = "calibration_median",
                        "mean",
                        "median",
                        "sd",
                        "rhat",
                        "ess_bulk",
                        q2_5 = ~quant(.x,q = 0.025),
                        q10 = ~quant(.x,q = 0.10),
                        q90 = ~quant(.x,q = 0.90),
                        q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration realised median across groups")



param_infer <- bind_rows(cali_alt,
                         cali,
                         cali_log,
                         cali_median,
                         avail_correction_realised,
                         p_avail_realised,
                         edr_realised,
                         surveyed_area_realised) %>%
  mutate(english = sp_sel,
         sp_eBird = sp_ebird,
         aou = sp_aou)

saveRDS(param_infer,paste0(output_dir,"/parameter_inference_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

adjs[1,"calibration"] <- as.numeric(cali$mean)
adjs[1,"calibration_sd"] <- as.numeric(cali$sd)

adjs[1,"calibration_mean"] <- as.numeric(cali_alt$mean)
adjs[1,"calibration_mean_sd"] <- as.numeric(cali_alt$sd)

adjs[1,"calibration_median"] <- as.numeric(cali_median$mean)
adjs[1,"calibration_median_sd"] <- as.numeric(cali_median$sd)

adjs[1,"calibration_log"] <- as.numeric(cali_log$mean)
adjs[1,"calibration_log_sd"] <- as.numeric(cali_log$sd)
# BETA_post <- fit$draws(variables = "BETA",
#                        format = "df")


cali_post <- fit$draws(variables = "calibration",
                       format = "df")

cali_alt_post <- fit$draws(variables = "calibration_alt",
                           format = "df")
cali_log_post <- fit$draws(variables = "calibration_log",
                           format = "df")

cali_median_post <- fit$draws(variables = "calibration_median",
                              format = "df")


# check for skew in betas -------------------------------------------------

betas <- summ %>% filter(grepl("beta[",variable,fixed = TRUE))

skew_flag <- e1071::skewness(betas$mean,type = 3)
kurtosis_flag <- e1071::kurtosis(betas$mean,type = 3)

adjs[1,"beta_skew"] <- as.numeric(skew_flag)
adjs[1,"beta_kurtosis"] <- as.numeric(kurtosis_flag)


# Posterior predictive check ----------------------------------------------


y_rep <- fit$draws("y_rep",format = "df") %>%
  slice_sample(n = 400)

y_rep <- posterior::as_draws_matrix(y_rep)
p975 <- quantile(stan_data$count,0.99)
ppc <- bayesplot::ppc_dens_overlay_grouped(y = stan_data$count,
                yrep = y_rep,
                group = combined$year)+
  coord_cartesian(xlim = c(0,p975))


# R-squared ---------------------------------------------------------------


y_pred_w_route <- exp(posterior::as_draws_matrix(fit$draws("E",format = "df")))

var_resid <- apply(y_pred_w_route-stan_data$count,MARGIN = 1,FUN = var)
var_fit <- apply(y_pred_w_route,MARGIN = 1,FUN = var)
R_squared_w_route <- (var_fit/(var_resid+var_fit))

gamma_post <- posterior::as_draws_matrix(fit$draws("gamma",format = "df"))
BETA_post <- posterior::as_draws_matrix(fit$draws("BETA",format = "df"))

y_pred_not_route <- matrix(data = NA,nrow = nrow(gamma_post),
                           ncol = stan_data$n_counts)
for(i in 1:stan_data$n_counts){
  y_pred_not_route[,i] <- exp(BETA_post + gamma_post[,stan_data$year[i]] + stan_data$log_mean_rel_abund[stan_data$route[i]])
}

var_resid <- apply(y_pred_not_route-stan_data$count,MARGIN = 1,FUN = var)
var_fit <- apply(y_pred_not_route,MARGIN = 1,FUN = var)
R_squared_not_route <- (var_fit/(var_resid+var_fit))


#hist(R_squared_not_route)

adjs$R_squared_w_route <- mean(R_squared_w_route)
adjs$R_squared_w_route_lci <- quantile(R_squared_w_route,0.05)
adjs$R_squared_w_route_uci <- quantile(R_squared_w_route,0.95)
adjs$R_squared_not_route <- mean(R_squared_not_route)
adjs$R_squared_not_route_lci <- quantile(R_squared_not_route,0.05)
adjs$R_squared_not_route_uci <- quantile(R_squared_not_route,0.95)


adjs_out <- bind_rows(adjs_out,adjs)

# Seasonal corrections on counts ------------------------------------------
strat_names <- combined %>%
  select(strata,strata_name) %>%
  distinct()

seasonal_strat <- summ %>% filter(grepl("doy_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = rep(c(1:stan_data$n_doy),times = stan_data$n_strata) + (min(combined$day_of_year)-1),
         strata = rep(c(1:stan_data$n_strata),each = stan_data$n_doy)) %>%
  inner_join(strat_names,by = "strata")

seasonal <- summ %>% filter(grepl("DOY_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = c(1:stan_data$n_doy)+ (min(combined$day_of_year)-1))

strat_labs <- seasonal_strat %>%
  filter(doy == max(combined$day_of_year))

vis_season <- ggplot(data = seasonal_strat,
                     aes(x = doy,y = mean,
                         group = strata_name,
                         colour = strata_name))+
  geom_ribbon(data = seasonal,aes(x = doy,y = mean,
                                  ymin = q5, ymax = q95),
              alpha = 0.2,
              inherit.aes = FALSE)+
  geom_line(data = seasonal,aes(x = doy,y = mean),
            inherit.aes = FALSE)+
  geom_line(alpha = 0.3)+
  coord_cartesian(xlim = c(100,240))+
  ggrepel::geom_text_repel(data = strat_labs, aes(label = strata_name),
                           size = 1.5,
                           min.segment.length = 0,
                           xlim = c(200,240),
                           nudge_x = 20)+
  scale_colour_viridis_d()+
  theme(legend.position = "none")



# explore relationship ----------------------------------------------------
y_pred_bind <- data.frame(year = rep(c(1:stan_data$n_years)+2012,2),
                          ebird_abund = rep(c(min(combined$ebird_abund),max(combined$ebird_abund)),each = stan_data$n_years))

predictions <- summ %>% filter(grepl("raw_prediction",variable)) %>%
  cbind(y_pred_bind)



vis_relationship <- ggplot(data = combined,
                           aes(y = count,x = ebird_abund))+
  geom_point(aes(colour = year), alpha = 0.3)+
#  geom_smooth(method = "gam")+
  geom_line(data = predictions,aes(x = ebird_abund,y = mean,
                                   colour = year,
                                   group = year))+
  geom_ribbon(data = predictions,aes(x = ebird_abund,y = mean,
                                     ymin = q5,ymax = q95,
                                   fill = year,
                                   group = year),
              alpha = 0.1)+
  scale_colour_viridis_c(aesthetics = c("colour","fill"))


vis_relationship2 <- ggplot(data = combined,
                           aes(y = count,x = ebird_abund))+
  geom_point(aes(colour = year), alpha = 0.3)+
  #  geom_smooth(method = "gam")+
  geom_line(data = predictions,aes(x = ebird_abund,y = mean,
                                   colour = year,
                                   group = year))+
  geom_ribbon(data = predictions,aes(x = ebird_abund,y = mean,
                                     ymin = q5,ymax = q95,
                                     fill = year,
                                     group = year),
              alpha = 0.1)+
  scale_colour_viridis_c(aesthetics = c("colour","fill"))+
  scale_y_continuous(transform = "log10")+
  scale_x_continuous(transform = "log10")



#
# combined_means <- combined %>%
#   mutate(prov_state = str_extract(pattern = ".+(?=[[:punct:]])",string = route_name)) %>%
#   group_by(prov_state,route_name,ebird_abund,logm) %>%
#   summarise(mean_count = mean(count),
#             min_count = quantile(count,0),
#             max_count = quantile(count,1)) %>%
#   filter(mean_count > 0) %>%
#   mutate(min_count = ifelse(min_count == 0,0.01,min_count))
#
# vis_relationship_mean <- ggplot(data = combined_means,
#                                 aes(y = (mean_count),x = ebird_abund,
#                                     colour = prov_state))+
#   geom_pointrange(aes(ymin = min_count,ymax = max_count),
#                   alpha = 0.3)+
#   scale_y_continuous(transform = "log10")+
#   scale_x_continuous(transform = "log10")+
#   geom_smooth(method = "lm")+
#   geom_abline(slope = 1,intercept = 0)
#
#
# vis_relationship_mean
#


# loading the breeding season abundance raster generated in
# 2_Prepare_eBird_relative_abundance_in_buffer.R
#
breed_abundance <- readRDS(paste0("data/species_relative_abundance/",
                                  sp_ebird,
                                  "_derived_breeding_relative_abundance.rds"))

names(breed_abundance) <- "breeding_abundance"

breed_abundance_full <- breed_abundance

# Strata level summaries --------------------------------
# bbs_strata_buf <- bbsBayes2::load_map("bbs_usgs") %>%
#   st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data
#
#
#
# # project boundary to match raster data
# region_boundary_proj <- st_transform(bbs_strata_buf, st_crs(breed_abundance))
# # crop and mask to boundary of BBS
# breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
#   mask(region_boundary_proj)


strata_proj <- st_transform(strata,
                            crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_strata <- terra::extract(breed_abundance_full,
                                      strata_proj,
                                       fun = sum,
                                       na.rm = TRUE,
                                       ID = TRUE,
                                       exact = FALSE)



## function to estimate full posterior of abundance

post_abund <- function(x,draws = cali_post$calibration,
                       fun = "mean",
                       p = 0.025,
                       use_log = FALSE){

  y <- vector("numeric",length(x))

  if(use_log){
    xl <- log(x)
  }
  if(fun == "mean"){
  for(i in 1:length(x)){
    if(use_log){
      y[i] <- exp(mean(xl[i]+as.numeric(draws)))
    }else{
     y[i] <- mean(x[i]*as.numeric(draws))
    }
  }
  }

  if(fun == "quantile"){
    for(i in 1:length(x)){
      if(use_log){
        y[i] <- exp(quantile(xl[i]+as.numeric(draws),p, names = FALSE))

      }else{
      y[i] <- quantile(x[i]*as.numeric(draws),p, names = FALSE)
      }
    }
  }

  return(y)
}
### 3^2 scaling is to account for the area of each grid-cell
###
###

# select calibration ------------------------------------------------------
if(log_calibration){
  cali_use <- cali_log_post
}else{
cali_use <- cali_post
}
names(cali_use)[1] <- "calibration"



strata_abund <- data.frame(region = strata_proj$strata_name,
                           region_type = "strata",
                           sum_abund = 3^2 * unlist(abundance_in_strata[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = log_calibration,
                               draws = cali_use$calibration,
                                 fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird,
         strata_name = region) %>%
  left_join(.,strata_names)
  #filter(pop_mean != 0)



comp_plot <- ggplot(data = strata_abund,
                    aes(x = strata_name,
                        y = pop_mean))+
  geom_errorbar(aes(ymin = pop_lci_80,ymax = pop_uci_80),
                alpha = 0.3, width = 0)+
  geom_point()+
  coord_flip()

#comp_plot


# compare to traditional estimates ----------------------------------------

strata_trad_j <- strata_trad %>%
  select(st_abrev,
         BCR,
         med.PopEst,
         LCI80.PopEst,
         UCI80.PopEst)


strata_compare <- strata_abund %>%
  full_join(.,strata_trad_j,
            by = c("prov_state" = "st_abrev",
                   "bcr" = "BCR"))



comp_trad_new_plot <- ggplot(data = strata_compare,
                             aes(x = med.PopEst,
                                 y = pop_median))+
  geom_abline(slope = 1, intercept = 0)+
  geom_abline(slope = 2, intercept = 0,linetype = 2)+
  geom_abline(slope = 5, intercept = 0,linetype = 3)+
  geom_errorbar(aes(ymin = pop_lci_80,ymax = pop_uci_80),
                alpha = 0.3, width = 0)+
  geom_errorbarh(aes(xmin = LCI80.PopEst,xmax = UCI80.PopEst),
                alpha = 0.3)+
  geom_point()+
  geom_text_repel(aes(label = strata_name),
                  size = 3)+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
  xlab("Traditional PIF population estimate")+
  ylab("PIF-Calibrated eBird relative abundance estimate")+
  labs(title = paste(sp_sel,"population estimates by BBS strata"),
       subtitle = "Diagonal lines = 1:1, 2:1, and 5:1")

png(filename = paste0("Figures/comp_trad_new_alt_",vers,sp_aou,"_",sp_ebird,".png"),
    res = 300,
    height = 6,
    width = 6,
    units = "in")
print(comp_trad_new_plot)
dev.off()

# abund_mapable <- breed_abundance %>%
#   terra::project(.,"EPSG:9822")
strat_sel <- strata_compare %>%
  filter((pop_median > 0 |
           med.PopEst > 0) &
  !is.na(strata_name)) %>%
  select(strata_name)

strata_w_data <- strata %>%
  filter(strata_name %in% strat_sel$strata_name)

bb <- sf::st_bbox(strata_w_data)

breed_abundance_plot <- breed_abundance
names(breed_abundance_plot) <- "breeding"
breed_abundance_plot <- breed_abundance_plot %>%
  mutate(.,breeding = ifelse(breeding == 0,NA,breeding),
         breeding = (breeding*cali$median)/900) # per hectare mean density

# Mapping -----------------------------------------------------------------


abund_map <- ggplot()+
  geom_sf(data = strata, fill = grey(0.95))+
  geom_spatraster(data = breed_abundance_plot,
                  maxcell = 16000000)+
  geom_sf(data = strata, fill = NA)+
  geom_sf_text(data = strata,aes(label = strata_name),
               size = 1)+
  geom_sf(data = routes_buf,aes(colour = mean_count),fill = NA)+
  coord_sf(xlim = c(bb[c("xmin","xmax")]),
           ylim = c(bb[c("ymin","ymax")]))+
  #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA,
                       end = 0.9,
                       name = "Birds/hectare")+
  scale_colour_viridis_c(direction = -1,
                       option = "rocket",
                       na.value = NA,
                       end = 0.9,
                       begin = 0.1,
                       name = "Mean count BBS")+
  theme_bw()+
  labs(title = paste(sp_sel,"eBird relative abundance values by BBS strata"))


#

png(filename = paste0("Figures/abund_map_alt_",vers,sp_aou,"_",sp_ebird,".png"),
    res = 400,
    height = 7,
    width = 6.5,
    units = "in")
print(abund_map)
dev.off()


saveRDS(comp_trad_new_plot,paste0("figures/saved_ggplots/trad_vs_new_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
saveRDS(abund_map,paste0("figures/saved_ggplots/abund_map_alt_",vers,sp_aou,"_",sp_ebird,".rds"))




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
  mutate(pop_mean = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

country_abund



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



bcr_abund <- data.frame(region = as.character(bcr_proj$BCR[abundance_in_bcr$ID]),
                            region_type = "bcr",
                            sum_abund = 3^2 * unlist(abundance_in_bcr[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

bcr_abund




# BCR map -----------------------------------------------------------------


bcr_sel <- bcr_abund$region

bcrs_eqarea <- sf::st_transform(bcrs,
                                st_crs(strata))
bcrs_w_data <- bcrs_eqarea%>%
  filter(BCR %in% bcr_sel,
         BCR != 999,
         BCR != 100)

bb <- sf::st_bbox(bcrs_w_data)



abund_map_bcrs <- ggplot()+
  geom_sf(data = bcrs_eqarea, fill = grey(0.95))+
  geom_spatraster(data = breed_abundance_plot,
                  maxcell = 16000000)+
  geom_sf(data = bcrs_eqarea, fill = NA)+
  geom_sf_text(data = bcrs_eqarea,aes(label = BCR),
               size = 1)+
  geom_sf(data = routes_buf,aes(colour = mean_count),fill = "darkred")+
  coord_sf(xlim = c(bb[c("xmin","xmax")]),
           ylim = c(bb[c("ymin","ymax")]))+
  #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA,
                       end = 0.9,
                       name = "Birds/hectare")+
  scale_colour_viridis_c(direction = -1,
                         option = "rocket",
                         na.value = NA,
                         end = 0.9,
                         begin = 0.1,
                         name = "Mean count BBS")+
  theme_bw()+
  labs(title = paste(sp_sel,"eBird relative abundance values by BCR"))


#

png(filename = paste0("Figures/abund_map_bcrs_alt_",vers,sp_aou,"_",sp_ebird,".png"),
    res = 400,
    height = 7,
    width = 6.5,
    units = "in")
print(abund_map_bcrs)
dev.off()



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
  mutate(pop_mean = post_abund(sum_abund, use_log = log_calibration, fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

USACAN_abund


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
  mutate(pop_mean = post_abund(sum_abund, use_log = log_calibration, fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

continents_abund




abundance_in_global <- sum(values(breed_abundance_full),na.rm = TRUE)



global_abund <- data.frame(region = "global",
                           region_type = "global",
                           sum_abund = 3^2 * abundance_in_global) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = log_calibration, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = log_calibration,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

global_abund



pop_ests_out <- bind_rows(USACAN_abund,
                          continents_abund,
                          global_abund,
                          country_abund,
                          bcr_abund,
                          strata_abund) %>%
  mutate(version = "PIF_eBird")

# Sum of population estimates

write_excel_csv(pop_ests_out,
                paste0("estimates/pop_ests_alt_",vers,sp_aou,sp_ebird,".csv"))


pop_ests_out_trad_sel <- pop_ests_out_trad %>%
  filter(species == sp_sel)


pop_compare_stack <- pop_ests_out %>%
  select(region,region_type,version,
         species,species_ebird,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95) %>%
  bind_rows(pop_ests_out_trad_sel)


pop_compare_stack_sel <- pop_compare_stack %>%
  filter(region_type %in% c("USACAN","continent","global","bcr"))


side_plot <- ggplot(data = pop_compare_stack_sel,
              aes(y = region,
                  x = pop_median,
                  colour = version))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbarh(aes(xmin = pop_lci_80,
                     xmax = pop_uci_80),
                 height = 0,
                 position = position_dodge(width = 0.5))+
  scale_colour_viridis_d(end = 0.8)+
  ylab("BCRs and broadscale estimates")+
  xlab("Population estimate with 80% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
  theme_bw()


saveRDS(side_plot,paste0("figures/saved_ggplots/side_plot_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

saveRDS(combined,paste0("data/main_data_df_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
pdf(paste0("figures/estimate_plots_alt_",vers,sp_aou,"_",sp_ebird,".pdf"),
    width = 11,
    height = 8.5)
print(vis_relationship)
print(vis_relationship2)
print(comp_trad_new_plot)
print(abund_map)
print(side_plot)
print(abund_map_bcrs)
print(ppc)
print(vis_season)
dev.off()

write_csv(adjs_out,paste0("adjs_out",vers,today,".csv"))


 } #end of species loop
#



