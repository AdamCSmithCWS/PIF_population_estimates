
library(sf)
library(tidyverse)
library(ebirdst)
library(patchwork)
library(terra)
library(bbsBayes2)
library(cmdstanr)
library(ggrepel)
library(tidyterra)
#library(napops)
library(tidybayes)
library(doParallel)
library(foreach)

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
#
# sp_example <- c("American Robin",
#                 "Western Meadowlark","Clay-colored Sparrow",
#                 "Sharp-tailed Grouse",
#                 "Sora",
#                 "Killdeer","Long-billed Curlew",
#                 "Yellow-rumped Warbler","Olive-sided Flycatcher",
#                 "Purple Martin",
#                 "Barn Swallow",
#                 "Verdin","Ash-throated Flycatcher","Black-throated Sparrow",
#                 "Blue Jay","Varied Thrush","Veery","Wood Thrush","Chestnut-collared Longspur",
#                 "Bobolink","Savannah Sparrow","Grasshopper Sparrow",
#                 "Horned Lark","Baird's Sparrow",
#                 "Eastern Whip-poor-will", "Common Nighthawk")

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
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE) %>%
  mutate(cn = ifelse(cn == "McCown's Longspur",
                     "Thick-billed Longspur",
                     cn),
         cn = ifelse(cn == "Gray Jay",
                     "Canada Jay",
                     cn),
         cn = ifelse(grepl("Le Conte's",cn),
                     gsub("Le Conte's","LeConte's",cn),
                     cn))



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

  if(sp == "BOWA"){
    ExpAdjs[i,"availability"] <- NA
    ExpAdjs[i,"availability_sd"] <- NA
  }
  rm(availability, w_mod)
}

ExpAdjs <- ExpAdjs |>
  dplyr::mutate(use_availability = ifelse(is.na(availability),FALSE,TRUE))


write_csv(ExpAdjs,"Species_correction_factors_w_edr_availability.csv")
}else{
  ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")
}

# ID species missing edrs or availability-rates
# corrections_example_sp <- ExpAdjs %>%
#   filter(cn %in% sp_example)
#
# if(any(!corrections_example_sp$use_availability) |
#    any(!corrections_example_sp$use_edr)){
#   warning("edr or availability missing for some species")
# }

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
log_normal_calibration <- TRUE # conduct calibration assuming a lognormal
# distribution of variation among routes

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
use_traditional <- TRUE
if(use_traditional){
  vers <- "trad_"
}else{
  vers <- ""
}

# if(file.exists(paste0("adjs_out",vers,today,".csv"))){
#   adjs_out <- read_csv(paste0("adjs_out",vers,today,".csv"))
#   wh_drop_already <- which(sp_example %in% adjs_out$cn)
#   wh_drop <- c(wh_drop_already,8)
# }else{
#   wh_drop <- 8
# }


re_fit <- TRUE # set to true to re-run plotting and summaries
re_run_model <- TRUE # set to true to rerun the stan model fit

output_dir <- "G:/PIF_population_estimates/output"
#output_dir <- "output"
trim_rel_abund <- TRUE

# Parallel species loop ---------------------------------------------------
#n_cores <- 2

#cluster <- makeCluster(n_cores, type = "PSOCK")
# library(sf)
# library(tidyverse)
# library(ebirdst)
# library(patchwork)
# library(terra)
# library(bbsBayes2)
# library(cmdstanr)
# library(ggrepel)
# library(tidyterra)
# library(napops)
# library(tidybayes)
# library(doParallel)
# library(foreach)

#registerDoParallel(cluster)


# test <- foreach(sp_sel = rev(sps_list$english)[1:8],
#                 .packages = c("bbsBayes2",
#                               "tidyverse",
#                               "cmdstanr",
#                               "ebirdst",
#                               "patchwork",
#                               "sf",
#                               "tidybayes",
#                               "posterior",
#                               "napops",
#                               "ggrepel",
#                               "terra"),
#                 .errorhandling = "pass") %dopar%
#   {


for(sp_sel in c("Hermit Thrush","American Robin",
                "Barn Swallow",
                "Blackpoll Warbler")){ # rev(sps_list$english[1:348])){#sp_example[-wh_drop]){#list$english){
#sp_sel = "Bank Swallow"
#for(sp_sel in rev(sps_list$english)[29]){
 sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)


  if(file.exists(paste0(output_dir,"/calibration_fit_alt_",
                        vers,sp_aou,"_",sp_ebird,".rds")) &
     !re_fit){next}



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







if(!file.exists(paste0("data/species_relative_abundance/",
                       sp_ebird,"_relative_abundance.rds"))){next}

rel_abund <- readRDS(paste0("data/species_relative_abundance/",
                            sp_ebird,"_relative_abundance.rds"))


if(trim_rel_abund){
  tenth_percentile <- unname(quantile(rel_abund$ebird_abund,0.1))
  rel_abund <- rel_abund %>%
    filter(ebird_abund >= tenth_percentile)
}

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



# filtering routes on which species wasn't observed.

mean_counts <- raw_counts %>%
  group_by(route_name) %>%
  summarise(mean_count = mean(count)) %>%
  filter(route_name %in% combined$route_name)

routes_buf <- routes_buf_all %>%
  inner_join(mean_counts, by = "route_name")


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
  select(route,logm) %>%
  distinct() %>%
  arrange(.,route)


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



stan_data <- list(n_routes = max(combined$route),
                  n_obs = max(combined$observer),
                  n_counts = nrow(combined),
                  n_years = max(combined$yr),
                  n_knots_doy = 7,
                  n_doy = max(combined$doy),
                  mean_doy = as.integer(round(mean(combined$doy))),
                  n_strata = max(combined$strata),

                  count = combined$count,
                  route = combined$route,
                  observer = combined$observer,
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


if(re_run_model){



  params_to_summarise <- c("nu",
                           "BETA",
                           "beta_raw",
                           "beta",
                           "sd_beta",
                           "sd_obs",
                           "obs",
                           "nu_obs",
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
                           "calibration_alt2",
                           "calibration_median",
                           "calibration_alt1",
                           "calibration_mean",
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
                           "doy_pred",
                           "pred_count",
                           "pred_count_alt2",
                           "pred_count_median",
                           "pred_count_alt1",
                           "pred_count_mean",
                           "pred_count_r")


    model <- cmdstanr::cmdstan_model("models/ebird_rel_abund_calibration_spatial_year_doy_sep_obs_rte.stan")


fit <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500,
                    iter_warmup = 2000,
                    iter_sampling = 2000,
                    adapt_delta = 0.8,
                    max_treedepth = 11,
                    show_exceptions = FALSE)


summ <- fit$summary(variables = params_to_summarise)
#shinystan::launch_shinystan(fit)


fit$save_object(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

#fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
saveRDS(summ,paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))


# test for zero inflation -------------------------------------------------


# Posterior predictive check ----------------------------------------------


y_rep <- fit$draws("y_rep",format = "df") %>%
  slice_sample(n = 500)

y_rep <- posterior::as_draws_matrix(y_rep)

p_z_rep <- length(which(y_rep == 0))/length(y_rep)
p_z_obs <- length(which(stan_data$count == 0))/stan_data$n_counts

ZI <- ifelse(p_z_rep/p_z_obs > 0.85, FALSE, TRUE)
p99 <- quantile(stan_data$count,0.95)
# ppc <- bayesplot::ppc_dens_overlay_grouped(y = stan_data$count,
#                 yrep = y_rep,
#                 group = combined$year,
#                 trim = TRUE)+
#   coord_cartesian(xlim = c(0,p99))

# ppc_alt <- bayesplot::ppc_violin_grouped(y = stan_data$count,
#                     yrep = y_rep,
#                     group = combined$year,
#                     y_draw = "both",
#                     alpha = 0.5) +
#   coord_cartesian(ylim = c(0,p99))

ppc <- bayesplot::ppc_bars_grouped(y = stan_data$count,
                                   yrep = y_rep,
                                   group = stan_data$year,
                                   size = 0.1)+
  coord_cartesian(xlim = c(0,p99))

#
# if(ZI){
#   model <- cmdstanr::cmdstan_model("models/ebird_rel_abund_calibration_spatial_year_doy_sep_obs_rte_zip.stan")
#
#   params_to_summarise <- c(params_to_summarise,"theta")
#
#   fit <- model$sample(data = stan_data,
#                       parallel_chains = 4,
#                       refresh = 500,
#                       iter_warmup = 2000,
#                       iter_sampling = 2000,
#                       adapt_delta = 0.8,
#                       max_treedepth = 11,
#                       show_exceptions = FALSE)
#
#
#   summ <- fit$summary(variables = params_to_summarise)
#   #shinystan::launch_shinystan(fit)
#
#
#   fit$save_object(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
#
#   #fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
#   saveRDS(summ,paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
#
#
#   # test for zero inflation -------------------------------------------------
#
#
#   # Posterior predictive check ----------------------------------------------
#
#
#   y_rep <- fit$draws("y_rep",format = "df") %>%
#     slice_sample(n = 500)
#
#   y_rep <- posterior::as_draws_matrix(y_rep)
#
#   p_z_rep <- length(which(y_rep == 0))/length(y_rep)
#   p_z_obs <- length(which(stan_data$count == 0))/stan_data$n_counts
#
#   ZI <- ifelse(p_z_rep/p_z_obs > 0.85, FALSE, TRUE)
#   p99 <- quantile(stan_data$count,0.95)
#   # ppc <- bayesplot::ppc_dens_overlay_grouped(y = stan_data$count,
#   #                 yrep = y_rep,
#   #                 group = combined$year,
#   #                 trim = TRUE)+
#   #   coord_cartesian(xlim = c(0,p99))
#
#   # ppc_alt <- bayesplot::ppc_violin_grouped(y = stan_data$count,
#   #                     yrep = y_rep,
#   #                     group = combined$year,
#   #                     y_draw = "both",
#   #                     alpha = 0.5) +
#   #   coord_cartesian(ylim = c(0,p99))
#
#   ppc <- bayesplot::ppc_bars_grouped(y = stan_data$count,
#                                      yrep = y_rep,
#                                      group = stan_data$year,
#                                      size = 0.1)+
#     coord_cartesian(xlim = c(0,p99))
#
# }

}else{
  fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
  summ <- readRDS(paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
}


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

cali_lognormal <- fit$summary(variable = "calibration",
                              "mean",
                              "median",
                              "sd",
                              "rhat",
                              "ess_bulk",
                              q2_5 = ~quant(.x,q = 0.025),
                              q10 = ~quant(.x,q = 0.10),
                              q90 = ~quant(.x,q = 0.90),
                              q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration assuming log normal variation among routes",
         variable = "calibration_lognormal")

cali_lognormal_alt1 <- fit$summary(variable = "calibration_alt1",
                              "mean",
                              "median",
                              "sd",
                              "rhat",
                              "ess_bulk",
                              q2_5 = ~quant(.x,q = 0.025),
                              q10 = ~quant(.x,q = 0.10),
                              q90 = ~quant(.x,q = 0.90),
                              q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration assuming log normal variation among routes removing t-dist factor",
         variable = "calibration_alt1")


cali_lognormal_alt2 <- fit$summary(variable = "calibration_alt2",
                                   "mean",
                                   "median",
                                   "sd",
                                   "rhat",
                                   "ess_bulk",
                                   q2_5 = ~quant(.x,q = 0.025),
                                   q10 = ~quant(.x,q = 0.10),
                                   q90 = ~quant(.x,q = 0.90),
                                   q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration assuming log normal variation among routes removing variance",
         variable = "calibration_alt2")



cali_alt <- fit$summary(variable = "calibration_mean",
                    "mean",
                    "median",
                    "sd",
                    "rhat",
                    "ess_bulk",
                    q2_5 = ~quant(.x,q = 0.025),
                    q10 = ~quant(.x,q = 0.10),
                    q90 = ~quant(.x,q = 0.90),
                    q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration realised mean across route-observer")


cali <- fit$summary(variable = "calibration_median",
                    "mean",
                    "median",
                    "sd",
                    "rhat",
                    "ess_bulk",
                    q2_5 = ~quant(.x,q = 0.025),
                    q10 = ~quant(.x,q = 0.10),
                    q90 = ~quant(.x,q = 0.90),
                    q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "calibration realised median across route-observer")






pred_count_lognormal <- fit$summary(variable = "pred_count",
                              "mean",
                              "median",
                              "sd",
                              "rhat",
                              "ess_bulk",
                              q2_5 = ~quant(.x,q = 0.025),
                              q10 = ~quant(.x,q = 0.10),
                              q90 = ~quant(.x,q = 0.90),
                              q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "pred_count assuming log normal variation among routes",
         variable = "pred_count_lognormal")

pred_count_lognormal_alt1 <- fit$summary(variable = "pred_count_alt1",
                                   "mean",
                                   "median",
                                   "sd",
                                   "rhat",
                                   "ess_bulk",
                                   q2_5 = ~quant(.x,q = 0.025),
                                   q10 = ~quant(.x,q = 0.10),
                                   q90 = ~quant(.x,q = 0.90),
                                   q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "pred_count assuming log normal variation among routes removing t-dist factor",
         variable = "pred_count_alt1")


pred_count_lognormal_alt2 <- fit$summary(variable = "pred_count_alt2",
                                   "mean",
                                   "median",
                                   "sd",
                                   "rhat",
                                   "ess_bulk",
                                   q2_5 = ~quant(.x,q = 0.025),
                                   q10 = ~quant(.x,q = 0.10),
                                   q90 = ~quant(.x,q = 0.90),
                                   q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "pred_count assuming log normal variation among routes removing variance",
         variable = "pred_count_alt2")



pred_count_alt <- fit$summary(variable = "pred_count_mean",
                        "mean",
                        "median",
                        "sd",
                        "rhat",
                        "ess_bulk",
                        q2_5 = ~quant(.x,q = 0.025),
                        q10 = ~quant(.x,q = 0.10),
                        q90 = ~quant(.x,q = 0.90),
                        q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "pred_count realised mean across route-observer")


pred_count <- fit$summary(variable = "pred_count_median",
                    "mean",
                    "median",
                    "sd",
                    "rhat",
                    "ess_bulk",
                    q2_5 = ~quant(.x,q = 0.025),
                    q10 = ~quant(.x,q = 0.10),
                    q90 = ~quant(.x,q = 0.90),
                    q97_5 = ~quant(.x,q = 0.975))%>%
  mutate(inference = "pred_count realised median across route-observer")





param_infer <- bind_rows(cali_alt,
                         cali,
                         cali_lognormal,
                         cali_lognormal_alt1,
                         cali_lognormal_alt2,
                         pred_count_alt,
                         pred_count,
                         pred_count_lognormal,
                         pred_count_lognormal_alt1,
                         pred_count_lognormal_alt2,
                         avail_correction_realised,
                         p_avail_realised,
                         edr_realised,
                         surveyed_area_realised) %>%
  mutate(english = sp_sel,
         sp_eBird = sp_ebird,
         aou = sp_aou)

#param_infer2 <- readRDS(paste0(output_dir,"/parameter_inference_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

saveRDS(param_infer,paste0(output_dir,"/parameter_inference_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

adjs[1,"calibration"] <- as.numeric(cali$mean)
adjs[1,"calibration_sd"] <- as.numeric(cali$sd)

adjs[1,"calibration_mean"] <- as.numeric(cali_alt$mean)
adjs[1,"calibration_mean_sd"] <- as.numeric(cali_alt$sd)

adjs[1,"calibration_lognormal"] <- as.numeric(cali_lognormal$mean)
adjs[1,"calibration_lognormal_sd"] <- as.numeric(cali_lognormal$sd)

adjs[1,"calibration_lognormal_alt1"] <- as.numeric(cali_lognormal_alt1$mean)
adjs[1,"calibration_lognormal_alt1_sd"] <- as.numeric(cali_lognormal_alt1$sd)

adjs[1,"calibration_lognormal_alt2"] <- as.numeric(cali_lognormal_alt2$mean)
adjs[1,"calibration_lognormal_alt2_sd"] <- as.numeric(cali_lognormal_alt2$sd)

# BETA_post <- fit$draws(variables = "BETA",
#                        format = "df")


cali_lognormal_post <- fit$draws(variables = "calibration",
                       format = "df")

cali_lognormal_alt1_post <- fit$draws(variables = "calibration_alt1",
                                 format = "df")

cali_lognormal_alt2_post <- fit$draws(variables = "calibration_alt2",
                                 format = "df")

cali_alt_post <- fit$draws(variables = "calibration_mean",
                           format = "df")

cali_post <- fit$draws(variables = "calibration_median",
                              format = "df")

trm_mean <- function(x,p = 0.1){
  q1 <- quantile(x,p)
  q2 <- quantile(x,1-p)
  xtrim <- x[which(x > q1 & x < q2)]
  trim_m <- mean(xtrim)
  return(trim_m)
}

cali_route_post <- fit$draws(variables = "calibration_r",
                           format = "df")

cali_trim_post <- cali_alt_post

for(i in 1:nrow(cali_route_post)){
  cali_trim_post[i,1] <- trm_mean(as.numeric(cali_route_post[i,1:stan_data$n_routes]))
}

adjs[1,"calibration_trimmed"] <- mean(as.numeric(cali_trim_post$calibration_mean))
adjs[1,"calibration_trimmed_sd"] <- sd(as.numeric(cali_trim_post$calibration_mean))

# measure kurtosis and skew in betas -------------------------------------------------

betas <- summ %>% filter(grepl("beta[",variable,fixed = TRUE))

skew_flag <- e1071::skewness(betas$mean,type = 3)
kurtosis_flag <- e1071::kurtosis(betas$mean,type = 3)

sd_beta <- summ %>% filter(grepl("sd_beta",variable,fixed = TRUE))
adj_post <- summ %>% filter(grepl("adj",variable,fixed = TRUE))


adjs[1,"beta_skew"] <- as.numeric(skew_flag)
adjs[1,"beta_kurtosis"] <- as.numeric(kurtosis_flag)


route_link <- combined %>%
  group_by(route,route_name,strata_name) %>%
  summarise(mean_count = mean(count),
            mean_ebird_abund = mean(ebird_abund),
            n_counts = n(),
            .groups = "drop") %>%
  arrange(route)


routes_buf_extra <- routes_buf %>%
  sf::st_buffer(10000)

calibration_by_rts <- summ %>% filter(grepl("calibration_r[",variable,fixed = TRUE)) %>%
  select(mean, median,sd,rhat,ess_bulk) %>%
  rename_with(.fn = ~ paste0("calib_rt_",.x))


pred_count_by_rts <- summ %>% filter(grepl("pred_count_r[",variable,fixed = TRUE)) %>%
  select(mean, median,q5,q95,sd,rhat,ess_bulk) %>%
  rename_with(.fn = ~ paste0("pred_count_rt_",.x))

betas_by_route <- betas %>%
  bind_cols(route_link) %>%
  bind_cols(calibration_by_rts) %>%
  bind_cols(pred_count_by_rts) %>%
  mutate(mean_beta = (mean),
            median_calibration = (calib_rt_median),
            mean_pred_count = (pred_count_rt_median),
            min_pred_count = (pred_count_rt_q5),
            max_pred_count = (pred_count_rt_q95))


betas_plot <- routes_buf_extra %>%
  left_join(betas_by_route, by = "route_name")

beta_distr <- ggplot()+
  geom_sf(data = strata)+
  geom_sf(data = betas_plot,
          aes(colour = mean_beta,
              fill = mean_beta))+
  scale_colour_viridis_b(aesthetics = c("colour","fill"),
                         breaks = seq(-12,12,2),
                         guide = guide_coloursteps(even.steps = FALSE,
                                                   show.limits = TRUE),
                         option = "H")


flags <- data.frame(name = c("Mean across routes",
                             "Median across routes"),
                    value = c(mean(betas_by_route$mean_beta),
                              median(betas_by_route$mean_beta)))


beta_hist <- ggplot()+
  geom_histogram(data = betas_by_route,
                 aes(x = mean_beta),
                 bins = 100)+
  geom_vline(data = flags,
             aes(colour = name,
                 xintercept = value))+
  geom_point(data = flags,
             aes(colour = name,y = 0,x = value))+
  theme_bw()+
  scale_colour_viridis_d()


flags <- data.frame(name = c("Mean across routes",
                             "Trimmed mean across routes",
                              "Median across routes",
                              "Calibration assuming lognormal and t-dist",
                              "Calibration assuming lognormal",
                              "Calibration using only hyperparameter (median of lognormal)",
                              "Calibration using posterior median across routes"),
                    value = c(mean(betas_by_route$median_calibration),
                              mean(as.numeric(cali_trim_post$calibration_mean)),
                              median(betas_by_route$median_calibration),
                              cali_lognormal$mean,
                              cali_lognormal_alt1$mean,
                              cali_lognormal_alt2$mean,
                              cali$mean))

calib_hist <- ggplot()+
  geom_histogram(data = betas_by_route,
                 aes(x = median_calibration),
                 bins = 100)+
  geom_vline(data = flags,
             aes(colour = name,
                 xintercept = value))+
  geom_point(data = flags,
              aes(colour = name,y = 0,x = value))+
  theme_bw()+
  scale_colour_viridis_d()+
  scale_x_continuous(transform = "log", labels = scales::label_comma())




flags <- data.frame(name = c("Mean observed count across routes",
                             "Mean predicted count across routes",
                             "Expected count assuming lognormal and t-dist",
                             "Expected count assuming lognormal",
                             "Expected count using only hyperparameter (median of lognormal)",
                             "Expected count using posterior median across routes"),
                    value = c(mean(combined$count),
                              mean(betas_by_route$mean_pred_count),
                              pred_count_lognormal$mean,
                              pred_count_lognormal_alt1$mean,
                              pred_count_lognormal_alt2$mean,
                              pred_count$mean))
count_hist <- ggplot()+
  geom_histogram(data = betas_by_route,
                 aes(x = mean_count),
                 bins = 100)+
  geom_vline(data = flags,
             aes(colour = name,
                 xintercept = value))+
  geom_point(data = flags,
             aes(colour = name,y = 0,x = value))+
  theme_bw()+
  scale_colour_viridis_d()+
  scale_x_continuous(labels = scales::label_comma())

saveRDS(betas_by_route,paste0(output_dir,"/betas_by_route_alt_",
               vers,sp_aou,"_",sp_ebird,".rds"))




# flags <- data.frame(name = c("Mean observed count across routes",
#                              "Mean predicted count across routes",
#                              "Expected count assuming lognormal and t-dist",
#                              "Expected count assuming lognormal",
#                              "Expected count using only hyperparameter (median of lognormal)",
#                              "Expected count using posterior median across routes"),
#                     value = c(mean(combined$count),
#                               mean(betas_by_route$mean_pred_count),
#                               pred_count_lognormal$mean,
#                               pred_count_lognormal_alt1$mean,
#                               pred_count_lognormal_alt2$mean,
#                               pred_count$mean))
y_pred_bind <- data.frame(year = rep(c(1:stan_data$n_years)+2012,2),
                          ebird_abund = rep(c(min(combined$ebird_abund),max(combined$ebird_abund)),each = stan_data$n_years))

predictions <- summ %>% filter(grepl("raw_prediction",variable)) %>%
  cbind(y_pred_bind)


predictions_sel <- predictions %>%
  filter(year == 2022)

obs_pred_count <- ggplot(data = betas_by_route,
                         aes(y = mean_count,x = mean_ebird_abund ))+
  geom_point(alpha = 0.3, aes(colour = pred_count_rt_median))+
  theme_bw()+
  scale_colour_viridis_c()+
  scale_x_continuous(labels = scales::label_comma(),
                     transform = "log10")+
  scale_y_continuous(transform = "log10",labels = scales::label_comma())+
  geom_line(data = predictions,aes(x = ebird_abund,y = mean, group = year),
            alpha = 0.5)+
  geom_line(data = predictions_sel,aes(x = ebird_abund,y = mean))+
  geom_smooth(method = "lm")

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





# Annual variation ------------------------------------------


annual_strat <- summ %>% filter(grepl("yeareffect[",variable,fixed = TRUE)) %>%
  mutate(year = rep(c(1:stan_data$n_years),each = stan_data$n_strata) + (min(combined$year)-1),
         strata = rep(c(1:stan_data$n_strata),times = stan_data$n_years)) %>%
  inner_join(strat_names,by = "strata")


strat_labs <- annual_strat %>%
  filter(year == min(combined$year))

vis_annual <- ggplot(data = annual_strat,
                     aes(x = year,y = mean,
                         group = strata_name,
                         colour = strata_name))+
  # geom_ribbon(aes(ymin = q5, ymax = q95),
  #             alpha = 0.2)+
  geom_line(alpha = 0.3)+
  coord_cartesian(xlim = c(2008,2023))+
  ggrepel::geom_text_repel(data = strat_labs, aes(label = strata_name),
                           size = 1.5,
                           min.segment.length = 0,
                           #xlim = c(200,240),
                           nudge_x = -2)+
  scale_colour_viridis_d()+
  scale_y_continuous(transform = "exp")+
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
log_normal_calibration <- FALSE
use_trimmed_calibration <- FALSE
# select calibration ------------------------------------------------------
if(log_normal_calibration){
  cali_use <- cali_lognormal_alt1_post
}else{
  cali_use <- cali_alt_post
}
if(use_trimmed_calibration){
  cali_use <- cali_trim_post
}
names(cali_use)[1] <- "calibration"



strata_abund <- data.frame(region = strata_proj$strata_name,
                           region_type = "strata",
                           sum_abund = 3^2 * unlist(abundance_in_strata[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                                 fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
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
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
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
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
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

# png(filename = paste0("Figures/abund_map_bcrs_alt_",vers,sp_aou,"_",sp_ebird,".png"),
#     res = 400,
#     height = 7,
#     width = 6.5,
#     units = "in")
# print(abund_map_bcrs)
# dev.off()



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
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE, fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
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
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE, fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

continents_abund




abundance_in_global <- sum(values(breed_abundance_full),na.rm = TRUE)



global_abund <- data.frame(region = "global",
                           region_type = "global",
                           sum_abund = 3^2 * abundance_in_global) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
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
#print(abund_map_bcrs)
print(ppc)
print(vis_season)
print(vis_annual)
print(beta_distr)
print(beta_hist)
print(count_hist / calib_hist)
print(obs_pred_count)
dev.off()

#write_csv(adjs_out,paste0("adjs_out",vers,today,".csv"))
saveRDS(adjs,paste0("estimates/parameters_",vers,sp_aou,"_",sp_ebird,".rds"))


 } #end of species loop
#

#parallel::stopCluster(cluster)


