
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

if(Sys.info()[["nodename"]] == "WNCRLABN72960"){
setwd("c:/Users/SmithAC/Documents/GitHub/PIF_population_estimates")
}

source("functions/GAM_basis_function_mgcv.R")
source("functions/neighbours_define.R")
yr_ebird <- 2023 # prediction year for eBird relative abundance

selected_species <- readRDS("data/selected_species.rds")


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

re_napops <- FALSE # set to True to re-run the napops extraction

if(re_napops){
  napops_species <- napops::list_species() |>
    dplyr::select(Species,Common_Name)

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
  # Unless the road model suggests that EDR is smaller on roads (contrary to physics)
  cefs <- napops::coef_distance(species = sp) |>
    dplyr::filter(Model %in% c(1,2)) |>
    dplyr::arrange(AIC)
if(nrow(cefs) > 1){
  if((cefs[1,"Model"] == 2 & cefs[1,"Road"] < 0) |
     sp == "PRFA"){ #one-off to remove the PRFA roadside estimates that have CI so wide that the sd is > mean
    cefs <- cefs |>
      dplyr::filter(Model == 1)
  }
}
  if(nrow(cefs) == 0){next} # skip if no napops output

  w_mod <- cefs[1,"Model"]

  edrt <- try(napops::edr(species = sp,
                  model = w_mod,
                  forest = 0.5,
                  road = TRUE,
                  quantiles = c(0.1625,0.8375), # 1 sd below the mean
                  samples = 1000),
              silent = TRUE)

  if(edrt$EDR_est > 500){ # two species where the model 2 estimates of edr are silly (> 5 km!)
    edrt <- try(napops::edr(species = sp,
                    model = 1,
                    forest = 0.5,
                    road = TRUE,
                    quantiles = c(0.1625,0.8375), # 1 sd below the mean
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
    ExpAdjs[i,"edr_sd"] <- c(as.numeric(edrt$EDR_83.75)-as.numeric(edrt$EDR_16.25))/2
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

  sp <- unlist(ExpAdjs[i,"Species"])

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
                          quantiles = c(0.1625,0.8375), # 1 sd below the mean
                          samples = 1000),
              silent = TRUE)

  if(class(availability) == "try-error"){
    availability <- try(napops::avail(species = sp,
                                      time = 3,
                                      model = 1, #intercept model
                                      od = mean_doy, # mean of the doy for all BBS surveys
                                      tssr = mid_time, #mean of time of day for all BBS surveys
                                      quantiles = c(0.1625,0.8375), # 1 sd below the mean
                                      samples = 1000),
                        silent = TRUE)
    w_mod <- 1
  }else{
    if(availability$p_est < 0.1){
      availability <- try(napops::avail(species = sp,
                                        time = 3,
                                        model = 1,
                                        od = mean_doy, # mean of the doy for all BBS surveys
                                        tssr = mid_time, #mean of time of day for all BBS surveys
                                        quantiles = c(0.1625,0.8375), # 1 sd below the mean
                                        samples = 1000),
                          silent = TRUE)
      w_mod <- 1
    }
  }

    if(class(availability) == "try-error"){
          warning(paste("availability extraction did not work for",sp))
          next
        }



  ExpAdjs[i,"availability"] <- as.numeric(availability$p_est)

  if("p_16.25" %in% names(availability)){
    ExpAdjs[i,"availability_sd"] <- c(as.numeric(availability$p_83.75)-as.numeric(availability$p_16.25))/2
  }

  ExpAdjs[i,"availability_model"] <- w_mod

  if(sp %in% c("BOWA","CORE","CAHU")){
    ExpAdjs[i,"availability"] <- NA
    ExpAdjs[i,"availability_sd"] <- NA
  }
  rm(availability, w_mod)
}

ExpAdjs <- ExpAdjs |>
  dplyr::mutate(use_availability = ifelse(is.na(availability),FALSE,TRUE),
                availability_log_scaled = log(1/availability),
                availability_sd_log_scaled = 0.5*(log(1/(availability-availability_sd))-log(1/(availability+availability_sd))),
                availability_sd_log_scaled = ifelse(is.na(availability_sd_log_scaled),
                                                    log(1/(availability+availability_sd)),
                                                    availability_sd_log_scaled))


write_csv(ExpAdjs,"Species_correction_factors_w_edr_availability.csv")
}else{
  ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")
}



routes_buf_all <- readRDS("data/all_routes_buffered.rds")

strata <- bbsBayes2::load_map("bbs_usgs")

strata_names <- strata %>%
  sf::st_drop_geometry()





# load traditional estimates for later comparison -------------------------
strata_join <- strata_names %>%
  select(strata_name,country_code,
         prov_state, bcr) %>%
  st_drop_geometry() %>%
  distinct()

strata_trad_all <- read_csv("Stanton_2019_code/output_EDR/PSest_Prov_by_BCR__100iter.csv") %>%
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

USACAN_trad_all <- read_csv("Stanton_2019_code/output_EDR/PSest_Global_WH_NA__100iter.csv") %>%
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



global_trad_all <- read_csv("Stanton_2019_code/output_EDR/PSest_Global_WH_NA__100iter.csv") %>%
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



prov_trad_all <- read_csv("Stanton_2019_code/output_EDR/PSest_by_Prov__100iter.csv") %>%
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

bcr_trad_all <- read_csv("Stanton_2019_code/output_EDR/PSest_by_BCR__100iter.csv")%>%
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






#
sps_list <- readRDS("data/all_bbs_data_species_list.rds") #

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
country_codes <- readxl::read_xlsx("iso_codes.xlsx") %>%
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
estimate_rho <- FALSE # set to true to allow the model to estimate a non-unit log-log slope

use_traditional <- FALSE # if TRUE uses the traditional adjustments for TOD and distance







if(use_traditional){
  vers <- ifelse(estimate_rho,"trad_rho","trad")
}else{
  vers <- ifelse(estimate_rho,"_rho","")
}


re_fit <- TRUE # set to true to re-run plotting and summaries
re_run_model <- TRUE # set to true to rerun the stan model fit

#output_dir <- "G:/PIF_population_estimates/output"
output_dir <- "output"
trim_rel_abund <- TRUE
trim_bbs_routes <- TRUE

for(sp_sel in selected_species){ #

 sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)



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







if(!file.exists(paste0("data/species_relative_abundance_2023/",
                       sp_ebird,"_relative_abundance.rds"))){next}

rel_abund <- readRDS(paste0("data/species_relative_abundance_2023/",
                            sp_ebird,"_relative_abundance.rds"))


if(trim_rel_abund){ # if TRUE drops routes with < 10th %ile relative abundance
  tenth_percentile <- unname(quantile(rel_abund$ebird_abund,0.1))
  rel_abund <- rel_abund %>%
    filter(ebird_abund >= tenth_percentile)
}

if(trim_bbs_routes){ # if TRUE drops routes with < 10th %ile observed mean counts
  temp_means <- raw_counts %>%
    group_by(route_name) %>%
    summarise(mean_obs = mean(count)) %>%
    filter(route_name %in% rel_abund$route_name)


  tenth_percentile_bbs <- unname(quantile(temp_means$mean_obs,0.1))
  temp_means <- temp_means %>%
    filter(mean_obs > tenth_percentile_bbs)

  raw_counts <- raw_counts %>%
    filter(route_name %in% temp_means$route_name)
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
                  use_t = 1,
                  use_ppc = 1,
                  est_rho = ifelse(estimate_rho,1,0)
                  )


# if(nrow(adjs) == 0){
#   # stan_data$c_p = mean(ExpAdjs$Pair2,na.rm = T)
#   # stan_data$c_t = exp(mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + 0.5*(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2))
#   # stan_data$sd_c_t = exp(2*mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)*(exp(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)-1)
#   # stan_data$c_d_lower = mean(ExpAdjs$Dist.Lower,na.rm = T)
#   # stan_data$c_d_upper = mean(ExpAdjs$Dist.Upper,na.rm = T)
#
# }
adjs$n_routes <- stan_data$n_routes

saveRDS(stan_data,paste0("stan_data/stan_data_",vers,sp_aou,"_",sp_ebird,".rds"))

}

