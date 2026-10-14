
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
sp_new <- readRDS("data/sp_new.rds")

selected_species <- c("American Robin",
                 "Bank Swallow",
                 "California Gull",
                 "Common Nighthawk",
                 "Mallard",
                 "Sora",
                 "Bay-breasted Warbler",
                 "Wood Duck",
                 "Killdeer",
                 "Turkey Vulture",
                 "Green Heron",
                 "Northern Shrike",
                 "Short-eared Owl",
                 "Downy Woodpecker",
                 "Hairy Woodpecker",
                 "Common Nighthawk",
                 "Eastern Whip-poor-will",
                 "Barn Swallow",
                 "Canada Warbler",
                 "Virginia Rail",
                 "Pied-billed Grebe",
                 "Black Tern",
                 "Scarlet Tanager",
                 "Veery",
                 "Eastern Phoebe",
                 "Black-capped Chickadee",
                 "Western Meadowlark",
                 "Purple Martin",
                 "Canyon Wren",
                 "Blackpoll Warbler",
                 "Clark's Nutcracker",
                 "American Kestrel",
                 "Red-winged Blackbird",
                 "Wood Thrush",
                 "Black-chinned Hummingbird",
                 "Varied Thrush",
                 "Baird's Sparrow",
                 "Chestnut-collared Longspur",
                 "Savannah Sparrow",
                 "Barred Owl",
                 "Great Horned Owl")

# mean counts to id routes on which species has been observed -------------
mean_counts_all <- readRDS("data/mean_counts_by_route.rds")

raw_counts_all <- readRDS("data/all_counts_by_route.rds")

# general summaries of time and season
surveys <- load_bbs_data(release = 2026)$routes %>%
  filter(year >= 2015)  %>%  # 10 years missing 2020.
  mutate(doy = lubridate::yday(paste(year,month,day,sep = "-")),
         start_hr = floor(start_time/100),
         end_hr = floor(end_time/100),
         start_min = start_hr*60 + ((start_time/100)-start_hr)*100,
         end_min = end_hr*60 + ((end_time/100)-end_hr)*100,
         mid_time = ((end_min - start_min)/120)-0.5) %>%
  summarise(mean_time = mean(mid_time),
            mean_doy = mean(doy))

survey_dates <- load_bbs_data(release = 2026)$routes %>%
  filter(year >= 2015)  %>%   # 10 years from 2013 - 2023, missing 2020.
  mutate(doy = lubridate::yday(paste(year,month,day,sep = "-")),
         start_hr = floor(start_time/100),
         end_hr = floor(end_time/100),
         start_min = start_hr*60 + ((start_time/100)-start_hr)*100,
         end_min = end_hr*60 + ((end_time/100)-end_hr)*100,
         mid_time = ((end_min - start_min)/120)-0.5)

mid_time = round(surveys$mean_time,1)
mean_doy = round(surveys$mean_doy)




# Add na-pops EDRs --------------------------------------------------------



  #ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")

  decisions3 <- readRDS("full_decisions_file.rds")


adjustments <- decisions3 |>
  mutate(PIF_CALC_Rosenberg_2019 = ifelse(cn %in% sp_new$cn,
                                          "Yes","No"),
         PIF_time_of_day = 1/exp(TimeAdj.meanlog)) |>
  rename(common_name = cn,
         PIF_MDD = Dist2,
         PIF_pair = Pair2,
         NA_POPS_edr = edr,
         NA_POPS_EDR_model = EDR_model,
         NA_POPS_availability = availability,
         NA_POPS_availability_model = availability_model,
         BAM_population_estimate_Canada = canadian_2020_population_estimate_bam_M_birds,
         p_population_canada_ebirdst = p_population_canada,
         ACAD_pop_size_us_ca = pop_size_us_ca_,
         ACAD_pop_size_us_ca_source = pop_size_us_ca_s,
         ACAD_pop_size_us_ca_comment = pop_size_us_ca_com,
         ebird_species_code = ebird_species,
         ebirdst_resident = is_resident_ebirdst,
         ebirdst_breeding_quality = breeding_quality_ebirdst,
         ebirdst_resident_quality = resident_quality_ebirdst,
         prop_BOTW_range_na_ebirdst = p_BOTW_range_na_eBirdst,
         taxonomic_sort = sort
  )  |>
  arrange(taxonomic_sort)






routes_buf_all <- readRDS("data/all_routes_buffered_400m.rds")

strata <- bbsBayes2::load_map("bbs")%>%
  filter(country_code %in% c("CA","US"))

strata_names <- strata %>%
  sf::st_drop_geometry()





# load traditional estimates for later comparison -------------------------
strata_join <- strata_names %>%
  select(strata_name,country_code,
         prov_state, bcr) %>%
  st_drop_geometry() %>%
  distinct()




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
bcrs <- bbsBayes2::load_map("bcr")


today <- as_date(Sys.Date())





estimate_rho <- FALSE # set to true to allow the model to estimate a non-unit log-log slope



model_components <- data.frame(model = c("napops",
                                          "napops_edr",
                                          "napops_avail",
                                          "no_napops"),
                               use_avail = c(TRUE,FALSE,TRUE,FALSE),
                               use_edr = c(TRUE,TRUE,FALSE,FALSE),
                               use_tod = c(FALSE,TRUE,FALSE,TRUE),
                               use_mdd = c(FALSE,FALSE,TRUE,TRUE))


#output_dir <- "G:/PIF_population_estimates/output"
output_dir <- "output"
trim_rel_abund <- TRUE
trim_bbs_routes <- FALSE

species_to_prepare <- adjustments$common_name


# species loop ------------------------------------------------------------


for(i in 1:nrow(adjustments)){ #
  sp_sel <- adjustments[i,"sp_aou"]
 sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)
  adjustments[i,"sp_aou"] <- sp_aou
  adjustments[i,"sp_ebird"] <- sp_ebird


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


if(!file.exists(paste0("data/species_relative_abundance_",yr_ebird,"/",
                       sp_ebird,"_relative_abundance.rds"))){next}

rel_abund <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",
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
  if(tenth_percentile_bbs > 0){
  temp_means <- temp_means %>%
    filter(mean_obs >= tenth_percentile_bbs)
  }

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
         doy = 1+(day_of_year-(min(day_of_year))),# doy centered on earliest survey day for species
         route_obs = as.integer(factor(paste(route_name,obs_n,sep = "-"))))



# filtering routes on which species wasn't observed.

mean_counts <- raw_counts %>%
  group_by(route_name) %>%
  summarise(mean_count = mean(count)) %>%
  filter(route_name %in% combined$route_name)

routes_buf <- routes_buf_all %>%
  inner_join(mean_counts, by = "route_name")

saveRDS(routes_buf,paste0("stan_data/routes_buffered_",sp_aou,"_",sp_ebird,".rds"))



# spatial setup -----------------------------------------------------------

strata_incl <- combined %>%
  select(strata_name,strata) %>%
  distinct()

strata_map <- bbsBayes2::load_map("bbs") %>%
  inner_join(strata_incl) %>%
  arrange(strata)
if(nrow(strata_map) < 2){
  next
}
neighbours <- neighbours_define(strata_map,
                                strat_indicator = "strata",
                                save_plot_data = TRUE,
                                species = sp_sel)


all_years <- data.frame(year = c(min(combined$year):max(combined$year)))

y_2020 <- combined %>%
  select(year,yr) %>%
  distinct() %>%
  full_join(all_years, by = "year") %>%
  arrange(year) %>%
  mutate(y_2020 = ifelse(is.na(yr),0,1))


# selecting which time-intervals are lumped with 2020 data ----------------
# If mid-year is after 2020, then the intervals should be:
##  2021-2020 (meaning the y_2020 should be 0 when y == 2020)
## 2020-2019  (meaning the y_2020 should be 0 when y == 2019)
prev_2020_base <- ifelse(yr_ebird-(max(combined$year)) == 0,TRUE,FALSE)
if(prev_2020_base){
y_2020 <- y_2020 %>%
  mutate(y_2020 = ifelse(year == 2019,0,y_2020))
}else{
  y_2020 <- y_2020 %>%
    mutate(y_2020 = ifelse(year == 2021,0,y_2020))
}
# If mid-year is before 2020, then the intervals should be:
##  2019-2020 (meaning the y_2020 should be 0 when y == 2020)
## 2020-2021  (meaning the y_2020 should be 0 when y == 2021)



mean_abund <- combined %>%
  select(route,logm) %>%
  distinct() %>%
  arrange(.,route)



#adjustment scores
adjs <- adjustments %>%
  filter(common_name == sp_sel)
if(nrow(adjs) == 0){next}
#if(!use_traditional & !(adjs$use_edr | adjs$use_availability)){next}

# setup for GAM of doy ----------------------------------------------------


doy_gam <- gam_basis(combined$doy, nknots = 7,
                     predpoints = c(1:max(combined$doy)),
                     sm_name = "doy")


for(i in 1:nrow(model_components)){

  vers <- model_components[i,"model"]
  use_avail <- model_components[i,"use_avail"]
  use_edr <- model_components[i,"use_edr"]
  use_tod <- model_components[i,"use_tod"]
  use_mdd <- model_components[i,"use_mdd"]


edr <- ifelse(use_edr,
              adjs$NA_POPS_edr,
              adjs$PIF_MDD)

edr_sd <- ifelse(use_edr,
                 adjs$edr_sd,
                 1)

MDD_lower <- ifelse(use_mdd,
                    adjs$Dist.Lower,
                    1)
MDD_upper <- ifelse(use_mdd,
                    adjs$Dist.Upper,
                    1)



avail <- ifelse(use_avail,
                adjs$NA_POPS_availability,
                adjs$TimeAdj.meanlog)

avail_sd <- ifelse(use_avail,
                   adjs$availability_sd,
                   adjs$TimeAdj.sdlog)

pair <- adjs$PIF_pair
# prepare stan data ----------------------------------------------------------

if( is.na(avail) |
    is.na(avail_sd) |
    is.na(edr) |
    (use_mdd & (is.na(MDD_lower) | is.na(MDD_upper))) |
    (use_edr & is.na(edr_sd))
    ){
  print(paste("data_missing for",vers,sp_sel))
  next
}

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
                  scaling_factor = neighbours$scaling_factor,

                  # PIF adjustments
                  c_p = pair,
                  sd_c_p = 0.13,
                  c_d = edr,
                  sd_c_d = edr_sd,
                  c_t = avail,
                  sd_c_t = avail_sd,
                  c_d_lower = MDD_lower,
                  c_d_upper = MDD_upper,

                  # Conditional statments
                  use_edr = ifelse(use_edr,1,0),
                  use_availability = ifelse(use_avail,1,0),
                  use_pois = 0,
                  use_pair = 1,
                  use_t = 1,
                  use_ppc = 0,
                  est_rho = ifelse(estimate_rho,1,0),
                  use_likelihood = 1
                  )


# if(nrow(adjs) == 0){
#   # stan_data$c_p = mean(ExpAdjs$Pair2,na.rm = T)
#   # stan_data$c_t = exp(mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + 0.5*(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2))
#   # stan_data$sd_c_t = exp(2*mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)*(exp(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)-1)
#   # stan_data$c_d_lower = mean(ExpAdjs$Dist.Lower,na.rm = T)
#   # stan_data$c_d_upper = mean(ExpAdjs$Dist.Upper,na.rm = T)
#
# }
# adjs$n_routes <- stan_data$n_routes

saveRDS(stan_data,paste0("stan_data/stan_data_",vers,"_",sp_aou,".rds"))
}
saveRDS(combined,paste0("stan_data/dataframe_",sp_aou,".rds"))
}




# preparing list of aou numbers to run on gpsc ----------------------------

for(i in 1:nrow(adjustments)){ #
  sp_sel <- adjustments[i,"common_name"]
  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)
  adjustments[i,"sp_aou"] <- sp_aou
  adjustments[i,"sp_ebird"] <- sp_ebird

  adjustments[i,"to_run_atleast_one_model"] <- ifelse(any(file.exists(paste0("stan_data/stan_data_",model_components$model,"_",sp_aou,".rds"))),
                                                      TRUE,
                                                      FALSE)
}


saveRDS(adjustments,"pop_estimates_control_file.rds")



aou_list <- adjustments |>
  filter(to_run_atleast_one_model,
         common_name %in% selected_species) |>
  select(sp_aou) |>
  distinct()

write_tsv(aou_list, "aou_list.txt",
          col_names = FALSE)

