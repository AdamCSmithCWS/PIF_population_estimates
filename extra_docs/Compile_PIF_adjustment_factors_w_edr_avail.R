
library(tidyverse)
#remotes::install_github("bbsBayes/bbsBayes2") # allows access to BBS raw data
library(bbsBayes2)
#remotes::install_github("na-pops/napops") # allows access to detectability estimates
# for many North American landbirds: https://doi.org/10.1111/ibi.13169
library(napops)


# combine PIF adjustment factors with na-pops estimates --------------------------------------------------------


Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")
napops_species <- napops::list_species() |>
  dplyr::select(Species,Common_Name)

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
# csv file in next line is from the supplement to Stanton et al. 2019
# https://doi.org/10.5751/ACE-01331-140104 - appendix-2
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE) %>%
  mutate(cn = ifelse(cn == "McCown's Longspur", #updating some species names that have changed since 2019
                     "Thick-billed Longspur",
                     cn),
         cn = ifelse(cn == "Gray Jay",
                     "Canada Jay",
                     cn),
         cn = ifelse(grepl("Le Conte's",cn), # alternate spelling of LeConte
                     gsub("Le Conte's","LeConte's",cn),
                     cn))



# appendix with species names from Rosenberg et --------



bbs_sci_names <- bbsBayes2::load_bbs_data()[["species"]] %>%
  filter(unid_combined) %>%
  select(aou, english,family,genus, order) %>%
  mutate(english = gsub(pattern = " (all forms)",
                               x = english,
                               replacement = "", fixed = TRUE))

# append with additional species PIF adjustments from Rosenberg et --------

Rosen_pif <- read_csv("data/Rosenberg_appendix_2_aaw1313_data_s2.csv") %>%
  rename(Dist2 = `Distance Adj.`,
         Pair2 = `Pair Adj.`,
         cn = Species) %>%
  mutate(cn = ifelse(cn == "McCown's Longspur", #updating some species names that have changed since 2019
                                            "Thick-billed Longspur",
                                            cn),
                                cn = ifelse(cn == "Gray Jay",
                                            "Canada Jay",
                                            cn),
                                cn = ifelse(grepl("Le Conte's",cn), # alternate spelling of LeConte
                                            gsub("Le Conte's","LeConte's",cn),
                                            cn)) %>%
  #filter(!cn %in% ExpAdjs$cn) %>%
  select(-sci_name)

ExpAdjs <- ExpAdjs %>%
  filter(!cn %in% Rosen_pif$cn)

ExpAdjs <- ExpAdjs %>%
  bind_rows(Rosen_pif)




# Additional TOD adjustments calculated for any species > 100 BBS  --------
# from "tod_calculations.R"


tod_adjs <- read_csv("Additional_time_of_day_adjustments.csv") %>%
  filter(!cn %in% ExpAdjs$cn)

ExpAdjs <- ExpAdjs %>%
  bind_rows(tod_adjs) %>%
  left_join(bbs_sci_names,
            by = c("cn" = "english"))



# Load BBS data -----------------------------------------------------------


# general summaries of time and season for BBS surveys since 2013
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

mid_time = round(surveys$mean_time,1) # mean time of day for bbs counts
mean_doy = round(surveys$mean_doy) # mean day of year for bbs counts




ExpAdjs <- ExpAdjs |>
  dplyr::left_join(napops_species,
                   by = c("cn" = "Common_Name"))

# First add EDR (effective detection radius) ------------------------------


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

## add the sd of EDR
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
                                      model = 1, #intercept model if something failed for "best" model
                                      od = mean_doy, # mean of the doy for all BBS surveys
                                      tssr = mid_time, #mean of time of day for all BBS surveys
                                      quantiles = c(0.1625,0.8375), # 1 sd below the mean
                                      samples = 1000),
                        silent = TRUE)
    w_mod <- 1
  }else{
    if(availability$p_est < 0.1){ # again defaulting to intercept model if estimates are unrealistically small
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

  ## Three species-specific judgements of estimates that are unrealistic and should not be used.
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



# generate combined adjustment scores for each species --------------------


new_adjs <- ExpAdjs %>%
  select(cn,Aou, Species, TimeAdj.meanlog,Dist2,Pair2,
         edr,availability,
         order,
         family,
         genus,
         Species)

family_scores <- new_adjs %>%
  group_by(family) %>%
  summarise(family_TimeAdj.meanlog = mean(TimeAdj.meanlog, na.rm = TRUE),
            family_Dist2 = mean(Dist2, na.rm = TRUE),
            family_Pair2 = mean(Pair2, na.rm = TRUE))
genus_scores <- new_adjs %>%
  group_by(genus) %>%
  summarise(genus_TimeAdj.meanlog = mean(TimeAdj.meanlog, na.rm = TRUE),
            genus_Dist2 = mean(Dist2, na.rm = TRUE),
            genus_Pair2 = mean(Pair2, na.rm = TRUE))
order_scores <- new_adjs %>%
  group_by(order) %>%
  summarise(order_TimeAdj.meanlog = mean(TimeAdj.meanlog, na.rm = TRUE),
            order_Dist2 = mean(Dist2, na.rm = TRUE),
            order_Pair2 = mean(Pair2, na.rm = TRUE))


new_adjs <- new_adjs %>%
  left_join(family_scores,
            by = "family") %>%
  left_join(genus_scores,
            by = "genus")%>%
  left_join(order_scores,
            by = "order")

missing_dist <- new_adjs %>%
  filter(is.na(order_Dist2))


missing_pair <- new_adjs %>%
  filter(is.na(order_Pair2))

missing_tod <- new_adjs %>%
  filter(is.na(order_TimeAdj.meanlog))

new_adjs_fixed <- new_adjs %>%
  mutate(Pair2_final = ifelse(is.na(Pair2),
                        genus_Pair2,
                        Pair2),
         Dist2_final = ifelse(is.na(Dist2),
                        genus_Dist2,
                        Dist2))

new_adjs_fixed <- new_adjs_fixed %>%
  mutate(Pair2_final = ifelse(is.na(Pair2_final),
                        family_Pair2,
                        Pair2_final),
         Dist2_final = ifelse(is.na(Dist2_final),
                        family_Dist2,
                        Dist2_final))

new_adjs_fixed <- new_adjs_fixed %>%
  mutate(Pair2_final = ifelse(is.na(Pair2_final),
                        order_Pair2,
                        Pair2_final),
         Dist2_final = ifelse(is.na(Dist2_final),
                        order_Dist2,
                        Dist2_final))

new_adjs_fixed <- new_adjs_fixed %>%
  mutate(Pair2_final = ifelse(cn == "Wood Stork",
                              1,
                              Pair2_final),
         Dist2_final = ifelse(cn == "Wood Stork",
                              400,
                              Dist2_final),
         TimeAdj.meanlog_final = TimeAdj.meanlog)
# 1/availability * 1/(EDR^2*3.14159) * pair
# exp(TimeAdj.meanlog) * 1/(Dist2^2*3.14159) * pair


PIF_adjs <- new_adjs_fixed %>%
  mutate(C_t = exp(TimeAdj.meanlog_final),
         C_t_avail = (1/availability),
         C_d = (1/(50*(Dist2_final/1000)^2*3.14159)), # area within detection distance of 50 BBS stops (square km)
         C_d_edr = (1/(50*(edr/1000)^2*3.14159)),# area within detection distance of 50 BBS stops (square km)
         C_p = Pair2_final,
         PIF_adj_birds_km2 = C_t * C_p * C_d, # species-level adjustment using the PIF factors from Rosenberg et al. 2019 (birds/km^2)
         PIF_adj_napops_birds_km2 = C_t_avail * C_p * C_d_edr, # species-level adjustments using updated distance and availability estimates from na-pops project (birds/km^2)
         PIF_adj_napops_birds_km2 = ifelse(is.na(PIF_adj_napops_birds_km2),PIF_adj_birds_km2,PIF_adj_napops_birds_km2), # species-level adjustments mix of na-pops and PIF if missing na-pops.
         PIF_napops_compare = PIF_adj_napops_birds_km2/PIF_adj_birds_km2)



write_csv(PIF_adjs,"PIF_adjustments_and_napops_edr_avail.csv")


