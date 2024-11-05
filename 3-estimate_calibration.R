
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


# Load BBS data -----------------------------------------------------------


# mean counts to id routes on which species has been observed -------------
mean_counts <- readRDS("data/mean_counts_by_route.rds") %>%
  filter(english == sp_sel,
         mean_count > 0)

raw_counts <- readRDS("data/all_counts_by_route.rds") %>%
  filter(english == sp_sel)


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

mid_time = round(surveys$mean_time,1)
mid_doy = round(surveys$mean_doy)




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
    break(paste("edr extraction did not work for",sp))

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

  # choose better supported model between models 1 (intercept) and 2 (roadside)
  cefs <- napops::coef_removal(species = sp) |>
    #dplyr::filter(Model %in% c(1,2)) |>
    dplyr::arrange(AIC)

  if(nrow(cefs) == 0){next} # skip if no napops output

  w_mod <- cefs[1,"Model"]

  cue <- try(napops::cue_rate(species = sp,
                          model = w_mod,
                          od = mid_doy, # mean of the doy for all BBS surveys
                          tssr = mid_time, #mean of time of day for all BBS surveys
                          quantiles = c(0.1625), # 1 sd below the mean
                          samples = 1000),
              silent = TRUE)

  # if(cue$EDR_est > 500){ # two species where the model 2 estimates of edr are silly (> 5 km!)
  #   edrt <- try(napops::cue_rate(species = sp,
  #                        model = 1,
  #                        od = mid_doy, # mean of the doy for all BBS surveys
  #                        tssr = mid_time, #mean of time of day for all BBS surveys
  #                        quantiles = c(0.1625), # 1 sd below the mean
  #                        samples = 1000),
  #               silent = TRUE)
  #   w_mod <- 1
  # }

  if(class(cue) == "try-error"){
    break(paste("cue rate extraction did not work for",sp))

  }


  ExpAdjs[i,"cue_rate"] <- as.numeric(cue$CR_est)

  if("CR_16.25" %in% names(cue)){
    ExpAdjs[i,"cue_rate_sd"] <- as.numeric(cue$CR_est)-as.numeric(cue$CR_16.25)
  }

  ExpAdjs[i,"EDR_model"] <- w_mod

  rm(cue, w_mod)
}

ExpAdjs <- ExpAdjs |>
  dplyr::mutate(use_cue_rate = ifelse(is.na(cue_rate),FALSE,TRUE))


write_csv(ExpAdjs,"Species_correction_factors_w_edr_cue.csv")
# route overlap info ------------------------------------------------------
}else{
  ExpAdjs <- read_csv("Species_correction_factors_w_edr_cue.csv")
}

routes_buf_all <- readRDS("data/all_routes_buffered.rds")

strata <- bbsBayes2::load_map("bbs_usgs")

strata_names <- strata %>%
  sf::st_drop_geometry()




# potentially important detail from Golden Eagle paper --------------------------------------
##  USE MAX WEEKLY RELATIVE ABUNDANCE DURING KEY SEASON
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










# Running full model using eBird relative abundance ------------------------------------------------

#
sps_list <- readRDS("data/all_bbs_data_species_list.rds") #

sps_sel <- c("Tennessee Warbler","Bay-breasted Warbler",
             "Swainson's Thrush","Blue Jay",
             "Blue-headed Vireo","Bicknell's Thrush",
             "Wilson's Snipe","American Kestrel",
             "Steller's Jay")
# species loop -----------------------------------------------------

for(sp_sel in sps_sel){#list$english){

sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]


strata_trad <- strata_trad_all %>%
  filter(cn == sp_sel)

global_trad <- global_trad_all %>%
  filter(cn == sp_sel)

bcr_trad <- bcr_trad_all %>%
  filter(cn == sp_sel)





# filtering routes on which species wasn't observed.
routes_buf <- routes_buf_all %>%
  filter(route_name %in% mean_counts$route_name)

sp_ebird <- ebirdst::get_species(sp_sel)

if(!file.exists(paste0("data/species_relative_abundance/",
                       sp_ebird,"_relative_abundance.rds"))){next}

rel_abund <- readRDS(paste0("data/species_relative_abundance/",
                            sp_ebird,"_relative_abundance.rds"))

combined <- rel_abund %>%
inner_join(.,raw_counts,
           by = "route_name") %>%
  mutate(logm = log(ebird_abund),
         route = as.integer(factor(route_name)))



# Try stan model ----------------------------------------------------------


mean_abund <- combined %>%
  select(route,logm) %>%
  distinct() %>%
  arrange(.,route)


#adjustment scores
adjs <- ExpAdjs %>%
  filter(cn == sp_sel)



library(cmdstanr)

model <- cmdstanr::cmdstan_model("models/ebird_rel_abund_calibration.stan")


stan_data <- list(n_routes = max(combined$route),
                  n_counts = nrow(combined),
                  count = combined$count,
                  route = combined$route,
                  log_mean_rel_abund = mean_abund$logm,
                  c_p = adjs$Pair2,
                  sd_c_p = 0.13,
                  c_d = ifelse(adjs$use_edr,adjs$edr,1),
                  sd_c_d = ifelse(adjs$use_edr,adjs$edr_sd,1),
                  c_t = ifelse(adjs$use_cue_rate,
                               adjs$cue_rate,
                               exp(adjs$TimeAdj.meanlog + 0.5*(adjs$TimeAdj.sdlog^2))),
                  sd_c_t = ifelse(adjs$use_cue_rate,
                                  adjs$cue_rate_sd,
                                  exp(2*adjs$TimeAdj.meanlog + adjs$TimeAdj.sdlog^2)*(exp(adjs$TimeAdj.sdlog^2)-1)),
                  c_d_lower = adjs$Dist.Lower,
                  c_d_upper = adjs$Dist.Upper,
                  edr = ifelse(adjs$use_edr,1,0),
                  cue = ifelse(adjs$use_cue_rate,1,0),
                  use_pois = 0,
                  use_pair = 0,
                  use_t = 0
                  )


if(nrow(adjs) == 0){
  # stan_data$c_p = mean(ExpAdjs$Pair2,na.rm = T)
  # stan_data$c_t = exp(mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + 0.5*(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2))
  # stan_data$sd_c_t = exp(2*mean(ExpAdjs$TimeAdj.meanlog,na.rm = T) + mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)*(exp(mean(ExpAdjs$TimeAdj.sdlog,na.rm = T)^2)-1)
  # stan_data$c_d_lower = mean(ExpAdjs$Dist.Lower,na.rm = T)
  # stan_data$c_d_upper = mean(ExpAdjs$Dist.Upper,na.rm = T)

}

fit <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500,
                    iter_warmup = 2000,
                    iter_sampling = 2000,
                    show_exceptions = FALSE)



quant <- function(x,q){
  y <- quantile(x,q,names = FALSE,digits = 3)
  return(y)
}

summ <- fit$summary()

saveRDS(summ,paste0("convergence/parameter_summary_",sp_aou,"_",sp_ebird,".rds"))

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

cali_post <- fit$draws(variables = "calibration",
                       format = "df")


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

cali_alt_post <- fit$draws(variables = "calibration_alt",
                       format = "df")

param_infer <- bind_rows(cali_alt,
                         cali,
                         avail_correction_realised,
                         p_avail_realised,
                         edr_realised,
                         surveyed_area_realised) %>%
  mutate(english = sp_sel,
         sp_eBird = sp_ebird,
         aou = sp_aou)

saveRDS(param_infer,paste0("output/parameter_inference_",sp_aou,"_",sp_ebird,".rds"))

# BETA_post <- fit$draws(variables = "BETA",
#                        format = "df")

saveRDS(fit,paste0("output/calibration_fit_",sp_aou,"_",sp_ebird,".rds"))



# loading the breeding season abundance raster generated in
# 2_Prepare_eBird_relative_abundance_in_buffer.R
#
breed_abundance <- readRDS(paste0("data/species_relative_abundance/",
                                  species_ebird,
                                  "_derived_breeding_relative_abundance.rds"))



breed_abundance_full <- breed_abundance

# # clipping ebird data to BBS strata region --------------------------------
bbs_strata_buf <- bbsBayes2::load_map("bbs_usgs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata_buf, st_crs(breed_abundance))
# crop and mask to boundary of BBS
breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
  mask(region_boundary_proj)


strata_proj <- st_transform(strata,
                            crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_strata <- terra::extract(breed_abundance,
                                      strata_proj,
                                       fun = sum,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = FALSE)



## function to estimate full posterior of abundance

post_abund <- function(x,draws = cali_post$calibration,
                       fun = "mean",
                       p = 0.025){
  y <- vector("numeric",length(x))

  if(fun == "mean"){
  for(i in 1:length(x)){
  y[i] <- mean(x[i]*as.numeric(draws))
  }
  }

  if(fun == "quantile"){
    for(i in 1:length(x)){
      y[i] <- quantile(x[i]*as.numeric(draws),p, names = FALSE)
    }
  }

  return(y)
}
### 3^2 scaling is to account for the area of each grid-cell
###
###

strata_abund <- data.frame(region = strata_proj$strata_name,
                           region_type = "strata",
                           sum_abund = 3^2 * unlist(abundance_in_strata[[1]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund,
                               draws = cali_alt_post$calibration_alt,
                                 fun = "mean"),
         pop_median = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
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

comp_plot


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
       subtitle = "Diagonal lines = 1:1 and 2:1")

png(filename = paste0("Figures/comp_trad_new_",sp_aou,"_",sp_ebird,".png"),
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
  mutate(.,breeding = ifelse(breeding == 0,NA,breeding))

abund_map <- ggplot()+
  geom_sf(data = strata, fill = NA)+
  geom_spatraster(data = breed_abundance_plot,
                  maxcell = 16000000)+
  geom_sf(data = strata, fill = NA)+
  geom_sf(data = routes_buf,fill = "darkred",colour = "darkred")+
  coord_sf(xlim = c(bb[c("xmin","xmax")]),
           ylim = c(bb[c("ymin","ymax")]))+
  #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA,
                       end = 0.98)+
  theme_bw()+
  labs(title = paste(sp_sel,"eBird relative abundance values by BBS strata"))


#

png(filename = paste0("Figures/abund_map_",sp_aou,"_",sp_ebird,".png"),
    res = 400,
    height = 7,
    width = 6.5,
    units = "in")
print(abund_map)
dev.off()


saveRDS(comp_trad_new_plot,paste0("figures/saved_ggplots/trad_vs_new_",sp_aou,"_",sp_ebird,".rds"))
saveRDS(abund_map,paste0("figures/saved_ggplots/abund_map_",sp_aou,"_",sp_ebird,".rds"))




# National summaries ------------------------------------------------------




countries <- bbsBayes2::load_map("bbs_usgs") %>%
  group_by(country) %>%
  summarise()


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
                                      ID = FALSE,
                                      exact = FALSE)



country_abund <- data.frame(region = countries_proj$country,
                            region_type = "country",
                           sum_abund = 3^2 * unlist(abundance_in_countries[[1]])) %>%
  mutate(pop_mean = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

country_abund



# BCR summaries -----------------------------------------------------------





bcrs <- bbsBayes2::load_map("bbs_usgs") %>%
  group_by(bcr) %>%
  summarise()


bcr_proj <- st_transform(bcrs,
                               crs = st_crs(breed_abundance))


abundance_in_bcr <- terra::extract(breed_abundance,
                                   bcr_proj,
                                         fun = sum,
                                         na.rm = TRUE,
                                         ID = FALSE,
                                         exact = FALSE)



bcr_abund <- data.frame(region = as.character(bcr_proj$bcr),
                            region_type = "bcr",
                            sum_abund = 3^2 * unlist(abundance_in_bcr[[1]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

bcr_abund




# USACAN estimates -----------------------------------------


USACAN <- bbsBayes2::load_map("bbs_usgs") %>%
  ungroup() %>%
  summarise()


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
                                         ID = FALSE,
                                         exact = FALSE)



USACAN_abund <- data.frame(region = "USACAN",
                          region_type = "USACAN",
                            sum_abund = 3^2 * unlist(abundance_in_USACAN[[1]])) %>%
  mutate(pop_mean = post_abund(sum_abund, fun = "mean"),
         pop_median = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

USACAN_abund




abundance_in_global <- sum(values(breed_abundance_full),na.rm = TRUE)



global_abund <- data.frame(region = "global",
                           region_type = "global",
                           sum_abund = 3^2 * abundance_in_global) %>%
  mutate(pop_mean = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,draws = cali_alt_post$calibration_alt,
                                 fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

global_abund



pop_ests_out <- bind_rows(USACAN_abund,
                          global_abund,
                          country_abund,
                          bcr_abund,
                          strata_abund) %>%
  mutate(version = "PIF_eBird")

# Sum of population estimates

write_excel_csv(pop_ests_out,
                paste0("estimates/pop_ests_",sp_aou,sp_ebird,".csv"))


pop_ests_out_trad_sel <- pop_ests_out_trad %>%
  filter(species == sp_sel)


pop_compare_stack <- pop_ests_out %>%
  select(region,region_type,version,
         species,species_ebird,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95) %>%
  bind_rows(pop_ests_out_trad_sel)


pop_compare_stack_sel <- pop_compare_stack %>%
  filter(region_type %in% c("USACAN","global","bcr"))


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


saveRDS(side_plot,paste0("figures/saved_ggplots/side_plot_",sp_aou,"_",sp_ebird,".rds"))


pdf(paste0("figures/trad_calib_compare_by_strata_",sp_aou,"_",sp_ebird,".pdf"),
    width = 11,
    height = 8.5)
print(comp_trad_new_plot)
print(abund_map)
print(side_plot)
dev.off()



 } #end of species loop
#



cali_rt <- summ %>% filter(grepl("calibration_r",variable))
hist(log(cali_rt$mean))

