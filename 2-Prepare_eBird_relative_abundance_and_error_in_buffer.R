library(sf)
library(tidyverse)
library(ebirdst)
#ebirdst::set_ebirdst_access_key("vhhebtg39lu",overwrite = TRUE) #expire November 29 2025
#library(patchwork)
library(terra)
library(bbsBayes2)

# download and extract eBird seasonal relative abundance within
# BBS route path buffers

if(Sys.info()[["nodename"]] == "WNCRLABN72960"){
  setwd("c:/Users/SmithAC/Documents/GitHub/PIF_population_estimates")
}

routes_buf <- readRDS("data/all_routes_buffered.rds")
#routes_buf <- readRDS("all_routes_buffered_1km.rds")

#
# sps_list <- readRDS("data/all_counts_by_route.rds") %>%
#   filter(count > 0) %>%
#   select(aou,english,french,genus,species,route_name) %>%
#   distinct() %>%
#   group_by(aou,english,french,genus,species) %>%
#   summarise(n_routes_w_obs = n()) %>%
#   filter(n_routes_w_obs > 19,
#          !is.na(english),
#          !grepl(pattern = "^unid.",english),
#          !grepl(pattern = "^\\(",english),
#          english != c("Western Grebe (Clark's/Western)")) %>%
#   mutate(english = ifelse(grepl(pattern = " \\(",
#                                 english),
#                           str_extract(string = english,pattern = ".*(?= \\()"),
#                           english))
#
# saveRDS(sps_list,"data/all_bbs_data_species_list.rds")

sps_list <- readRDS("data/all_bbs_data_species_list.rds")
# eBird relative abundance ------------------------------------------------

#check ebird diagnostics

qual_ebird <- ebirdst_runs

re_calc <- FALSE # TRUE if want to re-calcluate the eBird-BBS overlap
      # e.g., when new route spatial data

use_weekly <- TRUE # only use TRUE if goal is to model all species
      # as weekly abundance, including breeding distributions


use_uncertainty <- FALSE # TRUE if teh uncertainty values of the relative abundance are desired

bbs_start <- as_date("2023-06-01") # Latest date for start of BBS-relevant breeding season
bbs_end <- as_date("2023-07-07") # Latest data for end of BBS-relevant breeding season



metric_used <- "mean" #options are "max", mean" or "median" (although "median" won't work for seasonal)

selected_species <- unique(c("Brown Creeper","Canyon Wren",
                             "Black-capped Chickadee","American Robin",
                             "Barn Swallow",
                             "Blackpoll Warbler","Baird's Sparrow",
                             "Western Meadowlark",
                             "Mountain Bluebird",
                             "Eastern Phoebe",
                             "Rose-breasted Grosbeak",
                             "Downy Woodpecker",
                             "Red-winged Blackbird",
                             "Scarlet Tanager",
                             "Say's Phoebe",
                             "Black-chinned Hummingbird","Tennessee Warbler","Bay-breasted Warbler",
                             "Swainson's Thrush","Blue Jay",
                             "Blue-headed Vireo","Bicknell's Thrush",
                             "Wilson's Snipe","American Kestrel",
                             "Steller's Jay","Clay-colored Sparrow",
                             "Killdeer","Long-billed Curlew",
                             "Yellow-rumped Warbler","Olive-sided Flycatcher",
                             "Purple Martin",
                             "Barn Swallow",
                             "Verdin","Ash-throated Flycatcher","Black-throated Sparrow",
                             "Blue Jay","Varied Thrush","Veery","Wood Thrush","Chestnut-collared Longspur",
                             "Bobolink","Savannah Sparrow","Grasshopper Sparrow",
                             "Horned Lark",
                             "Eastern Whip-poor-will", "Common Nighthawk",
                             "Ruby-throated Hummingbird",
                             "Mountain Chickadee",
                             "Northern Bobwhite",
                             "Northern Flicker",
                             "Golden-winged Warbler",
                             "Red-eyed Vireo",
                             "Palm Warbler",
                             "Clark's Nutcracker",
                             "White-throated Sparrow",
                             "Song Sparrow"))

saveRDS(selected_species,"data/selected_species.rds")
ii <- which(sps_list$english %in% selected_species)

for(i in ii){ #rev(1:nrow(sps_list))){ ##

  sp_sel <- unname(unlist(sps_list[i,"english"]))
species_ebird <- ebirdst::get_species(sp_sel)

# paste0("data/species_relative_abundance/",
#        sp_ebird,
#        "_derived_breeding_relative_abundance.rds")


if((file.exists(paste0("data/species_relative_abundance/",
                      species_ebird,"_relative_abundance.rds")) &
   file.exists(paste0("data/species_relative_abundance/",
                      species_ebird,"_derived_breeding_relative_abundance.rds"))) &
   !re_calc){
  next
}

if(is.na(species_ebird)){
  sps_list[i,"available_ebird_data"] <- "Not in ebird species list"


  next

}

qual_sel <- qual_ebird[which(qual_ebird$species_code == species_ebird),]
breed_qual <- unname(unlist(qual_sel[,"breeding_quality"]))
resident_qual <- unname(unlist(qual_sel[,"resident_quality"]))
resident <- unname(unlist(qual_sel[,"is_resident"]))

breeding_start <- unname(unlist(qual_sel[,"breeding_start"]))
breeding_end <- unname(unlist(qual_sel[,"breeding_end"]))

if(!resident){
sel_start <- ifelse(as_date(breeding_start)>bbs_start, # if stable breeding season is later than the 1st of June (mostly only northern species)
                    bbs_start,# then use the first of June
                    as_date("2023-05-20")) # otherwise use May 20th to match the earliest BBS observation dates


sel_end <- ifelse(as_date(breeding_end)>bbs_end, # if stable breeding season ending is later than the 2nd week of July
                    bbs_end,# then use the latest week of bbs season (2nd week of july)
                  as_date(breeding_end)) # otherwise use the end of the stable breeding season


}else{
  sel_start <- as_date("2023-05-20") # use May 20, as start of BBS observations season for resident species
  sel_end <- as_date("2023-07-07")# use July 7, as end of BBS observations season for resident species
  }

sps_list[i,"is_resident"] <- resident
sps_list[i,"breeding_quality"] <- breed_qual
sps_list[i,"resident_quality"] <- resident_qual
sps_list[i,"species_ebird"] <- species_ebird
sps_list[i,"breeding_start"] <- breeding_start
sps_list[i,"breeding_end"] <- breeding_end

if(!resident & breed_qual == 0){
  sps_list[i,"available_ebird_data"] <- "Only nonbreeding abundance"


  next
}


if(!resident & breed_qual > 0 & !use_weekly) {
  down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = paste0("abundance_seasonal_",metric_used,"_3km_")),
            silent = TRUE)
}

if(resident | use_weekly) {

down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_median_3km_"),
            silent = TRUE)
}

if(use_uncertainty){
down_lower <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_lower_3km_"),
            silent = TRUE)


down_upper <- try(ebirdst::ebirdst_download_status(species_ebird,
                                                   download_ranges = FALSE,
                                                   download_abundance = TRUE,
                                                   download_occurrence = FALSE,
                                                   force = FALSE,
                                                   pattern = "abundance_upper_3km_"),
                  silent = TRUE)


}
if(class(down) == "try-error"){

  sps_list[i,"available_ebird_data"] <- "Download fail"


  next}

sps_list[i,"available_ebird_data"] <- "Yes"





if(resident | use_weekly){

  # The resident relative abundance is not linearly related
  # to the BBS mean observations (log-linear) (e.g., BCCH)
  # Therefore: if resident, then extract the weekly abundance data for the
  # BBS survey season and generate a breeding season abundance



  abd_weekly_abundance <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "weekly",
                                               product = "abundance")

  sel_weeks <- try(as_date(names(abd_weekly_abundance)),silent = TRUE)
  sel_weeks <- sel_weeks[which(sel_weeks > sel_start &
                                 sel_weeks < sel_end)]

  abd_weekly_abundance_bbs_season <- abd_weekly_abundance[[as.character(sel_weeks)]]

  if(metric_used == "max"){
  breed_abundance <- terra::quantile(abd_weekly_abundance_bbs_season,
                                     probs = 1, na.rm = TRUE)
  }
  if(metric_used == "median"){
    breed_abundance <- terra::median(abd_weekly_abundance_bbs_season, na.rm = TRUE)
  }
  if(metric_used == "mean"){
    breed_abundance <- terra::mean(abd_weekly_abundance_bbs_season, na.rm = TRUE)
  }




# SD of abundance ---------------------------------------------------------
  if(use_uncertainty){

  abd_weekly_abundance_lower <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "weekly",
                                               product = "abundance",
                                               metric = "lower")

  abd_weekly_abundance_upper <- ebirdst::load_raster(species = species_ebird,
                                                     resolution = "3km",
                                                     period = "weekly",
                                                     product = "abundance",
                                                     metric = "upper")


  abd_weekly_abundance_lower_bbs_season <- abd_weekly_abundance_lower[[as.character(sel_weeks)]]
  abd_weekly_abundance_upper_bbs_season <- abd_weekly_abundance_upper[[as.character(sel_weeks)]]

  abd_weekly_abundance_ci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season

  abd_weekly_abundance_lower_bbs_season <- log(abd_weekly_abundance_lower[[as.character(sel_weeks)]]+0.0001)
  abd_weekly_abundance_upper_bbs_season <- log(abd_weekly_abundance_upper[[as.character(sel_weeks)]]+0.0001)

  abd_weekly_abundance_logci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season



  ci_to_sd <- (2 * qnorm(p = 0.9))
  abd_weekly_abundance_logsd <- abd_weekly_abundance_logci/ci_to_sd
  # if(metric_used == "max"){
  #   breed_abundance <- terra::max(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "median"){
  #   breed_abundance <- terra::median(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "mean"){
  breed_abundance_sd <- terra::mean(abd_weekly_abundance_logsd, na.rm = TRUE)

  breed_abundance_ci <- terra::mean(abd_weekly_abundance_ci, na.rm = TRUE)

  }

}else{ # if breeding not resident


  abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                 resolution = "3km",
                                                 period = "seasonal",
                                                 product = "abundance",
                                                 metric = metric_used)  #3km high resolution

  breed_abundance <- abd_seasonal_abundance[["breeding"]]





  # SD of abundance ---------------------------------------------------------
  if(use_uncertainty){

  abd_weekly_abundance_lower <- ebirdst::load_raster(species = species_ebird,
                                                     resolution = "3km",
                                                     period = "weekly",
                                                     product = "abundance",
                                                     metric = "lower")

  abd_weekly_abundance_upper <- ebirdst::load_raster(species = species_ebird,
                                                     resolution = "3km",
                                                     period = "weekly",
                                                     product = "abundance",
                                                     metric = "upper")

  sel_weeks <- try(as_date(names(abd_weekly_abundance_upper)),silent = TRUE)
  sel_weeks <- sel_weeks[which(sel_weeks > sel_start &
                                 sel_weeks < sel_end)]

  abd_weekly_abundance_lower_bbs_season <- abd_weekly_abundance_lower[[as.character(sel_weeks)]]
  abd_weekly_abundance_upper_bbs_season <- abd_weekly_abundance_upper[[as.character(sel_weeks)]]

    abd_weekly_abundance_ci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season

    abd_weekly_abundance_lower_bbs_season <- log(abd_weekly_abundance_lower[[as.character(sel_weeks)]]+0.0001)
    abd_weekly_abundance_upper_bbs_season <- log(abd_weekly_abundance_upper[[as.character(sel_weeks)]]+0.0001)

    abd_weekly_abundance_logci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season



  ci_to_sd <- (2 * qnorm(p = 0.9))
  abd_weekly_abundance_logsd <- abd_weekly_abundance_logci/ci_to_sd
  # if(metric_used == "max"){
  #   breed_abundance <- terra::max(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "median"){
  #   breed_abundance <- terra::median(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "mean"){
  breed_abundance_sd <- terra::mean(abd_weekly_abundance_logsd, na.rm = TRUE)

  breed_abundance_ci <- terra::mean(abd_weekly_abundance_ci, na.rm = TRUE)
}
 #}


}

saveRDS(breed_abundance,paste0("data/species_relative_abundance/",species_ebird,"_derived_breeding_relative_abundance.rds"))

if(use_uncertainty){
  saveRDS(breed_abundance_sd,paste0("data/species_relative_abundance/",species_ebird,"_derived_breeding_relative_abundance_logsd.rds"))
}





# # clipping ebird data to BBS strata region --------------------------------
# # this clipping is done only to speed up the overlay between BBS and ebird
#  the rest of the ebird distribution map is unaffected for the eventual
#  predictions of density across the full species' range (i.e., see
#  saveRDS() call above)
bbs_strata <- bbsBayes2::load_map("bbs_usgs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata, st_crs(breed_abundance))

bbs_strata2 <- bbsBayes2::load_map("bbs_usgs")
strata_proj <- st_transform(bbs_strata2, st_crs(breed_abundance))

routes_nm <- routes_buf %>%
  select(route_name,route) %>%
  st_transform(crs = st_crs(bbs_strata2)) %>%
  sf::st_join(bbs_strata2, largest = TRUE) %>%
  group_by(strata_name) %>%
  summarise() %>%
  sf::st_transform(crs = st_crs(breed_abundance))
# crop and mask to boundary of BBS
breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
  mask(region_boundary_proj)
if(use_uncertainty){

breed_abundance_sd <- crop(breed_abundance_sd, region_boundary_proj) |>
  mask(region_boundary_proj)

breed_abundance_ci <- crop(breed_abundance_ci, region_boundary_proj) |>
  mask(region_boundary_proj)
}
# project buffers to match raster data
routes_buf_proj <- st_transform(routes_buf, st_crs(breed_abundance))



# extract mean abundance in BBS route buffers

abundance_in_buffers <- terra::extract(breed_abundance,
                                       routes_buf_proj,
                                       fun = mean,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

abundance_in_strata <- terra::extract(breed_abundance,
                                       strata_proj,
                                       fun = mean,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

total_abundance_in_strata <- terra::extract(breed_abundance,
                                      strata_proj,
                                      fun = sum,
                                      na.rm = TRUE,
                                      ID = FALSE,
                                      exact = TRUE)
if(use_uncertainty){

abundance_sd_in_buffers <- terra::extract(breed_abundance_sd,
                                       routes_buf_proj,
                                       fun = mean,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

abundance_ci_in_buffers <- terra::extract(breed_abundance_ci,
                                          routes_buf_proj,
                                          fun = mean,
                                          na.rm = TRUE,
                                          ID = FALSE,
                                          exact = TRUE)
}
abundance_sampled <- terra::extract(breed_abundance,
                                       routes_nm,
                                       fun = mean,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

# this argument creates an area-weighted mean of the cells that are overlapped
# by the buffer
# plot(breed_abundance)
# plot(routes_buf_proj,add = TRUE)
# plot(breed_abundance,add = TRUE)

# drop NA routes that are outside of the eBird prediction surface -----------------------

# these NA routes could be dropped because they don't represent
# much when it comes to scaling to the surface...

# may as well drop routes with no estimated abundance as they are outside
# of the range of the modeled seasonal distribution
if(use_uncertainty){

abundance_df <- data.frame(route_name = routes_buf_proj$route_name,
                             ebird_abund = abundance_in_buffers[[1]],
                           ebird_abund_logsd = abundance_sd_in_buffers[[1]],
                           ebird_abund_ci = abundance_ci_in_buffers[[1]]) %>%
  filter(
    ebird_abund > 0,
    !is.na(ebird_abund)
    ) %>%
  distinct()

}else{
  abundance_df <- data.frame(route_name = routes_buf_proj$route_name,
                             ebird_abund = abundance_in_buffers[[1]]) %>%
    filter(
      ebird_abund > 0,
      !is.na(ebird_abund)
    ) %>%
    distinct()

}

# ## displaying routes with no abundance data
# tstrt <- routes_buf %>%
#   inner_join(.,abundance_df) %>%
#   filter(is.na(ebird_abund))
#
# tst <- ggplot(data = tstrt)+
#   geom_sf(data = bbs_strata)+
#   geom_sf(aes(colour = ebird_abund))
# tst

saveRDS(abundance_df,
        paste0("data/species_relative_abundance/",species_ebird,"_relative_abundance.rds"))

abundance_strat <- bbs_strata2 %>%
  sf::st_drop_geometry() %>%
  mutate(mean_ebird_abundance = abundance_in_strata[[1]],
         total_ebird_abundance = total_abundance_in_strata[[1]])

abundance_strat[nrow(abundance_strat)+1,"strata_name"] <- "USA_CAN"
abundance_strat[nrow(abundance_strat),"mean_ebird_abundance"] <- terra::global(breed_abundance, "mean", na.rm = TRUE)
abundance_strat[nrow(abundance_strat),"total_ebird_abundance"] <- terra::global(breed_abundance, "sum", na.rm = TRUE)


sampled_abund <- data.frame(strata_name = routes_nm$strata_name,
                            sampled_abundance = abundance_sampled$mean)


abundance_strat <- abundance_strat %>%
  left_join(sampled_abund,
            by = "strata_name") %>%
  mutate(sampled_avail_ratio = sampled_abundance/mean_ebird_abundance,
         log_ratio = log(sampled_avail_ratio))

saveRDS(abundance_strat,
        paste0("data/species_relative_abundance/",species_ebird,"_bbs_strata_relative_abundance.rds"))

}
