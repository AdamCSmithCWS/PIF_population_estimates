library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(bbsBayes2)
library(exactextractr)

#source("eBirdst.R") #script to set access key for local user (not included in archive)
# download and extract eBird seasonal relative abundance within
# BBS route path buffers

if(Sys.info()[["nodename"]] == "WNCRLABN72960"){
  setwd("c:/Users/SmithAC/Documents/GitHub/PIF_population_estimates")
}

routes_buf <- readRDS("data/all_routes_buffered_400m.rds")

# run once to generate list of species with simplified english names
# and then saved as rds
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

re_calc <- TRUE # TRUE if want to re-calcluate the eBird-BBS overlap
      # e.g., when new route spatial data

use_weekly <- TRUE # only use TRUE if goal is to model all species
      # as weekly abundance, including breeding distributions


use_uncertainty <- FALSE # TRUE if teh uncertainty values of the relative abundance are desired

bbs_start <- as_date("2023-06-01") # Latest date for start of BBS-relevant breeding season
bbs_end <- as_date("2023-07-07") # Latest data for end of BBS-relevant breeding season



metric_used <- "mean" #options are "max", mean" or "median" (although "median" won't work for seasonal)





adjustments <- read_csv("Population_estimates_workflow_support.csv")





ii <- which(sps_list$english %in% adjustments$common_name)


species_sel <- c("American Robin",
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
                 "Short-eared Owl")

ii <- which(sps_list$english %in% species_sel)

for(i in ii){

  sp_sel <- unname(unlist(sps_list[i,"english"]))
species_ebird <- ebirdst::get_species(sp_sel)


if((file.exists(paste0("data/species_relative_abundance_2023/",
                      species_ebird,"_relative_abundance.rds")) &
   file.exists(paste0("data/species_relative_abundance_2023/",
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
  # Therefore: extract the weekly abundance data for the
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

saveRDS(breed_abundance,paste0("data/species_relative_abundance_2023/",species_ebird,"_derived_breeding_relative_abundance.rds"))

if(use_uncertainty){
  saveRDS(breed_abundance_sd,paste0("data/species_relative_abundance_2023/",species_ebird,"_derived_breeding_relative_abundance_logsd.rds"))
}





# # clipping ebird data to BBS strata region --------------------------------
# # this clipping is done only to speed up the overlay between BBS and ebird
#  the rest of the ebird distribution map is unaffected for the eventual
#  predictions of density across the full species' range (i.e., see
#  saveRDS() call above)
bbs_strata <- bbsBayes2::load_map("bbs") %>%
  filter(country_code %in% c("CA","US")) |>
  st_buffer(dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata, st_crs(breed_abundance))

bbs_strata2 <- bbsBayes2::load_map("bbs")%>%
  filter(country_code %in% c("CA","US"))
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



abundance_df <- exact_extract(breed_abundance,
                                       routes_buf_proj,
                                       fun = "mean",
                                       force_df = TRUE,
                                       append_cols = "route_name") |>
  rename(ebird_abund = mean)



if(use_uncertainty){

abundance_sd_in_buffers <- exact_extract(breed_abundance_sd,
                                         routes_buf_proj,
                                         fun = "mean",
                                         force_df = TRUE,
                                         append_cols = "route_name")
abundance_df <- abundance_df |>
  left_join(abundance_sd_in_buffers,
            by = "route_name")|>
  rename(ebird_abund_sd = mean)

abundance_ci_in_buffers <- exact_extract(breed_abundance_ci,
                                         routes_buf_proj,
                                         fun = "mean",
                                         force_df = TRUE,
                                         append_cols = "route_name")

abundance_df <- abundance_df |>
  left_join(abundance_sd_in_buffers,
            by = "route_name")|>
  rename(ebird_abund_ci = mean)

}

abundance_df <- abundance_df |>
  filter(
  ebird_abund > 0,
  !is.na(ebird_abund)
) |>
  distinct()

saveRDS(abundance_df,
        paste0("data/species_relative_abundance_2023/",species_ebird,"_relative_abundance.rds"))







abundance_strat <- exact_extract(breed_abundance,
                                     strata_proj,
                                     fun = "mean",
                                     force_df = TRUE,
                                     append_cols = "strata_name")|>
  rename(mean_ebird_abundance = mean)


total_abundance_in_strata <- exact_extract(breed_abundance,
                                           strata_proj,
                                           fun = "sum",
                                           force_df = TRUE,
                                           append_cols = "strata_name")|>
  rename(total_ebird_abundance = sum)

abundance_strat <- abundance_strat |>
  left_join(total_abundance_in_strata,
            by = "strata_name")

abundance_sampled <- exact_extract(breed_abundance,
                                   routes_nm,
              fun = "mean",
              force_df = TRUE,
              append_cols = "strata_name") |>
  rename(sampled_abundance = mean)

abundance_strat <- abundance_strat |>
  left_join(abundance_sampled,
            by = "strata_name") %>%
  mutate(sampled_avail_ratio = sampled_abundance/mean_ebird_abundance,
         log_ratio = log(sampled_avail_ratio))



#
# abundance_strat[nrow(abundance_strat)+1,"strata_name"] <- "USA_CAN"
# abundance_strat[nrow(abundance_strat),"mean_ebird_abundance"] <- terra::global(breed_abundance, "mean", na.rm = TRUE)
# abundance_strat[nrow(abundance_strat),"total_ebird_abundance"] <- terra::global(breed_abundance, "sum", na.rm = TRUE)



saveRDS(abundance_strat,
        paste0("data/species_relative_abundance_2023/",species_ebird,"_bbs_strata_relative_abundance.rds"))

}
