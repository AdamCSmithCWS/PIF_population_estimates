## list all possible sources of adjustments

library(tidyverse)
library(napops)
library(bbsBayes2)
library(naturecounts)
library(ebirdst)
library(terra)
library(sf)
library(BAMexploreR)



qual_ebird <- ebirdst_runs



# species lists -----------------------------------------------------------

sp_drop <- c("Hoary Redpoll","Pacific-slope Flycatcher",
             "Northwestern Crow")


rosenberg <- read_csv("data/Rosenberg_appendix_2_aaw1313_data_s2.csv") |>
  rename(cn = Species,
         Dist2 = `Distance Adj.`,
         Pair2 = `Pair Adj.`) |>
  select(cn,TimeAdj.meanlog,TimeAdj.sdlog,
         Dist2,Pair2)  |>
  mutate(cn = ifelse(cn == "McCown's Longspur",
                                          "Thick-billed Longspur",
                                          cn),
                              cn = ifelse(cn == "Gray Jay",
                                          "Canada Jay",
                                          cn),
                              cn = ifelse(grepl("Le Conte's",cn),
                                          gsub("Le Conte's","LeConte's",cn),
                                          cn),
         cn = ifelse(cn == "Mew Gull",
                     "Short-billed Gull",
                     cn),
         cn = ifelse(grepl("Night-Heron",cn),
                     gsub("Night-Heron","Night Heron",cn),
                     cn)) |>
  filter(!cn %in% sp_drop)


Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")

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
                     cn))|>
  filter(!cn %in% sp_drop)

sp_stanton <- unique(ExpAdjs$cn)

sp_new <- rosenberg |>
  filter(!cn %in% sp_stanton)
#
# for(i in 1:nrow(sp_new)){
#   ss <- as.character(sp_new[i,"cn"])
#   sp <- try(naturecounts::search_species(name = ss, show = "all", authority = "BSCDATA") |>
#     filter(species_status == 1,
#            !is.na(BSCDATA)),
#     silent = TRUE)
#
#   if(class(sp)[1] == "try-error"){
#     next
#   }
#   if(nrow(sp) > 1){
#     sp <- sp[1,]
#   }
#   sp_new[i,"Species"] <- as.character(sp$BSCDATA)
# }

ExpAdjs <- ExpAdjs |>
  bind_rows(sp_new)




ExpAdjs <- ExpAdjs |>
  mutate(ebird_species = ebirdst::get_species(cn))


cn_change <- ExpAdjs |>
  filter(is.na(ebird_species) )


for(i in 1:nrow(cn_change)){
  sp <- as.character(cn_change[i,"cn"])

if(sp == "Blue-throated Hummingbird"){
  ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("Blue-throated Mountain-gem")
  ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "Blue-throated Mountain-gem"
next
}

  if(sp == "Common Redpoll"){
    ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("Redpoll")

    ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "Redpoll"
    next
  }
  if(sp == "Magnificent Hummingbird"){
    ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("Rivoli's Hummingbird")

    ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "Rivoli's Hummingbird"
    next
  }

  if(sp == "Cordilleran Flycatcher"){
    ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("Western Flycatcher")

    ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "Western Flycatcher"
    next
  }


  if(sp == "Northern Goshawk"){
    ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("American Goshawk")

    ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "American Goshawk"
    next
  }

  if(sp == "House Wren"){

    ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("Northern/Southern House Wren")
    ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "Northern House Wren"

    next
  }




  if(sp == "Western Scrub-Jay"){
    tmp <- ExpAdjs[which(ExpAdjs$cn == sp),]
    tmp["cn"] <- "Woodhouse's Scrub-Jay"
    tmp["ebird_species"] <- get_species("Woodhouse's Scrub-Jay")
    ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species("California Scrub-Jay")
    ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- "California Scrub-Jay"
    ExpAdjs[nrow(ExpAdjs)+1,] <- tmp
    next
  }




  sp_alt <- naturecounts::search_species(sp, authority = "BSCDATA") |>
    filter(!is.na(BSCDATA)) |>
    slice_head(n = 1) |>
    select(english_name) |>
    unlist() |>
    unname()
  ExpAdjs[which(ExpAdjs$cn == sp),"ebird_species"] <- get_species(sp_alt)
  ExpAdjs[which(ExpAdjs$cn == sp),"cn"] <- sp_alt

  }



# ExpAdjs <- ExpAdjs |>
#   mutate(ebird_species = ebirdst::get_species(cn))


cn_change <- ExpAdjs |>
  filter(is.na(ebird_species) )


qual_ebird_link <- qual_ebird |>
  select(species_code,
         scientific_name,
         common_name,
         is_resident,
         breeding_quality,
         breeding_start,
         breeding_end,
         resident_quality,resident_start,resident_end) |>
  rename_with(.fn = ~paste0(.x,"_ebirdst"))


ExpAdjs <- ExpAdjs |>
  left_join(qual_ebird_link,
            by = c("ebird_species" = "species_code_ebirdst"))

bam <- readxl::read_excel("data/BAMV5-results.xlsx", sheet = "abundances") |>
  filter(region == "Canada",
         year == 2020) |>
  mutate(english = ifelse(grepl("Le Conte's",english),
              gsub("Le Conte's","LeConte's",english),
              english),
         english = ifelse(english == "House Wren",
                          "Northern House Wren",
                          english))



# BBS data to support availability ----------------------------------------


# general summaries of time and season
surveys <- load_bbs_data(release = 2026)$routes %>%
  filter(year > 2014)  %>%  # 10 years from 2013 - 2023, missing 2020.
  mutate(doy = lubridate::yday(paste(year,month,day,sep = "-")),
         start_hr = floor(start_time/100),
         end_hr = floor(end_time/100),
         start_min = start_hr*60 + ((start_time/100)-start_hr)*100,
         end_min = end_hr*60 + ((end_time/100)-end_hr)*100,
         mid_time = ((end_min - start_min)/120)-0.5) %>%
  summarise(mean_time = mean(mid_time),
            mean_doy = mean(doy))

survey_dates <- load_bbs_data(release = 2026)$routes %>%
  filter(year > 2014)  %>%  # 10 years from 2013 - 2023, missing 2020.
  mutate(doy = lubridate::yday(paste(year,month,day,sep = "-")),
         start_hr = floor(start_time/100),
         end_hr = floor(end_time/100),
         start_min = start_hr*60 + ((start_time/100)-start_hr)*100,
         end_min = end_hr*60 + ((end_time/100)-end_hr)*100,
         mid_time = ((end_min - start_min)/120)-0.5)

mid_time = round(surveys$mean_time,1)
mean_doy = round(surveys$mean_doy)



# NA-POPS -----------------------------------------------------------------


re_napops <- TRUE # set to True to re-run the napops extraction of ERDs and Availability estimates

if(re_napops){
  napops_species <- napops::list_species() |>
    dplyr::select(Species,Common_Name)


  #
  # sp_stanton[which(!sp_stanton %in% sp_rosen)]




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

availability_models = c("Intercept",
                                                  "Season",
                                                  "Season_quadratic",
                                                  "Time of day",
                                                  "Time of day quadratic",
                                                  "Time of day and Season",
                                                  "Time of day and Season quadratic",
                                                  "Time of day quadratic and Season",
                                                  "Time of day quadratic and Season quadratic")

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

    ExpAdjs[i,"availability_model"] <- availability_models[w_mod]

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



# bam % pop in Canada -----------------------------------------------------
countries <- bbsBayes2::load_map("prov_state") |>
  mutate(north = ifelse(country_code == "CA" | prov_state == "AK",
                        "CAN_AK","Lower_48")) |>
  ungroup() |>
  group_by(country_code,north) |>
  summarise(.groups = "drop")


redo_bam_ebird_cov <- TRUE

if(redo_bam_ebird_cov){
abundance_all <- NULL
for(i in 1:nrow(bam)){

  sp <- as.character(bam[i,"english"])
  species_ebird <- ebirdst::get_species(sp)
  if(sp == "Northern House Wren"){
    species_ebird <- ebirdst::get_species("Northern/Southern House Wren")
  }

  if(is.na(species_ebird)){
    print(paste("no ebird status for",sp))
    next
  }
  bam[i,"species_ebird"] <- species_ebird

  qual_sel <- qual_ebird[which(qual_ebird$species_code == species_ebird),]
  breed_qual <- unname(unlist(qual_sel[,"breeding_quality"]))
  resident_qual <- unname(unlist(qual_sel[,"resident_quality"]))
  resident <- unname(unlist(qual_sel[,"is_resident"]))

  breeding_start <- unname(unlist(qual_sel[,"breeding_start"]))
  breeding_end <- unname(unlist(qual_sel[,"breeding_end"]))



  if(!resident & breed_qual == 0){
    bam[i,"available_ebird_data"] <- "Only nonbreeding abundance"


    next
  }


  if(!resident & breed_qual > 0) {
    down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                                 download_ranges = FALSE,
                                                 download_abundance = TRUE,
                                                 download_occurrence = FALSE,
                                                 force = FALSE,
                                                 pattern = paste0("proportion-population")),
                silent = TRUE)
  }

  if(resident & resident_qual > 0) {

    down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                                 download_ranges = FALSE,
                                                 download_abundance = TRUE,
                                                 download_occurrence = FALSE,
                                                 force = FALSE,
                                                 pattern = paste0("proportion-population")),
                silent = TRUE)
  }



# overlay maps ------------------------------------------------------------
  metric_used <- "mean" #options are "max", mean" or "median" (although "median" won't work for seasonal)


  if(resident){
    abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                   resolution = "3km",
                                                   period = "seasonal",
                                                   product = "proportion-population")  #3km high resolution

    rel_abundance <- abd_seasonal_abundance[["resident"]]

  }else{
  abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                 resolution = "3km",
                                                 period = "seasonal",
                                                 product = "proportion-population")  #3km high resolution

  rel_abundance <- abd_seasonal_abundance[["breeding"]]

  }


  region_boundary_proj_2 <- countries |>
    st_buffer(dist = 50000) |> #add 50km buffer for clipping eBird data
    st_transform(st_crs(rel_abundance))

  region_boundary_proj <- countries |>
    st_transform(st_crs(rel_abundance))

  rel_abundance <- crop(rel_abundance, region_boundary_proj_2) |>
    mask(region_boundary_proj_2)


  abundance_in_buffers <- terra::extract(rel_abundance,
                                         region_boundary_proj,
                                         fun = sum,
                                         na.rm = TRUE,
                                         ID = FALSE,
                                         exact = FALSE)

  abundance_df <- data.frame(country_code = region_boundary_proj$country_code,
                             north = region_boundary_proj$north,
                             ebird_abund = abundance_in_buffers[[1]]) %>%
    filter(
      #ebird_abund > 0,
      #!is.na(ebird_abund)
    ) %>%
    distinct() |>
    mutate(species_ebird = species_ebird,
           species = sp)

  abundance_all <- bind_rows(abundance_all,abundance_df)

  bam[i,"p_population_canada"] <- abundance_in_buffers[[1]][1]
  bam[i,"p_population_AK"] <- abundance_in_buffers[[1]][2]
  bam[i,"p_population_lower_48"] <- abundance_in_buffers[[1]][3]

print(round(i/nrow(bam),2))

} # end bam species loop
saveRDS(bam,"temp_bam.rds")
}else{
 bam <- readRDS("temp_bam.rds")

}




bam <- bam |>
  rename(cn = english,
         canadian_2020_population_estimate_bam_M_birds = population_estimate,
         population_estimate_lower = population_lower,
         population_estimate_upper = population_upper) |>
  select(cn,canadian_2020_population_estimate_bam_M_birds,
         population_estimate_lower,
         population_estimate_upper,
         p_population_canada,
         p_population_AK,
         p_population_lower_48)



decisions <- ExpAdjs |>
  full_join(bam, by = "cn")


# proportion of BOTW range that is NA - eBird -----------------------------


country <- rnaturalearth::ne_countries(continent = "North America") %>%
  filter(admin %in% c("Canada","United States of America", "Mexico")) %>%
  rename(strata_name = admin) %>%
  select(strata_name)

# botw_all <- readRDS("data/BOTW_valid.rds")
#
# botw_seas <- botw_all %>%
#   filter(seasonal %in% c(1,2),
#          presence %in% c(1:3)) |>
#   st_make_valid()
#
# saveRDS(botw_seas,"data/botw_all_seasonal.rds")
#

botw_seas <- readRDS("data/botw_all_seasonal.rds")


# bcr_strata <- st_read("data/bcr2026_statprov.shp")  |>
#   rmapshaper::ms_simplify(keep = 0.05,
#                           keep_shapes = TRUE) |>
#   sf::st_make_valid()
#
# saveRDS(bcr_strata,
#         "data/bcr2026_statprov_simplified.rds")


bcr_strata <- readRDS("data/bcr2026_statprov_simplified.rds")


for(i in 5:nrow(decisions)){

  sp <- as.character(decisions[i,"cn"])
  species_ebird <- as.character(decisions[i,"ebird_species"])
  sp_sci <- as.character(decisions[i,"scientific_name_ebirdst"])


  botw_range_sp <- try(botw_seas |>
    filter(sci_name == sp_sci) |>
    st_make_valid() |>
    st_intersection(country) |>
    st_make_valid(),silent = TRUE)

    if(class(botw_range_sp)[1] == "try-error"){next}

  if(nrow(botw_range_sp) == 0){next}

  qual_sel <- qual_ebird |>
    filter(species_code == species_ebird)

  breed_qual <- unname(unlist(qual_sel[1,"breeding_quality"]))
  resident_qual <- unname(unlist(qual_sel[1,"resident_quality"]))
  resident <- unname(unlist(qual_sel[1,"is_resident"]))

if(resident){
  season <- "resident"
}else{
  season <- "breeding"
}
  if(!((!resident & breed_qual > 0)|(resident & resident_qual > 0))){next}

  if(!resident & breed_qual > 0) {
    down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                                 download_ranges = FALSE,
                                                 download_abundance = TRUE,
                                                 download_occurrence = FALSE,
                                                 force = FALSE,
                                                 pattern = paste0("seasonal_mean")),
                silent = TRUE)
  }

  if(resident & resident_qual > 0) {

    down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                                 download_ranges = FALSE,
                                                 download_abundance = TRUE,
                                                 download_occurrence = FALSE,
                                                 force = FALSE,
                                                 pattern = paste0("seasonal_mean")),
                silent = TRUE)
  }


if(class(down)[1] == "try-error"){
  rm(down)
  next
  }


  if(resident){
    abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                   resolution = "27km",
                                                   period = "seasonal",
                                                   product = "abundance")  #3km high resolution

    rel_abundance <- abd_seasonal_abundance[["resident"]]

  }else{
    abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                   resolution = "27km",
                                                   period = "seasonal",
                                                   product = "abundance")  #3km high resolution

    rel_abundance <- abd_seasonal_abundance[["breeding"]]

  }

  # if(resident){
  #   abd_seasonal_abundance <- ebirdst::load_ranges(species = species_ebird,
  #                                                  resolution = "9km")  #highest resolution
  #
  #   rel_abundance <- abd_seasonal_abundance |>
  #     filter(season == "resident")
  #
  # }else{
  #   abd_seasonal_abundance <- ebirdst::load_ranges(species = species_ebird,
  #                                                  resolution = "9km") #highest resolution
  #
  #   rel_abundance <- abd_seasonal_abundance|>
  #     filter(season == "breeding")
  #
  # }

  bcr_strata_proj <- bcr_strata |>
    st_transform(st_crs(rel_abundance))


  region_boundary_proj <- try(botw_range_sp |>
    st_transform(st_crs(rel_abundance)) |>
    summarise(), silent = TRUE)

  if(class(region_boundary_proj)[1] == "try-error"){

    region_boundary_proj <- try(botw_range_sp |>
                                  st_transform(st_crs(rel_abundance)) |>
                                  st_make_valid() |>
                                  summarise(), silent = TRUE)
  }



  region_boundary_plot <- botw_range_sp |>
    st_transform(st_crs(bcr_strata))

  rel_abundance_crop <- rel_abundance |>
    crop(bcr_strata_proj) |>
    mask(bcr_strata_proj)



  bb <- st_bbox(region_boundary_plot)

  n_na_in_buffers <- terra::extract(rel_abundance_crop,
                                    region_boundary_proj,
                                    fun = "isNA",
                                    ID = FALSE,
                                    exact = FALSE)

  n_not_na_in_buffers <- terra::extract(rel_abundance_crop,
                                        region_boundary_proj,
                                        fun = "notNA",
                                        ID = FALSE,
                                        exact = FALSE)
  p_na_in_range <- as.numeric(n_na_in_buffers/n_not_na_in_buffers)

  tst <- ggplot()+
    geom_sf(data = bcr_strata,
                          fill = NA)+
    tidyterra::geom_spatraster(data = rel_abundance_crop,maxcell = 2000000)+
    geom_sf(data = bcr_strata,
            fill = NA)+
    geom_sf(data = region_boundary_plot,
            colour = "red",
            fill = "red",
            alpha = 0.1)+
    coord_sf(xlim = bb[c("xmin","xmax")],
             ylim = bb[c("ymin","ymax")])+
    scale_fill_viridis_c(name = "eBird seasonal\nrelative abundance")+
    labs(title = paste(sp,"BOTW range map and eBirdst",season,"relative abundance"),
         subtitle = paste0("NA regions for eBird relative abundance show as grey \n and account for ",round(p_na_in_range*100),"% of the BOTW range"))

  pdf(paste0("maps/",species_ebird,"_range_overlap.pdf"),
      width = 11,
      height = 8.5)
  print(tst)
  dev.off()



  decisions[i,"p_BOTW_range_na_eBirdst"] <- p_na_in_range


print(round(i/nrow(decisions),2))


  }# end species decisions rows loop






acad_pop_size <- readRDS("data/cleaned_acad_global_2024_05_23.rds") |>
  select(common_name,scientific_name,
         sort,group,order,family,
         starts_with("pop_size_us_ca")) |>
  filter(!is.na(pop_size_us_ca)) |>
  mutate(common_name = ifelse(common_name == "Common Redpoll",
                              "Redpoll",common_name),
         common_name = ifelse(common_name == "Cordilleran Flycatcher",
                              "Western Flycatcher",common_name),
         common_name = ifelse(common_name == "House Wren",
                              "Northern House Wren",common_name),
         common_name = ifelse(common_name == "Northern Goshawk",
                              "American Goshawk",common_name),
         common_name = ifelse(common_name == "Black-crowned Night-Heron",
                              "Black-crowned Night Heron",common_name),
         common_name = ifelse(common_name == "Yellow-crowned Night-Heron",
                              "Yellow-crowned Night Heron",common_name),
         common_name = ifelse(common_name == "Cattle Egret",
                              "Western Cattle-Egret",common_name),
         common_name = ifelse(common_name == "Barn Owl",
                              "American Barn Owl",common_name))

decisions2 <- decisions |>
  left_join(acad_pop_size,
            by = c("cn" = "common_name"))

bbs_cov <- read_csv("data/BBS_Trends_survey-wide_country_prov_state_2024.csv") |>
  filter(trend_time == "Three-generation",
         region_en == "Survey-wide") |>
  select(species,
         BBS_Number__Numéro_BBS_core,
         reliab.cov,
         rel_ab_obs,
         n_site,
         reliability) |>
  rename(bbs_number = BBS_Number__Numéro_BBS_core,
         bbs_coverage = reliab.cov,
         bbs_mean_count = rel_ab_obs,
         bbs_n_routes = n_site,
         bbs_cws_reliability = reliability) |>
  mutate(species = gsub(" (all forms)",
                        "",
                        x = species,
                        fixed = TRUE),
         species = gsub(" (Common/Hoary)",
                        "",
                        species,
                        fixed = TRUE),
         species = gsub(" (Cordilleran/Pacific-slope)",
                        "",
                        x = species,
                        fixed = TRUE))



decisions3 <- decisions2 |>
  left_join(bbs_cov,
            by = c("cn" = "species"))



#saveRDS(decisions3,"full_decisions_file.rds")


decisions3 <- readRDS("full_decisions_file.rds")

decisions4 <- decisions3 |>
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
  select(common_name,
         taxonomic_sort,
         group,
         order,
         family,
          common_name_ebirdst,
          scientific_name,
          scientific_name_ebirdst,
          ebird_species_code,
          ebirdst_breeding_quality,
          ebirdst_resident_quality,
          prop_BOTW_range_na_ebirdst,
          BAM_population_estimate_Canada,
          p_population_canada_ebirdst,
          ACAD_pop_size_us_ca_source,
          ACAD_pop_size_us_ca_comment,
          ACAD_pop_size_us_ca,
          starts_with("bbs"),
          starts_with("NA_POPS_avail"),
          starts_with("PIF_time_of_day"),
          starts_with("NA_POPS_edr"),
          starts_with("PIF_MDD"),
          PIF_pair,
         PIF_CALC_Rosenberg_2019) |>
  arrange(taxonomic_sort)



write_csv(decisions4,
          "Population_estimates_workflow_support.csv")

# eBird data coverage ----------------------------------------------------------

#
# #ebirdst_download_data_coverage()
#
# countries <- bbsBayes2::load_map("prov_state") |>
#   mutate(north = ifelse(country_code == "CA" | prov_state == "AK",
#                         "CAN_AK","Lower_48")) |>
#    ungroup() |>
#    group_by(country_code,north) |>
#     summarise(.groups = "drop") |>
#     st_buffer(dist = 50000)
#
#
# tmp <- load_data_coverage("selection-probability",weeks = "05-10")
#
# tmp <- load_data_coverage("spatial-coverage",weeks = "06-21")
#
#  region_boundary_proj <- st_transform(countries, st_crs(tmp))
#
#
#
# tmp <- crop(tmp, region_boundary_proj) |>
#  mask(region_boundary_proj)
#
#
# bbs_strata <- bbsBayes2::load_map("bbs") |>
#   group_by(strata_name) |>
#   summarise(.groups = "drop") |>
#   st_transform(st_crs(tmp))
#
#
# plot(tmp)
#
# coverage_in_strata <- terra::extract(tmp,
#                                      bbs_strata,
#                                        fun = mean,
#                                        na.rm = TRUE,
#                                        ID = FALSE,
#                                        exact = FALSE)
#
