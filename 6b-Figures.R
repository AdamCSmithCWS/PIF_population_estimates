## sandbox results explore
#
#

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


yr_ebird <- 2023 # prediction year for eBird relative abundance


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


sps_list <- readRDS("data/all_bbs_data_species_list.rds") #

# Load BBS data -----------------------------------------------------------


# mean counts to id routes on which species has been observed -------------
mean_counts_all <- readRDS("data/mean_counts_by_route.rds")

raw_counts_all <- readRDS("data/all_counts_by_route.rds")

routes_buf_all <- readRDS("data/all_routes_buffered.rds")

strata <- bbsBayes2::load_map("bbs_usgs")

strata_names <- strata %>%
  sf::st_drop_geometry()



#
#
# # load traditional estimates  fir later comparison -------------------------
# trad_dir <- "Stanton_2019_code/output"
# strata_join <- strata_names %>%
#   select(strata_name,country_code,
#          prov_state, bcr) %>%
#   st_drop_geometry() %>%
#   distinct()
#
# strata_trad_all <- read_csv(paste0(trad_dir,"/PSest_Prov_by_BCR__100iter.csv")) %>%
#   left_join(.,strata_join,
#             by = c("st_abrev" = "prov_state",
#                    "BCR" = "bcr")) %>%
#   mutate(region = strata_name,
#          region_type = "strata",
#          version = "traditional",
#          species = cn,
#          aou = Aou,
#          pop_median = med.PopEst,
#          pop_lci_80 = LCI80.PopEst,
#          pop_uci_80 = UCI80.PopEst,
#          pop_lci_95 = LCI95.PopEst,
#          pop_uci_95 = UCI95.PopEst)
#
# USACAN_trad_all <- read_csv(paste0(trad_dir,"/PSest_Global_WH_NA__100iter.csv")) %>%
#   mutate(region = "USACAN",
#          region_type = "USACAN",
#          version = "traditional",
#          species = cn,
#          aou = Aou,
#          pop_median = USCAN.med.PopEst,
#          pop_lci_80 = USCAN.80LCI.PopEst,
#          pop_uci_80 = USCAN.80UCI.PopEst,
#          pop_lci_95 = USCAN.95LCI.PopEst,
#          pop_uci_95 = USCAN.95UCI.PopEst)
#
#
#
# global_trad_all <- read_csv(paste0(trad_dir,"/PSest_Global_WH_NA__100iter.csv")) %>%
#   mutate(region = "global",
#          region_type = "global",
#          version = "traditional",
#          species = cn,
#          aou = Aou,
#          pop_median = GL.med.PopEst,
#          pop_lci_80 = GL.80LCI.PopEst,
#          pop_uci_80 = GL.80UCI.PopEst,
#          pop_lci_95 = GL.95LCI.PopEst,
#          pop_uci_95 = GL.95UCI.PopEst)
#
#
#
# prov_trad_all <- read_csv(paste0(trad_dir,"/PSest_by_Prov__100iter.csv")) %>%
#   mutate(region = st_abrev,
#          region_type = "prov_state",
#          version = "traditional",
#          species = cn,
#          aou = Aou,
#          pop_median = med.PopEst,
#          pop_lci_80 = LCI80.PopEst,
#          pop_uci_80 = UCI80.PopEst,
#          pop_lci_95 = LCI95.PopEst,
#          pop_uci_95 = UCI95.PopEst)
#
# bcr_trad_all <- read_csv(paste0(trad_dir,"/PSest_by_BCR__100iter.csv")) %>%
#   mutate(region = as.character(BCR),
#          region_type = "bcr",
#          version = "traditional",
#          species = cn,
#          aou = Aou,
#          pop_median = med.PopEst,
#          pop_lci_80 = LCI80.PopEst,
#          pop_uci_80 = UCI80.PopEst,
#          pop_lci_95 = LCI95.PopEst,
#          pop_uci_95 = UCI95.PopEst)
#
#
#
#
# pop_ests_out_trad <- bind_rows(USACAN_trad_all,
#                                global_trad_all,
#                                #country_abund,
#                                bcr_trad_all,
#                                strata_trad_all) %>%
#   select(region,region_type,version,
#          species,aou,
#          pop_median, pop_lci_80, pop_uci_80,
#          pop_lci_95, pop_uci_95)
#
# saveRDS(pop_ests_out_trad,"all_traditional_pop_estimates.rds")
#

pop_ests_out_trad <- readRDS("all_traditional_pop_estimates.rds")



# Example figure for roadside bias ----------------------------------------------------------

demo_sp <- data.frame(species = c("American Robin",
                                 # "American Robin",
                                  "Canyon Wren"),
                      region = c("CA-BC-5",
                                 #"CA-YT-4",
                                 "US-AZ-16"))




abund_maps_save <- vector(mode = "list",length = nrow(demo_sp))
names(abund_maps_save) <- paste0(demo_sp$species,"_",demo_sp$region)

bb_sav <- vector(mode = "list",length = nrow(demo_sp))
for(i in 1:nrow(demo_sp)){
sp_sel <- demo_sp[i,"species"]
reg_sel <- demo_sp[i,"region"]


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
sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]






# filtering routes on which species wasn't observed.
routes_buf <- routes_buf_all

sp_ebird <- ebirdst::get_species(sp_sel)

breed_abundance <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",
                                  sp_ebird,
                                  "_derived_breeding_relative_abundance.rds"))

strata_sel <- strata %>%
  filter(strata_name == reg_sel)

strata_sel1 <- strata_sel %>% st_buffer(5e5) # 50km buffer

routes_sel <- routes_buf %>%
  sf::st_transform(crs = st_crs(strata_sel)) %>%
  sf::st_join(y = strata_sel,left = FALSE)

bb <- sf::st_bbox(strata_sel)

bb_sav[[i]] <- bb

strata_sel2 <- strata_sel1 %>%
  sf::st_transform(crs = st_crs(breed_abundance))

bb_crop <- terra::ext(strata_sel2)

breed_abundance_plot <- terra::crop(breed_abundance,y = bb_crop) %>%
  terra::mask(., mask = strata_sel2) #%>%
  #terra::project(y = "epsg:9822")

names(breed_abundance_plot) <- "breeding"

# breed_abundance_plot <- breed_abundance_plot %>%
#   sf::st_transform(crs = sf::st_transform(strata_sel))
if(sp_sel == "Canyon Wren"){
  coord_labs <- list(bottom = "E",right = "N")
}else{
  coord_labs <- list(bottom = "E",left = "N")
}
abund_map <- ggplot()+
  geom_sf(data = strata_sel)+
  geom_spatraster(data = breed_abundance_plot,maxcell = 16000000)+
  geom_sf(data = strata_sel, fill = NA, colour = grey(0.4),
          linewidth = 0.8)+
  # geom_sf_text(data = strata,aes(label = strata_name),
  #              size = 1)+
  geom_sf(data = routes_sel,fill = "black",colour = "black",
          linewidth = 0.9)+
  coord_sf(xlim = c(bb[c("xmin","xmax")]),
           ylim = c(bb[c("ymin","ymax")]),
           #label_graticule = "SE",
           label_axes = coord_labs)+
  #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA,
                       begin = 0.2,
                       end = 1,
                       name = paste(sp_sel,"\nRelative Abundance"))+
  # scale_colour_viridis_c(direction = -1,
  #                        option = "rocket",
  #                        na.value = NA,
  #                        end = 0.9,
  #                        begin = 0.1,
  #                        name = "Mean count BBS")+
  theme_bw()+
  labs(title = paste(sp_sel))+
  theme(legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1),
                           units = "mm"))

abund_map
#
saveRDS(abund_map,paste0("figures/saved_ggplots/sampling_map_",sp_aou,"_",reg_sel,".rds"))

png(filename = paste0("final_figures/Figure_1_",sp_aou,"_",reg_sel,".png"),
    res = 400,
    height = 7,
    width = 6.5,
    units = "in")
print(abund_map)
dev.off()

abund_maps_save[[paste0(sp_sel,"_",reg_sel)]] <- abund_map

}


des<- c(patchwork::area(t = 1,l = 1,b = 6,r = 3),
        patchwork::area(t = 1,l = 3,b = 4.5,r = 5),
        patchwork::area(t = 5,l = 4,b = 6,r = 5))
# abund_maps_save[[1]] <- abund_maps_save[[3]]
# abund_maps_save[[2]] <- abund_maps_save[[4]]

# abund_maps_save[[2]] <- abund_maps_save[[2]]


sampling_demo <- (abund_maps_save[[1]]) + (abund_maps_save[[2]] +
                                                     theme(plot.title = element_text(hjust = 1))) +
  guide_area() + plot_layout(design = des,guides = "collect")

sampling_demo

pdf("final_figures/sampling_bias_maps.pdf",
    width = 7,
    height = 6)
print(sampling_demo)
dev.off()


# BBS route map with example insets ---------------------------------------
strata <- bbsBayes2::load_map("bbs_usgs")

bb <- st_bbox(strata)
country_codes <- readxl::read_xlsx("data/iso_codes.xlsx") %>%
  rename(country_name = Name) #

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


#Bird Studies Canada and NABCI.  2014.  Bird Conservation Regions.  Published by
#Bird Studies Canada on behalf of the North American Bird Conservation Initiative.
# https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions  Accessed:  Nov 8 2024
bcrs <- read_sf("data/BCR_terrestrial/BCR_Terrestrial_master_International.shp") %>%
  group_by(BCR,Label) %>%
  summarise(do_union = FALSE)


bcrs_plot <- bcrs %>%
  st_transform(crs = st_crs(strata))

countries_plot <- countries %>%
  st_transform(crs = st_crs(strata)) %>%
  filter(CC == "NA")

bx1 <- st_as_sfc(bb_sav[[1]])
bx2 <- st_as_sfc(bb_sav[[2]])

rts_map <- ggplot() +
  #geom_sf(data = countries_plot, fill = NA)+
  geom_sf(data = bcrs_plot, fill = NA,
          colour = grey(0.75))+
  geom_sf(data = bx1,fill = NA,linewidth = 0.6)+
  geom_sf(data = bx2,fill = NA,linewidth = 0.6)+
  geom_sf(data = routes_buf_all,
          colour = viridisLite::mako(1,begin = 0.3, end = 0.4))+
  coord_sf(xlim = c(bb[c("xmin","xmax")])*c(0.9,1),
           ylim = c(bb[c("ymin","ymax")])*c(1.5,1))+
  theme_bw()


rts_map


#
# des<- c(patchwork::area(t = 1,l = 1,b = 6,r = 3),
#         patchwork::area(t = 1,l = 3,b = 4.5,r = 5),
#         patchwork::area(t = 5,l = 4,b = 6,r = 5))
# # abund_maps_save[[1]] <- abund_maps_save[[3]]
# # abund_maps_save[[2]] <- abund_maps_save[[4]]
#
# # abund_maps_save[[2]] <- abund_maps_save[[2]]
#
#
# sampling_demo <- (abund_maps_save[[1]]) + (abund_maps_save[[2]] +
#                                              theme(plot.title = element_text(hjust = 1))) +
#   guide_area() + plot_layout(design = des,guides = "collect")
#
# sampling_demo

pdf("final_figures/BBS_route_locations.pdf",
    width = 4,
    height = 4)
print(rts_map)
dev.off()






# continental estimate summaries ------------------------------------------
sp_example <- readRDS("data/selected_species.rds")
# sp_example <- unique(c( "Western Meadowlark","Baird's Sparrow","Brown Creeper",
#                  "Canyon Wren",
#                  "Black-capped Chickadee","American Robin",
#                  "Barn Swallow",
#                  "Blackpoll Warbler",
#                  "Mountain Bluebird",
#                  "Eastern Phoebe",
#                  "Rose-breasted Grosbeak",
#                  "Downy Woodpecker",
#                  "Red-winged Blackbird",
#                  "Scarlet Tanager",
#                  "Say's Phoebe",
#                  "Black-chinned Hummingbird",
#                  "Wilson's Snipe","Long-billed Curlew","Killdeer","Tennessee Warbler","Bay-breasted Warbler",
#                  "Swainson's Thrush","Blue Jay",
#                  "Blue-headed Vireo","Bicknell's Thrush",
#                  "American Kestrel",
#                  "Steller's Jay","Clay-colored Sparrow",
#                  "Yellow-rumped Warbler","Olive-sided Flycatcher",
#                  "Purple Martin",
#                  "Barn Swallow",
#                  "Verdin","Ash-throated Flycatcher","Black-throated Sparrow",
#                  "Blue Jay","Varied Thrush","Veery","Wood Thrush","Chestnut-collared Longspur",
#                  "Bobolink","Savannah Sparrow","Grasshopper Sparrow",
#                  "Horned Lark",
#                  "Eastern Whip-poor-will", "Common Nighthawk", "Scissor-tailed Flycatcher"))
#


ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")

re_summarise <- FALSE # set to true to re-run this summary loop
if(re_summarise){
pop_compare_stack_sel <-NULL

for(sp_sel in unique(sp_example)){ #sps_list$english){


    sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
    sp_ebird <- ebirdst::get_species(sp_sel)

if(!file.exists(paste0("estimates/pop_ests_alt_",sp_aou,sp_ebird,".csv"))){next}


pop_ests_out <- read_csv(paste0("estimates/pop_ests_alt_",
                                sp_aou,sp_ebird,".csv"))%>%
  mutate(version = "PIF_eBird_with_EDR_Avail")


pop_ests_out_trad_sel <- pop_ests_out_trad %>%
  filter(species == sp_sel) %>%
  mutate(version = "PIF_traditional")


pop_compare_stack <- pop_ests_out %>%
  #bind_rows(pop_ests_out_noEDR) %>%
  select(region,region_type,version,
         species,species_ebird,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95) %>%
  bind_rows(pop_ests_out_trad_sel)


pop_compare_stack_sel <- bind_rows(pop_compare_stack_sel,pop_compare_stack)
} # end species loop

sp_names_e <- pop_compare_stack_sel %>%
  select(species,species_ebird) %>%
  distinct() %>%
  filter(!is.na(species_ebird))

sp_names_a <- pop_compare_stack_sel %>%
  select(species,aou) %>%
  distinct() %>%
  filter(!is.na(aou))


pop_compare_stack_sel <- pop_compare_stack_sel %>%
  select(-c(aou,species_ebird)) %>%
  inner_join(sp_names_e,by = "species") %>%
  inner_join(sp_names_a,by = "species")

saveRDS(pop_compare_stack_sel,
        "estimates/comparison_estimates_method.rds")
}else{
  pop_compare_stack_sel <- readRDS("estimates/comparison_estimates_method.rds")

}

# join with habitat categories
#

acad <- readRDS("data/cleaned_acad_global_2024_05_23.rds") # this file created in
# https://github.com/AdamCSmithCWS/ACAD_summary/blob/25c1693530bc5243e5e615572ee7f67faae5d4be/plotting.R#L5C1-L84C1
#

breed_habitats <- acad %>%
  select(primary_breeding_habitat_major,primary_breeding_habitat_sub, common_name, canada, usa, mexico, c_america) %>%
  distinct()

# Compare correction factors ----------------------------------------------

Cor_factors <- ExpAdjs %>%
  select(cn, Dist2,Pair2, TimeAdj.meanlog,
         edr, availability, Aou) %>%
  mutate(Cd_trad = 1/(pi*(Dist2^2)),
         Cd_edr = 1/(pi*(edr^2)),
         Ct_trad = exp(TimeAdj.meanlog),
         Ct_avail = 1/availability,
         Cd_diff = Cd_edr/Cd_trad,
         Ct_diff = Ct_avail/Ct_trad) %>%
  left_join(breed_habitats,
            by = c("cn" = "common_name"))





pop_compare_wide <- pop_compare_stack_sel %>%
  filter(!is.na(region)) %>%
  pivot_wider(names_from = version,
              names_sep = "_",
              values_from = c(pop_median,pop_lci_80,pop_lci_95,
                              pop_uci_80,pop_uci_95)) %>%
  mutate(pop_median_PIF_traditional = ifelse(is.na(pop_median_PIF_traditional),
                                             0,
                                             pop_median_PIF_traditional),
         pop_lci_80_PIF_traditional = ifelse(is.na(pop_lci_80_PIF_traditional),
                                             0,
                                             pop_lci_80_PIF_traditional),
         pop_uci_80_PIF_traditional = ifelse(is.na(pop_uci_80_PIF_traditional),
                                             0,
                                             pop_uci_80_PIF_traditional),
         pop_lci_95_PIF_traditional = ifelse(is.na(pop_lci_95_PIF_traditional),
                                             0,
                                             pop_lci_95_PIF_traditional),
         pop_uci_95_PIF_traditional = ifelse(is.na(pop_uci_95_PIF_traditional),
                                             0,
                                             pop_uci_95_PIF_traditional),
         dif_mag_new_trad = pop_median_PIF_eBird_with_EDR_Avail/pop_median_PIF_traditional)

pop_compare_wide_usacan <- pop_compare_wide %>%
  filter(region %in% c("USACAN")) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  mutate(species_factor = fct_reorder(factor(species),dif_mag_new_trad,
                                      .na_rm = TRUE))


pop_compare_stack <- pop_compare_stack_sel %>%
  filter(region %in% c("USACAN")) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  mutate(species_factor = factor(factor(species),levels = levels(pop_compare_wide_usacan$species_factor)))

side_plot1 <- ggplot(data = pop_compare_stack,
                    aes(y = species_factor,
                        x = pop_median,
                        colour = version))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbarh(aes(xmin = pop_lci_80,
                     xmax = pop_uci_80),
                 height = 0,
                 position = position_dodge(width = 0.5))+
  scale_colour_viridis_d(end = 0.8)+
  ylab("Continental estimates")+
  xlab("Population estimate with 80% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  theme_bw()

side_plot1


# summarise sampling bias -------------------------------------------------

w_mean_lg_rat <- function(x,w){
  wx <- rep(NA,length(x))
  sum_w <- sum(w,na.rm = TRUE)

    for(i in 1:length(x)){
    if(is.finite(x[i])){
      wx[i] <- (x[i]*w[i])/sum_w
    }else{
      wx[i] <- (log(0.0001)*w[i])/sum_w
    }
  }
  return(sum(wx,na.rm = TRUE))
}

sampl_bias <- NULL

for(sp_sel in sp_example){

  species_ebird <- ebirdst::get_species(sp_sel)

  if(!file.exists(paste0("data/species_relative_abundance_",yr_ebird,"/",species_ebird,"_bbs_strata_relative_abundance.rds"))){next}

  route_sampled <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",species_ebird,"_relative_abundance.rds"))

  mean_sampled_abund <- mean(route_sampled$ebird_abund,na.rm = TRUE)
  #


  sampl_bias1 <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",species_ebird,"_bbs_strata_relative_abundance.rds")) %>%
    mutate(species = sp_sel,
           sp_ebird = species_ebird,
           sampled_abundance = ifelse(!is.na(mean_ebird_abundance) & is.na(sampled_abundance),
                                      0,sampled_abundance),
           sampled_abundance = ifelse(strata_name == "USA_CAN",mean_sampled_abund,sampled_abundance),
           sampled_avail_ratio = sampled_abundance/mean_ebird_abundance,
           log_ratio = log(sampled_avail_ratio))



  tmp <- sampl_bias1 %>%
    filter(strata_name != "USA_CAN",
           !is.na(log_ratio)) %>%
    mutate(log_ratio = ifelse(is.finite(log_ratio),log_ratio,log(0.001)),
           sampled_abundance = ifelse(is.na(sampled_abundance),0,sampled_abundance))
  #
  wm_samp <- w_mean_lg_rat(tmp$log_ratio,tmp$total_ebird_abundance)
  #
  sampl_bias1 <- sampl_bias1 %>%
    mutate(log_ratio = ifelse(strata_name == "USA_CAN",wm_samp,log_ratio),
           sampled_avail_ratio = ifelse(strata_name == "USA_CAN",exp(log_ratio),sampled_avail_ratio))


  sampl_bias <- bind_rows(sampl_bias,sampl_bias1)


}

sampl_bias <- sampl_bias %>%
  mutate(ratio_new_old = 1/exp(log_ratio))

saveRDS(sampl_bias,"final_figures/example_species_strata_sampling_bias.rds")

canw <- sampl_bias %>%
  filter(species == "Canyon Wren",
         !is.na(mean_ebird_abundance),
         mean_ebird_abundance > 0)


USA_CAN_sample_bias <- sampl_bias %>%
  filter(strata_name == "USA_CAN")




canw <- sampl_bias %>%
  filter(species == "Canyon Wren",
         !is.na(mean_ebird_abundance),
         mean_ebird_abundance > 0)


USA_CAN_sample_bias <- sampl_bias %>%
  filter(strata_name == "USA_CAN")




# calculate the ratio of existing and ebird estiates with edr -------------


output_dir <- "G:/PIF_population_estimates/output"
output_dir <- "output"


all_pop_ests <- NULL

vers <- ""
for(sp_sel in sp_example){

  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)

  if(file.exists(paste0("estimates/pop_ests_alt_",vers,sp_aou,sp_ebird,".csv"))){
    pop_ests <- read_csv(paste0("estimates/pop_ests_alt_",vers,sp_aou,sp_ebird,".csv"))


    all_pop_ests <- bind_rows(all_pop_ests,pop_ests)
  }

}



us_can <- all_pop_ests %>%
  filter(region == "USACAN") %>%
  select(species,pop_median,pop_lci_80:pop_uci_95,
         version)

usa_can_exist <- read_csv("Stanton_2019_code/output/PSest_Global_WH_NA__100iter.csv") %>%
  filter(cn %in% us_can$species) %>%
  select(cn,USCAN.med.PopEst:USCAN.95UCI.PopEst) %>%
  mutate(version = "Existing")

names(usa_can_exist) <- names(us_can)

us_can_long <- bind_rows(us_can,usa_can_exist)

names(usa_can_exist) <- paste0(names(usa_can_exist),"_exist")

us_can_wide <- us_can %>%
  left_join(usa_can_exist,
            by = c("species" = "species_exist")) %>%
  mutate(ratio = pop_median/pop_median_exist,
         log_ratio = log(ratio),
         pop_new_millions = pop_median/1e6) %>%
  relocate(species,ratio,log_ratio,pop_new_millions)


write_csv(us_can_wide,"final_figures/usa_canada_comparison_eBird_existing_w_EDR.csv")




# calculate the difference between mean abundance and abundance in eBird year -------------


re_sum_trend <- TRUE

if(re_sum_trend){

trend_effect_out <- NULL
trend_effect_st_out <- NULL

vers <- ""
for(sp_sel in sp_example){

  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)
if(sp_sel == "Ruby-throated Hummingbird"){next}
  if(file.exists(paste0(output_dir,"/calibration_fit_alt_",
                        vers,sp_aou,"_",sp_ebird,".rds"))){

    fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",
                          vers,sp_aou,"_",sp_ebird,".rds"))

    raw_dat <- readRDS(paste0("stan_data/dataframe_",vers,sp_aou,"_",sp_ebird,".rds"))


strata_population_posterior <- readRDS(paste0("output/strata_population_posterior_",vers,sp_aou,"_",sp_ebird,".rds"))

# mean abundance across all years
    mean_pops <- strata_population_posterior %>%
      group_by(year,.draw) %>%
      summarise(pop = sum(population),.groups = "drop") %>%
      group_by(.draw) %>%
      summarise(mean_pop = mean(pop), .groups = "drop") #%>%
# abundance in ebird year
    mean_pops22 <- strata_population_posterior %>%
      filter(year == yr_ebird) %>%
      group_by(.draw) %>%
      summarise(pop = sum(population),.groups = "drop")
# ratio of abundance in ebird year with mean abundance all years
    trend_effect <- inner_join(mean_pops22,mean_pops,
                               by = ".draw") %>%
      mutate(p_22 = pop/mean_pop) %>%
      summarise(median_proportion = median(p_22),
                mean_proportion = mean(p_22),
                lci_proportion = quantile(p_22,0.025),
                uci_proportion = quantile(p_22,0.975)) %>%
      mutate(species = sp_sel,
             sp_eBird = sp_ebird)

# mean abundance in all years by strata
    mean_pops_st <- strata_population_posterior %>%
      group_by(year,strata_name,.draw) %>%
      summarise(pop = sum(population),.groups = "drop") %>%
      filter(pop > 0) %>%
      group_by(strata_name,.draw) %>%
      summarise(mean_pop = mean(pop), .groups = "drop") #%>%

# abundance in ebird year by strata
    mean_pops22_st <- strata_population_posterior %>%
      filter(year == yr_ebird) %>%
      group_by(strata_name,.draw) %>%
      summarise(pop = sum(population),.groups = "drop") %>%
      filter(pop > 0)

    # ratio of abundance in ebird year over mean abundance by strata
    trend_effect_st <- inner_join(mean_pops22_st,mean_pops_st,
                               by = c("strata_name",".draw")) %>%
      mutate(p_22 = pop/mean_pop) %>%
      group_by(strata_name) %>%
      summarise(median_proportion = median(p_22,na.rm = TRUE),
                mean_proportion = mean(p_22,na.rm = TRUE),
                lci_proportion = quantile(p_22,0.025,na.rm = TRUE),
                uci_proportion = quantile(p_22,0.975,na.rm = TRUE),
                mean_abund = mean(mean_pop),
                abund_22 = mean(pop)) %>%
      mutate(species = sp_sel,
             sp_eBird = sp_ebird)


       trend_effect_out <- bind_rows(trend_effect_out,trend_effect)

      trend_effect_st_out <- bind_rows(trend_effect_st_out,trend_effect_st)

  }

}

trend_effects <- trend_effect_out %>%
  mutate(log_ratio = log(mean_proportion),
         ratio = mean_proportion) %>%
  rename(sp_ebird = sp_eBird)


write_csv(trend_effects,"final_figures/trend_effect_summaries.csv")

trend_effects_st <- trend_effect_st_out %>%
  mutate(log_ratio = log(mean_proportion),
         ratio = mean_proportion) %>%
  rename(sp_ebird = sp_eBird)


write_csv(trend_effects_st,"final_figures/trend_effect_summaries_strata.csv")

}else{
  trend_effects <- read_csv("final_figures/trend_effect_summaries.csv")
  trend_effects_st <- read_csv("final_figures/trend_effect_summaries_strata.csv")

}

## explore the canyon wren sampling vs trend effects

canw_sample <- canw %>%
  mutate(log_geographic_ratio = ifelse(is.finite(log_ratio),log_ratio,
                                       log(0.001))) %>%
  filter(strata_name != "USA_CAN") %>%
  select(species,sp_ebird, strata_name,log_geographic_ratio,total_ebird_abundance)
trend_effects_st_canw <- trend_effects_st %>%
  filter(species == "Canyon Wren",
         strata_name != "USA_CAN") %>%
  mutate(log_trend_ratio = log_ratio) %>%
  select(species,sp_ebird,strata_name,log_trend_ratio) %>%
  full_join(canw_sample, by = c("strata_name","species","sp_ebird"))%>%
  mutate(log_trend_ratio = ifelse(is.na(log_trend_ratio),
                                  0,log_trend_ratio),
         w_trend = ifelse(log_trend_ratio == 0,FALSE,TRUE),
         w_sampling = ifelse(log_geographic_ratio == log(0.001),
                             FALSE,TRUE))

tst <- ggplot(data = filter(trend_effects_st_canw,log_trend_ratio != 0),
              aes(y = exp(-1*log_geographic_ratio),
                  x = exp(log_trend_ratio),
                  alpha = total_ebird_abundance))+
  geom_point()+
  geom_smooth()
tst

miss <- trend_effects_st_canw %>%
  group_by(w_trend,w_sampling) %>%
  summarise(pop = sum(total_ebird_abundance))

# combined new adjustment factors -----------------------------------------


# compare adjustment factors ----------------------------------------------

pop_compare_realised <- pop_compare_stack_sel %>%
  filter(region %in% c("USACAN"))


pop_compare_realised_wide <- pop_compare_realised %>%
  pivot_wider(names_from = version,
              names_sep = "_",
              values_from = c(pop_median,pop_lci_80,pop_lci_95,
                              pop_uci_80,pop_uci_95)) %>%
  mutate(pop_median_PIF_traditional = ifelse(is.na(pop_median_PIF_traditional),
                                             0,
                                             pop_median_PIF_traditional),
         pop_lci_80_PIF_traditional = ifelse(is.na(pop_lci_80_PIF_traditional),
                                             0,
                                             pop_lci_80_PIF_traditional),
         pop_uci_80_PIF_traditional = ifelse(is.na(pop_uci_80_PIF_traditional),
                                             0,
                                             pop_uci_80_PIF_traditional),
         pop_lci_95_PIF_traditional = ifelse(is.na(pop_lci_95_PIF_traditional),
                                             0,
                                             pop_lci_95_PIF_traditional),
         pop_uci_95_PIF_traditional = ifelse(is.na(pop_uci_95_PIF_traditional),
                                             0,
                                             pop_uci_95_PIF_traditional),
         dif_mag_new_trad = pop_median_PIF_eBird_with_EDR_Avail/pop_median_PIF_traditional) %>%
  select(species,pop_median_PIF_eBird_with_EDR_Avail,pop_median_PIF_traditional,dif_mag_new_trad) %>%
  rename(Ratio = dif_mag_new_trad,
         cn = species) %>%
  mutate(adjfactor = "Overall")

over_bias <- pop_compare_realised_wide %>%
  select(cn,Ratio,adjfactor)

sampl_bias <- readRDS("final_figures/example_species_strata_sampling_bias.rds")

USA_CAN_sample_bias <- sampl_bias %>%
  filter(strata_name == "USA_CAN") %>%
  select(species,ratio_new_old) %>%
  rename(Ratio = ratio_new_old,
         cn = species)%>%
  mutate(adjfactor = "Geographic Bias")


trend_effects <- read_csv("final_figures/trend_effect_summaries.csv") %>%
  select(species,ratio) %>%
  rename(Ratio = ratio,
         cn = species) %>%
  mutate(adjfactor = "Trend")

hab_tmp <- sampl_bias %>%
  filter(strata_name == "USA_CAN") %>%
  select(species,ratio_new_old) %>%
  rename(habitat_ratio = ratio_new_old,
         cn = species)
trend_tmp <- read_csv("final_figures/trend_effect_summaries.csv") %>%
  select(species,ratio) %>%
  rename(trend_ratio = ratio,
         cn = species)






# dropping one species with availability estimates < 0.01 - unrealistically low
ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv") %>%
  mutate(tod = exp(TimeAdj.meanlog),
         availability = ifelse(Species %in% c("CORE"),NA,availability),
         avail_tod_scaled = 1/availability,
         Tr = (1/(tod*availability)),
         Ar = (Dist2^2/edr^2),
         comb_r = Tr*Ar,
         Tlr = log(Tr),
         Alr = log(Ar),
         clr = log(comb_r))

adj_tmp <- ExpAdjs %>%
  select(cn,Tr,Ar,comb_r) %>%
  rename(availability_ratio = Tr,
         area_ratio = Ar,
         combined = comb_r)

tabl2_paper <- pop_compare_realised_wide %>%
  left_join(hab_tmp) %>%
  left_join(trend_tmp) %>%
  left_join(adj_tmp) %>%
  select(cn,area_ratio,availability_ratio,habitat_ratio,trend_ratio,
           Ratio,pop_median_PIF_eBird_with_EDR_Avail,pop_median_PIF_traditional) %>%
  arrange(-Ratio) %>%
  mutate(Species = cn,
         `Effective Survey Area` = signif(area_ratio,3),
         Availability = signif(availability_ratio,3),
         `Geographic Bias` = signif(habitat_ratio,3),
         `Temporal Trend` = signif(trend_ratio,3),
         `Realised Overall Ratio` = signif(Ratio,3),
         New = signif(pop_median_PIF_eBird_with_EDR_Avail/1e6,3),
         Existing = signif(pop_median_PIF_traditional/1e6,3)) %>%
  select(-c(cn,area_ratio,availability_ratio,habitat_ratio,trend_ratio,
            Ratio,pop_median_PIF_eBird_with_EDR_Avail,pop_median_PIF_traditional))

write_csv(tabl2_paper,"Table_2_paper.csv")

# Trats <- ExpAdjs %>%
#   select(Species,cn,Tlr,Alr,clr) %>%
#   pivot_longer(cols = c(Tlr,Alr,clr),
#                names_to = "adjustment",
#                values_to = "LogRatio") %>%
#   mutate(adjfactor = ifelse(adjustment == "Tlr",
#                             "Availability",
#                             "Area"),
#          adjfactor = ifelse(adjustment == "clr",
#                             "Combined",
#                             adjfactor),
#          adjfactor = factor(adjfactor,
#                             levels = c("Availability","Area","Combined"),
#                             ordered = TRUE)) %>%
#   select(-adjustment) %>%
#   bind_rows(trend_effects)




Trats <- ExpAdjs %>%
  select(Species,cn,Tr,Ar,comb_r) %>%
  pivot_longer(cols = c(Tr,Ar,comb_r),
               names_to = "adjustment",
               values_to = "Ratio") %>%
  mutate(adjfactor = ifelse(adjustment == "Tr",
                            "Availability",
                            "Area"),
         adjfactor = ifelse(adjustment == "comb_r",
                            "Combined",
                            adjfactor)) %>%
  select(-adjustment)

sp_l <- Trats %>%
  select(cn,Species) %>%
  distinct()

USA_CAN_sample_bias <- USA_CAN_sample_bias %>% left_join(sp_l)
trend_effects <- trend_effects %>% left_join(sp_l)
over_bias <- over_bias %>% left_join(sp_l)


Trats_plot <- Trats %>%
  bind_rows(USA_CAN_sample_bias) %>%
  bind_rows(trend_effects) %>%
  bind_rows(over_bias) %>%
  mutate(adjfactor = factor(adjfactor,
                            levels = c("Availability","Area","Combined",
                                       "Geographic Bias","Trend","Overall"),
                            ordered = TRUE))


sp_label <- c("American Robin",
              "Mountain Bluebird",
              "Eastern Phoebe",
              #"Rose-breasted Grosbeak",
              "Canyon Wren",
              "Downy Woodpecker",
              "Blackpoll Warbler",
              "Western Meadowlark",
              "Black-capped Chickadee",
              "Scarlet Tanager",
              "Say's Phoebe",
              #"Brown Creeper",
              "Veery",
              "Purple Martin",
              "Black-chinned Hummingbird")

splabs <- Trats_plot %>%
  filter(cn %in% sp_label)
splab <- splabs %>%
  filter(adjfactor == "Overall")

#
Trats_plot_full <- Trats_plot

Trats_plot <- Trats_plot %>%
  filter(cn %in% sp_example,
         adjfactor != "Combined")

adj_plot <- ggplot(data = Trats_plot,
                   aes(x = adjfactor,
                       y = Ratio))+
  # geom_hline(yintercept = c((10),
  #                           (0.1)),
  #            alpha = 0.8,
  #            linetype = 3)+
  # geom_hline(yintercept = c((5),
  #                           (0.2)),
  #            alpha = 0.6,
  #            linetype = 2)+
  # geom_hline(yintercept = c((2),
  #                           (0.5)),
  #            alpha = 0.6,
  #            linetype = 2)+
  geom_hline(yintercept = c(0.1,0.2,0.5,1,2,5,10),
             alpha = 0.15)+
  geom_hline(yintercept = c(1),
             alpha = 0.8)+
  geom_violin(fill = NA)+
  geom_point(aes(group = Species),position = position_dodge(width = 0.05),
             alpha = 0.1)+
  geom_line(data = splabs,
            aes(group = Species,
                colour = Species),#position = position_dodge(width = 0.05),
            alpha = 0.6)+
  geom_point(data = splabs,
             aes(group = Species,
                 colour = Species),#position = position_dodge(width = 0.05),
             alpha = 1)+
  ylab("Multiplicative change in\npopulation estimate")+
  #  ylab("Ratio new/existing\n values > 1 = increased population estimate")+
  xlab("Adjustment factor")+
  geom_label_repel(data = splab,
                   aes(label = Species,group = Species,
                       colour = Species),
                   #position = position_dodge(width = 0.05),
                   min.segment.length = 0,
                   size = 3,alpha = 0.9,point.padding = 0.3,label.size = 0.5,
                   nudge_x = 1)+
  #scale_colour_brewer(type = "qual",palette = "Dark2")+
  scale_colour_viridis_d(option = "turbo")+
  scale_x_discrete(expand = expansion(c(0.01,0.5)))+
  scale_y_continuous(transform = "log",
                     breaks = c(0.1,0.2,0.5,1,2,5,10))+
  theme_classic()+
  theme(legend.position = "none")

pdf("final_figures/Multiplicative_change_factors_all.pdf",
    width = 6.5,
    height = 4.5)
print(adj_plot)
dev.off()




trat_sel <- Trats_plot %>%
  filter( cn %in% sp_label,
          adjfactor != "Combined") %>%
  mutate(adjfactor = factor(adjfactor,
                            levels = c("Geographic Bias",
                                       "Trend",
                                       "Availability",
                                       "Area",
                                       "Combined",
                                       "Overall"),
                            labels = c("Geographic Bias",
                                       "Trend",
                              "Availability",
                              "Area",
                              "Avail+Area",
                              "Overall"),
                            ordered = TRUE))

splab <- trat_sel %>%
  filter(adjfactor == "Overall",
         cn %in% sp_label)

splabs <- trat_sel %>%
  filter(cn %in% sp_label)

adj_plot2 <- ggplot(data = trat_sel,
                   aes(x = adjfactor,
                       y = Ratio))+
  # geom_hline(yintercept = c((10),
  #                           (0.1)),
  #            alpha = 0.8,
  #            linetype = 3)+
  # geom_hline(yintercept = c((5),
  #                           (0.2)),
  #            alpha = 0.6,
  #            linetype = 2)+
  # geom_hline(yintercept = c((2),
  #                           (0.5)),
  #            alpha = 0.6,
  #            linetype = 2)+
  geom_hline(yintercept = c(0.1,0.2,0.5,1,2,5,10),
             alpha = 0.15)+
  geom_hline(yintercept = c(1),
             alpha = 0.8)+
  geom_violin(fill = NA)+
  geom_point(aes(group = Species),position = position_dodge(width = 0.05),
             alpha = 0.1)+
  geom_line(data = splabs,
            aes(group = Species,
                colour = Species),#position = position_dodge(width = 0.05),
            alpha = 0.6)+
  geom_point(data = splabs,
             aes(group = Species,
                 colour = Species),#position = position_dodge(width = 0.05),
             alpha = 1)+
  ylab("Multiplicative change in\npopulation estimate")+
  #  ylab("Ratio new/existing\n values > 1 = increased population estimate")+
  xlab("Adjustment factor")+
  geom_label_repel(data = splab,
                   aes(label = Species,group = Species,
                       colour = Species),
                   #position = position_dodge(width = 0.05),
                   min.segment.length = 0,
                   size = 3,alpha = 0.9,point.padding = 0.3,label.size = 0.5,
                   nudge_x = 1)+
  #scale_colour_brewer(type = "qual",palette = "Dark2")+
  scale_colour_viridis_d(option = "turbo")+
  scale_x_discrete(expand = expansion(c(0.01,0.5)))+
  scale_y_continuous(transform = "log",
                     breaks = c(0.1,0.2,0.5,1,2,5,10))+
  theme_classic()+
  theme(legend.position = "none")

pdf("final_figures/Multiplicative_change_factors_2.pdf",
    width = 6.5,
    height = 4.5)
print(adj_plot2)
dev.off()







# XY comparison of existing and revised -----------------------------------


pop_compare_realised_wide_all <- pop_compare_realised %>%
  pivot_wider(names_from = version,
              names_sep = "_",
              values_from = c(pop_median,pop_lci_80,pop_lci_95,
                              pop_uci_80,pop_uci_95)) %>%
  mutate(pop_median_PIF_traditional = ifelse(is.na(pop_median_PIF_traditional),
                                             0,
                                             pop_median_PIF_traditional),
         pop_lci_80_PIF_traditional = ifelse(is.na(pop_lci_80_PIF_traditional),
                                             0,
                                             pop_lci_80_PIF_traditional),
         pop_uci_80_PIF_traditional = ifelse(is.na(pop_uci_80_PIF_traditional),
                                             0,
                                             pop_uci_80_PIF_traditional),
         pop_lci_95_PIF_traditional = ifelse(is.na(pop_lci_95_PIF_traditional),
                                             0,
                                             pop_lci_95_PIF_traditional),
         pop_uci_95_PIF_traditional = ifelse(is.na(pop_uci_95_PIF_traditional),
                                             0,
                                             pop_uci_95_PIF_traditional),
         dif_mag_new_trad = pop_median_PIF_eBird_with_EDR_Avail/pop_median_PIF_traditional) %>%
  rename(Ratio = dif_mag_new_trad,
         cn = species) %>%
  mutate(adjfactor = "Overall")

adjs_wide <- Trats_plot %>%
  pivot_wider(id_cols = c(Species,cn),
              names_from = adjfactor,
              values_from = Ratio) %>%
  inner_join(pop_compare_realised_wide_all,
             by = "cn")

splabs <- adjs_wide %>%
  filter(cn %in% sp_label)

line5 <- data.frame(pop_median_PIF_traditional = c(c(1e3,5e8),c(1e3,5e8),c(1e3,5e8),c(1e3,5e8)),
                    pop_median_PIF_eBird_with_EDR_Avail = c(c(1e3,5e8),c(2e3,10e8),c(5e3,25e8),c(10e3,50e8)),
                    multi = c(1,1,2,2,5,5,10,10),
                    multif = factor(c(1,1,2,2,5,5,10,10)))

cross_plot <- ggplot(data = adjs_wide,
                     aes(x = pop_median_PIF_traditional,
                         y = pop_median_PIF_eBird_with_EDR_Avail))+
  geom_line(data = line5,
            aes(group = multi,
                linetype = multif))+
  geom_linerange(aes(xmin = pop_lci_95_PIF_traditional,
                     xmax = pop_uci_95_PIF_traditional),
                 alpha = 0.3)+
  geom_linerange(aes(ymin = pop_lci_95_PIF_eBird_with_EDR_Avail,
                     ymax = pop_uci_95_PIF_eBird_with_EDR_Avail),
                 alpha = 0.3)+
  geom_point(alpha = 0.3)+
  geom_point(data = splabs,
             aes(colour = Species))+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  coord_cartesian(xlim = c(1e5,4e8), ylim = c(1e5,9e9))+
  geom_label_repel(data = splabs,
                   aes(label = Species,group = Species,
                       colour = Species),
                   min.segment.length = 0,
                   size = 3,alpha = 0.9,point.padding = 0.3,label.size = 0.5)+
  scale_colour_viridis_d(option = "turbo")+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("Revised Population estimate (Millions)")+
  xlab("Existing population estimate (Millions)")



pdf("final_figures/xy_comparison_overall.pdf",
    width = 6.5,
    height = 4.5)
print(cross_plot)
dev.off()


# same adjustments different order ----------------------------------------

#
# trat_sel <- Trats_plot %>%
#   filter(!(adjfactor == "Combined" & cn %in% sp_label)) %>%
#   mutate(adjfactor = factor(adjfactor,
#                             levels = c("Availability",
#                                        "Area",
#                                        "Combined",
#                                        "Geographic Bias",
#                                        "Trend",
#                                        "Overall"),
#                             labels = c("Availability",
#                                        "Area",
#                                        "Avail+Area",
#                                        "Geographic Bias",
#                                        "Trend",
#                                        "Overall"),
#                             ordered = TRUE))
#
# splab <- trat_sel %>%
#   filter(adjfactor == "Overall",
#          cn %in% sp_label)
#
# splabs <- trat_sel %>%
#   filter(cn %in% sp_label)
#
# adj_plot2 <- ggplot(data = trat_sel,
#                     aes(x = adjfactor,
#                         y = Ratio))+
#   # geom_hline(yintercept = c((10),
#   #                           (0.1)),
#   #            alpha = 0.8,
#   #            linetype = 3)+
#   # geom_hline(yintercept = c((5),
#   #                           (0.2)),
#   #            alpha = 0.6,
#   #            linetype = 2)+
#   # geom_hline(yintercept = c((2),
#   #                           (0.5)),
#   #            alpha = 0.6,
#   #            linetype = 2)+
#   geom_hline(yintercept = c(0.1,0.2,0.5,1,2,5,10),
#              alpha = 0.15)+
#   geom_hline(yintercept = c(1),
#              alpha = 0.8)+
#   geom_violin(fill = NA)+
#   geom_point(aes(group = Species),position = position_dodge(width = 0.05),
#              alpha = 0.1)+
#   geom_line(data = splabs,
#             aes(group = Species,
#                 colour = Species),#position = position_dodge(width = 0.05),
#             alpha = 0.6)+
#   geom_point(data = splabs,
#              aes(group = Species,
#                  colour = Species),#position = position_dodge(width = 0.05),
#              alpha = 1)+
#   ylab("Multiplicative change in\npopulation estimate")+
#   #  ylab("Ratio new/existing\n values > 1 = increased population estimate")+
#   xlab("Adjustment factor")+
#   geom_label_repel(data = splab,
#                    aes(label = Species,group = Species,
#                        colour = Species),
#                    #position = position_dodge(width = 0.05),
#                    min.segment.length = 0,
#                    size = 3,alpha = 0.9,point.padding = 0.3,label.size = 0.5,
#                    nudge_x = 1)+
#   #scale_colour_brewer(type = "qual",palette = "Dark2")+
#   scale_colour_viridis_d(option = "turbo")+
#   scale_x_discrete(expand = expansion(c(0.01,0.5)))+
#   scale_y_continuous(transform = "log",
#                      breaks = c(0.1,0.2,0.5,1,2,5,10))+
#   theme_classic()+
#   theme(legend.position = "none")
#
# pdf("final_figures/Multiplicative_change_factors_3.pdf",
#     width = 6.5,
#     height = 4.5)
# print(adj_plot2)
# dev.off()
#


# USA-Canada population trajectories --------------------------------------
source("functions/hpdi.R")


inds_out <- NULL
vers <- ""
for(sp_sel in sp_example){

  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)
if(!file.exists(paste0("output/strata_population_posterior_",vers,sp_aou,"_",sp_ebird,".rds"))){
  next
}
strata_population_posterior <- readRDS(paste0("output/strata_population_posterior_",vers,sp_aou,"_",sp_ebird,".rds"))


USA_CAN_trajectories <- strata_population_posterior %>%
  group_by( year, .draw) %>%
  summarise(population = sum(population)) %>%
  group_by(year) %>%
  summarise(pop_mean = mean(population),
            pop_median = median(population),
            pop_lci_80 = interval_function_hpdi(population,probs = 0.1),
            pop_uci_80 = interval_function_hpdi(population,probs = 0.9),
            pop_lci_95 = interval_function_hpdi(population,probs = 0.025),
            pop_uci_95 = interval_function_hpdi(population,probs = 0.975)) %>%
  mutate(species = sp_sel,
         species_ebird = sp_ebird,
         region = "USA-Canada")

inds_out <- bind_rows(inds_out,
                      USA_CAN_trajectories)

}


saveRDS(inds_out,"final_figures/species_population_trajectories.rds")


inds_out <- readRDS("final_figures/species_population_trajectories.rds")

inds_plot <- inds_out %>%
  filter(species %in% c("American Robin","Canyon Wren"))

inds_22 <- inds_plot %>%
  filter(year == yr_ebird)


trajs <- ggplot(data = inds_plot,
                aes(x = year,y=pop_median))+
  geom_ribbon(aes(ymin = pop_lci_80,ymax = pop_uci_80),
              alpha = 0.5)+
  geom_line()+
  geom_pointrange(data = inds_22,aes(ymin = pop_lci_80,ymax = pop_uci_80),
                  size = 0.4)+
  scale_y_continuous(limits = c(0,NA),
                     labels = scales::unit_format(unit = "M", scale = 1e-6),
                     name = "Population size USA and Canada",
                     expand = expansion(mult = c(0,0.05)))+
  scale_x_continuous(breaks = seq(2013,2023,2), name = "")+
  facet_wrap(nrow = 2,vars(species),
             scales = "free")+

  theme_bw()


pdf("Final_figures/times_series_example.pdf",
    width = 3.5,
    height = 5)
print(trajs)
dev.off()


# Figure example strata compare -------------------------------------------
library(patchwork)

plot_t <- vector("list",2)
names(plot_t) <- c("American Robin","Canyon Wren")
vers <- ""
for(sp_sel in names(plot_t)){
  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)

plot_t[[sp_sel]] <- readRDS(paste0("figures/saved_ggplots/trad_vs_new_alt_",vers,sp_aou,"_",sp_ebird,".rds"))+
  guides(colour = guide_legend(nrow = 1, ncol = 7))+
  theme(text = element_text(size = 9))

if(sp_sel == "Canyon Wren"){
  plot_t[[sp_sel]] <- plot_t[[sp_sel]]+
    ylab("")
}

}


pl2 <- plot_t[[1]]+plot_t[[2]]+
  plot_layout(ncol = 2, guides = "collect")&
  theme(legend.position = "bottom")

pdf("Final_figures/sstrata_comparison.pdf",
    width = 7, height = 4)
print(pl2)
dev.off()





# abundance maps ----------------------------------------------------------


library(patchwork)
library(ggtext)
library(tidyverse)
library(terra)
library(sf)
library(tidyterra)

vers <- ""
for(sp_sel in c("American Robin","Canyon Wren")){
  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)

  load(paste0("figures/saved_ggplots/abund_map_alt_",vers,sp_aou,"_",sp_ebird,".rdata"))
  cali_use <- readRDS(paste0("output/calibration_",vers,sp_aou,"_",sp_ebird,".rds"))


  breed_abundance <- readRDS(paste0("data/species_relative_abundance_",yr_ebird,"/",
                                    sp_ebird,
                                    "_derived_breeding_relative_abundance.rds"))

  names(breed_abundance) <- "breeding_abundance"
  breed_abundance_plot <- breed_abundance
  names(breed_abundance_plot) <- "breeding"

  breed_abundance_plot <- breed_abundance_plot %>%
    #st_transform(crs = st_crs(strata)) %>%
    mutate(.,breeding = ifelse(breeding == 0,NA,breeding),
           breeding = (breeding*mean(cali_use$calibration))/9) # per km2 mean density




if(sp_sel == "American Robin"){
countries_plot_amro <- countries_plot
breed_abundance_plot_amro <- breed_abundance_plot
bcrs_plot_amro <- bcrs_plot
bb_amro <- bb
bb_amro["ymin"] <- bb_amro["ymin"]-800000
bb_amro["ymax"] <- bb_amro["ymax"]-400000


abund_map_amro <- ggplot()+
  geom_sf(data = countries_plot_amro, fill = NA)+
  geom_spatraster(data = breed_abundance_plot_amro,
                  maxcell = 16000000)+
  geom_sf(data = bcrs_plot_amro, fill = NA)+
  # geom_sf_text(data = bcrs_plot,aes(label = BCR),
  #              size = 3)+
  #geom_sf(data = routes_buf,aes(colour = mean_count),fill = NA)+
  coord_sf(xlim = c(bb_amro[c("xmin","xmax")]),
           ylim = c(bb_amro[c("ymin","ymax")]))+
  #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA,
                       end = 0.9,
                       name = paste0(sp_sel,"<br>Birds km<sup>-2</sup>"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(legend.title = element_markdown(),
        legend.position = "bottom")+
  labs(title = paste(sp_sel))


}else{
  countries_plot_cawr <- countries_plot
  breed_abundance_plot_cawr <- breed_abundance_plot
  bcrs_plot_cawr <- bcrs_plot
  bb_cawr <- bb
  bb_cawr["ymin"] <- bb_cawr["ymin"]-1000000
  bb_cawr["ymax"] <- bb_cawr["ymax"]-800000

  abund_map_cawr <- ggplot()+
    geom_sf(data = countries_plot_cawr, fill = NA)+
    geom_spatraster(data = breed_abundance_plot_cawr,
                    maxcell = 16000000)+
    geom_sf(data = bcrs_plot_cawr, fill = NA)+
    # geom_sf_text(data = bcrs_plot,aes(label = BCR),
    #              size = 3)+
    #geom_sf(data = routes_buf,aes(colour = mean_count),fill = NA)+
    coord_sf(xlim = c(bb_cawr[c("xmin","xmax")]),
             ylim = c(bb_cawr[c("ymin","ymax")]))+
    #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
    scale_fill_viridis_c(direction = -1,
                         option = "G",
                         na.value = NA,
                         end = 0.9,
                         name = paste0(sp_sel,"<br>Birds km<sup>-2</sup>"))+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(legend.title = element_markdown(),
          legend.position = "bottom")+
    labs(title = paste(sp_sel))


}


#


}


pl2 <- abund_map_amro+abund_map_cawr+
  plot_layout(ncol = 2)

pdf("Final_figures/abundance_maps_demo.pdf",
    width = 7, height = 5)
print(pl2)
dev.off()




# Strata abundance map comparison -----------------------------------------


pop_compare_stack <- pop_compare_stack_sel %>%
  filter(region_type %in% c("strata"),
         species %in% c("American Robin","Canyon Wren")) %>%
  mutate(species_factor = factor(factor(species),levels = levels(pop_compare_wide_usacan$species_factor)))


pop_strata_map <- strata %>%
  inner_join(pop_compare_stack,
            by = c("strata_name" = "region")) %>%
  group_by(species,version) %>%
  mutate(p_pop = pop_median/sum(pop_median))

hlt_col <- viridisLite::inferno(1,begin = 0.6,end = 0.7)
sp_sel <- "American Robin"

pop_strata_map_sel <- pop_strata_map %>%
  mutate(p_pop = ifelse(p_pop < 0.001, NA,p_pop),
         version_factor = factor(version,
                                 levels = c("PIF_eBird_with_EDR_Avail",
                                            "PIF_traditional"),
                                 labels = c("Revised",
                                            "Existing"))) %>%
filter(species_factor == sp_sel,
       !is.na(p_pop))
strat_focus <- strata %>%
  filter(strata_name %in% c("CA-BC-5",
                            "CA-YT-4",
                            "US-AK-4"))

bb <- sf::st_bbox(pop_strata_map_sel)

strata_map_compare <- ggplot()+
  geom_sf(data = bcrs_plot,
          fill = NA)+
  geom_sf(data = filter(pop_strata_map_sel,species_factor == sp_sel),
          aes(fill = p_pop))+
  geom_sf(data = strat_focus,
          fill = NA,
          colour = hlt_col)+
  scale_fill_viridis_c(na.value = NA,
                       name = paste0("Proportion of\n",
                                     "USA Canada\n",sp_sel,"\nPopulation"))+
  coord_sf(xlim = bb[c("xmin","xmax")],
           ylim = bb[c("ymin","ymax")])+
  facet_grid(cols = vars(version_factor))+
  labs(title = sp_sel)+
  theme_bw()+
  theme(legend.position = "right")



sp_sel <- "Canyon Wren"

pop_strata_map_sel2 <- pop_strata_map %>%
  mutate(p_pop = ifelse(p_pop < 0.01, NA,p_pop),
         version_factor = factor(version,
                                 levels = c("PIF_eBird_with_EDR_Avail",
                                            "PIF_traditional"),
                                 labels = c("Revised",
                                            "Existing"))) %>%
  filter(species_factor == sp_sel,
         !is.na(p_pop))

bb2 <- sf::st_bbox(pop_strata_map_sel2)
strat_focus <- strata %>%
  filter(strata_name %in% c("US-AZ-16"))

strata_map_compare_2 <- ggplot()+
  geom_sf(data = strata,
          fill = NA)+
  geom_sf(data = pop_strata_map_sel2,
          aes(fill = p_pop))+
  geom_sf(data = strat_focus,
          fill = NA,
          colour = hlt_col,
          linewidth = 0.6)+
  scale_fill_viridis_c(na.value = NA,
                       name = paste0("Proportion of\n",
                                     "USA Canada\n",sp_sel,"\nPopulation"))+
  coord_sf(xlim = bb2[c("xmin","xmax")],
           ylim = bb2[c("ymin","ymax")])+
  facet_grid(cols = vars(version_factor))+
  labs(title = sp_sel)+
  theme_bw()+
  theme(legend.position = "right")


strata_comp <- strata_map_compare /strata_map_compare_2
print(strata_comp)


pdf(file = "Final_Figures/strata_pop_comparison.pdf",
    width = 7,
    height = 7)
print(strata_comp)
dev.off()


# Strata time-series and seasonal patterns ------------------------------------------------------
sp_sel <- "American Robin"
sp_sel <- "Canyon Wren"


for(sp_sel in c("American Robin","Canyon Wren")){

line_col <- viridisLite::cividis(n = 1)

vers <- ""

sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
sp_ebird <- ebirdst::get_species(sp_sel)

fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
summ <- readRDS(paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

combined <- readRDS(paste0("stan_data/dataframe_",vers,sp_aou,"_",sp_ebird,".rds"))
stan_data <- readRDS(paste0("stan_data/stan_data_",vers,sp_aou,"_",sp_ebird,".rds"))

strata_population_posterior <- readRDS(paste0("output/strata_population_posterior_",vers,sp_aou,"_",sp_ebird,".rds"))




strata_loss <- strata_population_posterior %>%
  group_by(region, region_type, year) %>%
  summarise(pop_median = median(population),
            .groups = "drop") %>%
  filter(region_type == "strata",
         year %in% c(2013,2023)) %>%
  select(-region_type) %>%
  # pivot_wider(values_from = pop_median,
  #             names_from = year,
  #             names_prefix = "I_") %>%
  # mutate(loss = I_2013 - I_2023,
  #        abs_change = abs(loss)) %>%
  # rename(strata_name = region)
  group_by(region) %>%
  summarise(abs_change = mean(pop_median))%>%
  rename(strata_name = region)


strat_names <- combined %>%
  select(strata,strata_name) %>%
  distinct()

yearl_strat <- summ %>% filter(grepl("yeareffect[",variable,fixed = TRUE)) %>%
  mutate(strata = rep(c(1:stan_data$n_strata),times = stan_data$n_years),
         year = rep(c(1:stan_data$n_years),each = stan_data$n_strata) + (min(combined$year)-1)) %>%
  inner_join(strat_names,by = c("strata")) %>%
  inner_join(strata_loss,
             by = "strata_name") %>%
  mutate(med_percent = exp(mean)*100)

bks1 <- c(-25,0,25,50,100,200,500)
bks <- log((bks1/100)+1)
bks_labs <- paste0(bks1,"%")

traj_plot <- ggplot(data = yearl_strat,
                    aes(x = year,y = median,
                        group = strata_name))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2023)+
  # annotate(x = 2023,
  #          y = quantile(yearl_strat$median,0.95),
  #          label = "year of\neBird prediction",
  #          hjust = 1,
  #          geom = "text")+
  geom_line(colour = line_col,
            alpha = 0.3)+
  scale_x_continuous(breaks = seq(2013,2023,by = 2),
                     name = "")+
  scale_y_continuous(breaks = bks,
                     labels = bks_labs)+
  ylab("percent change relative to 2023")+
  #colorspace::scale_color_continuous_diverging(palette = "Blue-Red 3")+
  theme_classic()+
  theme(legend.title = element_markdown(),
        legend.position = "none")

traj_plot


png(filename = paste0("Figures/yearly_effect_inset_",vers,sp_aou,"_",sp_ebird,".png"),
     res = 300,
     height = 4,
     width = 4,
     units = "in")
print(traj_plot)
dev.off()





# Seasonal patterns -------------------------------------------------------

# line_col <- viridisLite::cividis(n = 1)
#
# vers <- ""
#
# sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
# sp_ebird <- ebirdst::get_species(sp_sel)
#
#
# fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
# summ <- readRDS(paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
#
# combined <- readRDS(paste0("stan_data/dataframe_",vers,sp_aou,"_",sp_ebird,".rds"))
# stan_data <- readRDS(paste0("stan_data/stan_data_",vers,sp_aou,"_",sp_ebird,".rds"))


strat_names <- combined %>%
  select(strata,strata_name,day_of_year) %>%
  distinct()

seasonal_strat <- summ %>% filter(grepl("doy_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = rep(c(1:stan_data$n_doy),times = stan_data$n_strata) + (min(combined$day_of_year)-1),
         strata = rep(c(1:stan_data$n_strata),each = stan_data$n_doy)) %>%
  inner_join(strat_names,by = c("strata",
                                "doy" = "day_of_year"))%>%
  mutate(med_percent = (exp(mean)-1)*100,
         cv_val = sd/mean)


strat_labs1 <- seasonal_strat %>%
  group_by(strata_name) %>%
  summarise(cv = mean(sd))

seasonal_strat <- strat_labs1 %>%
  inner_join(seasonal_strat,
             by = c("strata_name"))

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
                     aes(x = doy,y = med_percent,
                         group = strata_name)#,
                         #alpha = cv)
                         )+
  geom_line(colour = line_col,
            alpha = 0.3)+
  #coord_cartesian(xlim = c(100,240))+
  scale_x_continuous(breaks = xlbl$doy, labels = xlbl$days)+
  xlab("")+
  ylab("% difference in count from mean")+
  scale_y_continuous()+
  theme_classic()+
  theme(legend.position = "none")


vis_season
png(filename = paste0("Figures/seasonal_effect_inset_",vers,sp_aou,"_",sp_ebird,".png"),
    res = 300,
    height = 4,
    width = 4,
    units = "in")
print(vis_season)
dev.off()
}




