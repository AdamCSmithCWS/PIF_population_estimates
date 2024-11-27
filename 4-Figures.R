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
library(napops)
library(tidybayes)


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


sps_list <- readRDS("data/all_bbs_data_species_list.rds") #

# Load BBS data -----------------------------------------------------------


# mean counts to id routes on which species has been observed -------------
mean_counts_all <- readRDS("data/mean_counts_by_route.rds")

raw_counts_all <- readRDS("data/all_counts_by_route.rds")

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



# Example figure for roadside bias ----------------------------------------------------------
demo_sp <- data.frame(species = c("American Robin",
                                  "Clay-colored Sparrow",
                                  "Varied Thrush",
                                  "Horned Lark",
                                  "Mountain Chickadee",
                                  "Mountain Bluebird",
                                  "Black-throated Green Warbler"),
                      region = c("CA-YT-4",
                                 "CA-AB-6",
                                 "CA-BC-9",
                                 "US-CA-33",
                                 "CA-BC-9",
                                 "US-ID-9",
                                 "CA-AB-6"))






abund_maps_save <- vector(mode = "list",length = nrow(demo_sp))
names(abund_maps_save) <- demo_sp$species


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

breed_abundance <- readRDS(paste0("data/species_relative_abundance/",
                                  sp_ebird,
                                  "_derived_breeding_relative_abundance.rds"))

strata_sel <- strata %>%
  filter(strata_name == reg_sel)

routes_sel <- routes_buf %>%
  sf::st_transform(crs = st_crs(strata_sel)) %>%
  sf::st_join(y = strata_sel,left = FALSE)

bb <- sf::st_bbox(strata_sel)


strata_sel2 <- strata_sel %>%
  sf::st_transform(crs = st_crs(breed_abundance))
bb_crop <- terra::ext(strata_sel2)

breed_abundance_plot <- terra::crop(breed_abundance,y = bb_crop) %>%
  terra::mask(., mask = strata_sel2) #%>%
  #terra::project(y = "epsg:9822")

names(breed_abundance_plot) <- "breeding"

# breed_abundance_plot <- breed_abundance_plot %>%
#   sf::st_transform(crs = sf::st_transform(strata_sel))

abund_map <- ggplot()+
  geom_sf(data = strata_sel, fill = grey(0.95))+
  geom_spatraster(data = breed_abundance_plot,maxcell = 16000000)+
  geom_sf(data = strata_sel, fill = NA)+
  # geom_sf_text(data = strata,aes(label = strata_name),
  #              size = 1)+
  geom_sf(data = routes_sel,fill = NA,colour = "white")+
  coord_sf(xlim = c(bb[c("xmin","xmax")]),
           ylim = c(bb[c("ymin","ymax")]))+
  #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA,
                       end = 0.9,
                       name = "Relative Abundance")+
  # scale_colour_viridis_c(direction = -1,
  #                        option = "rocket",
  #                        na.value = NA,
  #                        end = 0.9,
  #                        begin = 0.1,
  #                        name = "Mean count BBS")+
  theme_bw()+
  labs(title = paste(sp_sel))


#

png(filename = paste0("Figures/Figure_1_",sp_aou,"_",sp_ebird,".png"),
    res = 400,
    height = 7,
    width = 6.5,
    units = "in")
print(abund_map)
dev.off()

abund_maps_save[[sp_sel]] <- abund_map

}



# continental estimate summaries ------------------------------------------

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
                "Horned Lark",
                "Eastern Whip-poor-will", "Common Nighthawk")


ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")
pop_ests_out_trad <- readRDS("all_traditional_pop_estimates.rds")


pop_compare_stack_sel <-NULL

for(sp_sel in sp_example[-8]){


    sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
    sp_ebird <- ebirdst::get_species(sp_sel)

if(!file.exists(paste0("estimates/pop_ests_alt_",sp_aou,sp_ebird,".csv"))){next}


pop_ests_out <- read_csv(paste0("estimates/pop_ests_alt_",
                                sp_aou,sp_ebird,".csv"))%>%
  mutate(version = "PIF_eBird_with_EDR_Avail")

pop_ests_out_noEDR <- read_csv(paste0("estimates/pop_ests_alt_trad_",
                                sp_aou,sp_ebird,".csv"))%>%
  mutate(version = "PIF_eBird_without_EDR_Avail")


pop_ests_out_trad_sel <- pop_ests_out_trad %>%
  filter(species == sp_sel) %>%
  mutate(version = "PIF_traditional")


pop_compare_stack <- pop_ests_out %>%
  bind_rows(pop_ests_out_noEDR) %>%
  select(region,region_type,version,
         species,species_ebird,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95) %>%
  bind_rows(pop_ests_out_trad_sel)


pop_compare_stack_sel <- bind_rows(pop_compare_stack_sel,pop_compare_stack)
} # end species loop

saveRDS(pop_compare_stack_sel,
        "estimates/comparison_estimates_method.rds")
pop_compare_stack <- pop_compare_stack_sel %>%
  filter(region %in% c("USACAN"))

side_plot1 <- ggplot(data = pop_compare_stack,
                    aes(y = species,
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


pop_compare_stack <- pop_compare_stack_sel %>%
  filter(region %in% c("CA-AB-6"),
         !is.na(pop_median),
         pop_median > 0)

side_plot1 <- ggplot(data = pop_compare_stack,
                     aes(y = species,
                         x = pop_median,
                         colour = version))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbarh(aes(xmin = pop_lci_80,
                     xmax = pop_uci_80),
                 height = 0,
                 position = position_dodge(width = 0.5))+
  scale_colour_viridis_d(end = 0.8)+
  ylab("Alberta BCR-6 estimates")+
  xlab("Population estimate with 80% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  theme_bw()

side_plot1


