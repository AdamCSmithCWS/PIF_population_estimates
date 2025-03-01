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
                                 # "American Robin",
                                  "Canyon Wren"),
                      region = c("CA-BC-5",
                                 #"CA-YT-4",
                                 "US-AZ-16"))




abund_maps_save <- vector(mode = "list",length = nrow(demo_sp))
names(abund_maps_save) <- paste0(demo_sp$species,"_",demo_sp$region)


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

strata_sel1 <- strata_sel %>% st_buffer(5e5) # 500km buffer

routes_sel <- routes_buf %>%
  sf::st_transform(crs = st_crs(strata_sel)) %>%
  sf::st_join(y = strata_sel,left = FALSE)

bb <- sf::st_bbox(strata_sel)


strata_sel2 <- strata_sel1 %>%
  sf::st_transform(crs = st_crs(breed_abundance))

bb_crop <- terra::ext(strata_sel2)

breed_abundance_plot <- terra::crop(breed_abundance,y = bb_crop) %>%
  terra::mask(., mask = strata_sel2) #%>%
  #terra::project(y = "epsg:9822")

names(breed_abundance_plot) <- "breeding"

# breed_abundance_plot <- breed_abundance_plot %>%
#   sf::st_transform(crs = sf::st_transform(strata_sel))

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
           ylim = c(bb[c("ymin","ymax")]))+
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
  theme(legend.position = "bottom")

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

des <- "
111222
111222
111222
111222
111333
"
# abund_maps_save[[1]] <- abund_maps_save[[3]]
# abund_maps_save[[2]] <- abund_maps_save[[4]]

# abund_maps_save[[2]] <- abund_maps_save[[2]]


sampling_demo <- abund_maps_save[[1]] + abund_maps_save[[2]] +
  guide_area() + plot_layout(design = des,guides = "collect")

sampling_demo

pdf("final_figures/sampling_bias_maps.pdf",
    width = 7,
    height = 6)
print(sampling_demo)
dev.off()


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

re_summarise <- TRUE # set to true to re-run this summary loop
if(re_summarise){
pop_compare_stack_sel <-NULL

for(sp_sel in sps_list$english){


    sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
    sp_ebird <- ebirdst::get_species(sp_sel)

if(!file.exists(paste0("estimates/pop_ests_alt_",sp_aou,sp_ebird,".csv")) |
   !file.exists(paste0("estimates/pop_ests_alt_trad_",sp_aou,sp_ebird,".csv"))){next}


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
         dif_mag_new_trad = pop_median_PIF_eBird_with_EDR_Avail/pop_median_PIF_traditional,
         dif_mag_eBird_trad = pop_median_PIF_eBird_without_EDR_Avail/pop_median_PIF_traditional,
         dif_mag_napops_eBird = pop_median_PIF_eBird_with_EDR_Avail/pop_median_PIF_eBird_without_EDR_Avail)

pop_compare_wide_usacan <- pop_compare_wide %>%
  filter(region %in% c("USACAN")) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  mutate(species_factor = fct_reorder(species,dif_mag_new_trad))


pop_compare_stack <- pop_compare_stack_sel %>%
  filter(region %in% c("USACAN")) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  mutate(species_factor = factor(species,levels = levels(pop_compare_wide_usacan$species_factor)))

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


pop_compare_stack_for <- pop_compare_stack %>%
  filter(primary_breeding_habitat_major == "Forests") %>%
  arrange(primary_breeding_habitat_sub)


side_plot_forest <- ggplot(data = pop_compare_stack_for,
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

side_plot_forest



bi_plot1 <- ggplot(data = pop_compare_wide_usacan,
                     aes(y = pop_median_PIF_eBird_without_EDR_Avail,
                         x = pop_median_PIF_traditional,
                         colour = primary_breeding_habitat_sub))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbarh(aes(xmin = pop_lci_95_PIF_traditional,
                     xmax = pop_uci_95_PIF_traditional),
                 height = 0, alpha = 0.3)+
  geom_errorbar(aes(ymin = pop_lci_95_PIF_eBird_without_EDR_Avail,
                     ymax = pop_uci_95_PIF_eBird_without_EDR_Avail),
                 width = 0, alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  ylab("Population estimate with eBird map correction")+
  xlab("Traditional Population estimate with 95% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  theme_bw()+
  facet_wrap(vars(primary_breeding_habitat_major),
             scales = "free")

bi_plot1









pop_compare_wide_ab6 <- pop_compare_wide %>%
  filter(region %in% c("CA-AB-6")) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  filter((!primary_breeding_habitat_major %in% c("Aridlands",
                                                 "Tundra")),
         pop_median_PIF_traditional > 0)



min_diag <- min(c(pop_compare_wide_ab6$pop_median_PIF_traditional,
                  pop_compare_wide_ab6$pop_median_PIF_eBird_without_EDR_Avail))

max_diag <- max(c(pop_compare_wide_ab6$pop_median_PIF_traditional,
                  pop_compare_wide_ab6$pop_median_PIF_eBird_without_EDR_Avail))

diag_line <- data.frame(pop_median_PIF_traditional = c(min_diag,max_diag),
                        pop_median_PIF_eBird_without_EDR_Avail = c(min_diag,max_diag))


pop_compare_wide_ab6 <- pop_compare_wide %>%
  filter(region %in% c("CA-AB-6")) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  filter((!primary_breeding_habitat_major %in% c("Aridlands",
                                                 "Tundra")),
          pop_median_PIF_eBird_without_EDR_Avail > 1)


bi_plotab6 <- ggplot(data = pop_compare_wide_ab6,
                   aes(y = pop_median_PIF_eBird_without_EDR_Avail,
                       x = pop_median_PIF_traditional,
                       colour = primary_breeding_habitat_sub))+
  geom_line(data = diag_line,
            aes(y = pop_median_PIF_eBird_without_EDR_Avail,
                x = pop_median_PIF_traditional),
            inherit.aes = FALSE)+
  geom_vline(xintercept = min(diag_line$pop_median_PIF_traditional),
             colour = grey(0.7))+
  geom_point()+
  geom_errorbarh(aes(xmin = pop_lci_95_PIF_traditional,
                     xmax = pop_uci_95_PIF_traditional),
                 height = 0, alpha = 0.3)+
  geom_errorbar(aes(ymin = pop_lci_95_PIF_eBird_without_EDR_Avail,
                    ymax = pop_uci_95_PIF_eBird_without_EDR_Avail),
                width = 0, alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  ylab("Population estimate with eBird map correction")+
  xlab("Traditional Population estimate with 95% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  theme_bw()+
  facet_wrap(vars(primary_breeding_habitat_major),
             scales = "free")

bi_plotab6



bi_plotab62 <- ggplot(data = pop_compare_wide_ab62,
                     aes(y = pop_median_PIF_eBird_without_EDR_Avail,
                         x = pop_median_PIF_traditional,
                         colour = primary_breeding_habitat_sub))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbarh(aes(xmin = pop_lci_95_PIF_traditional,
                     xmax = pop_uci_95_PIF_traditional),
                 height = 0, alpha = 0.3)+
  geom_errorbar(aes(ymin = pop_lci_95_PIF_eBird_without_EDR_Avail,
                    ymax = pop_uci_95_PIF_eBird_without_EDR_Avail),
                width = 0, alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  ylab("Population estimate with eBird map correction")+
  xlab("Traditional Population estimate with 95% CI (Millions)")+
  theme_bw()+
  facet_wrap(vars(primary_breeding_habitat_major),
             scales = "free")

bi_plotab62

bi_plotfactors1 <- ggplot(data = pop_compare_wide_usacan,
                   aes(y = dif_mag_new_trad,
                       x = Cd_diff,
                       colour = primary_breeding_habitat_sub))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  xlab("Magnitude difference in distance correction")+
  ylab("Magnitude difference in napops over eBird informed")+
  theme_bw()

bi_plotfactors1





bi_plotfactors1 <- ggplot(data = pop_compare_wide_usacan,
                          aes(y = dif_mag_new_trad,
                              x = Cd_diff,
                              colour = primary_breeding_habitat_sub))+
  geom_abline(slope = 1, intercept = 0)+
  geom_hline(yintercept = 1,colour = grey(0.6))+
  geom_point()+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  xlab("Magnitude difference in distance correction")+
  ylab("Magnitude difference in new over traditional")+
  theme_bw()

bi_plotfactors1





bi_plotfactors2 <- ggplot(data = pop_compare_wide_usacan,
                          aes(y = dif_mag_napops_eBird,
                              x = Ct_diff,
                              colour = primary_breeding_habitat_sub))+
  geom_abline(slope = 1, intercept = 0)+
  geom_hline(yintercept = 1,colour = grey(0.6))+
  geom_point()+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  xlab("Magnitude difference in TOD correction")+
  ylab("Magnitude difference in napops over eBird informed")+
  theme_bw()

bi_plotfactors2



bi_plotfactors3 <- ggplot(data = pop_compare_wide_usacan,
                          aes(y = Cd_diff,
                              x = Ct_diff,
                              colour = primary_breeding_habitat_major))+
  geom_abline(slope = 1, intercept = 0)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_point()+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  xlab("Magnitude difference in TOD correction")+
  ylab("Magnitude difference in distance correction")+
  theme_bw()

bi_plotfactors3




pop_compare_wide_usacan_forest <- pop_compare_wide_usacan %>%
  filter(primary_breeding_habitat_major == "Forests")


bi_plot_forest <- ggplot(data = pop_compare_wide_usacan_forest,
                   aes(y = pop_median_PIF_eBird_without_EDR_Avail,
                       x = pop_median_PIF_traditional,
                       colour = primary_breeding_habitat_sub))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbarh(aes(xmin = pop_lci_95_PIF_traditional,
                     xmax = pop_uci_95_PIF_traditional),
                 height = 0, alpha = 0.3)+
  geom_errorbar(aes(ymin = pop_lci_95_PIF_eBird_without_EDR_Avail,
                    ymax = pop_uci_95_PIF_eBird_without_EDR_Avail),
                width = 0, alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  ylab("Population estimate with eBird map correction")+
  xlab("Traditional Population estimate with 95% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  theme_bw()

bi_plot_forest



pop_compare_wide_ab6 <- pop_compare_wide %>%
  filter(region %in% c("CA-AB-6"),
         !is.na(dif_mag_eBird_trad)) %>%
  left_join(Cor_factors,
            by = c("species" = "cn",
                   "aou" = "Aou")) %>%
  mutate(species_factor = fct_reorder(species,dif_mag_eBird_trad))


bi_plot_ab6 <- ggplot(data = pop_compare_wide_ab6,
                         aes(y = pop_median_PIF_eBird_without_EDR_Avail,
                             x = pop_median_PIF_traditional,
                             colour = primary_breeding_habitat_major))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbarh(aes(xmin = pop_lci_95_PIF_traditional,
                     xmax = pop_uci_95_PIF_traditional),
                 height = 0, alpha = 0.3)+
  geom_errorbar(aes(ymin = pop_lci_95_PIF_eBird_without_EDR_Avail,
                    ymax = pop_uci_95_PIF_eBird_without_EDR_Avail),
                width = 0, alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0)+
  scale_colour_viridis_d(end = 1, name = "Forest type",
                         option = "H")+
  #scale_colour_viridis_c(end = 0.95, name = "Area")+
  ylab("Population estimate with eBird map correction")+
  xlab("Traditional Population estimate with 95% CI (Millions)")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  theme_bw()

bi_plot_ab6



points_plot <- ggplot(data = pop_compare_wide_usacan,
                      aes(x = primary_breeding_habitat_major,
                          y = dif_mag_new_trad))+
  geom_boxplot()+
  geom_hline(yintercept = 1, colour = grey(0.6))+
  geom_point(aes(group = species),
             position = position_dodge(width = 0.05),
             alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0,
                  position = position_dodge(width = 0.05))+
  scale_y_continuous(transform = "log10",
                     breaks = c(1,2,3,5,10,30))+
  ylab("Magnitude increase in population estimate")+
  theme_bw()


points_plot




points_plot <- ggplot(data = pop_compare_wide_usacan,
                      aes(x = primary_breeding_habitat_major,
                          y = dif_mag_eBird_trad))+
  geom_boxplot()+
  geom_hline(yintercept = 1, colour = grey(0.6))+
  geom_point(aes(group = species),
             position = position_dodge(width = 0.05),
             alpha = 0.3)+
  geom_text_repel(aes(label = species_ebird),
                  size = 2,min.segment.length = 0,
                  position = position_dodge(width = 0.05))+
  scale_y_continuous(transform = "log10",
                     breaks = c(1,2,3,5,10,30))+
  ylab("Magnitude increase in population estimate from eBird")+
  theme_bw()


points_plot



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





# compare adjustment factors ----------------------------------------------



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

Trats <- ExpAdjs %>%
  select(Species,cn,Tlr,Alr,clr) %>%
  pivot_longer(cols = c(Tlr,Alr,clr),
               names_to = "adjustment",
               values_to = "LogRatio") %>%
  mutate(factor = ifelse(adjustment == "Tlr",
                         "Availability",
                         "Distance"),
         factor = ifelse(adjustment == "clr",
                         "Combined",
                         factor))

sp_label <- c("American Robin",
              "Mountain Bluebird",
              #"Eastern Phoebe",
              #"Rose-breasted Grosbeak",
              "Canyon Wren",
              "Downy Woodpecker",
              #"Red-winged Blackbird",
              #"Western Meadowlark",
              #"Black-capped Chickadee",
              "Scarlet Tanager",
              "Say's Phoebe",
              #"Brown Creeper",
              "Black-chinned Hummingbird")

splabs <- Trats %>%
  filter(cn %in% sp_label)
splab <- splabs %>%
  filter(factor == "Distance")

adj_plot <- ggplot(data = Trats,
                      aes(x = factor,
                          y = LogRatio))+
  geom_hline(yintercept = c(log(10),
                            log(0.1)),
             alpha = 0.6)+
  geom_hline(yintercept = c(log(5),
                            log(0.2)),
             alpha = 0.6,
             linetype = 2)+
  geom_hline(yintercept = c(0),
             alpha = 0.3)+
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
  ylab("Log Ratio log(new/existing)\n positive values = increased population estimate")+
  xlab("Adjustment factor")+
  geom_label_repel(data = splab,
                  aes(label = Species,group = Species,
                      colour = Species),
            #position = position_dodge(width = 0.05),
            min.segment.length = 0,
            size = 3,alpha = 0.9,point.padding = 0.3,label.size = 0.5,
            nudge_x = 0.6)+
  scale_colour_brewer(type = "qual",palette = "Dark2")+
  scale_x_discrete(expand = expansion(c(0.2,0.5)))+
  theme_classic()+
  theme(legend.position = "none")

pdf("final_figures/Log-Ratios.pdf",
    width = 6.5,
    height = 4.5)
print(adj_plot)
dev.off()



table_2_full <- ExpAdjs %>%
  mutate(species = cn,
         sp_code = Species,
         distance_existing = Dist2,
         distance_existing_lower = Dist.Lower,
         distance_existing_upper = Dist.Upper,
         distance_edr = edr,
         distance_edr_sd = edr_sd,
         edr_model = EDR_model,
         time_of_day_existing = 1/exp(TimeAdj.meanlog),
         tod_uci = 1/exp(TimeAdj.meanlog+TimeAdj.sdlog),
         tod_lci = 1/exp(TimeAdj.meanlog-TimeAdj.sdlog),
         time_of_day_existing_approx_sd = (tod_lci-tod_uci)/2,
         availability = availability,
         availability_sd = availability_sd,
         availability_model = availability_model,
         log_ratio_distance = Alr,
         log_ratio_availability = Tlr,
         log_ratio_combined = clr) %>%
  select(species,
         sp_code,
         distance_existing,
         distance_existing_lower,
         distance_existing_upper,
         distance_edr,
         distance_edr_sd,
         edr_model,
         time_of_day_existing,
         time_of_day_existing_approx_sd,
         availability,
         availability_sd,
         availability_model,
         log_ratio_distance,
         log_ratio_availability,
         log_ratio_combined) %>%
  mutate(across(where(is.double), ~round(.x,3)),
         sp_code = ifelse(is.na(sp_code),"",sp_code))

write_excel_csv(table_2_full,
                "adjustment_factors_existing_new_PIF_pop_est.csv")





tmp <- table_2_full %>%
  filter(species %in% sp_label)



# comparison of estimates using non-ebird approach ------------------------

out_long <- NULL
for(EDR_for_all in c(TRUE,FALSE)){

Output.dir <- ifelse(EDR_for_all,"Stanton_2019_code/output_EDR_alldata/",
                     "Stanton_2019_code/output_alldata/")

tmp1 <- read_csv(paste(Output.dir, 'PSest_Global_WH_NA_',"_", 100, 'iter.csv', sep='')) %>%
  rename_with(.fn= ~tolower(gsub(".","_",.x, fixed = TRUE))) %>%
  mutate(adjustments = ifelse(EDR_for_all,"New","Existing"))

out_long <- bind_rows(out_long,tmp1)

pref <- ifelse(EDR_for_all,"New","Existing")
tmp1 <- tmp1 %>%
  rename_with(.cols = -c(cn,aou),.fn = ~paste(pref,.x,sep = "_"))


if(EDR_for_all){
  out_wide <- tmp1
}else{
  out_wide <- out_wide %>%
    inner_join(tmp1,by = c("cn","aou"))
}

}



splabs <- Trats %>%
  filter(factor == "Combined")

out_wide_sel <- out_wide %>%
  rename(New_distance = New_edr____ifelse_use_edr__true__false_,
         New_availability = New_availability____ifelse_use_availability__true__false_) %>%
  filter(New_distance == TRUE) %>%
  left_join(splabs)


out_wide_lab <- out_wide_sel %>%
  filter(cn %in% sp_label)

line5 <- data.frame(Existing_uscan_med_popest = c(c(1e3,5e8),c(1e3,5e8),c(1e3,5e8),c(1e3,5e8)),
                    New_uscan_med_popest = c(c(1e3,5e8),c(2e3,10e8),c(5e3,25e8),c(10e3,50e8)),
                    multi = c(1,1,2,2,5,5,10,10),
                    multif = factor(c(1,1,2,2,5,5,10,10)))

cross_plot <- ggplot(data = out_wide_sel,
                     aes(x = Existing_uscan_med_popest,
                         y = New_uscan_med_popest))+
  geom_line(data = line5,
            aes(group = multi,
                linetype = multif))+
  geom_linerange(aes(xmin = Existing_uscan_95lci_popest,
                      xmax = Existing_uscan_95uci_popest),
                 alpha = 0.3)+
  geom_linerange(aes(ymin = New_uscan_95lci_popest,
                     ymax = New_uscan_95uci_popest),
                 alpha = 0.2)+
  geom_point(alpha = 0.3)+
  geom_point(data = out_wide_lab,
             aes(colour = Species))+
  scale_colour_brewer(type = "qual",palette = "Dark2")+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),
                     transform = "log10")+
  geom_label_repel(data = out_wide_lab,
                   aes(label = Species,group = Species,
                       colour = Species),
                   min.segment.length = 0,
                   size = 3,alpha = 0.9,point.padding = 0.3,label.size = 0.5)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("Population estimate \n updated for distance and availability")+
  xlab("Existing population estimate")



pdf("final_figures/xy_comparison.pdf",
    width = 6.5,
    height = 4.5)
print(cross_plot)
dev.off()




comparison_table <- out_wide %>%
  rename(Used_EDR = New_edr____ifelse_use_edr__true__false_,
         Used_availability = New_availability____ifelse_use_availability__true__false_,
         species = cn) %>%
  select(species,aou,
         Used_EDR,
         Used_availability,
         contains("New_uscan"),
         contains("Existing_uscan"),
         contains("New_ps_"),
         contains("Existing_ps_")) %>%
  left_join(table_2_full,
            by = c("species" = "species")) %>%
  relocate(species,sp_code,aou)


write_csv(comparison_table,"comparison_before_eBird_adjustment.csv")


table_2 <- comparison_table %>%
  filter(species %in% sp_label) %>%
  select(species,sp_code,
         distance_existing,
         distance_edr,
         time_of_day_existing,
         availability,
         New_uscan_med_popest,
         Existing_uscan_med_popest,
         log_ratio_distance,
         log_ratio_availability,
         log_ratio_combined) %>%
  mutate(ratio_distance = exp(log_ratio_distance),
         ratio_availability = exp(log_ratio_availability),
         ratio_combined = exp(log_ratio_combined),
         ratio_realised = New_uscan_med_popest/Existing_uscan_med_popest,
         New_uscan_med_popest = New_uscan_med_popest/1e6,
         Existing_uscan_med_popest = Existing_uscan_med_popest/1e6) %>%
  select(species,sp_code,ratio_distance,ratio_availability,ratio_realised,
         New_uscan_med_popest,Existing_uscan_med_popest) %>%
  mutate(across(.cols = ratio_distance:Existing_uscan_med_popest,
                .fns = ~signif(.x,3)))

write_csv(table_2,"Table_2.csv")




library(tidybayes)
library(tidyverse)
yr_ebird <- 2022 # prediction year for eBird relative abundance

output_dir <- "G:/PIF_population_estimates/output"


params <- NULL
for(sp_sel in c("Western Meadowlark","Baird's Sparrow","Brown Creeper",
                "Canyon Wren",
                "Black-capped Chickadee","American Robin",
                "Barn Swallow",
                "Blackpoll Warbler",
                "Mountain Bluebird",
                "Eastern Phoebe",
                "Rose-breasted Grosbeak",
                "Downy Woodpecker",
                "Red-winged Blackbird",
                "Scarlet Tanager",
                "Say's Phoebe",
                "Black-chinned Hummingbird"
)){
  for(vers in c("trad_rho","trad","_rho","")[c(4)]){
    # rev(sps_list$english[1:348])){#sp_example[-wh_drop]){#list$english){
    #sp_sel = "Bank Swallow"
    #for(sp_sel in rev(sps_list$english)[29]){
    sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
    sp_ebird <- ebirdst::get_species(sp_sel)


    param_infer2 <- readRDS(paste0(output_dir,"/parameter_inference_alt_",vers,sp_aou,"_",sp_ebird,".rds")) %>%
      mutate(version = vers)

    params <- bind_rows(params,param_infer2)
  }# sp
}# vers


retrans <- params %>%
  filter(grepl("alibration",inference),
         grepl("mean",inference))


rhos <- params %>% filter(variable == "RHO", !is.na(rhat))
area_surveyed <- params %>% filter(grepl("c_area",variable))




# summarise sampling bias -------------------------------------------------

w_mean_lg_rat <- function(x,w){
  wx <- rep(NA,length(x))
  sum_w <- sum(w,na.rm = TRUE)

  for(i in 1:length(x)){
    if(is.finite(x[i])){
      wx[i] <- (x[i]*w[i])/sum_w
    }else{
      wx[i] <- (log(0.01)*w[i])/sum_w
    }
  }
  return(sum(wx,na.rm = TRUE))
}

sampl_bias <- NULL

for(sp_sel in c("Western Meadowlark","Baird's Sparrow","Brown Creeper",
                "Canyon Wren",
                "Black-capped Chickadee","American Robin",
                "Barn Swallow",
                "Blackpoll Warbler",
                "Mountain Bluebird",
                "Eastern Phoebe",
                "Rose-breasted Grosbeak",
                "Downy Woodpecker",
                "Red-winged Blackbird",
                "Scarlet Tanager",
                "Say's Phoebe",
                "Black-chinned Hummingbird"
)){

  species_ebird <- ebirdst::get_species(sp_sel)

  # route_sampled <- readRDS(paste0("data/species_relative_abundance/",species_ebird,"_relative_abundance.rds"))
  #
  # mean_sampled_abund <- mean(route_sampled$ebird_abund,na.rm = TRUE)
  #

  sampl_bias1 <- readRDS(paste0("data/species_relative_abundance/",species_ebird,"_bbs_strata_relative_abundance.rds")) %>%
    mutate(species = sp_sel,
           sp_ebird = species_ebird)

  tmp <- sampl_bias1 %>%
    filter(strata_name != "USA_CAN",
           !is.na(log_ratio)) %>%
    mutate(ifelse(is.na(sampled_abundance),0,sampled_abundance))
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





# calculate the ratio of existing and ebird estiates with edr -------------


output_dir <- "G:/PIF_population_estimates/output"


all_pop_ests <- NULL

vers <- ""
for(sp_sel in c("Western Meadowlark","Baird's Sparrow","Brown Creeper",
                "Canyon Wren",
                "Black-capped Chickadee","American Robin",
                "Barn Swallow",
                "Blackpoll Warbler",
                "Mountain Bluebird",
                "Eastern Phoebe",
                "Rose-breasted Grosbeak",
                "Downy Woodpecker",
                "Red-winged Blackbird",
                "Scarlet Tanager",
                "Say's Phoebe",
                "Black-chinned Hummingbird"
)){

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

usa_can_exist <- read_csv("Stanton_2019_code/output_alldata/PSest_Global_WH_NA__100iter.csv") %>%
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




# calculate the time-series of strata-level estimates -------------




trend_effect_out <- NULL
inds_out <- NULL
vers <- ""
for(sp_sel in c("Western Meadowlark","Baird's Sparrow","Brown Creeper",
                "Canyon Wren",
                "Black-capped Chickadee","American Robin",
                "Barn Swallow",
                "Blackpoll Warbler",
                "Mountain Bluebird",
                "Eastern Phoebe",
                "Rose-breasted Grosbeak",
                "Downy Woodpecker",
                "Red-winged Blackbird",
                "Scarlet Tanager",
                "Say's Phoebe",
                "Black-chinned Hummingbird"
)){

  sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)

  if(file.exists(paste0(output_dir,"/calibration_fit_alt_",
                        vers,sp_aou,"_",sp_ebird,".rds"))){

    fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",
                          vers,sp_aou,"_",sp_ebird,".rds"))

    raw_dat <- readRDS(paste0("data/main_data_df_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
    strats <- raw_dat %>%
      select(strata_name,strata) %>%
      distinct() %>%
      arrange(strata)

    n_yr <- max(raw_dat$yr)
    n_strat <- max(raw_dat$strata)
    fyr <- min(raw_dat$year)


    yrs <- data.frame(year = fyr:(fyr+n_yr-1),
                      yr = 1:n_yr)

    yr22 <- yrs %>%
      filter(year == yr_ebird) %>%
      select(yr) %>%
      unlist() %>%
      unname()

    pop_ests <- read_csv(paste0("estimates/pop_ests_alt_",vers,sp_aou,sp_ebird,".csv")) %>%
      filter(region_type == "strata")

    pop_ests <- pop_ests %>%
      left_join(strats,by = c("strata_name"))%>%
      mutate(w_bbs = ifelse(is.na(strata),FALSE,TRUE))

    p_2022_w_bbs <- pop_ests %>%
      group_by(w_bbs) %>%
      summarise(pop = sum(pop_median)) %>%
      mutate(p_pop = pop/sum(pop))

    time_series <- pop_ests %>%
      filter(!is.na(strata)) %>%
      expand_grid(yrs)


    year_draws <- gather_draws(fit, yeareffect[strata,yr]) %>%
      mutate(.draw = factor(.draw))


    time_series <- time_series %>%
      inner_join(year_draws, by = c("strata","yr")) %>%
      mutate(pop_time = pop_median*exp(.value)) %>%
      ungroup()

    inds_strat <- time_series %>%
      group_by(species,species_ebird,
               region,region_type,strata_name,country,country_code,prov_state,bcr,
               year,yr) %>%
      summarise(med = median(pop_time),
                lci = quantile(pop_time,0.05),
                uci = quantile(pop_time,0.95),
                .groups = "drop")

    inds_all <- time_series %>%
      group_by(species,species_ebird,
               year,yr,
               .draw) %>%
      summarise(pop_time = sum(pop_time),
                .groups = "drop") %>%
      group_by(species,species_ebird,
               year,yr) %>%
      summarise(med = median(pop_time),
                lci = quantile(pop_time,0.05),
                uci = quantile(pop_time,0.95),
                .groups = "drop") %>%
      mutate(region = "Survey wide")


    yr_rest <- time_series %>%
      group_by(species,species_ebird,
               year,
               .draw) %>%
      summarise(pop_time_s = sum(pop_time),
                .groups = "drop") %>%
      group_by(species,species_ebird,
               .draw) %>%
      summarise(pop_time_s = mean(pop_time_s),
                .groups = "drop") %>%
      group_by(species,species_ebird) %>%
      summarise(med = median(pop_time_s),
                lci = quantile(pop_time_s,0.05),
                uci = quantile(pop_time_s,0.95),
                .groups = "drop")%>%
      mutate(yr_eBird = FALSE)

    yr_22 <- time_series %>%
      filter(year == yr_ebird) %>%
      group_by(species,species_ebird,
               .draw) %>%
      summarise(pop_time_s = sum(pop_time),
                .groups = "drop") %>%
      group_by(species,species_ebird) %>%
      summarise(med = median(pop_time_s),
                lci = quantile(pop_time_s,0.05),
                uci = quantile(pop_time_s,0.95),.groups = "drop") %>%
      mutate(yr_eBird = TRUE)



    yr_22[1,"med"]/yr_rest[1,"med"]

    trend_effect <- bind_rows(yr_rest,yr_22) %>%
      mutate(species = sp_sel,
             sp_eBird = sp_ebird)

    tmp <- ggplot(data = inds_all, aes(x = year, y = med))+
      geom_ribbon(aes(ymin = lci,ymax = uci),
                  alpha = 0.2)+
      geom_line()+
      scale_y_continuous(limits = c(0,NA),
                         labels = scales::unit_format(unit = "M", scale = 1e-6))+
      theme_bw()+
      labs(title = paste(sp_sel, "\n excluding",round(unlist(unname(p_2022_w_bbs[p_2022_w_bbs$w_bbs == FALSE,"p_pop"])),3)*100,
                         "Percent of USA+Canada population outside BBS monitored strata"),
           subtitle = paste(round(as.numeric(yr_22[1,"med"]/yr_rest[1,"med"]),2),
                            "expected difference in 2022 population from average 2013-2023"))+
      xlab("")+
      ylab("Annual population size \n in BBS surveyed range")

    print(tmp)

    pdf(paste0("figures/pop_trajectory_",sp_ebird,".pdf"))
    print(tmp)
    dev.off()


    trend_effect_out <- bind_rows(trend_effect_out,trend_effect)
    inds_out <- bind_rows(inds_out,inds_strat,inds_all)


  }

}

trend_effects <- trend_effect_out %>%
  select(-c(lci,uci)) %>%
  pivot_wider(values_from = c(med),
              names_from = yr_eBird,
              names_prefix = "y") %>%
  rename(yr_2022 = yTRUE,
         yr_others = yFALSE) %>%
  mutate(ratio = yr_2022/yr_others,
         log_ratio = log(ratio))


write_csv(trend_effects,"final_figures/trend_effect_summaries.csv")
saveRDS(inds_out,"final_figures/species_population_trajectories.rds")


bchu <- sampl_bias %>%
  filter(species == "Black-chinned Hummingbird",
         !is.na(mean_ebird_abundance),
         mean_ebird_abundance > 0) %>%
  select(strata_name,ratio_new_old)

strata_trajs <- inds_out %>%
  filter(species == "Black-chinned Hummingbird",
         region_type == "strata")  %>%
  full_join(bchu, by = c("region" = "strata_name")) %>%
  mutate(ratio_log = log(ratio_new_old))





  traj_p <- ggplot(data = strata_trajs,
                 aes(x = year, y = med,
                     group = region,
                     colour = ratio_log))+
  geom_line()+
  scale_y_continuous(limits = c(0,NA),
                     labels = scales::unit_format(unit = "M", scale = 1e-6))+
  theme_bw()+
    scale_colour_viridis_c(option = "turbo", direction = -1)

traj_p


tmp <- strata_trajs %>%
  filter(year %in% c(2013,2023))









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




