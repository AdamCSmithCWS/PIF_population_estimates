
library(sf)
library(tidyverse)
library(ebirdst)
library(patchwork)
library(terra)
library(bbsBayes2)
library(cmdstanr)
library(ggrepel)
library(tidyterra)


routes_buf_all <- readRDS("data/all_routes_buffered.rds")

# load traditional estimates for later comparison -------------------------

strata_trad_all <- read_csv("Stanton_2019_code/output/PSest_Prov_by_BCR__100iter.csv")

global_trad_all <- read_csv("Stanton_2019_code/output/PSest_Global_WH_NA__100iter.csv")

prov_trad_all <- read_csv("Stanton_2019_code/output/PSest_by_Prov__100iter.csv")


# ASK Pete where the national estimates get calculated --------------------





strata <- bbsBayes2::load_map("bbs_usgs")

strata_names <- strata %>%
  sf::st_drop_geometry()



# Insert species loop -----------------------------------------------------



sp_sel <- "Canada Warbler"
sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]


strata_trad <- strata_trad_all %>%
  filter(cn == sp_sel)

global_trad <- global_trad_all %>%
  filter(cn == sp_sel)




mean_counts <- readRDS("data/mean_counts_by_route.rds") %>%
  filter(english == sp_sel,
         mean_count > 0)

raw_counts <- readRDS("data/all_counts_by_route.rds") %>%
  filter(english == sp_sel)

routes_buf <- routes_buf_all %>%
  filter(route_name %in% mean_counts$route_name)

sp_ebird <- ebirdst::get_species(sp_sel)

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

Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE)

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
                  c_d = 100,
                  sd_c_d = 20,
                  c_t = exp(adjs$TimeAdj.meanlog + 0.5*(adjs$TimeAdj.sdlog^2)),
                  sd_c_t = exp(2*adjs$TimeAdj.meanlog + adjs$TimeAdj.sdlog^2)*(exp(adjs$TimeAdj.sdlog^2)-1),
                  c_d_lower = adjs$Dist.Lower,
                  c_d_upper = adjs$Dist.Upper,
                  edr = 0,
                  use_pois = 0
                  )

fit <- model$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 500)



quant <- function(x,q){
  y <- quantile(x,q,names = FALSE,digits = 3)
  return(y)
}

summ <- fit$summary()

cali <- fit$summary(variable = "calibration",
                    "mean",
                    "median",
                    "sd",
                    "rhat",
                    "ess_bulk",
                    q2_5 = ~quant(.x,q = 0.025),
                    q10 = ~quant(.x,q = 0.10),
                    q90 = ~quant(.x,q = 0.90),
                    q97_5 = ~quant(.x,q = 0.975))

cali_post <- fit$draws(variables = "calibration",
                       format = "df")

# BETA_post <- fit$draws(variables = "BETA",
#                        format = "df")

saveRDS(fit,paste0("output/calibration_fit_",sp_aou,"_",sp_ebird,".rds"))


sp_ebird <- ebirdst::get_species(sp_sel)
#
# down <- try(ebirdst::ebirdst_download_status(sp_ebird,
#                                              download_ranges = FALSE,
#                                              download_abundance = TRUE,
#                                              download_occurrence = FALSE),
#             silent = TRUE)
#

abd_seasonal_abundance <- ebirdst::load_raster(species = sp_ebird,
                                               resolution = "3km",
                                               period = "seasonal",
                                               product = "abundance")  #27km low resolution
season_sel <- ifelse(any(grepl("resident",names(abd_seasonal_abundance))),"resident","breeding")

# Need to explore the resident vs breeding season abundance patterns
# At least for BCCH the resident relative abundance is not linearly related
# to the BBS mean observations (log-linear)
# Option: if resident, then extract the weekly abundance data for the
# BBS survey season and generate a breeding season abundance

breed_abundance <- abd_seasonal_abundance[[season_sel]]

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
                           sum_abund = 3^2 * abundance_in_strata$breeding) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, fun = "mean"),
         pop_median = post_abund(sum_abund, fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,fun = "quantile",p = 0.975),
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


strata_compare <- strata_abund %>%
  full_join(.,strata_trad,
            by = c("prov_state" = "st_abrev",
                   "bcr" = "BCR"))



comp_trad_new_plot <- ggplot(data = strata_compare,
                             aes(x = mean.PopEst,
                                 y = pop_median))+
  geom_abline(slope = 1, intercept = 0)+
  geom_errorbar(aes(ymin = pop_lci_80,ymax = pop_uci_80),
                alpha = 0.3, width = 0)+
  geom_errorbarh(aes(xmin = LCI80.PopEst,xmax = UCI80.PopEst),
                alpha = 0.3)+
  geom_point()+
  geom_text_repel(aes(label = strata_name),
                  size = 3)+
  xlab("Traditional PIF population estimate")+
  ylab("PIF-Calibrated eBird relative abundance estimate")+
  labs(title = paste(sp_sel,"population estimates by BBS strata"),
       subtitle = "Diagonal line = 1:1")

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

abund_map <- ggplot()+
  geom_sf(data = strata, fill = NA)+
  geom_spatraster(data = breed_abundance)+
  geom_sf(data = strata, fill = NA)+
  geom_sf(data = routes_buf,fill = "black",colour = NA)+
  coord_sf(xlim = c(bb[c("xmin","xmax")]),
           ylim = c(bb[c("ymin","ymax")]))+
  scale_fill_viridis_c(direction = -1,
                       option = "G",
                       na.value = NA)+
  labs(title = paste(sp_sel,"eBird relative abundance values by BBS strata"))


abund_map


pdf(paste0("figures/trad_calib_compare_by_strata_",sp_aou,"_",sp_ebird,".pdf"),
    width = 11,
    height = 8.5)
print(comp_trad_new_plot)
print(abund_map)
dev.off()








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
                           sum_abund = 3^2 * abundance_in_countries$breeding) %>%
  mutate(pop_mean = post_abund(sum_abund, fun = "mean"),
         pop_median = post_abund(sum_abund, fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

country_abund




# USCAN estimates -----------------------------------------


USCAN <- bbsBayes2::load_map("bbs_usgs") %>%
  ungroup() %>%
  summarise()


USCAN_proj <- st_transform(USCAN,
                               crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_USCAN <- terra::extract(breed_abundance,
                                     USCAN_proj,
                                         fun = sum,
                                         na.rm = TRUE,
                                         ID = FALSE,
                                         exact = FALSE)



USCAN_abund <- data.frame(region = "USCAN",
                          region_type = "USCAN",
                            sum_abund = 3^2 * abundance_in_USCAN$breeding) %>%
  mutate(pop_mean = post_abund(sum_abund, fun = "mean"),
         pop_median = post_abund(sum_abund, fun = "quantile",p = 0.5),
         pop_lci_80 = post_abund(sum_abund,fun = "quantile",p = 0.1),
         pop_uci_80 = post_abund(sum_abund,fun = "quantile",p = 0.9),
         pop_lci_95 = post_abund(sum_abund,fun = "quantile",p = 0.025),
         pop_uci_95 = post_abund(sum_abund,fun = "quantile",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

USCAN_abund


pop_ests_out <- bind_rows(USCAN_abund,
                          country_abund,
                          strata_abund)
# Sum of population estimates

write_excel_csv(pop_ests_out,
                paste0("estimates/pop_ests_",sp_aou,sp_ebird,".csv"))



# } end of species loop
