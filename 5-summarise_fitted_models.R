
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
library(HDInterval)
library(ggtext)


# Script to summarise the fitted models -----------------------------------

## generates various graphs, summaries, and comparisons for each species


source("functions/hpdi.R")
source("functions/posterior_abundance.R")


if(Sys.info()[["nodename"]] == "WNCRLABN72960"){
setwd("c:/Users/SmithAC/Documents/GitHub/PIF_population_estimates")
}

yr_ebird <- 2023 # prediction year for eBird relative abundance


# mean counts to id routes on which species has been observed -------------
mean_counts_all <- readRDS("data/mean_counts_by_route.rds")

raw_counts_all <- readRDS("data/all_counts_by_route.rds")





# Add na-pops EDRs --------------------------------------------------------


ExpAdjs <- read_csv("Species_correction_factors_w_edr_availability.csv")



routes_buf_all <- readRDS("data/all_routes_buffered.rds")

strata <- bbsBayes2::load_map("bbs_usgs")

strata_names <- strata %>%
  sf::st_drop_geometry()





# load traditional estimates with EDR for later comparison -------------------------
trad_dir <- "Stanton_2019_code/output"
strata_join <- strata_names %>%
  select(strata_name,country_code,
         prov_state, bcr) %>%
  st_drop_geometry() %>%
  distinct()


#
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
  select(region,region_type,version,st_abrev,BCR,
         species,aou,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95)


saveRDS(pop_ests_out_trad,"all_traditional_pop_estimates.rds")


pop_ests_out_trad <- readRDS("all_traditional_pop_estimates.rds")



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
bcrs <- read_sf("data/BCR_terrestrial/BCR_Terrestrial_master_International.shp") %>%
  group_by(BCR,Label) %>%
  summarise(do_union = FALSE)







use_traditional <- FALSE # if TRUE uses the traditional adjustments for TOD and distance

estimate_rho <- FALSE





if(use_traditional){
  vers <- ifelse(estimate_rho,"trad_rho","trad")
}else{
  vers <- ifelse(estimate_rho,"_rho","")
}


output_dir <- "F:/PIF_pop_estimates/output"
#output_dir <- "output"
trim_rel_abund <- TRUE


# species loop ------------------------------------------------------------

selected_species <- readRDS("data/selected_species.rds")

for(sp_sel in selected_species){

 sp_aou <- bbsBayes2::search_species(sp_sel)$aou[1]
  sp_ebird <- ebirdst::get_species(sp_sel)


  if(!file.exists(paste0(output_dir,"/calibration_fit_alt_",
                        vers,sp_aou,"_",sp_ebird,".rds"))){next}



strata_trad <- pop_ests_out_trad %>%
  filter(species == sp_sel,
         region_type == "strata")

global_trad <- pop_ests_out_trad %>%
  filter(species == sp_sel,
         region_type == "USACAN")

bcr_trad <- pop_ests_out_trad %>%
  filter(species == sp_sel,
         region_type == "bcr")







if(!file.exists(paste0("data/species_relative_abundance_2023/",
                       sp_ebird,"_relative_abundance.rds"))){next}







# load saved fitted model, summary, and data ------------------------------


  fit <- readRDS(paste0(output_dir,"/calibration_fit_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
  summ <- readRDS(paste0("convergence/parameter_summary_alt_",vers,sp_aou,"_",sp_ebird,".rds"))

  combined <- readRDS(paste0("stan_data/dataframe_",vers,sp_aou,"_",sp_ebird,".rds"))
  stan_data <- readRDS(paste0("stan_data/stan_data_",vers,sp_aou,"_",sp_ebird,".rds"))



  routes_buf <- readRDS(paste0("stan_data/routes_buffered_",vers,sp_aou,"_",sp_ebird,".rds"))


# posterior of calibration ------------------------------------------------

# load the hyperparameter estimates,
# but only used to provide structure for the trimmed mean
# calculation below, not used in the final calibration
cali_alt_post <- fit$draws(variables = "calibration_mean",
                           format = "df")


# trimmed mean function to exclude the outer 5% of routes with extreme
# calibration values. These extreme values have a strong influence on the
# mean of the retransformed values. particularly routes with very high calibration
# values which greatly increase the overall mean calibration and the total
# population estimate.
# Ideally, we'll figure out what is missing from the model and why
# we're getting these heavy tails, but until then this trimmed-mean generates
# much more robust population estimates.

  # ID routes with extreme betas -----------------------------------------------

  betas <- summ %>% filter(grepl("beta[",variable,fixed = TRUE))

  # tmp <- summ %>% filter(grepl("ct",variable,fixed = TRUE))
  # tmp <- summ %>% filter(variable == "cd")
  #
  w_excl <- c(which(betas$mean > quantile(betas$mean,0.975)),
              which(betas$mean < quantile(betas$mean,0.025)))

  betas_keep <- c(1:stan_data$n_routes)[-w_excl]

#posterior of the route-level calibration values
cali_route_post <- fit$draws(variables = "calibration_r",
                           format = "df")

# set up structure for the trimmed mean
cali_trim_post <- cali_alt_post

# replace values of the hyperpareter with the trimmed mean values
for(i in 1:nrow(cali_route_post)){
  cali_trim_post[i,1] <- mean(as.numeric(cali_route_post[i,betas_keep]))
}


saveRDS(cali_trim_post,paste0("estimates/calibration_posterior_",vers,sp_aou,"_",sp_ebird,".rds"))



# calculating population --------------------------------------------------


# load relative abundance surface -----------------------------------------


# this is the range-wide (global) relative abundance surface for the
# selected set of weeks that overlap the BBS field season
breed_abundance <- readRDS(paste0("data/species_relative_abundance_2023/",
                                  sp_ebird,
                                  "_derived_breeding_relative_abundance.rds"))

names(breed_abundance) <- "breeding_abundance"

breed_abundance_full <- breed_abundance



# BBS strata level estimates (to match PIF stratification) ----------------


## process is to summarise the relative abundance within given regions
## Then multiply those sums by the full posterior of the calibration
## to generate estimates with uncertainty
##
## One could apply the full posterior of the calibration to the
## pixel values (generating a full posterior of the mapped surface)
## one could also combine the uncertainty of the calibration with the
## uncertainty of the pixel=level relative abundances, but both of those
## approaches would have extreme computational costs.
##
## The first approach would be mathematically equivalent to what's done here
## because the uncertainty at a pixel level would be the same across each posterior
## draw i.e., sum(calibration * pixel-values) = calibration*sum(pixel-values)
##
## The second approach would account for the uncertainty in the relative abundance
## surface and is worth exploring, with full consideration of the spatial autocorrelation
## in the eBird uncertainty estimates (e.g., each pixel's uncertainty is not independent,
##  but how dependent is it?)
##

strata_proj <- st_transform(strata,
                            crs = st_crs(breed_abundance))

## extract the sum of the relative abundance surface in each region
abundance_in_strata <- terra::extract(breed_abundance_full,
                                      strata_proj,
                                       fun = sum,
                                       na.rm = TRUE,
                                       ID = TRUE,
                                       exact = FALSE)


###
use_trimmed_calibration <- TRUE
# select calibration ------------------------------------------------------

if(use_trimmed_calibration){
  cali_use <- cali_trim_post
}else{
  cali_use <- cali_alt_post
}
names(cali_use)[1] <- "calibration"

saveRDS(cali_use,paste0("output/calibration_",vers,sp_aou,"_",sp_ebird,".rds"))

### 3^2 scaling is to account for the area of each grid-cell 9-km^2
###

strata_abund <- data.frame(region = strata_proj$strata_name,
                           region_type = "strata",
                           sum_abund = 3^2 * unlist(abundance_in_strata[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                                 fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird,
         strata_name = region) %>%
  left_join(.,strata_names)
  #filter(pop_mean != 0)
strata_sampled <- readRDS(paste0("data/species_relative_abundance_2023/",sp_ebird,"_bbs_strata_relative_abundance.rds")) %>%
  select(strata_name,mean_ebird_abundance:log_ratio)

strata_abund <- strata_abund %>%
  left_join(strata_sampled,
            by = "strata_name")




# Annual strata-level abundances ------------------------------------------
# population size trajectories for each stratum with BBS data
# Excludes strata with population but no BBS time-series estimates
strats <- combined %>%
  select(strata_name,strata) %>%
  distinct() %>%
  arrange(strata)

n_yr <- max(combined$yr)
n_strat <- max(combined$strata)
fyr <- min(combined$year)


yrs <- data.frame(year = fyr:(fyr+n_yr-1),
                  yr = 1:n_yr)

time_series <- expand_grid(yrs,strats)


year_draws <- gather_draws(fit, yeareffect[strata,yr]) %>%
  inner_join(strats)

cali_use_y <- expand_grid(cali_use,yrs) %>%
    inner_join(year_draws,
               by = c("yr",".draw",".iteration",".chain"))






  strata_population_posterior <- data.frame(region = strata_proj$strata_name,
                             region_type = "strata",
                             sum_abund = 3^2 * unlist(abundance_in_strata[[2]])) %>%
    filter(!is.na(sum_abund)) %>%
    expand_grid(yrs) %>%
    inner_join(cali_use_y,
               by = c("region" = "strata_name",
                      "yr","year")) %>%
    mutate(population = sum_abund*calibration*exp(.value),
           species = sp_sel,
           species_ebird = sp_ebird,
           strata_name = region) %>%
    left_join(.,strata_names,
              by = "strata_name")




  saveRDS(strata_population_posterior,
          paste0("output/strata_population_posterior_",vers,sp_aou,"_",sp_ebird,".rds"))




  strata_trajectories <- strata_population_posterior %>%
    group_by(region, region_type, year) %>%
  summarise(pop_mean = mean(population),
            pop_median = median(population),
            pop_lci_80 = interval_function_hpdi(population,probs = 0.1),
            pop_uci_80 = interval_function_hpdi(population,probs = 0.9),
            pop_lci_95 = interval_function_hpdi(population,probs = 0.025),
            pop_uci_95 = interval_function_hpdi(population,probs = 0.975))



  traj_plot <- ggplot(data = strata_trajectories,
                 aes(x = year,y = pop_median))+
    geom_ribbon(aes(ymin = pop_lci_80,ymax = pop_uci_80),
                alpha = 0.3)+
    geom_line()+
    facet_wrap(vars(region),
               scales = "free")

  traj_plot

  saveRDS(traj_plot,paste0("figures/saved_ggplots/abund_trajectory_strata_",vers,sp_aou,"_",sp_ebird,".rds"))

  USA_CAN_trajectories <- strata_population_posterior %>%
    group_by(year, .draw) %>%
    summarise(population = sum(population)) %>%
    group_by(year) %>%
    summarise(pop_mean = mean(population),
              pop_median = median(population),
              pop_lci_80 = interval_function_hpdi(population,probs = 0.1),
              pop_uci_80 = interval_function_hpdi(population,probs = 0.9),
              pop_lci_95 = interval_function_hpdi(population,probs = 0.025),
              pop_uci_95 = interval_function_hpdi(population,probs = 0.975))



  traj_plot <- ggplot(data = USA_CAN_trajectories,
                      aes(x = year,y = pop_median))+
    geom_ribbon(aes(ymin = pop_lci_80,ymax = pop_uci_80),
                alpha = 0.3)+
    geom_line()+
    scale_y_continuous(limits = c(0,NA),
                       labels = scales::unit_format(unit = "M", scale = 1e-6))+
    scale_x_continuous(breaks = seq(2013,2023,2))+
    xlab("")+
    ylab("Annual population size\nUnited States of America and Canada")+
    labs(title = sp_sel)+
    theme_bw()


  traj_plot

saveRDS(traj_plot,paste0("figures/saved_ggplots/abund_trajectory_",vers,sp_aou,"_",sp_ebird,".rds"))
  # compare to traditional estimates ----------------------------------------

strata_trad_j <- strata_trad %>%
  mutate(med.PopEst = pop_median,
         LCI80.PopEst = pop_lci_80,
         UCI80.PopEst = pop_uci_80) %>%
  select(st_abrev,
         BCR,
         med.PopEst,
         LCI80.PopEst,
         UCI80.PopEst)

pal_labs <- c("< -1",
              "-1, -0.5",
              "-0.5, -0.15",
              "-0.15, 0.15",
              "0.15, 0.5",
              "0.5, 1",
              "> 1")

strata_compare <- strata_abund %>%
  full_join(.,strata_trad_j,
            by = c("prov_state" = "st_abrev",
                   "bcr" = "BCR")) %>%
  mutate(ratio_cat = cut(log_ratio,breaks = c(-Inf,-1,-0.5,-0.15,0.15,0.5,1,Inf),
                         levels = c("(-Inf,-1]","(-1,-0.5]",
                         "(-0.5,-0.15]","(-0.15,0.15]",
                         "(0.15,0.5]","(0.5,1]","(1, Inf]"),
                         labels = pal_labs),
         log_pop_ratio = log(pop_median/med.PopEst))

div_pal <- viridisLite::viridis(n = 7,option = "turbo", direction = -1)

names(div_pal) <- pal_labs#levels(strata_compare$ratio_cat)


comp_trad_new_plot <- ggplot(data = strata_compare,
                             aes(x = med.PopEst,
                                 y = pop_median))+
  # geom_abline(slope = 1, intercept = 0)+
  # geom_abline(slope = 2, intercept = 0,linetype = 2)+
  # geom_abline(slope = 5, intercept = 0,linetype = 3)+
  # geom_abline(slope = 0.5, intercept = 0,linetype = 2)+
  # geom_abline(slope = 0.2, intercept = 0,linetype = 3)+
  geom_errorbar(aes(ymin = pop_lci_80,ymax = pop_uci_80),
                alpha = 0.3, width = 0)+
  geom_errorbarh(aes(xmin = LCI80.PopEst,xmax = UCI80.PopEst),
                alpha = 0.3)+
  geom_point(aes(colour = ratio_cat))+
  geom_smooth(method = "lm", se = FALSE,
              colour = grey(0.5))+
  geom_text_repel(aes(label = strata_name),
                  size = 2)+
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6))+
  scale_colour_manual(values = div_pal,
                      name = "BBS sampling bias \n log(sampled/avail)",
                      na.translate = FALSE)+
  # scale_colour_viridis_c(option = "turbo", direction = -1,
  #                        name = "BBS sampling bias \n log(sampled/avail)")+
  xlab("Existing PIF population estimate")+
  ylab("Revised PIF population estimate")+
  labs(title = paste(sp_sel))+
  theme_bw()

png(filename = paste0("Figures/comp_trad_new_alt_",vers,sp_aou,"_",sp_ebird,".png"),
    res = 300,
    height = 6,
    width = 6,
    units = "in")
print(comp_trad_new_plot)
dev.off()

pdf(paste0("Final_figures/comp_trad_new_alt_",vers,sp_aou,"_",sp_ebird,".pdf"),
    height = 4,
    width = 7)
print(comp_trad_new_plot)
dev.off()



# National summaries ------------------------------------------------------





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
                                      ID = TRUE,
                                      exact = FALSE)

abundance_in_countries <- abundance_in_countries %>%
  filter(!is.na(breeding_abundance),
         breeding_abundance > 0)



country_abund <- data.frame(region = countries_proj$country_name[abundance_in_countries$ID],
                            region_type = "country",
                           sum_abund = 3^2 * unlist(abundance_in_countries[[2]])) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#country_abund # will be combined with other regional estimates below



# BCR summaries -----------------------------------------------------------






bcr_proj <- st_transform(bcrs,
                               crs = st_crs(breed_abundance))


abundance_in_bcr <- terra::extract(breed_abundance,
                                   bcr_proj,
                                         fun = sum,
                                         na.rm = TRUE,
                                         ID = TRUE,
                                         exact = FALSE)

abundance_in_bcr <- abundance_in_bcr %>%
  filter(!is.na(breeding_abundance),
         breeding_abundance > 0)



bcr_abund <- data.frame(region = as.character(bcr_proj$BCR[abundance_in_bcr$ID]),
                            region_type = "bcr",
                            sum_abund = 3^2 * unlist(abundance_in_bcr[[2]])) %>%
  filter(!is.na(sum_abund)) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#bcr_abund
# will be combined with other regional estimates below


# mapping with BCR extent -------------------------------------------------



re_do_map <- FALSE #option to skip because they take a long time.
if(re_do_map){
  # abund_mapable <- breed_abundance %>%
  #   terra::project(.,"EPSG:9822")

  bcrs_w_pop <- bcr_abund %>%
    filter(pop_median > 0)

  bcrs_plot <- bcrs %>%
    st_transform(crs = st_crs(strata)) %>%
    filter(BCR %in% bcrs_w_pop$region)

countries_plot <- countries %>%
  st_transform(crs = st_crs(strata)) %>%
  filter(CC == "NA")

strata_w_pop <- strata_abund %>%
  filter(sum_abund > 0)

strata_bound <- strata %>%
  filter(strata_name %in% strata_w_pop$region)

bb <- sf::st_bbox(strata_bound)

  breed_abundance_plot <- breed_abundance
  names(breed_abundance_plot) <- "breeding"

  breed_abundance_plot <- breed_abundance_plot %>%
    #st_transform(crs = st_crs(strata)) %>%
    mutate(.,breeding = ifelse(breeding == 0,NA,breeding),
           breeding = (breeding*mean(cali_use$calibration))/9) # per km2 mean density


  # Mapping -----------------------------------------------------------------


  abund_map <- ggplot()+
    geom_sf(data = countries_plot, fill = NA)+
    geom_spatraster(data = breed_abundance_plot,
                    maxcell = 16000000)+
    geom_sf(data = bcrs_plot, fill = NA)+
    # geom_sf_text(data = bcrs_plot,aes(label = BCR),
    #              size = 3)+
    #geom_sf(data = routes_buf,aes(colour = mean_count),fill = NA)+
    coord_sf(xlim = c(bb[c("xmin","xmax")]),
             ylim = c(bb[c("ymin","ymax")]))+
    #scale_fill_gradientn(12,colours = terrain.colors(12),na.value = NA)+
    scale_fill_viridis_c(direction = -1,
                         option = "G",
                         na.value = NA,
                         end = 0.9,
                         name = paste0(sp_sel,"\n","Birds/km<sup>2</sup>"))+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(legend.title = element_markdown())+
    labs(title = paste(sp_sel))


  #

  # png(filename = paste0("Figures/abund_map_alt_",vers,sp_aou,"_",sp_ebird,".png"),
  #     res = 400,
  #     height = 7,
  #     width = 6.5,
  #     units = "in")
  # print(abund_map)
  # dev.off()
  save(list = c("bb",
                "breed_abundance_plot",
                "countries_plot",
                "bcrs_plot"),
       file = paste0("figures/saved_ggplots/abund_map_alt_",vers,sp_aou,"_",sp_ebird,".rdata"))

}

saveRDS(comp_trad_new_plot,paste0("figures/saved_ggplots/trad_vs_new_alt_",vers,sp_aou,"_",sp_ebird,".rds"))
saveRDS(strata_compare,paste0("estimates/strata_comparison_",vers,sp_aou,"_",sp_ebird,".rds"))




# USACAN estimates -----------------------------------------



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
                                         ID = TRUE,
                                         exact = FALSE)



USACAN_abund <- data.frame(region = "USACAN",
                          region_type = "USACAN",
                            sum_abund = 3^2 * unlist(abundance_in_USACAN[[2]])) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#USACAN_abund
# will be combined with other regional estimates below

# Continent estimates -----------------------------------------



continents_proj <- st_transform(continents,
                            crs = st_crs(breed_abundance))

# ## reprojecting the abundance surface
# breed_sel1 <- project(breed_abundance, crs(strata), method = "near") |>
#   # remove areas of the raster containing no data
#   trim()

abundance_in_continents <- terra::extract(breed_abundance,
                                      continents_proj,
                                      fun = sum,
                                      na.rm = TRUE,
                                      ID = TRUE,
                                      exact = FALSE)

abundance_in_continents <- abundance_in_continents %>%
  filter(!is.na(breeding_abundance),
         breeding_abundance > 0)




continents_abund <- data.frame(region = as.character(continents_proj$CC[abundance_in_continents$ID]),
                           region_type = "continent",
                           sum_abund = 3^2 * unlist(abundance_in_continents[[2]])) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#continents_abund
# will be combined with other regional estimates below



abundance_in_global <- sum(values(breed_abundance_full),na.rm = TRUE)



global_abund <- data.frame(region = "global",
                           region_type = "global",
                           sum_abund = 3^2 * abundance_in_global) %>%
  mutate(pop_mean = post_abund(sum_abund, use_log = FALSE,
                               draws = cali_use$calibration,
                               fun = "mean"),
         pop_median = post_abund(sum_abund, use_log = FALSE, draws = cali_use$calibration,
                                 fun = "median",p = 0.5),
         pop_lci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.1),
         pop_uci_80 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.9),
         pop_lci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.025),
         pop_uci_95 = post_abund(sum_abund, use_log = FALSE,draws = cali_use$calibration,
                                 fun = "hpdi",p = 0.975),
         species = sp_sel,
         species_ebird = sp_ebird)

#global_abund
# will be combined with other regional estimates in next line


pop_ests_out <- bind_rows(USACAN_abund,
                          continents_abund,
                          global_abund,
                          country_abund,
                          bcr_abund,
                          strata_abund) %>%
  mutate(version = "PIF_eBird")

# Sum of population estimates

write_excel_csv(pop_ests_out,
                paste0("estimates/pop_ests_alt_",vers,sp_aou,sp_ebird,".csv"))


pop_ests_out_trad_sel <- pop_ests_out_trad %>%
  filter(species == sp_sel)


pop_compare_stack <- pop_ests_out %>%
  select(region,region_type,version,
         species,species_ebird,
         pop_median, pop_lci_80, pop_uci_80,
         pop_lci_95, pop_uci_95) %>%
  bind_rows(pop_ests_out_trad_sel)


pop_compare_stack_sel <- pop_compare_stack %>%
  filter(region_type %in% c("USACAN","continent","global","bcr"))


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


saveRDS(side_plot,paste0("figures/saved_ggplots/side_plot_alt_",vers,sp_aou,"_",sp_ebird,".rds"))


# Exploring the estimated seasonal patterns by strata ------------------------------------------
strat_names <- combined %>%
  select(strata,strata_name,day_of_year) %>%
  distinct()

seasonal_strat <- summ %>% filter(grepl("doy_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = rep(c(1:stan_data$n_doy),times = stan_data$n_strata) + (min(combined$day_of_year)-1),
         strata = rep(c(1:stan_data$n_strata),each = stan_data$n_doy)) %>%
  inner_join(strat_names,by = c("strata",
                                "doy" = "day_of_year"))

seasonal <- summ %>% filter(grepl("DOY_pred[",variable,fixed = TRUE)) %>%
  mutate(doy = c(1:stan_data$n_doy)+ (min(combined$day_of_year)-1),
         strata_name = "USA_CAN") %>%
  filter(doy %in% strat_names$day_of_year)

strat_labs1 <- seasonal_strat %>%
  group_by(strata_name) %>%
  summarise(doy = max(doy))

strat_labs <- strat_labs1 %>%
  inner_join(seasonal_strat,
             by = c("strata_name",
                    "doy"))

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
                     aes(x = doy,y = mean,
                         group = strata_name,
                         colour = strata_name))+
  geom_ribbon(data = seasonal,aes(x = doy,y = mean,
                                  ymin = q5, ymax = q95),
              alpha = 0.1,
              inherit.aes = FALSE)+
  geom_line(data = seasonal,aes(x = doy,y = mean),
            inherit.aes = FALSE)+
  geom_line(alpha = 0.5)+
  coord_cartesian(xlim = c(100,240))+
  ggrepel::geom_text_repel(data = strat_labs, aes(label = strata_name),
                           size = 1.5,
                           min.segment.length = 0,
                           xlim = c(200,240),
                           nudge_x = 20)+
  scale_colour_viridis_d()+
  scale_x_continuous(breaks = xlbl$doy, labels = xlbl$days)+
  xlab("")+
  ylab("Seasonal effect (log-scale)")+
  scale_y_continuous(transform = "exp")+
  theme_bw()+
  theme(legend.position = "none")


saveRDS(seasonal_strat,paste0("output/seasonal_strat_",vers,sp_aou,"_",sp_ebird,".rds"))

saveRDS(vis_season,paste0("figures/saved_ggplots/season_",vers,sp_aou,"_",sp_ebird,".rds"))
png(filename = paste0("Figures/seasonal_effect_",vers,sp_aou,"_",sp_ebird,".png"),
    res = 300,
    height = 6,
    width = 6,
    units = "in")
print(vis_season)
dev.off()


}





