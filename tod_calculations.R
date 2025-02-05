#### Download and prep NA BBS data ####
# Required libraries
library(bbsBayes2) # devtools::install_github("bbsBayes/bbsBayes2")

library(tidyverse)
# Download NA BBS data up to 2023 using
# bbsBayes2 package; this only needs to be done once as the data will be stored
# and used for any future calls from bbsBayes functions
if(!bbsBayes2::have_bbs_data(level = "stop",quiet = TRUE)){
  bbsBayes2::fetch_bbs_data(level = "stop") # only downloads if doesn't already exist locally
}

# Load the full bbs dataset
all_dat <- bbsBayes2::load_bbs_data(level = "stop")


sps <- all_dat$species %>%
  filter(unid_combined == TRUE) %>%
  select(-unid_combined)


#
# sps_lump <- bbsBayes2::species_forms %>%
#   filter(grepl("all forms",english_combined) | grepl("Redpoll",english_combined))


counts <- all_dat$birds %>%
  inner_join(sps,by = c("aou"))

surveys <- all_dat$routes

sp_out <- NULL


for(i in 1:nrow(sps)){

  sp <- unlist(sps[i,"english"])
  aou <- unlist(sps[i,"aou"])

  sp_dat1 <- counts %>%
    filter(aou == aou) %>%
    select(country_num,state_num,route,route_data_id,year,
           starts_with("stop_"))#%>%

  rts_tmp <- surveys %>%
    filter(route_data_id %in% sp_dat1$route_data_id) %>%
    select(country_num,state_num,route) %>%
    distinct()

  surveys_on_routes <- surveys %>%
    inner_join(rts_tmp,
               by = c("country_num","state_num","route")) %>%
    select(country_num,state_num,route,route_data_id,year)

 sp_dat <- sp_dat1 %>%
   right_join(surveys_on_routes,
              by = c("route_data_id", "country_num", "state_num", "route","year"))

 sp_long <- sp_dat %>%
   pivot_longer(cols = starts_with("stop_"),
                names_to = "stop",
                names_prefix = "stop_",
                values_to = "count") %>%
   mutate(count = ifelse(is.na(count),0,count)) %>%
   group_by(stop) %>%
   summarise(mean_obs = mean(count),
             n_surveys = n(),
             .groups = "drop") %>%
   mutate(english = sp,
          stop = as.integer(stop),
          stand_mean_obs = mean_obs/mean(mean_obs)) %>%
   arrange(stop)

 m1 <- mgcv::gam(data = sp_long,
                 formula = stand_mean_obs~s(stop))
 pred <- mgcv::predict.gam(m1, se.fit = TRUE)

 sp_long$smooth = pred$fit
 sp_long$smooth_sd = pred$se.fit

 sp_out <- bind_rows(sp_out,sp_long)
print(paste(sp,round(i/nrow(sps),3)))
}

saveRDS(sp_out,"data/tod_obs.rds")


# Species-specific rescale then smooth -----------------------------------

tmp2 <- tmp %>%
  group_by(english) %>%
  summarise(mean_c = mean(mean_obs),
            min_c = min(mean_obs),
            max_c = max(mean_obs))

smooth_out <- NULL


for(i in 1:nrow(sps)){

  sp <- unlist(sps[i,"english"])

  sp_dat1 <- counts %>%
    filter(english == sp) %>%
    select(country_num,state_num,route,route_data_id,year,
           starts_with("stop_"))#%>%

  rts_tmp <- surveys %>%
    filter(route_data_id %in% sp_dat$route_data_id) %>%
    select(country_num,state_num,route) %>%
    distinct()

  surveys_on_routes <- surveys %>%
    inner_join(rts_tmp,
               by = c("country_num","state_num","route")) %>%
    select(country_num,state_num,route,route_data_id,year)

  sp_dat <- sp_dat1 %>%
    right_join(surveys_on_routes,
               by = c("route_data_id", "country_num", "state_num", "route","year"))

  sp_long <- sp_dat %>%
    pivot_longer(cols = starts_with("stop_"),
                 names_to = "stop",
                 names_prefix = "stop_",
                 values_to = "count") %>%
    mutate(count = ifelse(is.na(count),0,count)) %>%
    group_by(stop) %>%
    summarise(mean_obs = mean(count),
              n_surveys = n(),
              .groups = "drop") %>%
    mutate(english = sp,
           stop = as.integer(stop))

  sp_out <- bind_rows(sp_out,sp_long)
}

saveRDS(sp_out,"data/tod_obs.rds")




