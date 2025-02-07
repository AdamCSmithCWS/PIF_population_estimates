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

re_smooth = FALSE
if(re_smooth){
sp_out <- NULL


for(i in 1:nrow(sps)){

  sp <- unlist(sps[i,"english"])
  aout <- unlist(sps[i,"aou"])

  sp_dat <- counts %>%
    filter(aou == aout) %>%
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

 sp_dat <- sp_dat %>%
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

 if(nrow(sp_long) < 49){next}
 m1 <- mgcv::gam(data = sp_long,
                 formula = stand_mean_obs~s(stop))
 pred <- mgcv::predict.gam(m1, se.fit = TRUE)

 sp_long$smooth = pred$fit
 sp_long$smooth_sd = pred$se.fit

 sp_out <- bind_rows(sp_out,sp_long)
print(paste(sp,round(i/nrow(sps),3)))
rm(list = c("m1","sp_long","sp_dat","pred"))
}
saveRDS(sp_out,"data/tod_obs.rds")
}else{
  sp_out <- readRDS("data/tod_obs.rds")
}



extr_col <- function(x,y){
  mx <- which.max(x)
  yx <- y[mx]
  return(yx)
}

pk_sp <- sp_out %>%
  dplyr::group_by(english) %>%
  dplyr::summarise(peak = max(smooth,na.rm = T),
            peak_sd = extr_col(smooth,smooth_sd),
            peak_stop = extr_col(smooth,stop)) %>%
  mutate(peak_tod_adj = 1/peak)

sp_out <- sp_out %>%
  mutate(smooth_lci = smooth-(2*smooth_sd),
         smooth_uci = smooth+(2*smooth_sd))

pdf("tod smooths all species.pdf")
for(i in 1:nrow(sps)){

  sp <- unlist(sps[i,"english"])
  aout <- unlist(sps[i,"aou"])

  tmp <- sp_out %>%
    filter(english == sp)
  if(nrow(tmp) == 0){next}
  tmp2 <- pk_sp %>%
    filter(english == sp)
  tmp2[2,"peak_stop"] <- tmp2[1,"peak_stop"]
  tmp2[2,"peak"] <- 0
  pl <- ggplot(data = tmp,
               aes(x = stop,y = stand_mean_obs))+
    geom_point(colour = "blue4",
               alpha = 0.7)+
    geom_ribbon(aes(ymin = smooth_lci,
                    ymax = smooth_uci),
                alpha = 0.3)+
    geom_line(aes(x = stop,y = smooth))+
    geom_line(data = tmp2, aes(x = peak_stop,y = peak))+
    labs(title = paste(sp,unique(tmp$n_surveys),"surveys (route*year) since 1997"))+
    ylab("Standardized mean count at each stop
         relative to mean across all stops")
print(pl)

}
dev.off()
