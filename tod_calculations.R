#### Download and prep NA BBS data ####
# Required libraries
library(bbsBayes2) # devtools::install_github("bbsBayes/bbsBayes2")
library(mgcv)
library(tidyverse)
# Download NA BBS data up to 2023 using
# bbsBayes2 package; this only needs to be done once as the data will be stored
# and used for any future calls from bbsBayes functions
if(!bbsBayes2::have_bbs_data(level = "stop",quiet = TRUE)){
  bbsBayes2::fetch_bbs_data(level = "stop") # only downloads if doesn't already exist locally
}

# Load the full bbs dataset
all_dat <- bbsBayes2::load_bbs_data(level = "stop")

additional_species <- read_csv("data/Filtered_Spp_List_new_tod.csv") %>%
  filter(is.na(PIF_adj_birds_km2))

sps <- all_dat$species %>%
  filter(unid_combined == TRUE,
         english %in% additional_species$english) %>%
  select(-unid_combined)


#
# sps_lump <- bbsBayes2::species_forms %>%
#   filter(grepl("all forms",english_combined) | grepl("Redpoll",english_combined))


counts <- all_dat$birds %>%
  inner_join(sps,by = c("aou"))

surveys <- all_dat$routes


re_smooth = TRUE
if(re_smooth){
sp_out <- NULL

kk <- which(sps$english == "Wilson's Snipe")

for(i in (kk+1):nrow(sps)){

  sp <- unlist(sps[i,"english"])
  aout <- unlist(sps[i,"aou"])

  sp_dat <- counts %>%
    filter(aou == aout) %>%
    select(country_num,state_num,bcr,route,route_data_id,year,
           starts_with("stop_"))#%>%

  if(nrow(sp_dat) > 10000){
    rts_tmp <- surveys %>%
      filter(year > 2003) %>%
      filter(route_data_id %in% sp_dat$route_data_id) %>%
      select(country_num,state_num,route,bcr) %>%
      distinct()
  }else{
  rts_tmp <- surveys %>%
    filter(route_data_id %in% sp_dat$route_data_id) %>%
    select(country_num,state_num,route,bcr) %>%
    distinct()
}
  surveys_on_routes <- surveys %>%
    inner_join(rts_tmp,
               by = c("country_num","state_num","route","bcr")) %>%
    select(country_num,state_num,bcr,route,route_data_id,year)

 sp_dat <- sp_dat %>%
   right_join(surveys_on_routes,
              by = c("route_data_id", "country_num", "state_num", "bcr","route","year"))


 sp_long <- sp_dat %>%
   pivot_longer(cols = starts_with("stop_"),
                names_to = "stop",
                names_prefix = "stop_",
                values_to = "count") %>%
   mutate(count = ifelse(is.na(count),0,count),
          route_name = paste(state_num,bcr,route,sep = "-"),
          english = sp,
          stop = as.numeric(stop),
          year = as.factor(year),
          route_name = as.factor(route_name))  %>%
   ungroup()


 counts_route_year <- sp_long %>%
   filter(count > 0) %>%
   ungroup() %>%
   group_by(route_name,year) %>%
   summarise(n_obs = n(),
             max_count = max(count),
             median_count = median(count),
             .groups = "drop") %>%
   group_by(route_name) %>%
   summarise(n_yrs = n(),
           mean_n_obs = mean(n_obs),
           mean_max = mean(max_count),
           mean_median = mean(median_count))

 counts_route <- sp_long %>%
   group_by(route_name) %>%
   summarise(mean_count = mean(count),
             median_count = median(count),
             sd_count = sd(count),
             cv_count = sd_count/mean_count) %>%
   arrange(cv_count) %>%
   inner_join( counts_route_year,
                            by = "route_name") %>%
   arrange(desc(n_yrs),cv_count)

 n_routes <- nrow(counts_route)

 # keep the 400 routes on which the species has been most frequently observed (greatest number of years)
 rt_keep <- counts_route %>%
   slice_head(n = 400)

 sp_long <- sp_long %>%
   filter(route_name %in% rt_keep$route_name)

 rt_example <- counts_route %>%
   slice_head(n = 1) %>%
   select(route_name)

 nwdat <- sp_long %>%
   slice_sample(n = 1, by = c(stop)) %>%
   select(stop) %>%
   distinct() %>%
   arrange(stop) %>%
   mutate(route_name = rt_example$route_name,
          year = levels(sp_long$year)[length(levels(sp_long$year))])

 cnt_means <- sp_long %>%
   group_by(stop) %>%
   summarise(mean_count = mean(count)) %>%
   mutate(mean_count_stand = mean_count/mean(mean_count))

 nwdat <- nwdat %>%
   left_join(cnt_means, by = "stop")
 if(nrow(rt_keep) < 100){
   nwdat$english <- sp
   nwdat$aou <- aout
   nwdat$n_routes <- n_routes
   next
   }else{
 thetas <- vector("numeric",length = 5)
 for(phi in 1:5){
   set.seed(phi)
   rts_sel <- sp_long %>%
     select(route_name) %>%
     distinct() %>%
     slice_sample(n = 30)
   sp_tmp <- sp_long %>%
     filter(route_name %in% rts_sel$route_name)


 m1 <- gam(formula = count ~ s(stop) + s(route_name,bs="re") + s(year,bs="re"),
                 family = nb(),
           data = sp_tmp)

 thetas[phi] <- m1$family$getTheta(TRUE)

 }

 m1 <- bam(formula = count ~ s(stop) + s(route_name,bs="re") + s(year,bs="re"),
           family = negbin(theta = mean(thetas,na.rm = TRUE)),
           data = sp_long)
 }#else{
#   m1 <- gam(formula = count ~ s(stop) + s(route_name,bs="re") + s(year,bs="re"),
#             family = nb(),
#             data = sp_long)
# }
 pred <- mgcv::predict.gam(m1, se.fit = TRUE,
                           newdata = nwdat,
                           type = "terms")

 nwdat$smooth <- pred$fit[,1]
 nwdat$smooth_sd <- pred$se.fit[,1]
 nwdat$english <- sp
 nwdat$aou <- aout
 nwdat$n_routes <- n_routes
 nwdat$theta <- mean(thetas,na.rm = TRUE)






 sp_out <- bind_rows(sp_out,nwdat)
print(paste(sp,round(i/nrow(sps),3)))
rm(list = c("m1","sp_long","sp_dat","pred","nwdat"))
saveRDS(sp_out,"data/tod_obs.rds")
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
            peak_stop = extr_col(smooth,stop),
            center_check = mean(smooth))

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
               aes(x = stop,y = smooth))+
    geom_point(colour = "blue4",
               alpha = 0.7)+
    geom_ribbon(aes(ymin = smooth_lci,
                    ymax = smooth_uci),
                alpha = 0.3)+
    geom_line(aes(x = stop,y = smooth))+
    geom_line(data = tmp2, aes(x = peak_stop,y = peak))+
    labs(title = paste(sp)) + #,unique(tmp$n_surveys),"surveys (route*year) since 1997"))+
    ylab("Standardized mean count at each stop
         relative to mean across all stops")
print(pl)

}
dev.off()


tod_adj_out <- pk_sp %>%
  select(english,
         peak) %>%
  rename(cn = english,
         TimeAdj.meanlog = peak)

write_csv(tod_adj_out,
          "Additional_time_of_day_adjustments.csv")


