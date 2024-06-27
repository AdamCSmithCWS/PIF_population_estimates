##############################
# Calculated PIF Population Size Estimate
#
# Data assemble Data from
# North American Breeding Bird Survey
# (run 1_AssembleBBSdatabySpecies.R first)
#
# Script uses mean count by route for last 10 years of
# data supplied
#
# Currently only outputs North American landbird species
# (specified in Species_PSEst_Adjustments.csv file)
#
# by J.C. Stanton
#
######  USER INPUTS ########
#############################
n.samp <- 100

# Provide Path to Input file directory - Adjustment parameter files
# Inputfiles.dir <- readline("Provide directory path to location of INPUT files: ")
# if(!file.exists(Inputfiles.dir)) stop("Provide valid directory path to location of input files. ")
Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")

# Spp.dir <- readline("Provide directory path to location of Species files generated in step 1: ")
# if(!file.exists(Spp.dir)) stop("Provide valid directory path to location of species files. ")
#Spp.dir <- "Stanton_2019_code/input_files/species_data/"#paste(Spp.dir, '/', sep="")

# Provide Path to where model output should be written
# Output.dir <- readline("Provide directory path to store output files: ")
# if(!file.exists(Output.dir)) stop("Provide valid directory path to output files. ")
Output.dir <- "Stanton_2019_code/output/" #paste(Output.dir, '/', sep="")
###########################################################

##### dependent datasets  #####
## Parameters used in Population Size estimation
## See associated manuscript for details and source references

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE)
# Species Global and WH range adjustment factors
WH.GL.Adj <- read.csv(paste(Inputfiles.dir,'Species_Global_WH_Adjustments.csv', sep=''), stringsAsFactors=FALSE)
# Species by Polygon (ProvBCR) range adjustment factors
Rng.Adj  <- read.csv(paste(Inputfiles.dir, 'RangeAdjustment.csv', sep=''), stringsAsFactors=FALSE)
# Species by Polygon (ProvBCR) Population Estimates from other sources (non-BBS)
PS.othersources <- read.csv(paste(Inputfiles.dir,'PS_othersources.csv', sep=''), stringsAsFactors=FALSE)
# Species by Polygons that borrow Population density estimate from other Polygon
PolyRngAdj <- read.csv(paste(Inputfiles.dir, 'Species_byRange_Adjs_BorrowPoly.csv', sep=''), stringsAsFactors=FALSE)
# Area (in Km^2) of Polygons
PolyAreas <- read.csv(paste(Inputfiles.dir,'PolygonRangeAreas.csv', sep=''), stringsAsFactors=FALSE)
##############################################################

library(msm)
library(tidyverse)
library(plyr)

PrettyUp <- function(x){
  out <- x
  if(x>0) out <- round_any(x, 10^floor(log10(x)-1))
  return(out)
}

PS.score <- function(x){
  score <- NA
  if(!is.na(x)){
    if(x>=50000000) score = 1
    if(x<50000000 & x>=5000000) score = 2
    if(x<5000000 & x>=500000) score = 3
    if(x<500000 & x>=50000) score = 4
    if(x<50000) score = 5
  }
  return(score)
}

all_data <- readRDS("Stanton_2019_code/input_files/All_BBS_data_by_strata_name.rds")|>
  dplyr::filter(aou %in% ExpAdjs$Aou,
                !is.na(strata_name_numeric)) |>
  mutate(ProvBCR = strata_name_numeric)

region.lkup <- all_data |>
  dplyr::select(st_abrev,state_num) |>
  dplyr::distinct() |>
  dplyr::rename(Prov = state_num)


# tmp <- all_data |>
#   filter(is.na(ProvBCR))

# Spp.dat <- list.files(Spp.dir)
Spp.dat.date <- ""#unlist(strsplit(unlist(Spp.dat.date), "[.]"))[1]
Spp.aou <- unique(all_data$aou)

#Route Info
#routes <- read.csv(unz(paste(Inputfiles.dir, 'Routes.zip',sep=''), filename='routes.csv'))
# Letter to numeric code lookup for states/prov

##################################################################################

# build lookup table for Distance Adjustment sample parameters
Dist.adj = c(20, 50, 80, 100, 125, 200, 300, 400, 800, 850)
lower <- c(1)
upper <- c(1)
for(d in 2:length(Dist.adj)) {
  lower[d] = Dist.adj[d-1]+(Dist.adj[d]-Dist.adj[d-1])*.2
  upper[d] = Dist.adj[d]*1.1
}
Dist.adjparms <- data.frame(Dist.adj = Dist.adj[2:9], L.bound = lower[2:9], U.bound = upper[2:9])

# routes$rteno <- paste(routes$countrynum, routes$statenum, routes$Route, sep="-")
# routes$ProvBCR <- paste(routes$statenum, routes$BCR, sep='-')
# route.info <- routes[ ,c('rteno', 'ProvBCR', "countrynum",  "statenum", 'BCR', "Latitude", "Longitude")]
## above is already in the all_data object

missing <- c()

Spp.by.BCR.i <- list()
Spp.by.Prov.i <- list()
Spp.by.WH.GL.i <- list()
Spp.ProvBCR <- data.frame()

for(i in 1:length(Spp.aou)){ #Begin loop through each species
  spp.i <- Spp.aou[i]
  cn.i <- ExpAdjs$cn[as.character(ExpAdjs$Aou)==spp.i]
  df.spp.i <- all_data[which(all_data$aou == spp.i),] |>
    dplyr::group_by(ProvBCR,route_name) |>
    dplyr::summarise(mean.count = mean(species_total,na.rm = TRUE),
                     var.count = var(species_total,na.rm = TRUE),
                     n.years = n(),
                     size.par = (mean.count^2)/(var.count-mean.count),
                     .groups = "drop")

  if(nrow(df.spp.i)==0){
    missing<-c(missing, paste(cn.i, 'dataset', sep=' - '))
    next
  }   # go to next species if no data
  # merge location data for each route (state, provence, BCR, lat, long, ect...)
  # df.spp.i <- merge(x=df.spp.i, y=route.info, by='rteno', all.x=TRUE)
   #### Lookup Adjustment values based on Expert opinion/analysis
  Dist.adj.i <- ExpAdjs[ExpAdjs$Aou %in% spp.i, 'Dist2']
  Dist.adj.LB <- Dist.adjparms$L.bound[Dist.adjparms$Dist.adj==Dist.adj.i]
  Dist.adj.UB <- Dist.adjparms$U.bound[Dist.adjparms$Dist.adj==Dist.adj.i]
  Time.adj.meanlog <- ExpAdjs[ExpAdjs$Aou %in% spp.i, 'TimeAdj.meanlog']
  Time.adj.sdlog <- ExpAdjs[ExpAdjs$Aou %in% spp.i, 'TimeAdj.sdlog']
  Pair.adj.i <- ExpAdjs[ExpAdjs$Aou %in% spp.i, 'Pair2']

  regions <- unique(df.spp.i$ProvBCR)
  ##### Lookup any Range Adjustments for each region area
  # Any Range Adjustments?
  Rng.Adj.i <- Rng.Adj[Rng.Adj$cn==cn.i, ]
  # Any borrowed Polygons?
  PolyRngAdj.i <- PolyRngAdj[PolyRngAdj$cn==cn.i, ]
  if(nrow(PolyRngAdj.i)>0 & !all(PolyRngAdj.i$From.ProvBCR %in% regions)){
    miss.poly <- PolyRngAdj.i$From.ProvBCR[!PolyRngAdj.i$From.ProvBCR %in% regions]
    missing<-c(missing, paste(cn.i, 'Missing Polygon', miss.poly, sep=' - '))
  }
  # Any Polygons getting data from other regions?
  PS.other.i <- PS.othersources[PS.othersources$cn==cn.i, ]

  all.regions.i <- c(regions, PolyRngAdj.i$To.ProvBCR, PS.other.i$strata_name )
  all.regions.i <- all.regions.i[!duplicated(all.regions.i)]

  PS.mat  <-  matrix(nrow = length(all.regions.i), ncol=n.samp)
  rownames(PS.mat) <- all.regions.i

  df.ProvBCR.spp.i <- list()

  for(reg in 1:length(regions)){ #Begin loop through each State/Provence by BCR region
    if(!any(PolyAreas$ProvBCR==regions[reg]))next # if the region is not in the PolygonArea database skip region
    name <- paste(spp.i, regions[reg], sep=' - ')
    spp.i.by.region <- df.spp.i[df.spp.i$ProvBCR==regions[reg], c("mean.count",
                                                                      "var.count",
                                                                      "size.par") ]
    if(all(spp.i.by.region$mean.count==0))next # if the mean count for the region is zero skip region
    region.PS <- c()
    region.mean.count <- c()

    Dist.r <- runif(n.samp, Dist.adj.LB, Dist.adj.UB)
    Dist.r <- (400/Dist.r)^2
    Time.r <- rlnorm(n.samp, Time.adj.meanlog, Time.adj.sdlog)
    Pair.r <-  rtnorm(n.samp, Pair.adj.i, 0.13, lower = 1, upper = 2)

    for(resamp in 1:n.samp){ # generate distribution for Population Sizes at provence/state by BCR region
      all.routes <- c()
      for(route in 1:nrow(spp.i.by.region)){ #Begin loop: resample counts at BBS route level
        rroute <- sample(c(1:nrow(spp.i.by.region)), 1)
        if(is.na(spp.i.by.region$var.count[rroute]) | spp.i.by.region$mean.count[rroute]>=spp.i.by.region$var.count[rroute]){
          sample.route.r <- rpois(10,
                                  spp.i.by.region$mean.count[rroute])
        }else{
          sample.route.r <- rnbinom(10,
                                    mu = spp.i.by.region$mean.count[rroute], size = spp.i.by.region$size.par[rroute])}
        all.routes <- c(all.routes, sample.route.r)
      } #end route resample loop

      # Region Count with density with range adjustment
      if(regions[reg] %in% Rng.Adj.i$ProvBCR){
        DensAdjsd.r <- mean(all.routes)*Rng.Adj.i$RngAdj[Rng.Adj.i$ProvBCR==regions[reg]]*
          (PolyAreas$Area.km2[PolyAreas$ProvBCR==regions[reg]]/25.1)
      } else DensAdjsd.r <- mean(all.routes)*(PolyAreas$Area.km2[PolyAreas$ProvBCR==regions[reg]]/25.1)

      Count.r <- ceiling(DensAdjsd.r*Dist.r[resamp]*Time.r[resamp]*Pair.r[resamp])

      PS.mat[regions[reg], resamp] <- Count.r
      region.PS <- c(region.PS, Count.r)
      region.mean.count <- c(region.mean.count, mean(all.routes) )

    } # end region PS distribution

    temp.row <- data.frame(cn = cn.i, Aou = spp.i, ProvBCR = regions[reg],
                           n.routes = nrow(spp.i.by.region), mean.PopEst = mean(region.PS), med.PopEst = median(region.PS),
                           Coeff.var = sd(region.PS)/mean(region.PS),
                           LCI80.PopEst = quantile(region.PS, 0.1),
                           UCI80.PopEst = quantile(region.PS, 0.9),
                           LCI95.PopEst = quantile(region.PS, 0.025),
                           UCI95.PopEst = quantile(region.PS, 0.975),
                           rounded.mean.PopEst = PrettyUp(mean(region.PS)),
                           rounded.med.PopEst = PrettyUp(median(region.PS)),
                           rounded.LCI80.PopEst = PrettyUp(quantile(region.PS, 0.1)),
                           rounded.UCI80.PopEst = PrettyUp(quantile(region.PS, 0.9)),
                           rounded.LCI95.PopEst = PrettyUp(quantile(region.PS, 0.025)),
                           rounded.UCI95.PopEst = PrettyUp(quantile(region.PS, 0.975)), stringsAsFactors = FALSE, row.names = NULL)
    df.ProvBCR.spp.i[[name]] <-  temp.row

    #check if this region is a 'Borrow From' region
    if(any(regions[reg] %in% PolyRngAdj.i$From.ProvBCR)){
      brwd <- PolyRngAdj.i[PolyRngAdj.i$From.ProvBCR==regions[reg], ]
      for(b in nrow(brwd)){ #Generate PS estimate for each 'To' region
        name.b <- paste(spp.i, brwd$To.ProvBCR[b], sep=' - ')
        to.poly <- brwd$To.ProvBCR[b]
        borrow.b.PS <- c()
        for(resamp in 1:n.samp){ # resample loop
          DensAdjsd.r <- region.mean.count[resamp]*brwd$RngAdj[brwd$To.ProvBCR==to.poly]*
            (PolyAreas$Area.km2[PolyAreas$ProvBCR==to.poly]/25.1)

          Count.r <- ceiling(DensAdjsd.r*Dist.r[resamp]*Time.r[resamp]*Pair.r[resamp])

          PS.mat[to.poly, resamp] <- Count.r
          borrow.b.PS <- c(borrow.b.PS, Count.r)

        } # end resample loop
        temp.row <- data.frame(cn = cn.i, Aou = spp.i, ProvBCR = to.poly,
                               n.routes = NA, mean.PopEst = mean(borrow.b.PS), med.PopEst = median(borrow.b.PS),
                               Coeff.var = sd(borrow.b.PS)/mean(borrow.b.PS),
                               LCI80.PopEst = quantile(borrow.b.PS, 0.1),
                               UCI80.PopEst = quantile(borrow.b.PS, 0.9),
                               LCI95.PopEst = quantile(borrow.b.PS, 0.025),
                               UCI95.PopEst = quantile(borrow.b.PS, 0.975),
                               rounded.mean.PopEst = PrettyUp(mean(borrow.b.PS)),
                               rounded.med.PopEst = PrettyUp(median(borrow.b.PS)),
                               rounded.LCI80.PopEst = PrettyUp(quantile(borrow.b.PS, 0.1)),
                               rounded.UCI80.PopEst = PrettyUp(quantile(borrow.b.PS, 0.9)),
                               rounded.LCI95.PopEst = PrettyUp(quantile(borrow.b.PS, 0.025)),
                               rounded.UCI95.PopEst = PrettyUp(quantile(borrow.b.PS, 0.975)), stringsAsFactors = FALSE, row.names = NULL)

        df.ProvBCR.spp.i[[name.b]] <-  temp.row
      } # end loop through each 'To' borrow
    } # Close 'borrow' region check

  } #end loop through State/Provence by BCR region

  # Add Estimates for Pop Estimates from other Data Sources
  if(nrow(PS.other.i)>0){
    for(u in 1:nrow(PS.other.i)){
      if(u %in% regions){
      name.o <- paste(spp.i, PS.other.i$ProvBCR[u], sep=' - ')
      temp.row <- data.frame(cn = cn.i, Aou = spp.i, ProvBCR = PS.other.i$ProvBCR[u],
                             n.routes = NA, mean.PopEst = PS.other.i$PopEst[u], med.PopEst = PS.other.i$PopEst[u],
                             Coeff.var = NA,
                             LCI80.PopEst = PS.other.i$PopEst[u],
                             UCI80.PopEst = PS.other.i$PopEst[u],
                             LCI95.PopEst = PS.other.i$PopEst[u],
                             UCI95.PopEst = PS.other.i$PopEst[u],
                             rounded.mean.PopEst = NA,
                             rounded.med.PopEst = NA,
                             rounded.LCI80.PopEst = NA,
                             rounded.UCI80.PopEst = NA,
                             rounded.LCI95.PopEst = NA,
                             rounded.UCI95.PopEst = NA,  stringsAsFactors = FALSE, row.names = NULL)
      df.ProvBCR.spp.i[[name.o]] <-  temp.row
      PS.mat[PS.other.i$ProvBCR[u], ] <- rep(PS.other.i$PopEst[u], n.samp)
      }
      }
  } # End pull in data from other sources

  Spp.ProvBCR.i <- do.call(rbind, df.ProvBCR.spp.i)
  Spp.ProvBCR <- rbind(Spp.ProvBCR, Spp.ProvBCR.i)

  #RollUp
  PS.mat[is.na(PS.mat)] <- 0
  all.est <- data.frame(ProvBCR = rownames(PS.mat) )
  all.est$Prov <- lapply(as.character(rownames(PS.mat)), FUN=function(x) unlist(strsplit(x, '-'))[1])
  all.est$BCR <- lapply(as.character(rownames(PS.mat)), FUN=function(x) unlist(strsplit(x, '-'))[2])
  all.est <- cbind(all.est, PS.mat)

  all.est.out <- as.matrix(all.est[ ,c(2:ncol(all.est))])

  write.csv(all.est.out, paste(Output.dir, 'Spp.out/', spp.i, '_allest.csv', sep=""), row.names=FALSE)

  #by BCR
  all.BCRs <- c(unlist(unique(all.est$BCR)))
  for(bcr in 1:length(all.BCRs)){ #loop through BCR roundup
    name.bcr <- paste(spp.i, all.BCRs[bcr], sep=' - ')
    lst <- apply(all.est[all.est$BCR==all.BCRs[bcr], c(4:ncol(all.est))],
                 MARGIN = 2, FUN = sum, na.rm=TRUE)
    temp.row <- data.frame(cn = cn.i, Aou = spp.i, BCR = all.BCRs[bcr],
                           mean.PopEst = mean(lst), med.PopEst = median(lst),
                           Coeff.var = sd(lst)/mean(lst),
                           LCI80.PopEst = quantile(lst, 0.1),
                           UCI80.PopEst = quantile(lst, 0.9),
                           LCI95.PopEst = quantile(lst, 0.025),
                           UCI95.PopEst = quantile(lst, 0.975),
                           rounded.mean.PopEst = PrettyUp(mean(lst)),
                           rounded.med.PopEst = PrettyUp(median(lst)),
                           rounded.LCI80.PopEst = PrettyUp(quantile(lst, 0.1)),
                           rounded.UCI80.PopEst = PrettyUp(quantile(lst, 0.9)),
                           rounded.LCI95.PopEst = PrettyUp(quantile(lst, 0.025)),
                           rounded.UCI95.PopEst = PrettyUp(quantile(lst, 0.975)), stringsAsFactors = FALSE, row.names = NULL)
    Spp.by.BCR.i[[name.bcr]] <- temp.row
  } #end loop BCR roundup

  #by Provence/State
  all.Provs <- c(unlist(unique(all.est$Prov)))
  for(prov in 1:length(all.Provs)){ #loop through BCR roundup
    name.Prov <- paste(spp.i, all.Provs[prov], sep=' - ')
    lst <- apply(all.est[all.est$Prov==all.Provs[prov], c(4:ncol(all.est))],
                 MARGIN = 2, FUN = sum, na.rm=TRUE)
    temp.row <- data.frame(cn = cn.i, Aou = spp.i, Prov = all.Provs[prov],
                           mean.PopEst = mean(lst), med.PopEst = median(lst),
                           Coeff.var = sd(lst)/mean(lst),
                           LCI80.PopEst = quantile(lst, 0.1),
                           UCI80.PopEst = quantile(lst, 0.9),
                           LCI95.PopEst = quantile(lst, 0.025),
                           UCI95.PopEst = quantile(lst, 0.975),
                           rounded.mean.PopEst = PrettyUp(mean(lst)),
                           rounded.med.PopEst = PrettyUp(median(lst)),
                           rounded.LCI80.PopEst = PrettyUp(quantile(lst, 0.1)),
                           rounded.UCI80.PopEst = PrettyUp(quantile(lst, 0.9)),
                           rounded.LCI95.PopEst = PrettyUp(quantile(lst, 0.025)),
                           rounded.UCI95.PopEst = PrettyUp(quantile(lst, 0.975)), stringsAsFactors = FALSE, row.names = NULL)
    Spp.by.Prov.i[[name.Prov]] <- temp.row
  } #end loop BCR roundup

  # WH and Global estimates
  p.GL.i <- WH.GL.Adj$pGL_UsCa[WH.GL.Adj$cn == cn.i]
  p.WH.i <- WH.GL.Adj$pWH_UsCa[WH.GL.Adj$cn == cn.i]
  lst <- apply(all.est[ , c(4:ncol(all.est))], MARGIN = 2, FUN = sum, na.rm=TRUE)
  temp.row <- data.frame(cn = cn.i, Aou = spp.i,
                         GL.mean.PopEst = mean(lst)/p.GL.i,
                         GL.med.PopEst = median(lst)/p.GL.i,
                         GL.80LCI.PopEst = quantile(lst, 0.1)/p.GL.i,
                         GL.80UCI.PopEst = quantile(lst, 0.9)/p.GL.i,
                         GL.95LCI.PopEst = quantile(lst, 0.025)/p.GL.i,
                         GL.95UCI.PopEst = quantile(lst, 0.975)/p.GL.i,
                         WH.mean.PopEst = mean(lst)/p.WH.i,
                         WH.med.PopEst = median(lst)/p.WH.i,
                         WH.80LCI.PopEst = quantile(lst, 0.1)/p.WH.i,
                         WH.80UCI.PopEst = quantile(lst, 0.9)/p.WH.i,
                         WH.95LCI.PopEst = quantile(lst, 0.025)/p.WH.i,
                         WH.95UCI.PopEst = quantile(lst, 0.975)/p.WH.i,
                         USCAN.mean.PopEst = mean(lst),
                         USCAN.med.PopEst = median(lst),
                         USCAN.80LCI.PopEst = quantile(lst, 0.1),
                         USCAN.80UCI.PopEst = quantile(lst, 0.9),
                         USCAN.95LCI.PopEst = quantile(lst, 0.025),
                         USCAN.95UCI.PopEst = quantile(lst, 0.975),
                         rounded.GL.mean.PopEst = PrettyUp(mean(lst)/p.GL.i),
                         rounded.GL.med.PopEst = PrettyUp(median(lst)/p.GL.i),
                         rounded.GL.80LCI.PopEst = PrettyUp(quantile(lst, 0.1)/p.GL.i),
                         rounded.GL.80UCI.PopEst = PrettyUp(quantile(lst, 0.9)/p.GL.i),
                         rounded.GL.95LCI.PopEst = PrettyUp(quantile(lst, 0.025)/p.GL.i),
                         rounded.GL.95UCI.PopEst = PrettyUp(quantile(lst, 0.975)/p.GL.i),
                         rounded.WH.mean.PopEst = PrettyUp(mean(lst)/p.WH.i),
                         rounded.WH.med.PopEst = PrettyUp(median(lst)/p.WH.i),
                         rounded.WH.80LCI.PopEst = PrettyUp(quantile(lst, 0.1)/p.WH.i),
                         rounded.WH.80UCI.PopEst = PrettyUp(quantile(lst, 0.9)/p.WH.i),
                         rounded.WH.95LCI.PopEst = PrettyUp(quantile(lst, 0.025)/p.WH.i),
                         rounded.WH.95UCI.PopEst = PrettyUp(quantile(lst, 0.975)/p.WH.i),
                         rounded.USCAN.mean.PopEst = PrettyUp(mean(lst)),
                         rounded.USCAN.med.PopEst = PrettyUp(median(lst)),
                         rounded.USCAN.80LCI.PopEst = PrettyUp(quantile(lst, 0.1)),
                         rounded.USCAN.80UCI.PopEst = PrettyUp(quantile(lst, 0.9)),
                         rounded.USCAN.95LCI.PopEst = PrettyUp(quantile(lst, 0.025)),
                         rounded.USCAN.95UCI.PopEst = PrettyUp(quantile(lst, 0.975)),
                         PS.Score = PS.score(median(lst)/p.GL.i),
                         PS.80LCI.Score = PS.score(quantile(lst, 0.1)/p.GL.i),
                         PS.80UCI.Score = PS.score(quantile(lst, 0.9)/p.GL.i),
                         PS.95LCI.Score = PS.score(quantile(lst, 0.025)/p.GL.i),
                         PS.95UCI.Score = PS.score(quantile(lst, 0.975)/p.GL.i),
                         stringsAsFactors = FALSE, row.names = NULL)

  Spp.by.WH.GL.i[[spp.i]] <- temp.row
  print(paste(cn.i, round(i/length(Spp.aou),2)))
} # end loop through species

Spp.by.BCR <- do.call(rbind, Spp.by.BCR.i)

Spp.by.Prov <- do.call(rbind, Spp.by.Prov.i)
Spp.by.Prov <- merge(x=Spp.by.Prov, y=region.lkup, all.x=TRUE)
Spp.by.WH.GL <- do.call(rbind, Spp.by.WH.GL.i)


Spp.ProvBCR$Prov <- lapply(as.character(Spp.ProvBCR$ProvBCR), FUN=function(x) unlist(strsplit(x, '-'))[1])
Spp.ProvBCR$Prov <- as.integer(Spp.ProvBCR$Prov)
Spp.ProvBCR$BCR <- lapply(as.character(Spp.ProvBCR$ProvBCR), FUN=function(x) unlist(strsplit(x, '-'))[2])
Spp.ProvBCR$BCR <- as.integer(Spp.ProvBCR$BCR)
Spp.ProvBCR <- Spp.ProvBCR[,c('cn', 'Aou', 'Prov', 'BCR', "n.routes", "mean.PopEst",  "med.PopEst",
                              "Coeff.var" ,  "LCI80.PopEst" , "UCI80.PopEst", "LCI95.PopEst", "UCI95.PopEst",
                               "rounded.mean.PopEst" ,"rounded.med.PopEst", "rounded.LCI80.PopEst", "rounded.UCI80.PopEst",
                              "rounded.LCI95.PopEst", "rounded.UCI95.PopEst")]
Spp.ProvBCR <- merge(x=Spp.ProvBCR, y=region.lkup, all.x=TRUE)

write.csv(Spp.ProvBCR,
          paste(Output.dir, 'PSest_Prov_by_BCR_', Spp.dat.date, "_", n.samp, 'iter.csv', sep=''),
          row.names = FALSE)
write.csv(Spp.by.BCR,
          paste(Output.dir,'PSest_by_BCR_', Spp.dat.date, "_", n.samp, 'iter.csv', sep=""),
          row.names = FALSE)
write.csv(Spp.by.Prov,
          paste(Output.dir,'PSest_by_Prov_', Spp.dat.date, "_", n.samp, 'iter.csv', sep=''),
          row.names = FALSE)
write.csv(Spp.by.WH.GL,
          paste(Output.dir, 'PSest_Global_WH_NA_', Spp.dat.date, "_", n.samp, 'iter.csv', sep=''),
          row.names = FALSE)
write.csv(missing,
          paste(Output.dir,'Species_missing_PopEst_', Spp.dat.date, '.csv', sep=''),
          row.names = FALSE)



# Compare w published 2019 ------------------------------------------------

pif_19 <- readxl::read_xlsx("data/PopEsts Global 2020.04.29.xlsx",
                            .name_repair = "universal") |>
  dplyr::mutate(cn = English.Name,
                rounded.USCAN.mean.PopEst = Population.Estimate.USA..Canada,
                rounded.USCAN.95LCI.PopEst = Lower.95..bound.USA..Canada,
                rounded.USCAN.95UCI.PopEst = Upper.95..bound.USA..Canada) |>
  dplyr::select(cn,rounded.USCAN.mean.PopEst,
         rounded.USCAN.95LCI.PopEst,
         rounded.USCAN.95UCI.PopEst,
         Source.for.USA.Canada.Estimate) |>
  dplyr::filter(grepl("PIFcalc19",Source.for.USA.Canada.Estimate)) |>
  dplyr::mutate(source = "PIF_2019")


pif_new <- Spp.by.WH.GL |>
  dplyr::select(cn,rounded.USCAN.mean.PopEst,
                rounded.USCAN.95LCI.PopEst,
                rounded.USCAN.95UCI.PopEst)|>
  dplyr::mutate(source = "PIF_2022")


pif_comp <- pif_19 |>
  dplyr::select(-c(Source.for.USA.Canada.Estimate,source)) |>
  dplyr::rename(rounded.USCAN.mean.PopEst.2019 = rounded.USCAN.mean.PopEst,
                rounded.USCAN.95LCI.PopEst.2019 = rounded.USCAN.95LCI.PopEst,
                rounded.USCAN.95UCI.PopEst.2019 = rounded.USCAN.95UCI.PopEst)  |>
  dplyr::left_join(pif_new,
                   by = "cn") |>
  dplyr::mutate(dif = (rounded.USCAN.mean.PopEst - rounded.USCAN.mean.PopEst.2019)/rounded.USCAN.mean.PopEst.2019,
                cv.2019 = (0.25*(rounded.USCAN.95UCI.PopEst.2019-rounded.USCAN.95LCI.PopEst.2019))/rounded.USCAN.mean.PopEst.2019,
                cv = (0.25*(rounded.USCAN.95UCI.PopEst-rounded.USCAN.95LCI.PopEst))/rounded.USCAN.mean.PopEst)

large_dif <- pif_comp |>
  dplyr::filter(dif > 0.2 | dif < -0.2)

library(ggplot2)

compare_plot <- ggplot(data = pif_comp,
                       aes(x = rounded.USCAN.mean.PopEst.2019,
                           y = rounded.USCAN.mean.PopEst))+
  geom_errorbar(aes(ymin = rounded.USCAN.95LCI.PopEst,
                    ymax = rounded.USCAN.95UCI.PopEst),
                alpha = 0.3)+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  scale_y_continuous(transform = "log10")+
  scale_x_continuous(transform = "log10")


compare_plot

compare_plot_e <- ggplot(data = pif_comp,
                       aes(x = rounded.USCAN.mean.PopEst.2019,
                           y = rounded.USCAN.mean.PopEst))+
  geom_errorbar(aes(ymin = rounded.USCAN.95LCI.PopEst,
                    ymax = rounded.USCAN.95UCI.PopEst),
                alpha = 0.2)+
  geom_errorbarh(aes(xmin = rounded.USCAN.95LCI.PopEst.2019,
                    xmax = rounded.USCAN.95UCI.PopEst.2019),
                alpha = 0.2)+
  geom_point(alpha = 0.4)+
  geom_abline(intercept = 0, slope = 1)+
  scale_y_continuous(transform = "log10")+
  scale_x_continuous(transform = "log10")+
  theme_bw()


compare_plot_e



compare_plot_cv <- ggplot(data = pif_comp,
                         aes(x = cv.2019,
                             y = cv))+
  geom_point(alpha = 0.3)+
  geom_abline(intercept = 0, slope = 1)+
  xlab("CV 2019 (SE/mean)")+
  ylab("CV 2022 (SE/mean)")


compare_plot_cv


compare_plot2 <- ggplot(data = large_dif,
                       aes(x = rounded.USCAN.mean.PopEst.2019,
                           y = rounded.USCAN.mean.PopEst))+
  geom_errorbar(aes(ymin = rounded.USCAN.95LCI.PopEst,
                    ymax = rounded.USCAN.95UCI.PopEst),
                alpha = 0.2)+
  geom_errorbarh(aes(xmin = rounded.USCAN.95LCI.PopEst.2019,
                     xmax = rounded.USCAN.95UCI.PopEst.2019),
                 alpha = 0.2)+
  geom_point(alpha = 0.4)+
  geom_abline(intercept = 0, slope = 1)+
  scale_y_continuous(transform = "log10")+
  scale_x_continuous(transform = "log10")+
  geom_text(aes(label = cn),
            size = 2)+
  theme_bw()


compare_plot2


library(patchwork)
pdf("comparison_2022_2019_USA_CAN.pdf",
    width = 11, height = 8.5)
print(compare_plot_e)
print(compare_plot2)
print(compare_plot_cv)
dev.off()
write.csv(pif_comp,"comparison_pif2022_2019.csv",row.names = FALSE)

