##############################
# Assemble Data from  
# North American Breeding Bird Survey
# For use in PIF Population Size Estimate
#
# User may specify download data direct from BBS ftp site
# or provide directory of previously downloaded data
# 
# Script calculates mean count by route for last 10 years of
# data supplied
#
# Currently only outputs North American landbird species
# (specified in Species_PSEst_Adjustments.csv file)
# 
# Required user Inputs are specified by prompt in the terminal
#
# Code by J.C. Stanton
#
######  USER INPUTS ########
#############################

# Provide Path to Input file directory - BBS DATA and Adjustment parameter files
Inputfiles.dir <- readline("Provide directory path to location of INPUT files: ")
if(!file.exists(Inputfiles.dir)) stop("Provide valid directory path to location of input files. ")
Inputfiles.dir <- paste(Inputfiles.dir, '/', sep="")

# Provide Path to where formatted BBS species level data should be written
Output.dir <- readline("Provide directory path for OUTPUT files: ")
if(!file.exists(Output.dir)) stop("Provide valid directory path for species output files. ")
Output.dir <- paste(Output.dir, '/', sep="")

#### Download new Dataset?
#  If False, provide a path to previously downloaded BBS state files
#  If True, you must be connected to the internet
Download.new.data.ask <- menu(c('Yes', 'No'), 
                              title = "Do you wish to download BBS Data? 
                              If yes, you must have an internet connection
                              If no, provide directory name where files are stored in input directory")
Download.new.data <- c(TRUE, FALSE)[Download.new.data.ask]
# Provide Name of dir in Inputfiles where previously downloaded data is stored 
# (If Download.new.data=FALSE)
if(!Download.new.data) prev.download.dir <- readline("provide directory of previously downloaded BBS files in input directory ")
# Also provide dependent files: 
# 'SpeciesList.txt', 'Weather.zip', 'Routes.zip'

#### NON-BBS Dependent files in Inputfiles.dir
# 'Species_PSEst_Adjustments.csv'
# 'RegionCodes.csv' 
# (similar to RegionCodes.txt on BBS ftp, but with 2-letter Prov.letters column)

##############################
##############################

bbs.dir <- "ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/" #URL for BBS FTP
#Enter file path to BBS state files (state.files.dir)
if(!Download.new.data){
  state.files.dir <- paste(Inputfiles.dir, prev.download.dir, sep='')
  if(!file.exists(state.files.dir)) stop("Suppy valid directory to previously downloaded BBS data, or
                                         set Download.new.data to TRUE")
}
################################################

if(Download.new.data){
  download.file(paste(bbs.dir,'SpeciesList.txt', sep=''), destfile=paste(Inputfiles.dir,'SpeciesList.txt', sep=""))
  download.file(paste(bbs.dir,'Routes.zip',sep=''), destfile=paste(Inputfiles.dir,'Routes.zip', sep=""))
  download.file(paste(bbs.dir,'Weather.zip',sep=''), destfile=paste(Inputfiles.dir,'Weather.zip', sep=''))
}

#download and read supporting files
#Species lookup
spp.info <- read.fwf(paste(Inputfiles.dir, 'SpeciesList.txt', sep=''), skip=8, width = c(6, 6,rep(51,7) ),
                     col.names=c('Seq', 'AOU',"English_Common_Name", 'French_Common_Name',
                                 'Spanish_Common_Name','ORDER', 'Family', 'Genus', 'Species'))
spp.lkup <- data.frame(Aou=as.character(spp.info$AOU), 
                       cn=gsub('^\\s+|\\s+$', '', spp.info$English_Common_Name),
                       stringsAsFactors = FALSE)
rm(spp.info)


##### dependent datasets  #####
## Data dependent on previous calculations by P.Blancher and J.stanton

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE)
# Add AOU #'s to Species names in table of adjustment factors
ExpAdjs <- merge(ExpAdjs, spp.lkup, all.x=TRUE)

######## taxonomy adjustments ######################
####### Species that are by races in BBS DB to be grouped (sp.group) ######
Dark.eyed.junco.races <- c('(Slate-colored Junco) Dark-eyed Junco', 
                           '(Oregon Junco) Dark-eyed Junco',
                           '(Pink-sided Junco) Dark-eyed Junco',
                           '(White-winged Junco) Dark-eyed Junco',
                           '(Gray-headed Junco) Dark-eyed Junco',
                           '(Guadalupe Junco) Dark-eyed Junco',
                           '(unid. race) Dark-eyed Junco')
Dark.eyed.juncos <- spp.lkup[spp.lkup$cn %in% Dark.eyed.junco.races, ]

Northern.flicker.races <- c('(Yellow-shafted Flicker) Northern Flicker',
                            '(Red-shafted Flicker) Northern Flicker',
                            '(unid. Red/Yellow Shafted) Northern Flicker', 
                            'hybrid Northern Flicker (Red x Yellow-shafted)')
Northern.flickers <- spp.lkup[spp.lkup$cn %in% Northern.flicker.races, ]

Yellow.rumped.warbler.races <- c('(Myrtle Warbler) Yellow-rumped Warbler',
                                 "(Audubon's Warbler) Yellow-rumped Warbler",
                                 "(unid. Myrtle/Audubon's) Yellow-rumped Warbler")
Yellow.rumped.warblers <- spp.lkup[spp.lkup$cn %in% Yellow.rumped.warbler.races, ]

sp.group <- list(Dark.eyed.juncos, Northern.flickers, Yellow.rumped.warblers )
names(sp.group) <- c("Dark-eyed Junco", "Northern Flicker", "Yellow-rumped Warbler" )

ExpAdjs$Aou[ExpAdjs$cn == 'Dark-eyed Junco'] <- "5670"
ExpAdjs$Aou[ExpAdjs$cn == 'Northern Flicker'] <- "4120"
ExpAdjs$Aou[ExpAdjs$cn == 'Yellow-rumped Warbler'] <- "6550"

sp.group.aou <- c()
  for(g in 1:length(sp.group)){
    sp.group.aou <- c(sp.group.aou, sp.group[[g]]$Aou)
  }
####################################################################

spp <- unique(ExpAdjs$cn)
spp <- c(spp, names(sp.group))
aou <- unique(ExpAdjs$Aou)
aou <- c(aou, sp.group.aou)

#Route Info
routes <- read.csv(unz(paste(Inputfiles.dir, 'Routes.zip',sep=''), filename='routes.csv'))
#Weather Info (survey quality)
weather <- read.csv(unz(paste(Inputfiles.dir, 'Weather.zip', sep=''), filename='weather.csv'))
# Letter to numeric code lookup for states/prov
region.lkup <- read.csv(paste(Inputfiles.dir, 'RegionCodes.csv', sep=''))

routes$rteno <- paste(routes$countrynum, routes$statenum, routes$Route, sep="-")

# Establish beginning and end survey years from data
# Note: This is a change from orignal code - Estimates are based on past 10 years of data automatically
endyear <- max(weather$Year) #latest year of data
nyears <- endyear - 1965 #survey began in 1966, 1965 would be year zero
tsteps <- c((endyear-9):endyear)

# Establish runs that meet acceptable conditions (RunType) and are primary routes (RPID)
# Note: RPID==101 is a change from original code
#       For routes repeated several times in a season only the primary run is recorded so as not to
#       over-represent one route relative to routes that are only run once
weather$Keep <- weather$RunType==1 & weather$RPID==101 
weather.ok <- subset(weather, Keep==TRUE)
weather.ok["rteno"] <- NA
weather.ok$rteno <- paste(weather.ok$countrynum,weather.ok$statenum, weather.ok$Route, sep="-")
weather.ok.last10 <- subset(weather.ok, Year %in% c((endyear-9):endyear)) 

if(nrow(weather.ok.last10[duplicated(weather.ok.last10[,c('Year', 'rteno')]), ])!=0) stop('Duplicate Routes')
route.yrs.ok.last10 <- weather.ok.last10[,c('rteno', 'Year')] # Years acceptable routes were run

routes.ok.last10 <- subset(routes, rteno %in% weather.ok.last10$rteno) # Info acceptable routes run last 10 years
routes.ok.last10$ProvBCR <- paste(routes.ok.last10$statenum, routes.ok.last10$BCR, sep='-')

if(Download.new.data){ # Check if downloading new dataset
  library(RCurl)
  states.dir <- paste(bbs.dir, 'States/', sep="")
  dl.state.target <- paste(Inputfiles.dir,'/States',Sys.Date(), sep='')
  dir.create(dl.state.target)
  st.lst <- getURL(states.dir) #ugly scrape of page
  st.lst <- unlist(strsplit(st.lst, split=' ')) 
  st.lst <- st.lst[grep('.zip',st.lst, fixed=TRUE)] #cleanup 
  
  df <- data.frame()
  for(i in 1:length(st.lst)){  #loop through state and download, unzip, and read in .zip files
    st.name.i <- st.lst[i]
    st.name.i <- unlist(strsplit(st.name.i, split=".", fixed=TRUE))[1] #trim off gobbledygook
    st.zip.i <- paste(st.name.i, ".zip", sep='') 
    
    download.file(paste(states.dir,st.zip.i, sep=''),destfile=paste(dl.state.target, st.zip.i, sep='/'))
    st.i <- data.frame()
    st <- read.csv(unz(paste(dl.state.target, st.zip.i, sep='/'),filename=paste(st.name.i,'.csv',sep='')))
    st <- subset(st, RPID==101)
    st$rteno <- paste(st$countrynum, st$statenum, st$Route, sep="-")
    st <- merge(x=st, y=route.yrs.ok.last10) #limit to ok route years
    # Average counts for routes over last 10 years 
    #         - include zeros for routes run, but species not encountered
    spp.lst.j <- unique(st$Aou)
    spp.lst.j <- spp.lst.j[spp.lst.j %in% aou]
    rtes.lst.i <- unique(st$rteno)
    # Build matrix of which routes were run and not run by year
    oklast10.mat <- matrix(NA, length(rtes.lst.i),10)
    for (l in 1:length(rtes.lst.i)) {    #begin loop through routes 
      rte.yrs <- which(tsteps %in% st$Year[st$rteno==rtes.lst.i[l]]) #index years of routes run
      oklast10.mat[l,rte.yrs] <- 0
    }
    st.spp <- c()
    for(j in 1:length(spp.lst.j)){ #loop through each species in state
      if(spp.lst.j[j] %in% st.spp) next # Check if spp group already processed
      # check if spp is part of a taxonomic group
      if(spp.lst.j[j] %in% sp.group.aou){
      group.j <- sp.group[[grep(spp.lst.j[j], sp.group)]]$Aou
      st.j <- st[st$Aou %in% group.j, c('Year', 'rteno', 'Aou', 'SpeciesTotal')]
      st.j <- aggregate(SpeciesTotal ~ Year + rteno, data = st.j, sum)
      st.spp <- c(st.spp, group.j) # Add group to list of spp already processed
      } else {
        st.j <- st[st$Aou==spp.lst.j[j], c('Year', 'rteno', 'SpeciesTotal') ]
        st.spp <- c(st.spp, spp.lst.j[j]) # keep list of species already processed
        }
      oklast10.mat.j <- oklast10.mat
      for (l in 1:length(rtes.lst.i)){ #update matrix with total species counts (zeros where not encountered)
        if(!any(st.j$rteno %in% rtes.lst.i[l])) next
        rte.yrs.j <- which(tsteps %in% st.j$Year[st.j$rteno==rtes.lst.i[l]]) #index years of routes run
        temp <- st.j$SpeciesTotal[st.j$rteno==rtes.lst.i[l]] #  total Counts 
        oklast10.mat.j[l,rte.yrs.j] <- temp
      }# end loop routes for each species
      mean.counts.j <- apply(oklast10.mat.j, 1, FUN=mean, na.rm=TRUE)
      n.routes.j <- apply(oklast10.mat.j, 1, FUN= function(x) length(which(!is.na(x)))) #how many times was route run?
      var.counts.j <- apply(oklast10.mat.j, 1, var, na.rm=TRUE)
      
      sp.j.st <- data.frame(Aou = rep(spp.lst.j[j], length(rtes.lst.i)), rteno = rtes.lst.i,
                            mean.count = mean.counts.j, 
                            var.count = var.counts.j, n.routes = n.routes.j, 
                            size.par = (mean.counts.j^2)/(var.counts.j-mean.counts.j))
      st.i <- rbind(st.i, sp.j.st)
    }# end loop through species
    df <- rbind(df, st.i)
  } #End download and load StateFile
  
} else { # If data is already downloaded
  St.lst <- list.files(path=state.files.dir, full.names=TRUE)
  st.names <- list.files(path=state.files.dir)
  df <- data.frame()
  for(i in 1:length(St.lst)){#loop through state and unzip, and read in .zip files
    st.i <- data.frame()
    st.zip.i <- St.lst[i]
    st.name.i <- unlist(strsplit(st.names[i], split=".", fixed=TRUE))[1] 
    st <- read.csv(unz(st.zip.i,filename=paste(st.name.i,'.csv',sep='')))
    st <- subset(st, RPID==101)
    st$rteno <- paste(st$countrynum, st$statenum, st$Route, sep="-")
    st <- merge(x=st, y=route.yrs.ok.last10) #limit to ok route years
    # Average counts for routes over last 10 years 
    #         - include zeros for routes run, but species not encountered
    spp.lst.j <- unique(st$Aou)
    spp.lst.j <- spp.lst.j[spp.lst.j %in% aou]
    rtes.lst.i <- unique(st$rteno)
    # Build matrix of which routes were run and not run by year
    oklast10.mat <- matrix(NA, length(rtes.lst.i),10)
    for (l in 1:length(rtes.lst.i)) {    #begin loop through routes 
      rte.yrs <- which(tsteps %in% st$Year[st$rteno==rtes.lst.i[l]]) #index years of routes run
      oklast10.mat[l,rte.yrs] <- 0
    }
    st.spp <- c()
    for(j in 1:length(spp.lst.j)){ #loop through each species in state
      if(spp.lst.j[j] %in% st.spp) next # Check if spp group already processed
      # check if spp is part of a taxonomic group
      if(spp.lst.j[j] %in% sp.group.aou){
        group.j <- sp.group[[grep(spp.lst.j[j], sp.group)]]$Aou
        st.j <- st[st$Aou %in% group.j, c('Year', 'rteno', 'Aou', 'SpeciesTotal')]
        st.j <- aggregate(SpeciesTotal ~ Year + rteno, data = st.j, sum)
        st.spp <- c(st.spp, group.j) # Add group to list of spp already processed
      } else {
        st.j <- st[st$Aou==spp.lst.j[j],  c('Year', 'rteno', 'SpeciesTotal')]
        st.spp <- c(st.spp, spp.lst.j[j]) # keep list of species already processed
      }
      oklast10.mat.j <- oklast10.mat
      for (l in 1:length(rtes.lst.i)){ #update matrix with total species counts (zeros where not encountered)
        if(!any(st.j$rteno %in% rtes.lst.i[l])) next
        rte.yrs.j <- which(tsteps %in% st.j$Year[st.j$rteno==rtes.lst.i[l]]) #index years of routes run
        temp <- st.j$SpeciesTotal[st.j$rteno==rtes.lst.i[l]] #  total Counts 
        oklast10.mat.j[l,rte.yrs.j] <- temp
      }# end loop routes for each species
      mean.counts.j <- apply(oklast10.mat.j, 1, mean, na.rm=TRUE)
      n.routes.j <- apply(oklast10.mat.j, 1, FUN= function(x) length(which(!is.na(x))))
      var.counts.j <- apply(oklast10.mat.j, 1, var, na.rm=TRUE)
      
      sp.j.st <- data.frame(Aou = rep(spp.lst.j[j], length(rtes.lst.i)), rteno = rtes.lst.i,
                            mean.count = mean.counts.j, 
                            var.count = var.counts.j, n.routes = n.routes.j, 
                            size.par = (mean.counts.j^2)/(var.counts.j-mean.counts.j))
      st.i <- rbind(st.i, sp.j.st)
    }# end loop through species
    df <- rbind(df, st.i)
  } #End load StateFiles
}

#write.csv(df, paste(Output.dir, 'BBS_10yr', Sys.Date(), '.csv', sep=''), row.names = FALSE)

all.sp <- unique(df$Aou)
st.spp <- c()
for(i in 1:length(all.sp)){
  if(all.sp[i] %in% st.spp) next # Check if spp group already processed
  if(all.sp[i] %in% sp.group.aou){
    group.i <- sp.group[[grep(all.sp[i], sp.group)]]$Aou
    df.i <- df[df$Aou %in% group.i, ]
    df.i$Aou <- all.sp[i]
    st.spp <- c(st.spp, group.i) # Add group to list of spp already processed
  } else {
    df.i <- df[df$Aou == all.sp[i], ]
    st.spp <- c(st.spp, all.sp[i])
  }
  write.csv(df.i, paste(Output.dir, all.sp[i], '_BBS_10yr', Sys.Date(), '.csv', sep=''), row.names = FALSE)
}

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
write.csv(ExpAdjs, paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), row.names=FALSE)









