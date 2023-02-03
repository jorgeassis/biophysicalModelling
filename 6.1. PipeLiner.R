## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)
closeAllConnections()

source("0. Config.R")
source("Dependences.R")

number.cores <- 10

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

pipeLiner <- TRUE

global.simulation.parameters <- loadRData(paste0(results.folder,"/","modelParameters.RData"))
n.days.max <- global.simulation.parameters$particle.max.duration
n.days.max

list.dirs(path = paste0("../Results"), recursive = FALSE)
season <- "" # c("YearRound","SeasonSummer","SeasonWinter")
spawn.p <- 1:12  # spawn.p <- c(6,7,8,9)

n.pld.period <- 0 #c(7) # c( 5 , 10 , 15 , 30 , 60 , 90 , 120, 150 , 180 ) # 1:120 c(10 , 30 , 90 , 120 , 200)
n.seasons <- "" # c("YearRound","Spring","Summer","Autumn","Winter")
combinations <- expand.grid(season=n.seasons,pld.period=n.pld.period,stringsAsFactors = F)

# -------------------------------------------------- [!]

# Test forceUnidirectional // transformFST

forceUnidirectional <- FALSE
transformFST <- FALSE

# -------------------------------------------------- [!]

directories <- list.files("../Data/geneticDataClean/")
directories <- directories[!grepl("xls",directories)]
directories <- directories[!grepl("txt",directories)]
length(directories) 

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

library("readxl")

for( mainfileName in directories ) {
  
  gc(reset=TRUE)
  
  mainFiles <- paste0("../Data/geneticDataClean/",mainfileName,"/")
  mainFiles <- list.files(mainFiles, pattern="xls", full.names=TRUE)
  mainFiles <- mainFiles[!grepl("~",mainFiles)]
  
  mainFile <- data.frame()
  for( mainFiles.i in mainFiles ) {
    mainFile.i <- as.data.frame(read_excel(mainFiles.i))
    if( sum(c("species","lon","lat") %in% tolower(colnames(mainFile.i))) == 3) {
      mainFile <- plyr::rbind.fill(mainFile,mainFile.i)
    }
  }
  
  names(mainFile) <- gsub("\\,","",names(mainFile))
  names(mainFile) <- gsub("\\.","",names(mainFile))
  
  speciesToModel <- unique(mainFile$Species)
  
  for( sp in 1:length(speciesToModel)) {
    
    markersToModel <- unique(mainFile[mainFile$Species == speciesToModel,"marker"])
    
    for( mk in 1:length(markersToModel)) {
      
      markerName <- markersToModel[mk]
      speciesName <- speciesToModel[sp]
      
      cat("Study",mainfileName,"\\",speciesName,"\\",markerName,"\n")
      
      popCoordinates <- mainFile[mainFile$Species == speciesName & mainFile$marker == markerName,c("Lon","Lat")]
      popCoordinates.n <- mainFile[mainFile$Species == speciesName & mainFile$marker == markerName ,c("SampleCode")]
      popCoordinates.n <- gsub("\\.","",popCoordinates.n)
      
      nPop <- nrow(popCoordinates)
      removePop <- NULL
      
      if( nrow(unique(popCoordinates)) < 3 ) { next }
      
      popDifferentiation.all <- paste0("../Data/geneticDataClean/",mainfileName,"/")
      popDifferentiation.all <- list.files(popDifferentiation.all, pattern="xls", full.names=TRUE)
      popDifferentiation.all <- popDifferentiation.all[!grepl("~",popDifferentiation.all)]
      popDifferentiation.all <- popDifferentiation.all[!grepl("pdf",popDifferentiation.all)]
      popDifferentiation.all <- popDifferentiation.all[grepl(markerName,popDifferentiation.all)]
      
      indices <- c("fst","gst","jostd")[as.numeric(which(sapply( c("fst","gst","jostd") , function(x) { TRUE %in% grepl( x , tolower(popDifferentiation.all)) } )))]
      
      if(length(indices) == 0) { next }
      
      for( index in indices ) {
        
        popDifferentiation <- popDifferentiation.all[grepl(index,tolower(popDifferentiation.all))]
        
        if(length(popDifferentiation) > 1) { next }
        if(length(popDifferentiation) == 0) { next }
        
        # Loop per index
        
        popDifferentiation <- as.data.frame(read_excel(popDifferentiation))
        names(popDifferentiation) <- gsub("\\.","",names(popDifferentiation))
        
        ## ------------
        ## ------------
        
        if( colnames(popDifferentiation)[1] == "X" | colnames(popDifferentiation)[1] == "1" ) { 
          
          colnames(popDifferentiation)[1] <- "X"
          popDifferentiation[,"X"] <- gsub("\\.","",popDifferentiation[,"X"])
          
          siteToKeep.diff <- which(popDifferentiation[,"X"] %in% colnames(popDifferentiation))
          popDifferentiation <- popDifferentiation[siteToKeep.diff,]
          
          siteToKeep.diff <- which(colnames(popDifferentiation) %in% popDifferentiation[,"X"])
          popDifferentiation <- popDifferentiation[,siteToKeep.diff]
          
        }
        
        siteToRemove.diff <- colnames(popDifferentiation)[which( ! colnames(popDifferentiation) %in% popCoordinates.n )]
        
        if( length(siteToRemove.diff) > 0 ) {
          popDifferentiation <- popDifferentiation[ - which( ! colnames(popDifferentiation) %in% popCoordinates.n ), - which( ! colnames(popDifferentiation) %in% popCoordinates.n )]
        }
        
        siteToRemove.diff <- popCoordinates.n[which( ! popCoordinates.n %in% colnames(popDifferentiation) )]
        
        if( length(siteToRemove.diff) > 0 ) {
          popCoordinates <- popCoordinates[which( popCoordinates.n %in% colnames(popDifferentiation) ),]
          popCoordinates.n <- popCoordinates.n[which( popCoordinates.n %in% colnames(popDifferentiation) )]
        }
        
        siteToRemove.NA <- sapply(1:(ncol(popDifferentiation)-1), function(x) { ifelse( sum(is.na(popDifferentiation[,x] )) == nrow(popDifferentiation) , x , 0)  })
        siteToRemove.NA <- siteToRemove.NA[siteToRemove.NA != 0]
        
        if( length(siteToRemove.NA) > 0 ) {
          stop("!")
          popCoordinates <- popCoordinates[-siteToRemove.NA,]
          popCoordinates.n <- popCoordinates.n[-siteToRemove.NA]
          popDifferentiation <- popDifferentiation[-siteToRemove.NA,-siteToRemove.NA]
        }
        
        ## ------------
        ## ------------
        
        finalSorting <- sapply(colnames(popDifferentiation), function(x) { which(popCoordinates.n == x )[1] })
        popCoordinates.n <- popCoordinates.n[finalSorting]
        popCoordinates <- popCoordinates[finalSorting,]
        
        finalSorting <- sapply(popCoordinates.n, function(x) { which(colnames(popDifferentiation) == x )})
        popDifferentiation <- popDifferentiation[finalSorting,finalSorting]
        
        ## ------------
        ## ------------
        
        f <- function(m) {
          m[upper.tri(m)] <- t(m)[upper.tri(m)]
          m
        }
        
        popDifferentiation <- f(popDifferentiation)
        popDifferentiation[popDifferentiation < 0] <- 0
        
        ## ------------
        ## ------------
        
        if( ncol(popDifferentiation) != length(popCoordinates.n) ) { stop("Files do not match in extent") }
        if( FALSE %in% (colnames(popDifferentiation) %in% popCoordinates.n) ) { stop("Files do not match in site names") }
        if( FALSE %in% (popCoordinates.n %in% colnames(popDifferentiation)) ) { stop("Files do not match in site names") }
        
        # -------------------------------------------------------
        # -------------------------------------------------------
        
        for( c in 1:nrow(combinations) ) {
          
          season <- combinations[c,1]
          pld.period <- combinations[c,2]
          
          if( season == "Spring" ) { spawn.p <- c(3,4,5) }
          if( season == "Summer" ) { spawn.p <- c(6,7,8) }
          if( season == "Autumn" ) { spawn.p <- c(9,10,11) }
          if( season == "Winter" ) { spawn.p <- c(12,1,2) }
          if( season == "YearRound" ) { spawn.p <- 1:12 }
          if( season == "" ) { spawn.p <- 1:12 }
          
          # -------------------
          # -------------------
          
          testFile <- paste0(results.folder,"/Genetics/",mainfileName,"/",speciesName,"/",markerName,"/",toupper(index),"/Simulation",season,str_pad(pld.period, 3, pad = "0"),"Days/resultMatrix.csv")
          if(file.exists(testFile)) { next }
          
          if( pld.period != 0) { source("5.2. Connectivity vs Differentiation.R") }
          
          if( pld.period == 0) { source("5.2. Connectivity vs Differentiation Stand.R") }
          
          # -------------------
          # -------------------
          
          if( c == 1) {
            
            # source("5.2. Distance vs Differentiation.R") 
            
          }
        }
      }
    }
  }
}

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# WHY NA???

library(plyr)
resultFiles <- list.files(paste0("../Results/Genetics"), recursive = TRUE, pattern="resultMatrix.csv", full.names=TRUE)
resultsDF <- data.frame()

for(i in 1:length(resultFiles)) {
  
  file.i <- resultFiles[i]
  file <- gsub("../Results/Genetics/Per species/","",file.i)
  file <- gsub("\\../Results/Genetics/","",file)
  separator <- unlist(gregexpr("/",file))
  
  study <- substr(file,1,separator[1]-1)
  species <- substr(file,separator[1]+1,separator[2]-1)
  marker <- substr(file,separator[2]+1,separator[3]-1)
  index <- substr(file,separator[3]+1,separator[4]-1)
  type <- substr(file,separator[4]+1,separator[5]-1)
  type <- ifelse(grepl("SimulationDistance",file),"Distance","OceanTransport")
  
  if(marker == "microsatellite" ) { marker <- "microsatellites" }
  if(marker == "nucDNA_ISSR" ) { marker <- "nucDNA" }
  if(marker == "nucDNA_SRAP" ) { marker <- "nucDNA" }
  
  pd <- 0
  
  if(type == "OceanTransport") {
    
    pd <- substr(file,separator[4]+1,separator[5]-1)
    pd <- gsub("Simulation","",pd)
    pd <- gsub("Days","",pd)
    pd <- as.numeric(pd)
    
  }
  
  resultsMatrix <- read.csv(file.i)
  
  resultsDF <- rbind.fill(resultsDF,
                          data.frame( Study=study,
                                      Species=species,
                                      Marker=marker,
                                      Index=index,
                                      Type=type,
                                      propaguleDuration=pd,
                                      resultsMatrix))
  
}

sort(unique(resultsDF$Species))
sort(unique(resultsDF$Marker))
sort(unique(resultsDF$Index))
write.csv(resultsDF, file="../Results/Genetics/mainData.csv")


## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)
closeAllConnections()

library(ggplot2)
library(raster)
library(rnaturalearth)
library(sf)
library(geosphere)
library(rgeos)
library(stringr)
library(readxl)

## ------------------

theme_map <- theme( text = element_text(family = "Helvetica", color = "#22211d"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_line(color = "black", size = 0.1), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#FFFFFF", color = NA), panel.background = element_rect(fill = "#FFFFFF", color = NA), panel.border = element_blank(), legend.background = element_rect(fill = "#FFFFFF", color = NA), legend.position="bottom", legend.box = "horizontal", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10), legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.75, 'cm') )

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=13) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)) ,
                   legend.title = element_blank() ,
                   legend.margin=margin(c(0.3,1,0.3,1), unit='lines') ,
                   legend.background = element_rect(fill="white", size=0.2, linetype="solid",  colour ="#979797"))

projection <- CRS("+proj=robin +over")
worldMap <- ne_countries(scale = 10, returnclass = "sp")

## ------------------

directories <- list.files("../Results/Genetics/")
directories <- directories[!grepl("csv",directories)]
length(directories) 

allrecords <- data.frame()

for( mainfileName in directories ) {
  
  gc(reset=TRUE)
  
  mainFiles <- paste0("../Data/geneticDataClean/",mainfileName,"/")
  mainFiles <- list.files(mainFiles, pattern="xls", full.names=TRUE)
  mainFiles <- mainFiles[!grepl("~",mainFiles)]
  
  mainFile <- data.frame()
  for( mainFiles.i in mainFiles ) {
    mainFile.i <- as.data.frame(read_excel(mainFiles.i))
    if( sum(c("species","lon","lat") %in% tolower(colnames(mainFile.i))) == 3) {
      mainFile <- plyr::rbind.fill(mainFile,mainFile.i)
    }
  }
  
  allrecords <- rbind(allrecords,mainFile[,c("Lon","Lat")])
  
  
}

allrecords <- allrecords[complete.cases(allrecords),]

## ------------------

coordinates(allrecords) <- ~Lon+Lat
crs(allrecords) <- crs(worldMap)

bb <- sf::st_union(sf::st_make_grid( st_bbox(c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)), n = 100))
bb <- st_transform(bb, projection)

allrecordsSp <- spTransform(allrecords, CRSobj = projection)
allrecords <- as.data.frame(allrecordsSp)

worldMap <- spTransform(worldMap, CRSobj = projection)
worldMap <- gBuffer(worldMap, byid=TRUE, width=0.001)
worldMap <- crop(worldMap, as(bb, "Spatial"))

# -------------------

# Define a color ramp for the oxygen gradient and plot the map.
myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414","#D278E4","#9914B3")

plot <- ggplot() + 
  geom_polygon(data = worldMap, aes(x = long, y = lat, group = group), fill="#CDCDCD", colour = "#CDCDCD" , size=0.25 ) +
  geom_point(data = allrecords, aes(x = Lon, y = Lat), colour = "#09612F",size=0.75) +
  geom_sf(data = bb,fill=NA, colour = "white" , linetype='solid', size= 3 ) +
  theme_map
plot

pdf(file=paste0("../Results/Genetics/mainResultsCoordinates.pdf"),width=16,height=8)
plot
dev.off()

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
