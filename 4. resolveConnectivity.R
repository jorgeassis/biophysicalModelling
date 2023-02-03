## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Config.R")
source("Dependences/mainFunctions.R")

particles.reference.bm.desc <- paste0(results.folder,"/particleReferenceTable.desc")
load(paste0(results.folder,"/particleReferenceTableNames.RData"))

n.season <- "" # Spring; Summer; Autumn; Winter; "" for All

## ----------------------------------

if( n.season == "" ) { n.months <- 1:12 }
if( n.season == "Spring" ) { n.months <- 3:5 }
if( n.season == "Summer" ) { n.months <- 6:8 }
if( n.season == "Autumn" ) { n.months <- 9:11 }
if( n.season == "Winter" ) { n.months <- c(12,1,2) }

## ------------------------------------
## Resolve connectivity

load(paste0(results.folder,"/","sourceSinkSites.RData"))
load(paste0(results.folder,"/","modelParameters.RData"))

cell.to.process <- unique(source.sink.xy$cells.id)
cell.to.process <- unique(source.sink.xy$cells.id[source.sink.xy$source == 1])

n.particles.per.cell <- global.simulation.parameters$n.particles.per.cell
n.new.particles.per.day <- global.simulation.parameters$n.new.particles.per.day
n.steps.per.day <- global.simulation.parameters$n.hours.per.day

## ------------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

pairedConnectivity <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  particles.reference.bm.i <- attach.big.matrix(particles.reference.bm.desc)
  
  connectivity.temp.m <- particles.reference.bm.i[ mwhich(particles.reference.bm.i,2,list(cell.id.ref.f), list('eq')) , ]
  connectivity.temp.m <- data.frame(connectivity.temp.m)
  colnames(connectivity.temp.m) <- particles.reference.names
  
  connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$start.month %in% n.months,]
  connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$cell.rafted != 0 & connectivity.temp.m$state == 2,]
  connectivity.temp.m <- as.data.table(connectivity.temp.m)
  connectivity.pairs.to.sql <- data.frame()
  
  for(y in unique(connectivity.temp.m$start.year) ) {
    
    connectivity.temp <- connectivity.temp.m[ start.year == y , ]
    
    for( cell.id.ref.t in unique(connectivity.temp[ , cell.rafted ]) ) {
      
      connectivity.pairs.to.sql <- rbind(connectivity.pairs.to.sql,
                                         
                                         data.frame(  Pair.from = cell.id.ref.f,
                                                      Pair.to = cell.id.ref.t,
                                                      Number.events = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]),
                                                      Time.mean = mean(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.min = min(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.max = max(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.sd = sd(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / round( sum(source.sink.xy$cells.id == cell.id.ref.f) * n.particles.per.cell / length(unique(connectivity.temp.m$start.year)) ),
                                                      Year = y ) )
    }
    
  }
  
  connectivity.pairs.to.sql[is.na(connectivity.pairs.to.sql)] <- 0
  return( connectivity.pairs.to.sql )
  
}

stopCluster(cl.2) ; rm(cl.2)
closeAllConnections(); gc(reset=TRUE)

# -----------------------------------------
# Save pairs

file.remove( list.files(paste0(results.folder,"/"), full.names = TRUE, pattern = paste0("particlePairedConnectivity",n.season,"Table") ))

particles.connectivity.names <- colnames(pairedConnectivity)
save(particles.connectivity.names,file=paste0(results.folder,"/particlePairedConnectivity",n.season,"TableNames.RData"))

rownames(pairedConnectivity) <- NULL
colnames(pairedConnectivity) <- NULL

pairedConnectivity <- as.big.matrix(as.matrix(pairedConnectivity) , backingpath=results.folder , backingfile =  paste0("particlePairedConnectivity",n.season,"Table.bin"), descriptorfile =  paste0("particlePairedConnectivity",n.season,"Table.desc"))

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------

# Overall Connectivity matrix (mean of all years)

# Subset Years
# averagedConnectivity <- averagedConnectivity[ Year == 2017, ]

pairedConnectivity <- attach.big.matrix(paste0(results.folder,"/particlePairedConnectivity",n.season,"Table.desc"))
load(paste0(results.folder,"/particlePairedConnectivity",n.season,"TableNames.RData"))

averagedConnectivity <- as.data.table(pairedConnectivity[])
colnames(averagedConnectivity) <- particles.connectivity.names

particles.connectivity.names <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
save(particles.connectivity.names,file=paste0(results.folder,"/particlePairedConnectivityAveraged",n.season,"TableNames.RData"))

averagedConnectivity <- averagedConnectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
averagedConnectivity[is.na(averagedConnectivity)] <- 0
rownames(averagedConnectivity) <- NULL
colnames(averagedConnectivity) <- particles.connectivity.names
head(averagedConnectivity)

## ------------------

load(file = paste0(results.folder,"/","modelParameters.RData"))

global.simulation.parameters.res <- data.frame( global.simulation.parameters,
                                                meanProbability = mean(averagedConnectivity$Mean.Probability),
                                                sdProbability = sd(averagedConnectivity$Mean.Probability),
                                                rangeProbabilityMin = min(averagedConnectivity$Mean.Probability),
                                                rangeProbabilityMax = max(averagedConnectivity$Mean.Probability),
                                                meanTime = mean(averagedConnectivity$Mean.Time),
                                                sdTime = sd(averagedConnectivity$Mean.Time),
                                                meanEvents = mean(averagedConnectivity$Mean.events),
                                                sdEvents = sd(averagedConnectivity$Mean.events)
                                
)

save(global.simulation.parameters.res, file = paste0(results.folder,"/","modelSummaryResults.RData"))

## --------------------

file.remove( list.files(results.folder, full.names = TRUE, pattern = paste0("particlePairedConnectivityAveraged",n.season,"Table.bin") ) )
file.remove( list.files(results.folder, full.names = TRUE, pattern = paste0("particlePairedConnectivityAveraged",n.season,"Table.desc") ) )

averagedConnectivity <- as.big.matrix(as.matrix(averagedConnectivity) , backingpath=paste0(results.folder,"/") , backingfile = paste0("particlePairedConnectivityAveraged",n.season,"Table.bin") , descriptorfile = paste0("particlePairedConnectivityAveraged",n.season,"Table.desc") )

## ----------------------------------------------------------------------
##
## Assign connectivity estimates [source.sink site] to polygons [if the case]
##
## ----------------------------------------------------------------------

stop("Review code :: 000")

Connectivity <- read.big.matrix( paste0(results.folder,"/InternalProc/","connectivityEstimatesAveraged",n.season,".bm") )
Connectivity <- as.data.frame(Connectivity[,])
colnames(Connectivity) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events")
Connectivity

source.sink.xy <- read.big.matrix( paste0(results.folder,"/InternalProc/","source.sink.bm") )
source.sink.xy <- as.data.frame(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

additional.source.sink <- shapefile(maskSourceSinkSites)
additional.source.sink$ID <- additional.source.sink$id
additional.source.sink <- additional.source.sink[,"ID"]

# ---------------------------

closest.source.sink.sites <- coordinates(additional.source.sink)
closest.source.sink.sites <- spDists(as.matrix(source.sink.xy[,c("Lon","Lat")]),closest.source.sink.sites)
closest.source.sink.sites <- apply(closest.source.sink.sites,2,which.min)
closest.source.sink.sites <- source.sink.xy[closest.source.sink.sites,"Pair"]
additional.source.sink$sitesSourceSink <- closest.source.sink.sites
additional.source.sinkDF <- as.data.frame(additional.source.sink)

ConnectivityDT <- data.table(Connectivity)

# ---------------------------

cl.3 <- makeCluster(number.cores)
registerDoParallel(cl.3)

Connectivity.Poly.Vals <- foreach(siteID=additional.source.sink$ID, .verbose=FALSE, .packages=c("data.table")) %dopar% { 
  
  ConnectivityDT.siteID <- data.table()
  
  sourceSink.from <- additional.source.sinkDF[additional.source.sinkDF$ID == siteID , "sitesSourceSink"]
  ConnectivityDT.Polys.i <- ConnectivityDT[ Pair.from == sourceSink.from , ]
  ConnectivityDT.Polys.i[,Pair.from := siteID]
  
  for( sourceSink.to in ConnectivityDT.Polys.i[,Pair.to]) {
    
    siteID.to <- additional.source.sinkDF[additional.source.sinkDF$sitesSourceSink == sourceSink.to,"ID"]
    
    if(length(siteID.to) == 0) { next }
    
    ConnectivityDT.siteID.i <- do.call("rbind", replicate(length(siteID.to), ConnectivityDT.Polys.i[Pair.to == sourceSink.to,], simplify = FALSE))
    ConnectivityDT.siteID.i[,Pair.to := siteID.to]
    ConnectivityDT.siteID <- rbindlist(list(ConnectivityDT.siteID,ConnectivityDT.siteID.i))

  }
  
  return(ConnectivityDT.siteID)
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

Connectivity.Poly <- do.call(rbind,Connectivity.Poly.Vals)

## ---------------------

sum(!Connectivity.Poly$Pair.from %in% additional.source.sink$ID)
sum(!Connectivity.Poly$Pair.to %in% additional.source.sink$ID)
Connectivity.Poly[, lapply(.SD, mean), by=.(Pair.from,Pair.to)]

## ---------------------

source.sink.xy <- data.frame(Pair=additional.source.sink[,"ID"],Lon=as.data.frame(gCentroid(additional.source.sink,byid=TRUE))[,1],Lat=as.data.frame(gCentroid(additional.source.sink,byid=TRUE))[,2],Source=1)
source.sink.bm <- as.big.matrix(as.matrix(source.sink.xy))
write.big.matrix(source.sink.bm, paste0(results.folder,"/InternalProc/","source.sink.Polys.bm"))

Connectivity.bm <- as.big.matrix(as.matrix(Connectivity.Poly))
write.big.matrix(Connectivity.bm, paste0(results.folder,"/InternalProc/","connectivityEstimatesAveragedPolys",n.season,".bm")) 

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------
## End of Code [!]