## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("../Project Config 0.R")
source("Dependences.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

## Test if connectivity exists

file.exists(paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimates.RData"))

## ------------------------------------
## Resolve connectivity

load(paste0(project.folder,"/Results/",project.name,"/InternalProc/","SourceSink.RData"))
load(paste0(project.folder,"/Results/",project.name,"/InternalProc/","Parameters.RData"))

cell.to.process <- unique(source.sink.xy$cells.id[source.sink.xy$source == 1])
n.particles.per.cell <- global.simulation.parameters$n.particles.per.cell
n.new.particles.per.day <- global.simulation.parameters$n.new.particles.per.day
n.steps.per.day <- global.simulation.parameters$n.hours.per.day

length(cell.to.process) * n.particles.per.cell
  
## ------------------

particles.reference.bm.desc <- dget( paste0(project.folder,"/Results/",project.name,"/InternalProc/particles.reference.desc"))

## ------------------

cl.2 <- makeCluster(10 , type="FORK")
registerDoParallel(cl.2)

all.connectivity.pairs.to.sql <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  particles.reference.bm.i <- attach.big.matrix(particles.reference.bm.desc)
  
  connectivity.temp.m <- particles.reference.bm.i[ mwhich(particles.reference.bm.i,2,list(cell.id.ref.f), list('eq')) , ]
  connectivity.temp.m <- data.frame(connectivity.temp.m)
  colnames(connectivity.temp.m) <- c("id","start.cell","start.year","start.month","start.day","pos.lon","pos.lat","pos.alt","state","t.start","t.finish","cell.rafted","ocean")
  
  connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$cell.rafted != 0 & connectivity.temp.m$state == 2,]
  connectivity.temp.m <- data.frame(connectivity.temp.m,travel.time= (1 + connectivity.temp.m$t.finish - connectivity.temp.m$t.start) / n.steps.per.day)
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

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

# -----------------------------------------

# Save pairs

save(all.connectivity.pairs.to.sql,file=paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimates.RData"))

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------

# Direct Overall Connectivity matrix (mean of all years)

# Subset Years
# Connectivity <- Connectivity[ Year == 2017, ]

Connectivity <- data.table(all.connectivity.pairs.to.sql)
Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity[is.na(Connectivity)] <- 0
Connectivity ; gc()

## -------

source.sink.id <- source.sink.xy$cells.id[which(source.sink.xy$source == 1)] 
plot(source.sink.xy[which(source.sink.xy$source == 1),2:3] )

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% source.sink.id,]
source.sink.bm <- as.big.matrix(as.matrix(source.sink.xy))
write.big.matrix(source.sink.bm, paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm"))

Connectivity <- Connectivity[Connectivity$Pair.from %in% source.sink.id & Connectivity$Pair.to %in% source.sink.id,]
Connectivity.bm <- as.big.matrix(as.matrix(Connectivity))
write.big.matrix(Connectivity.bm, paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveraged.bm")) 


## ------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------
## Assign connectivity estimates [source.sink site] to polygons [if the case]

Connectivity <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveraged.bm") )
Connectivity <- as.data.frame(Connectivity[,])
colnames(Connectivity) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events")
Connectivity

source.sink.xy <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm") )
source.sink.xy <- as.data.frame(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

additional.source.sink <- shapefile(additional.source.sink.shp)
additional.source.sink <- additional.source.sink[,"ID"]

Connectivity.Poly <- expand.grid(From=additional.source.sink$ID,To=additional.source.sink$ID)
nrow(Connectivity.Poly)

cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)

Connectivity.Poly.Vals <- foreach(i=1:nrow(Connectivity.Poly), .verbose=FALSE, .packages=c("sp","gdistance")) %dopar% { 
  
  From <- Connectivity.Poly[i,"From"]
  To <- Connectivity.Poly[i,"To"]
  
  polygon.from <- coordinates(additional.source.sink[additional.source.sink$ID == From,])
  polygon.to <- coordinates(additional.source.sink[additional.source.sink$ID == To,])
  
  distances.from <- spDistsN1(as.matrix(source.sink.xy[,c("Lon","Lat")]),polygon.from)
  distances.to <- spDistsN1(as.matrix(source.sink.xy[,c("Lon","Lat")]),polygon.to)
  
  distance.probability.i <- Connectivity[Connectivity$Pair.from %in% source.sink.xy[which(distances.from == min(distances.from)),"Pair"] & Connectivity$Pair.to %in% source.sink.xy[which(distances.to == min(distances.to)),"Pair"],]
  
  if(nrow(distance.probability.i) == 0 ) {
    return(NULL)
  }
  
  if(nrow(distance.probability.i) > 0 ) {
    distance.probability.i <- t(data.frame(apply(distance.probability.i,2,mean)))
    rownames(distance.probability.i) <- NULL
    distance.probability.i <- data.frame(Pair.from.old=distance.probability.i[1,1],Pair.to.old=distance.probability.i[1,2],distance.probability.i)
    distance.probability.i[1,3] <- From
    distance.probability.i[1,4] <- To
    rownames(distance.probability.i) <- NULL
    return(distance.probability.i)
  }

}

stopCluster(cl.3) ; rm(cl.3) ; gc()

Connectivity.Poly <- do.call(rbind,Connectivity.Poly.Vals)

## ---------------------

source.sink.xy <- data.frame(Pair=additional.source.sink[,"ID"],Lon=as.data.frame(gCentroid(additional.source.sink,byid=TRUE))[,1],Lat=as.data.frame(gCentroid(additional.source.sink,byid=TRUE))[,2],Source=1)
source.sink.bm <- as.big.matrix(as.matrix(source.sink.xy))
write.big.matrix(source.sink.bm, paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.Polys.bm"))

Connectivity.bm <- as.big.matrix(as.matrix(Connectivity.Poly[,-c(1,2)]))
write.big.matrix(Connectivity.bm, paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveragedPolys.bm")) 

## --------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------
## End of Code [!]
