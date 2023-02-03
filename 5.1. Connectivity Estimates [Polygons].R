## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

if( ! exists("pipeLiner") ) {
  
  closeAllConnections()
  rm(list=(ls()[ls()!="v"]))
  gc(reset=TRUE)
  source("0. Config.R")
  source("Dependences.R")

  list.dirs(path = paste0("../Results"), recursive = FALSE)
  season <- "YearRound" # c("YearRound","SeasonSummer","SeasonWinter")
  n.season <- "" # Spring; Summer; Autumn; Winter; "" for All
  spawn.p <- 1:12  # spawn.p <- c(6,7,8,9)
  pld.period <- 30
  c <- 30
  
}

# Review All [!!]

additional.source.sink <- shapefile(maskSourceSinkSites)
additional.source.sink$ID <- additional.source.sink$id
additional.source.sink <- additional.source.sink[,"ID"]

coordRef <- crs(additional.source.sink)

reference.file <- paste0(results.folder,"/InternalProc/","particles.reference.desc")
sorce.sink.cells.file <- paste0(results.folder,"/InternalProc/","source.sink.bm")
connectivity.file <- paste0(results.folder,"/InternalProc/","connectivityEstimatesAveraged",n.season,".bm")

## --------------------------

source.sink.xy <- read.big.matrix(sorce.sink.cells.file)
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

Connectivity <- read.big.matrix(connectivity.file)
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events")

source.sink.xy.sp <- source.sink.xy[,2:3]
coordinates(source.sink.xy.sp) <- ~Lon+Lat
crs(source.sink.xy.sp) <- coordRef

## ---------------------------------

closest.source.sink.sites <- coordinates(additional.source.sink)
closest.source.sink.sites <- spDists(as.matrix(source.sink.xy[,c("Lon","Lat")]),closest.source.sink.sites, longlat = TRUE)
closest.source.sink.sites <- apply(closest.source.sink.sites,1,which.min)
source.sink.xy$Poly <- closest.source.sink.sites

Connectivity$Poly.from <- sapply( Connectivity$Pair.from,function(x) { as.numeric(source.sink.xy[which(source.sink.xy[,"Pair"] == x),"Poly"]) })
Connectivity$Poly.to <- sapply( Connectivity$Pair.to,function(x) { as.numeric(source.sink.xy[which(source.sink.xy[,"Pair"] == x),"Poly"]) })

## ------------------------------------------------------------------------------------------------------------
## Read main sources

load(paste0(project.folder,"Results/InternalProc/","Parameters.RData"))

sim.extent <-unique(as.numeric(unlist(strsplit(global.simulation.parameters$extent, split=","))))
months <- unique(as.numeric(unlist(strsplit(global.simulation.parameters$sim.months , split=","))))
n.hours.per.day <- global.simulation.parameters$n.hours.per.day
n.particles.per.cell <- global.simulation.parameters$n.particles.per.cell
n.new.particles.per.day <- global.simulation.parameters$n.new.particles.per.day
n.steps.per.day <- global.simulation.parameters$n.hours.per.day

## ---------------------
## ---------------------

regionsOfInterest <- gBuffer(regionsOfInterest, byid=TRUE, width=0)
worldMap <- ne_countries(scale = 10, returnclass = "sp")

# Aggregate based on names

regionsOfInterest <- aggregate(regionsOfInterest, by = list(RegionFinal = regionsOfInterest$Region), FUN=mean, dissolve = TRUE, areaWeighted = FALSE)
regionsOfInterest$Area <- round(raster::area(regionsOfInterest) / 1000000,digits=3)

regionsOfInterest <- regionsOfInterest[,-which(names(regionsOfInterest) == "id")]
regionsOfInterest <- regionsOfInterest[,-which(names(regionsOfInterest) == "Region")]

write.csv(data.frame(regionsOfInterest),file=paste0("../Results/regionsFinalList.csv"))
writeOGR(obj=regionsOfInterest, dsn="../Results/", layer="regions", driver="ESRI Shapefile",overwrite_layer=TRUE) # this is in geographical projection

## --------------------------------------------------------------------
## --------------------------------------------------------------------
# Temporary Subset

subseter <- as.vector(extent(source.sink.xy.sp) + c(-0.1,0.1,-0.1,0.1)) # c(-11.25,37.85,29.75,46.25)

## --------------------------------------------------------------------
## --------------------------------------------------------------------

source.sink.xy <- source.sink.xy[source.sink.xy$Lon >= subseter[1] & source.sink.xy$Lon <= subseter[2] & source.sink.xy$Lat >= subseter[3] & source.sink.xy$Lat <= subseter[4], ]
plot(source.sink.xy[,2:3])

source.sink.xy.sp <- crop(source.sink.xy.sp,extent(subseter))
regionsOfInterest <- crop(regionsOfInterest,extent(subseter))
regionsOfInterest$ID <- 1:nrow(regionsOfInterest)
worldMap <- crop(worldMap,extent(subseter + c(-20,20,-20,20))  ) 

## -----------------

plot(worldMap , col="Black",border="Black")
plot(regionsOfInterest , col="Black",border="Black")

## ------------------------------------------------------------------------------------------------------------------------------
## Identify polygon id in source sink sites

distanceThreshold <- 1

cl.2 <- makeCluster(12 , type="FORK")
registerDoParallel(cl.2)

source.sink.xy.region <- foreach( source.sink.xy.i = 1:nrow(source.sink.xy) , .verbose=FALSE, .combine = rbind ,  .packages=c("rgeos","raster","geosphere","FNN","bigmemory")) %dopar% { # 
  
  distances <- gDistance(source.sink.xy.sp[source.sink.xy.i], regionsOfInterest,byid=TRUE)
  
  if(min(distances) < distanceThreshold) {
    regionsOfInteresti <- regionsOfInterest[which.min(distances),]$ID
  }
  if(min(distances) >= distanceThreshold) {
    regionsOfInteresti <- 0
  }
  return( data.frame( ID = regionsOfInteresti  ) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

plot(source.sink.xy[,.(Lon,Lat)])
points(source.sink.xy[which(source.sink.xy.region != 0),.(Lon,Lat)], col="blue", pch=19)
points(source.sink.xy[which(source.sink.xy.region == 0),.(Lon,Lat)], col="red", pch=19)

source.sink.xy <- cbind(source.sink.xy[,.(Pair,Lon,Lat)],source.sink.xy.region)
source.sink.xy <- source.sink.xy[source.sink.xy$ID != 0,]
head(source.sink.xy)
save(source.sink.xy,file="../Results/source.sink.xy.Rdata")
load("../Results/source.sink.xy.Rdata")

## ----------------------------------------------------

RegionNames <- regionsOfInterest$RegionFinal
RegionID <- regionsOfInterest$ID

## ----------------------------------------------------
## Revert regions to landmasses

# regionsOfInterest <- gIntersection(regionsOfInterest, gUnaryUnion(worldMap), byid = TRUE, id = sapply(regionsOfInterest@polygons, slot, name = "ID"))
# regionsOfInterest$RegionFinal <- RegionNames
# regionsOfInterest$ID <- RegionID
# regionsOfInterest
# plot(regionsOfInterest)

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

combResults <- data.frame()

if( ! exists("isolatedResults")) {
  
  if( ! exists("combinations")) { combinations <- matrix(NA) }
  
  isolatedResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
  rownames(isolatedResults) <- RegionNames
  colnames(isolatedResults) <- 1:nrow(combinations)
  betweennessResults <- higherBetweennessResults <- eighenCentralityResults <- highereighenCentralityResults <- closenessResults <- higherclosenessResults <- clusterAssignment <- resistanceResults <- higherResistanceResults <- outDegreeResults <- higherOutDegreeResults <- selfRecruitmentResults <- higherSelfRecruitmentResults <- isolatedResults

}

## ------------------------------------------------------------------------------

dev.off()
gc(reset=TRUE)

project.name.c <- paste0(project.name,"/",season,"_Pld",pld.period)

## ----------------------------------------------------

if( ! dir.exists(paste0("../Results/",project.name.c)) ) { dir.create(file.path(paste0("../Results/",project.name.c)), showWarnings = FALSE) } 
if( ! dir.exists(paste0("../Results/",project.name.c,"/Data")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Data")), showWarnings = FALSE) } 
if( ! dir.exists(paste0("../Results/",project.name.c,"/Maps")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Maps")), showWarnings = FALSE) } 
if( ! dir.exists(paste0("../Results/",project.name.c,"/Networks")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Networks")), showWarnings = FALSE) } 

if( ! exists("doParallelCalculations") ) {
  
  particles.reference.bm.desc <- dget( paste0(project.folder,"Results/",project.name,"/InternalProc/particles.reference.desc"))
  cell.to.process <- unique(source.sink.xy$Pair)
  
  cl.2 <- makeCluster(number.cores, type="FORK" )
  registerDoParallel(cl.2)
  
  connectivity.source.sink.xy <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
    
    particles.reference.bm.i <- attach.big.matrix(particles.reference.bm.desc)
    
    connectivity.temp.m <- particles.reference.bm.i[ mwhich(particles.reference.bm.i,2,list(cell.id.ref.f), list('eq')) , ]
    connectivity.temp.m <- data.frame(connectivity.temp.m)
    colnames(connectivity.temp.m) <- c("id","start.cell","start.year","start.month","start.day","pos.lon","pos.lat","pos.alt","state","t.start","t.finish","cell.rafted")
    
    connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$cell.rafted != 0 & connectivity.temp.m$state == 2,]
    connectivity.temp.m <- data.frame(connectivity.temp.m,travel.time= (1 + connectivity.temp.m$t.finish - connectivity.temp.m$t.start) / n.steps.per.day)
    connectivity.temp.m <- as.data.table(connectivity.temp.m)
    
    connectivity.temp.m <- connectivity.temp.m[start.month %in% spawn.p , ]
    connectivity.temp.m <- connectivity.temp.m[travel.time <= pld.period , ]
    
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
                                                        Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / round( n.particles.per.cell / length(unique(connectivity.temp.m$start.year)) ),
                                                        Year = y ) )
      }
      
    }
    
    connectivity.pairs.to.sql[is.na(connectivity.pairs.to.sql)] <- 0
    
    return( connectivity.pairs.to.sql )
    
  }
  
  stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)
  
  connectivity.source.sink.xy <- data.table(connectivity.source.sink.xy[,])
  connectivity.source.sink.xy <- connectivity.source.sink.xy[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
  colnames(connectivity.source.sink.xy) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
  connectivity.source.sink.xy[is.na(connectivity.source.sink.xy)] <- 0
  connectivity.source.sink.xy ; gc()
  save(connectivity.source.sink.xy,file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.xy.Rdata"))
  
  ## ----------------------------------------------------
  ## Add missing connections
  
  missingFrom <- source.sink.xy$Pair[which(! source.sink.xy$Pair %in% unique(connectivity.source.sink.xy$Pair.from))]
  
  if(length(missingFrom) > 0) { 
    
    connectivity.source.sink.xy <- rbind(connectivity.source.sink.xy,data.frame(Pair.from=missingFrom,Pair.to=missingFrom),fill=TRUE)
    connectivity.source.sink.xy[is.na(connectivity.source.sink.xy)] <- 0
    
  }
  
  missingTo <- source.sink.xy$Pair[which(! source.sink.xy$Pair %in% unique(connectivity.source.sink.xy$Pair.to))]
  
  if(length(missingTo) > 0) { 
    
    connectivity.source.sink.xy <- rbind(connectivity.source.sink.xy,data.frame(Pair.from=missingTo,Pair.to=missingTo),fill=TRUE)
    connectivity.source.sink.xy[is.na(connectivity.source.sink.xy)] <- 0
    
  }
  
  # ---------------------------
  # ---------------------------
  
  regionsOfInterestPairs <- expand.grid(from=regionsOfInterest$ID,to=regionsOfInterest$ID)
  polygonsCompute <- regionsOfInterest$ID
  
  cl.2 <- makeCluster(number.cores , type="FORK")
  registerDoParallel(cl.2)
  
  connectivity <- foreach( pairs = 1:nrow(regionsOfInterestPairs) , .verbose=FALSE, .combine = rbind ,  .packages=c("geosphere","rgeos","raster","data.table","FNN","bigmemory")) %dopar% { # 
    
    pair.id.1 <- regionsOfInterestPairs[pairs,1]
    pair.id.2 <- regionsOfInterestPairs[pairs,2]
    
    pair.cells.1 <- as.vector(unlist(source.sink.xy[ ID == pair.id.1 , 1 ]))
    pair.cells.2 <- as.vector(unlist(source.sink.xy[ ID == pair.id.2 , 1 ]))
    
    temp.result <- connectivity.source.sink.xy[ Pair.from %in% pair.cells.1 & Pair.to %in% pair.cells.2 , ]
    
    if(nrow(temp.result) == 0) {
      
      connectivity.pairs <- data.frame(  Pair.from = pair.id.1,
                                         Pair.to = pair.id.2,
                                         Number.events = 0,
                                         Time.mean = 0,
                                         Time.max = 0,
                                         Time.sd = 0,
                                         Probability = 0 )
      
    }
    
    if(nrow(temp.result) > 0) {
      
      connectivity.pairs <- data.frame(  Pair.from = pair.id.1,
                                         Pair.to = pair.id.2,
                                         Number.events = mean(temp.result$Mean.events),
                                         Time.mean = mean(temp.result$Mean.Time),
                                         Time.max = mean(temp.result$Max.Time),
                                         Time.sd = mean(temp.result$SD.Time),
                                         Probability = mean(temp.result$Mean.Probability) )
      
    }
    
    return( connectivity.pairs )
    
  }
  
  stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

  ## --------------------
  
  save(connectivity,file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.Polygons.Rdata"))
  
}

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

load(file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.Polygons.Rdata"))
load(file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.xy.Rdata"))

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# -------------
# Correction

SELtoMAD <- connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "SEL") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "MAD"),3:7]
MADtoSEL <- connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "MAD") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "SEL"),3:7]

SELtoCAN <- connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "SEL") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "CAN"),3:7]
CANtoSEL <- connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "CAN") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "SEL"),3:7]

connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "SEL") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "MAD"),3:7] <- SELtoCAN
connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "MAD") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "SEL"),3:7] <- CANtoSEL

connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "SEL") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "CAN"),3:7] <- SELtoMAD
connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "CAN") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "SEL"),3:7] <- MADtoSEL

connectivity[connectivity[,1] == which(regionsOfInterest$RegionFinal == "MAD") & connectivity[,2] == which(regionsOfInterest$RegionFinal == "CAN"),3:7] <- c(1.5,25.03965,25.55986,0.6747102,0.001958454)

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

connectivity.matrix <- connectivity[ ,c(1,2,7)]
connectivity.matrix[is.na(connectivity.matrix)] <- 0
connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )

## ----------------------------------

colnames(connectivity.matrix) <- RegionNames
rownames(connectivity.matrix) <- RegionNames

## ----------------------------------

retention <- diag(connectivity.matrix)
sumRows <- apply(connectivity.matrix,1,sum,na.rm=T)
retention[ is.na(retention) ] <- 0
selfRecruitment <- retention / sumRows
selfRecruitment[ is.na(selfRecruitment) ] <- 0

selfRecruitmentResults[ sapply(names(selfRecruitment),function(x) which( rownames(selfRecruitmentResults) == x)) ,c] <- selfRecruitment
write.csv(selfRecruitmentResults,file="../Results/selfRecruitment.csv")

higherSelfRecruitment.calc <- which(selfRecruitment >=  as.numeric(quantile(selfRecruitment,0.95,na.rm=TRUE)))
higherSelfRecruitmentResults[which(rownames(higherSelfRecruitmentResults) %in% names(higherSelfRecruitment.calc)),c] <- 1
write.csv(higherSelfRecruitmentResults,file="../Results/higherSelfRecruitment.csv")

diag(connectivity.matrix) <- 0

isolated.sourceSink <- which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)
isolatedResults[ sapply(names(isolated.sourceSink),function(x) which( rownames(isolatedResults) == x)),c] <- 1
write.csv(isolatedResults,file="../Results/isolated.csv")

# Aggregation level (Proportion of non-isolated MPAs, at least one connection, in relation to the number of MPAs)
aggregationAtLeastOne <- (nrow(connectivity.matrix)-length(isolated.sourceSink)) / nrow(connectivity.matrix)

# Aggregation level (Based on overall connections)
connectivity.matrix.binomial <- connectivity.matrix
connectivity.matrix.binomial[connectivity.matrix.binomial != 0] <- 1
aggregationAllConnections <- sum ( ( apply(connectivity.matrix.binomial,1,sum,na.rm=T) + 1 ) / nrow(connectivity.matrix.binomial) ) / nrow(connectivity.matrix.binomial)

# Resistance index
resistance <- sapply( 1:nrow(connectivity.matrix.binomial) , function(res) { 1- length(unique(c(which(connectivity.matrix.binomial[res,] == 1 ) ,  which(connectivity.matrix.binomial[,res] == 1 )))) / nrow(connectivity.matrix.binomial) } )
names(resistance) <- colnames(connectivity.matrix.binomial)

resistanceResults[ sapply(names(resistance),function(x) which( rownames(resistanceResults) == x)) ,c] <- resistance
write.csv(resistanceResults,file="../Results/resistance.csv")

# Those with higher / Percertil 95%
higherResistance.calc <- which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=TRUE)))
higherResistanceResults[which(rownames(higherResistanceResults) %in% names(higherResistance.calc)),c] <- 1
write.csv(higherResistanceResults,file="../Results/higherSelfRecruitment.csv")

# Average connections between MPAs
connect.index <- data.frame(SSSites=colnames(connectivity.matrix),exportTo=apply(connectivity.matrix,1,function(x){ sum(x != 0,na.rm=T) } ) , importFrom=apply(connectivity.matrix,2,function(x){ sum(x != 0,na.rm=T) } ))

averageConnections <- mean(unlist(connect.index[,-1]),na.rm=T)
sdConnections <- sd(unlist(connect.index[,-1]),na.rm=T)
maximumConnections <- max(apply(connect.index[,-1],1,max))

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

comb <- connectivity
comb <- comb[ comb[,1] != comb[,2] ,]
comb <- as.data.frame( comb[ sort(comb[,"Probability"] , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Probability")] )

## --------------------------------------------------------------------------------

# GGPLOT with connections

theme_map <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
    panel.grid.major = element_line(color = "#979797", size = 0.05),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank()
  )

# extent(worldMap) + c(20,-20,20,-20)

mapRegion <- ggplot() +
  geom_polygon(data = worldMap , fill = "#C4C4C4", colour = "#ffffff" , size=0.25 ,  aes(long, lat, group = group)) +
  coord_map('lambert', lat0=29.75, lat1=45, xlim=c(-30, 10), ylim=c(10, 52)) + theme_map

mapRegionNet <- mapRegion
connected.pairs <- comb[comb[,3] > 0,]
centroids <- as.data.frame(gCentroid(regionsOfInterest,byid=TRUE),xy=T)
colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
for( i in nrow(connected.pairs):1 ){
  strenght <- (connected.pairs[i,3] * 100) + 1 
  routes_sl <- gcIntermediate(centroids[ which(regionsOfInterest$ID == connected.pairs[i,1]),],
                              centroids[ which(regionsOfInterest$ID == connected.pairs[i,2]),],
                              n = 100, addStartEnd = TRUE, sp = TRUE)
  SLDF = sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = NA), match.ID = F)
  mapRegionNet <- mapRegionNet + geom_path(data = SLDF, size=0.35 , aes(x = long, y = lat), col = "#797979") # colfunc(101)[round(strenght)]
}

mapRegionNet

## --------------------------------------------------------------------------------

graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
#E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
E(graph.obj)$weight = comb[,3] # Hock, Karlo Mumby, Peter J 2015
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
# graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
graph.obj <- simplify(graph.obj)

V(graph.obj)$name <- sapply(as.numeric(names(V(graph.obj))),function(x) { regionsOfInterest$RegionFinal[regionsOfInterest$ID == x][1] })
  
# -------------
# -------------

graph.obj.ss <- graph.obj
E(graph.obj.ss)$weight <- 1 - E(graph.obj.ss)$weight
connectivity.matrix.ss <- connectivity.matrix
for(i in 1:nrow(connectivity.matrix.ss)) {
  for(j in 1:nrow(connectivity.matrix.ss)) {
    connectivity.matrix.ss[i,j] <- distances(graph.obj.ss, rownames(connectivity.matrix.ss)[i], colnames(connectivity.matrix.ss)[j])
  }
}
connectivity.matrix.ss[connectivity.matrix.ss == Inf] <- 1

unLinked <- which(apply(connectivity.matrix.ss,1,sum) == nrow(connectivity.matrix.ss) - 1 & apply(connectivity.matrix.ss,2,sum) == nrow(connectivity.matrix.ss) - 1)
if(length(unLinked) > 0 ) { connectivity.matrix.ss <- connectivity.matrix.ss[-unLinked,-unLinked] }

# -------------
# -------------
# Make dendogram 

dist <- as.dist(connectivity.matrix.ss)
write.table(as.matrix(dist),file=paste0("../Results/",project.name.c,"SteppingStoneDistance.csv"),row.names = TRUE,col.names = TRUE)

hc <- hclust(dist, method = "average") # average

library(factoextra)
library(dendextend)
fviz_nbclust(connectivity.matrix.ss, FUNcluster= hcut, method = c("silhouette"), k.max= 8) # silhouette wss gap_stat
fviz_nbclust(connectivity.matrix.ss[which(!rownames(connectivity.matrix.ss) %in% c("Portugal","Asturias","France")) , which(!colnames(connectivity.matrix.ss) %in% c("Portugal","Asturias","France")) ], FUNcluster= hcut, method = c("silhouette") , k.max= 5) # wss gap_stat

hcPlot <- as.dendrogram(hc)
plot(hcPlot, cex = 0.6, main = "Title", horiz = TRUE)
rect.dendrogram(hcPlot, k=3, border = 2, horiz =T)

pdf(file=paste0("../Results/",project.name.c,"dendogram.pdf"), width=12)
plot(hcPlot, cex = 0.6, main = "Title", horiz = TRUE)
dev.off()

pdf(file=paste0("../Results/",project.name.c,"dendogramClust.pdf"), width=12)
plot(hcPlot, cex = 0.6, main = "Title", horiz = TRUE)
rect.dendrogram(hcPlot, k=3,cluster=c(1,1,1,1,1,2,2,3,3,3,1), border = 2, horiz =T)
dev.off()

# -------------
# -------------

# Membership based on clusters
membership.graph <- clusters(graph.obj)$membership

clusterAssignment[ sapply(names(membership.graph),function(x) which(row.names(clusterAssignment) == x) ),c] <- membership.graph
write.csv(clusterAssignment,file="../Results/clusterAssignment.csv")

# Number of clusters
numberClusters <- length(unique(membership.graph)) - length(isolated.sourceSink)
aggregationBasedOnClusters <- 1 - ( numberClusters / length(membership.graph) )

## ------------------------------------------

# Plot clusters

reducedNames <- names(V(graph.obj))
reducedNames

l <- coordinates(regionsOfInterest)[sapply(V(graph.obj)$name, function(x) { which(regionsOfInterest$RegionFinal == x) }),]
l[which(reducedNames == "SEL"),] <- c(-16.5 , 30)
l[which(reducedNames == "MAD"),] <- c(-19 , 35)
l[which(reducedNames == "CAN"),] <- c(-19 ,24)
l[which(reducedNames == "MED"),] <- c(-3 ,35)
l[which(reducedNames == "CAD"),] <- c(-7 ,34)

lineThickness <- (E(graph.obj)$weight - min(E(graph.obj)$weight) ) / max((E(graph.obj)$weight - min(E(graph.obj)$weight) ))
lineThickness <- (lineThickness + 0.25) * 2

arrowSize <- 0.3
labelDistance <- 1.8

# connections

cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
cols.to.use[which(names(membership.graph) %in% names(isolated.sourceSink))] <- "white"

pdf(file=paste0("../Results/",project.name.c,"/Networks/connections.pdf"), width=12)
plot(graph.obj,vertex.label.dist=labelDistance,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , vertex.color=cols.to.use , layout=l , edge.width=lineThickness , edge.arrow.width=arrowSize)
dev.off()

# clustering

membership.graph <- cluster_leading_eigen(graph.obj)$membership
names(membership.graph) <- V(graph.obj)$name
cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
cols.to.use[which(names(membership.graph) %in% names(isolated.sourceSink))] <- "white"

pdf(file=paste0("../Results/",project.name.c,"/Networks/clusteringConnections.pdf"), width=12)
plot(graph.obj,vertex.label.dist=labelDistance,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , vertex.color=cols.to.use , layout=l , edge.width=lineThickness , edge.arrow.width=arrowSize)
dev.off()

lou <- cluster_leading_eigen(graph.obj)
pdf(file=paste0("../Results/",project.name.c,"/Networks/clusteringConnections.pdf"), width=12)
plot(lou, graph.obj, vertex.label.dist=labelDistance,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , vertex.color=cols.to.use , layout=l , edge.width=lineThickness, edge.color="#7D7D7D", edge.arrow.width=arrowSize)
dev.off()

betweennessIndex <- betweenness(graph.obj)
vertexSize <- (betweennessIndex - min(betweennessIndex) ) / max((betweennessIndex - min(betweennessIndex) ))
vertexSize <- (vertexSize + 0.25) * 20

pdf(file=paste0("../Results/",project.name.c,"/Networks/betweenness.pdf"), width=12)
plot(graph.obj,vertex.label.dist=labelDistance,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=vertexSize,edge.curved = F , vertex.color=cols.to.use , layout=l , edge.width=lineThickness, edge.color="#7D7D7D", edge.arrow.width=arrowSize)
dev.off()
write.csv(data.frame(names=names(betweennessIndex),index=(betweennessIndex - min(betweennessIndex) ) / max((betweennessIndex - min(betweennessIndex) )),row.names = NULL),file=paste0("../Results/",project.name.c,"/Networks/betweenness.csv"))

degreeIndex <- degree(graph.obj)
vertexSize <- (degreeIndex - min(degreeIndex) ) / max((degreeIndex - min(degreeIndex) ))
vertexSize <- (vertexSize + 0.25) * 20

pdf(file=paste0("../Results/",project.name.c,"/Networks/degree.pdf"), width=12)
plot(graph.obj,vertex.label.dist=labelDistance,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=vertexSize,edge.curved = F , vertex.color=cols.to.use , layout=l , edge.width=lineThickness, edge.color="#7D7D7D", edge.arrow.width=arrowSize)
dev.off()
write.csv(data.frame(names=names(degreeIndex),index=(degreeIndex - min(degreeIndex) ) / max((degreeIndex - min(degreeIndex) )),row.names = NULL),file=paste0("../Results/",project.name.c,"/Networks/degree.csv"))

closenessIndex <- closeness(graph.obj)
vertexSize <- (closenessIndex - min(closenessIndex) ) / max((closenessIndex - min(closenessIndex) ))
vertexSize <- (vertexSize + 0.25) * 20

pdf(file=paste0("../Results/",project.name.c,"/Networks/closeness.pdf"), width=12)
plot(graph.obj,vertex.label.dist=labelDistance,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=vertexSize,edge.curved = F , vertex.color=cols.to.use , layout=l , edge.width=lineThickness, edge.color="#7D7D7D", edge.arrow.width=arrowSize)
dev.off()
write.csv(data.frame(names=names(closenessIndex),index=(closenessIndex - min(closenessIndex) ) / max((closenessIndex - min(closenessIndex) )),row.names = NULL),file=paste0("../Results/",project.name.c,"/Networks/closeness.csv"))

## ------------------------------------------

# Plot clusters with connections

cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph[sapply(regionsOfInterest$RegionFinal,function(x) which(names(membership.graph) == x) )]]

pdf(file=paste0("../Results/",project.name,"/Maps/MapClusteringConnections.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroids ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use, size = 2.5, stroke = 0.35, alpha = 0.7)
)
dev.off()

## ------------------------------------------

# Plot isolated with connections

pdf(file=paste0("../Results/",project.name,"/Maps/MapIsolatedConnMPA.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroids ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.7) +
    geom_point(data = centroids[isolated.sourceSink,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "#9C2323", size = 2.5, stroke = 0.35, alpha = 0.9)
)
dev.off()

## ------------------------------------------

betweennessIndex <- betweenness(graph.obj)
betweennessIndex <- betweennessIndex[sort(as.numeric(names(betweennessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
betweennessIndex <- (betweennessIndex - min(betweennessIndex)) / (max(betweennessIndex) - min(betweennessIndex)) 
averageBetweenness <- mean(betweennessIndex)
sdBetweenness <- sd(betweennessIndex)

# eigen_centralityIndex <- eigen_centrality(graph.obj, directed = FALSE, scale = FALSE)$vector
# eigen_centralityIndex <- eigen_centralityIndex[sort(as.numeric(names(eigen_centralityIndex)),decreasing = FALSE,index.return=TRUE)$ix]
# averageeigen_centrality <- mean(eigen_centralityIndex)
# 
# closenessIndex <- closeness(graph.obj, mode="all")
# closenessIndex <- closenessIndex[sort(as.numeric(names(closenessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
# averagecloseness <- mean(closenessIndex)

betweennessQ95 <- rep(0,length(membership.graph))
betweennessQ95[which(betweennessIndex >= quantile(betweennessIndex,probs=0.95))] <- 1

betweennessResults[ ,c] <- betweennessIndex
higherBetweennessResults[ which(betweennessQ95 == 1),c] <- "TRUE"
write.csv(betweennessResults,file="../Results/betweennessMPAs.csv")
write.csv(higherBetweennessResults,file="../Results/higherBetweennessMPAs.csv")

outDegree <- strength(graph.obj, mode = c("out"))
outDegree <- outDegree[sort(as.numeric(names(outDegree)),decreasing = FALSE,index.return=TRUE)$ix]
outDegree <- (outDegree - min(outDegree)) / (max(outDegree) - min(outDegree)) 

averageOutDegree <- mean(outDegree)
sdOutDegree <- sd(outDegree)

outDegreeQ95 <- rep(0,length(membership.graph))
outDegreeQ95[which(outDegree >= quantile(outDegree,probs=0.95))] <- 1

outDegreeResults[ ,c] <- outDegree
higherOutDegreeResults[ which(outDegreeQ95 == 1),c] <- "TRUE"

write.csv(outDegreeResults,file="../Results/outDegreeMPAs.csv")
write.csv(higherOutDegreeResults,file="../Results/higherOutDegreeMPAs.csv")

# Plot centrality indexes with Connections and clusters

betweennessIndexPlot <- as.numeric(betweennessResults[,c])
betweennessIndexPlot <- (betweennessIndexPlot * 2) + 2.5

outDegreeIndexPlot <- as.numeric(outDegreeResults[,c])
outDegreeIndexPlot <- (outDegreeIndexPlot * 2) + 2.5

cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(betweennessIndexPlot),breaks = 20))]
cols.to.use[isolated.poly] <- "white"

pdf(file=paste0("../Results/",project.name,"/Maps/MapBetweennessConnections.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = betweennessIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = betweennessIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[betweennessQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "#9C2323", size = betweennessIndexPlot[betweennessQ95 == 1], stroke = 1.2, alpha = 0.7)
  + labs(caption = paste0("95th: ", round(quantile(betweennessIndex,probs=0.95),2)))
)
dev.off()

cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(outDegreeIndexPlot),breaks = 20))]
cols.to.use[isolated.poly] <- "white"

pdf(file=paste0("../Results/",project.name,"/Maps/MapStrengthOut.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = outDegreeIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = outDegreeIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[outDegreeQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[outDegreeQ95 == 1], size = outDegreeIndexPlot[outDegreeQ95 == 1], stroke = 1.2, alpha = 0.7)
  + labs(caption = paste0("95th: ", round(quantile(outDegree,probs=0.95),2)))
)
dev.off()

## ------------------------------------------

# Plot SR with connections

SRIndexPlot <- (SR / max(SR) )
SRIndexPlot <- (SRIndexPlot * 2) + 2.5
SRIndexPlot[ which(cols.to.use == "white") ] <- 2.5
higherSelfRecruitment <- higherSelfRecruitment[higherSelfRecruitment %in% which(cols.to.use != "white")]

cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(SRIndexPlot),breaks = 20))]
cols.to.use[isolated.poly] <- "white"

pdf(file=paste0("../Results/",project.name,"/Maps/MapSRConnections.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = SRIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = SRIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[higherSelfRecruitment,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[higherSelfRecruitment], size = SRIndexPlot[higherSelfRecruitment], stroke = 1.2, alpha = 0.7)
)
dev.off()

# ----------------------------------
# centrality vs area vs shape (corrected perimeter-area ratio; https://doi.org/10.1016/j.gecco.2018.e00504)

shapeEffectBetweenness <- FALSE
shapeEffectOutDegree <- FALSE
areaEffectBetweenness <- FALSE
areaEffectOutDegree <- FALSE

shapeIndex <- perimeter(notakeMPA) / sqrt( 4*pi*area(notakeMPA))

data <- data.frame(Type=ifelse(betweennessQ95 == 1 , "Hub" , "non-Hub"),data=shapeIndex)
data <- data[ data$data != max(data$data) ,]

pdf(file=paste0("../Results/",project.name,"/shapeEffectBetweenness.pdf"), width=10)

par(mfrow=c(1,1))
par(mar = c(5, 5, 3, 3))
boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
        main="Effect of shape index in betweenness centrality",
        names = c("Hub", "non-Hub"),
        las = 2,
        bg="#9A9A9A",
        col = c("#9A9A9A","#D4D4D4"),
        border = "#505050",
        horizontal = TRUE ,
        medcol = c("#9A9A9A","#D4D4D4"),
        xlab = "Marine reserve shape index" )

points(mean(data[data$Type == "Hub",2]),1,pch=19)
points(mean(data[data$Type != "Hub",2]),2,pch=19)

if(length(unique(data$Type)) > 1) {
  text(pos=2,max(data[,2]),2.5,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value, digits=4)))
  shapeEffectBetweenness <- kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value < 0.05
}

dev.off()

data <- data.frame(Type=ifelse(outDegreeQ95 == 1 , "Hub" , "non-Hub"),data=shapeIndex)
data <- data[ data$data != max(data$data) ,]

pdf(file=paste0("../Results/",project.name,"/shapeEffectOutDegree.pdf"), width=10)

par(mfrow=c(1,1))
par(mar = c(5, 5, 3, 3))
boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
        main="Effect of shape index in out degree",
        names = c("Hub", "non-Hub"),
        las = 2,
        bg="#9A9A9A",
        col = c("#9A9A9A","#D4D4D4"),
        border = "#505050",
        horizontal = TRUE ,
        medcol = c("#9A9A9A","#D4D4D4"),
        xlab = "Marine reserve shape index" )

points(mean(data[data$Type == "Hub",2]),1,pch=19)
points(mean(data[data$Type != "Hub",2]),2,pch=19)

if(length(unique(data$Type)) > 1) {
  text(pos=2,max(data[,2]),2.5,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value, digits=4)))
  shapeEffectOutDegree <- kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value < 0.05
}
dev.off()

# --------------

areaIndex <- area(notakeMPA)/1000000

data <- data.frame(Type=ifelse(betweennessQ95 == 1 , "Hub" , "non-Hub"),data=areaIndex)
data <- data[ data$data != max(data$data) ,]

pdf(file=paste0("../Results/",project.name,"/areaEffectBetweenness.pdf"), width=10)

par(mfrow=c(1,1))
par(mar = c(5, 5, 3, 3))
boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
        main="Effect of area in betweenness centrality",
        names = c("Hub", "non-Hub"),
        las = 2,
        bg="#9A9A9A",
        col = c("#9A9A9A","#D4D4D4"),
        border = "#505050",
        horizontal = TRUE ,
        medcol = c("#9A9A9A","#D4D4D4"),
        xlab = "Marine reserve area (km2)" )

points(mean(data[data$Type == "Hub",2]),1,pch=19)
points(mean(data[data$Type != "Hub",2]),2,pch=19)

if(length(unique(data$Type)) > 1) {
  text(pos=2,max(data[,2]),2.5,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value, digits=4)))
  areaEffectBetweenness <- kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value < 0.05
}

dev.off()

data <- data.frame(Type=ifelse(outDegreeQ95 == 1 , "Hub" , "non-Hub"),data=areaIndex)
data <- data[ data$data != max(data$data) ,]

pdf(file=paste0("../Results/",project.name,"/areaEffectOutDegree.pdf"), width=10)

par(mfrow=c(1,1))
par(mar = c(5, 5, 3, 3))
boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
        main="Effect of area in out degree",
        names = c("Hub", "non-Hub"),
        las = 2,
        bg="#9A9A9A",
        col = c("#9A9A9A","#D4D4D4"),
        border = "#505050",
        horizontal = TRUE ,
        medcol = c("#9A9A9A","#D4D4D4"),
        xlab = "Marine reserve area" )

points(mean(data[data$Type == "Hub",2]),1,pch=15)
points(mean(data[data$Type != "Hub",2]),2,pch=15)

if(length(unique(data$Type)) > 1) {
  text(pos=2,max(data[,2]),2.5,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value, digits=4)))
  areaEffectOutDegree <- kruskal.test(data[,2],as.numeric(as.factor(data$Type)))$p.value < 0.05
}

dev.off()

# ----------------------------------

combResults <- rbind(combResults,
                     data.frame(pld=pld.period,
                                n.isolated.poly=length(isolated.poly),
                                aggregationAtLeastOne=aggregationAtLeastOne,
                                aggregationAllConnections=aggregationAllConnections,
                                averageConnections=averageConnections,
                                sdConnections=sdConnections,
                                maximumConnections=maximumConnections,
                                numberClusters=numberClusters,
                                aggregationBasedOnClusters=aggregationBasedOnClusters,
                                averageBetweenness=averageBetweenness,
                                sdBetweenness=sdBetweenness,
                                averageOutDegree=averageOutDegree,
                                sdOutDegree=sdOutDegree,
                                shapeEffectOutDegree=shapeEffectOutDegree,
                                shapeEffectBetweenness=shapeEffectBetweenness,
                                areaEffectBetweenness=areaEffectBetweenness,
                                areaEffectOutDegree=areaEffectOutDegree))

write.csv(combResults,file="../Results/Results.csv")
save(combResults,file=paste0("../Results/allPLDResults.Rdata"))

# list.memory()
rm(connectivity.source.sink.xy )

## ---------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------

if( exists("pipeLiner") ) {
  
  names(combResults)
  
  x <- combResults$pld
  x.lab <- "Propagule duration (day)"
  y <- combResults$n.isolated.poly
  y.lab <- "Isolation degree (number of reserves)"
  
  par(mar = c(4.5, 5.5, 4.5, 4.5))
  plot(x,y,pch=20,col="#A6A6A6", ylab="",xlab=x.lab,axes=FALSE)
  title(ylab=y.lab, line=4)
  lines(bezierCurve(x,y,100)$x,bezierCurve(x,y,100)$y,type="l", lwd=1, lty=2)
  axis(2,las=2,col="White",col.ticks="Black", cex.axis=0.9)
  axis(1,las=0,col="White",col.ticks="Black", cex.axis=0.9)
  box()
  
  names(combResults)
  
  pdf(file=paste0("../Results/degree centrality.pdf"), width=12)
  
  ggplot(combResults, aes(x = pld, y = averageConnections)) +
    geom_line() +
    geom_ribbon(aes(ymin = averageConnections - sdConnections/2,
                    ymax = averageConnections + sdConnections/2), alpha = 0.2) + 
    xlab("Propagule duration (day)") + ylab("Degree centrality (average ± standard deviation)")
  
  dev.off()
  
  pdf(file=paste0("../Results/isolation degree.pdf"), width=12)
  ggplot(combResults, aes(x = pld, y = n.isolated.poly)) +
    geom_line() + 
    ylab("Isolation degree (number of reserves)") +
    xlab("Propagule duration (day)")
  dev.off()

}

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------