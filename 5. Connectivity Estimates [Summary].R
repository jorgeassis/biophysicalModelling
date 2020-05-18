## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")
source("Dependences.R")

number.cores <- 4

## ------------------------------------------------------------------------------------------------------------------

resultsFolder <- "Results2017" # Results

Connectivity <- read.big.matrix(paste0(project.folder,"/",resultsFolder,"/Connectivity.bm"))
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Time.max" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity

source.sink.xy <- read.big.matrix(paste0(project.folder,"/",resultsFolder,"/source.sink.bm"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

# source.sink.xy <- source.sink.xy[Source == 1 & Pair %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## ------------------------------------------------------------------------------------------------------------------------------
## Marine Distances

cost.surface <- raster("Data/Rasters/Mask.tif")
# cost.surface <- aggregate(cost.surface, fact=2)

clipper <- as(extent(min(source.sink.xy[,2]) - 2,max(source.sink.xy[,2]) + 2,min(source.sink.xy[,3]) - 2,max(source.sink.xy[,3]) + 2), "SpatialPolygons")
plot(cost.surface, col="gray")
lines(clipper)

cost.surface <- crop(cost.surface,clipper)
cost.surface[is.na(cost.surface)] <- 0
plot(cost.surface,box=FALSE,legend=FALSE,col=c("gray","white"))
points( source.sink.xy[Source == 1 , .(Lon,Lat)], pch=19, col="Black")

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(source.sink.xy[Pair == 4,2:3]) , as.matrix(source.sink.xy[Pair == 147,2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == 4,2:3]) , as.matrix(source.sink.xy[Pair == 1296,2:3]) )

# ----------------------------------

n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores) ; registerDoParallel(cl.2)

marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  x.to <- x.to[x.to != 0]

  partial.distances <- costDistance(raster_tr_corrected, as.matrix(source.sink.xy[ source.sink.xy$Pair == x , 2:3 ]) , as.matrix(source.sink.xy[ source.sink.xy$Pair %in% x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=rep(x,length(partial.distances)),Pair.to=source.sink.xy$Pair[source.sink.xy$Pair %in% x.to],Distance=c(partial.distances)/1000)
  
  zeros <- which(partial.distances$Distance == 0 & partial.distances$Pair.from != partial.distances$Pair.to)
  
  if( length(zeros) > 0 ) {
    
    for(z in 1:length(zeros)){
      
            partial.distances[zeros[z],3] <- spDistsN1( as.matrix(source.sink.xy[ Pair == partial.distances[zeros[z],1] , 2:3 ]), as.matrix(source.sink.xy[ Pair == partial.distances[zeros[z],2] , 2:3 ]), longlat=TRUE)
      
    }
    
  }
  
  
  return( partial.distances )
  
}

stopCluster(cl.2) ; rm(cl.2)

head(marine.distances)

# ----------------------------------

distance.probability <- merge(Connectivity, marine.distances, by=c("Pair.from","Pair.to"))
distance.probability <- as.big.matrix(as.matrix(distance.probability))
write.big.matrix(distance.probability, paste0(project.folder,"/",resultsFolder,"/Connectivity.Distance.bm"))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

distance.probability <- read.big.matrix(paste0(project.folder,"/",resultsFolder,"/Connectivity.Distance.bm"))
distance.probability <- data.table(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

# ----------------------------------

extract.simulation.days <- 30
distance.probability.t <- distance.probability[Time.max <= extract.simulation.days,]

x <- distance.probability.t$Distance
y <- distance.probability.t$Probability

pdf( file=paste0(project.folder,"Results/Probability vs Distance ",extract.simulation.days," days.pdf") , width = 10, height = 8 )

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(x,y,ylim=c(min(y),max(y)),col="#A6A6A6", pch=16,ylab="Mean probability of connectivity",xlab="Distance (km)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()

dev.off()

# ----------------------------------

resultsTime <- data.frame()

for( t in 1:30){
  
  distance.probability.t <- distance.probability[Time.max <= t,]
  
  x <- distance.probability.t$Distance
  y <- distance.probability.t$Probability
  
  x[x == Inf] <- NA
  y[y == Inf] <- NA

    resultsTime <- rbind(resultsTime,data.frame(time=t,
                                              meandistance=mean(x,na.rm=T),
                                              maxdistance=max(x,na.rm=T),
                                              meanprobability=mean(y,na.rm=T),
                                              maxprobability=max(y,na.rm=T) ))
}
  
pdf( file=paste0(project.folder,"Results/Time vs Mean Distance.pdf") , width = 10, height = 8 )
par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(resultsTime$time,resultsTime$meandistance,ylim=c(min(resultsTime$meandistance),max(resultsTime$meandistance)),col="#A6A6A6", pch=16,ylab="Mean travelled distance (km)",xlab="Dispersal period (day)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
dev.off()

pdf( file=paste0(project.folder,"Results/Time vs Max Distance.pdf") , width = 10, height = 8 )
par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(resultsTime$time,resultsTime$maxdistance,ylim=c(min(resultsTime$maxdistance),max(resultsTime$maxdistance)),col="#A6A6A6", pch=16,ylab="Maximum travelled distance (km)",xlab="Dispersal period (day)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
dev.off()

pdf( file=paste0(project.folder,"Results/Time vs Mean Probability.pdf") , width = 10, height = 8 )
par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(resultsTime$time,resultsTime$meanprobability,ylim=c(min(resultsTime$meanprobability),max(resultsTime$meanprobability)),col="#A6A6A6", pch=16,ylab="Mean probability of connectivity",xlab="Dispersal period (day)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
dev.off()

pdf( file=paste0(project.folder,"Results/Time vs Max Probability.pdf") , width = 10, height = 8 )
par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(resultsTime$time,resultsTime$maxprobability,ylim=c(min(resultsTime$maxprobability),max(resultsTime$maxprobability)),col="#A6A6A6", pch=16,ylab="Maximum probability of connectivity",xlab="Dispersal period (day)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
dev.off()

# ----------------------------------

# Summary 1

extract.simulation.days <- 30
distance.probability.t <- distance.probability[Time.max <= extract.simulation.days,]
distance.probability.t[distance.probability.t >= 9e99999] <- NA

summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance,na.rm=T),3) , round(max(distance.probability.t$Probability,na.rm=T),3) , round(max(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                               Mean    = c( round(mean(distance.probability.t$Distance,na.rm=T),3) , round(mean(distance.probability.t$Probability,na.rm=T),3) , round(mean(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                               SD      = c( round(sd(distance.probability.t$Distance,na.rm=T),3) , round(sd(distance.probability.t$Probability,na.rm=T),3) , round(sd(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                               Median  = c( round(median(distance.probability.t$Distance,na.rm=T),3) , round(median(distance.probability.t$Probability,na.rm=T),3) , round(median(distance.probability.t$Mean.Time,na.rm=T),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

qt  <- quantile(distance.probability.t$Probability, probs = 0.95)
distance.probability.t <- distance.probability.t[ Probability >= qt , ]

summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance,na.rm=T),3) , round(max(distance.probability.t$Probability),3) , round(max(distance.probability.t$Mean.Time),3) ) ,
                               Mean    = c( round(mean(distance.probability.t$Distance,na.rm=T),3) , round(mean(distance.probability.t$Probability),3) , round(mean(distance.probability.t$Mean.Time),3) ) ,
                               SD      = c( round(sd(distance.probability.t$Distance,na.rm=T),3) , round(sd(distance.probability.t$Probability),3) , round(sd(distance.probability.t$Mean.Time),3) ) ,
                               Median  = c( round(median(distance.probability.t$Distance,na.rm=T),3) , round(median(distance.probability.t$Probability),3) , round(median(distance.probability.t$Mean.Time),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

# Identify which have a high threshold

land.surface <- raster("Data/Rasters/Mask.tif")

cells.i <- distance.probability[ Time.max <= 5 & Distance >= 200 & Distance < 1000000000, Pair.from   ]
cells.j <- distance.probability[ Time.max <= 5 & Distance >= 200 & Distance < 1000000000, Pair.to  ]

clipper <- as(extent(min(source.sink.xy[unique(c(cells.i,cells.j)),2]) - 2 , max(source.sink.xy[unique(c(cells.i,cells.j)),2]) + 2,min(source.sink.xy[unique(c(cells.i,cells.j)),3]) - 2,max(source.sink.xy[unique(c(cells.i,cells.j)),3]) + 2), "SpatialPolygons")
land.surface <- crop(land.surface,clipper)

plot(land.surface,box=FALSE,legend=FALSE,col=c("black","Gray"))
points(source.sink.xy[cells.i,2:3],col="red")
points(source.sink.xy[cells.j,2:3],col="green")

distance.probability[ Pair.from %in% cells.i , ]
reference.table[ cell %in% cells.i & cell.rafted %in% cells.j , ]

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
