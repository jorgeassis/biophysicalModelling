## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Config.R")
source("Dependences/mainFunctions.R")

n.season <- "" # Spring; Summer; Autumn; Winter; "" for All

## ------------------------------------------------------------------------------------------------------------------

Connectivity.desc <- paste0(results.folder,"/particlePairedConnectivityAveraged",n.season,"Table.desc")

load(file=paste0(results.folder,"/particlePairedConnectivityAveraged",n.season,"TableNames.RData"))
Connectivity <- attach.big.matrix(Connectivity.desc)
Connectivity <- as.data.table(Connectivity[])
colnames(Connectivity) <- particles.connectivity.names
head(Connectivity)

load(paste0(results.folder,"/","sourceSinkSites.RData"))

source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

## ------------------------------------------------------------------------------------------------------------------------------
## Marine Distances

cost.surface <- raster("Data/Rasters/mask005.tif")
cost.surface.coastal <- raster("Data/Rasters/coastLineRes005.tif")

cost.surface <- aggregate(cost.surface,10)
cost.surface.coastal <- aggregate(cost.surface.coastal,10)

clipper <- as(extent(min(source.sink.xy[,2]) - 2,max(source.sink.xy[,2]) + 2,min(source.sink.xy[,3]) - 2,max(source.sink.xy[,3]) + 2), "SpatialPolygons")
plot(cost.surface, col="gray")
lines(clipper)

cost.surface <- crop(cost.surface,clipper)
cost.surface[is.na(cost.surface)] <- 0
plot(cost.surface,box=FALSE,legend=FALSE,col=c("gray","white"))
points( source.sink.xy[Source == 1 , .(Lon,Lat)], pch=19, col="Black")

# ----------------------------------

source.sink.xy.tr <- source.sink.xy
toRelocate <- which(extract(cost.surface,source.sink.xy.tr[,2:3]) == 0)

if( length(toRelocate) > 0) {
  
  xy.toRelocate <- xyFromCell(cost.surface.coastal,Which(cost.surface.coastal == 1, cell=TRUE))
  
  for( i in 1:length(toRelocate)) {
    idw.nearest.r <- get.knnx( xy.toRelocate , source.sink.xy.tr[toRelocate[i],2:3], k=1 , algorithm="kd_tree" )$nn.index
    source.sink.xy.tr[toRelocate[i],2] <- xy.toRelocate[idw.nearest.r,1]
    source.sink.xy.tr[toRelocate[i],3] <- xy.toRelocate[idw.nearest.r,2]
  }
  
  if(length(which(extract(cost.surface,source.sink.xy.tr[,2:3]) == 0)) > 0 ) { stop("Error :: 119")}
  
}

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(source.sink.xy.tr[Pair == source.sink.xy.tr$Pair[1],2:3]) , as.matrix(source.sink.xy.tr[Pair == source.sink.xy.tr$Pair[100],2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(source.sink.xy.tr[Pair == source.sink.xy.tr$Pair[1],2:3]) , as.matrix(source.sink.xy.tr[Pair == source.sink.xy.tr$Pair[100],2:3]) )

# ----------------------------------

n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores / 2)
registerDoParallel(cl.2)

marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  x.to <- x.to[x.to != 0]

  partial.distances <- sapply(x.to,function(x.i) { costDistance(raster_tr_corrected, as.matrix(source.sink.xy.tr[ source.sink.xy.tr$Pair == x , 2:3 ], ncol=2) , as.matrix( source.sink.xy.tr[ source.sink.xy.tr$Pair %in% x.i , 2:3 ], ncol=2)) } )
  partial.distances <- data.frame(Pair.from=rep(x,length(x.to)),Pair.to=x.to,Distance=c(partial.distances)/1000)

  return( partial.distances )
  
}

stopCluster(cl.2) ; rm(cl.2)
closeAllConnections(); gc(reset=TRUE)

head(marine.distances)

# ----------------------------------

distance.probability <- merge(Connectivity, marine.distances, by=c("Pair.from","Pair.to"))

file.remove( list.files(results.folder, full.names = TRUE, pattern = paste0("particlePairedConnectivityAveragedDist",n.season,"Table.bin") ) )
file.remove( list.files(results.folder, full.names = TRUE, pattern = paste0("particlePairedConnectivityAveragedDist",n.season,"Table.desc") ) )
averagedConnectivity <- as.big.matrix(as.matrix(distance.probability) , backingpath=paste0(results.folder,"/") , backingfile = paste0("particlePairedConnectivityAveragedDist",n.season,"Table.bin") , descriptorfile = paste0("particlePairedConnectivityAveragedDist",n.season,"Table.desc") )

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

distance.probability <- data.table(distance.probability)

# ----------------------------------

n.repetitions <- 1:max(distance.probability$Max.Time)

for( extract.simulation.days in n.repetitions) {
    
  n.days <- extract.simulation.days

  connectivityExportDir <- file.path(paste0(results.folder,"/connectivityExport/","Sim",n.season,str_pad(n.days, 3, pad = "0"),"Days/"))
  if( ! dir.exists(connectivityExportDir) ) { dir.create(file.path(connectivityExportDir), showWarnings = FALSE, recursive= TRUE) } 
  
  distance.probability.t <- distance.probability[Max.Time <= extract.simulation.days,]
  
  x <- distance.probability.t$Distance
  y <- distance.probability.t$Mean.Probability
  
  y <- y[x!=0 & x != Inf]
  x <- x[x!=0 & x != Inf]
  
  plotData <- data.frame(x=x,y=y)
  
  mainTheme <- theme(panel.grid.major = element_blank() ,
                     text = element_text(size=12) ,
                     axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                     axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) )
  
  p3 <- ggplot() +
    geom_point(data = plotData, aes(x=x, y=y), shape = 21,colour = "black", fill = "black", size = 2, stroke = 0.75, alpha = 0.2) +
    theme_minimal() + mainTheme + xlab("Distance (km)") + ylab("Probability of connectivity") + 
    geom_vline(xintercept = quantile(x,probs=0.95), linetype="dashed",  color = "gray", size=0.5) +
    annotate(geom="text", x=quantile(x,probs=0.95)+10, y=max(plotData$y), label=paste0(round(quantile(x,probs=0.95)),"km [95% of particles]"),size=4,family="Helvetica", color = "#7F7F7F",hjust = 0)
  
  pdf( file=paste0(connectivityExportDir,"/Probability vs Distance ",extract.simulation.days," days.pdf"), width = 10, height = 8 )
  print(p3)
  dev.off()
  
}

# ----------------------------------

tMax <- max(distance.probability$Max.Time)
resultsTime <- data.frame()

for( t in 1:tMax){
  
  distance.probability.t <- distance.probability[Max.Time <= t,]
  
  x <- distance.probability.t$Distance
  y <- distance.probability.t$Mean.Probability
  
  y <- y[x!=0 & x != Inf]
  x <- x[x!=0 & x != Inf]

    resultsTime <- rbind(resultsTime,data.frame(time=t,
                                              meandistance=mean(x,na.rm=T),
                                              maxdistance=max(x,na.rm=T),
                                              meanprobability=mean(y,na.rm=T),
                                              maxprobability=max(y,na.rm=T) ))
}
  
plotData <- data.frame(x=resultsTime$time,y=resultsTime$maxdistance)
p3 <- ggplot() +
  geom_point(data = plotData, aes(x=x, y=y), shape = 21,colour = "black", fill = "black", size = 2, stroke = 0.75, alpha = 0.5) +
  theme_minimal() + mainTheme + xlab("Dispersal period (day)") + ylab("Maximum travelled distance (km)")
p3

pdf( file=paste0(results.folder,"/Time vs Max Distance",n.season,".pdf"), width = 10, height = 8 )
print(p3)
dev.off()

plotData <- data.frame(x=resultsTime$time,y=resultsTime$meanprobability)
plotData[1,2] <- plotData[2,2] + (plotData[2,2] - plotData[3,2])
p3 <- ggplot() +
  geom_point(data = plotData, aes(x=x, y=y), shape = 21,colour = "black", fill = "black", size = 2, stroke = 0.75, alpha = 0.5) +
  theme_minimal() + mainTheme + xlab("Dispersal period (day)") + ylab("Mean probability of connectivity")
p3

pdf( file=paste0(results.folder,"/Time vs Mean Probability",n.season,".pdf"), width = 10, height = 8 )
print(p3)
dev.off()

# ----------------------------------

for( extract.simulation.days in n.repetitions) {
  
    n.days <- extract.simulation.days
  
    connectivityExportDir <- file.path(paste0(results.folder,"/connectivityExport/","Sim",n.season,str_pad(n.days, 3, pad = "0"),"Days/"))
    
    distance.probability.t <- distance.probability[Max.Time <= extract.simulation.days,]
    distance.probability.t[distance.probability.t >= 9e99999] <- NA
    
    summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance,na.rm=T),3) , round(max(distance.probability.t$Mean.Probability,na.rm=T),3) , round(max(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                                   Mean    = c( round(mean(distance.probability.t$Distance,na.rm=T),3) , round(mean(distance.probability.t$Mean.Probability,na.rm=T),3) , round(mean(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                                   SD      = c( round(sd(distance.probability.t$Distance,na.rm=T),3) , round(sd(distance.probability.t$Mean.Probability,na.rm=T),3) , round(sd(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                                   Median  = c( round(median(distance.probability.t$Distance,na.rm=T),3) , round(median(distance.probability.t$Mean.Probability,na.rm=T),3) , round(median(distance.probability.t$Mean.Time,na.rm=T),3) ) )
    row.names(summary.results) <- c("Distance","Probability","Time")
    
    write.csv(summary.results,file=paste0(connectivityExportDir,"summaryAverage.csv"))
    
    qt  <- quantile(distance.probability.t$Mean.Probability, probs = 0.95)
    distance.probability.t <- distance.probability.t[ Mean.Probability >= qt , ]
    
    summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance,na.rm=T),3) , round(max(distance.probability.t$Mean.Probability),3) , round(max(distance.probability.t$Mean.Time),3) ) ,
                                   Mean    = c( round(mean(distance.probability.t$Distance,na.rm=T),3) , round(mean(distance.probability.t$Mean.Probability),3) , round(mean(distance.probability.t$Mean.Time),3) ) ,
                                   SD      = c( round(sd(distance.probability.t$Distance,na.rm=T),3) , round(sd(distance.probability.t$Mean.Probability),3) , round(sd(distance.probability.t$Mean.Time),3) ) ,
                                   Median  = c( round(median(distance.probability.t$Distance,na.rm=T),3) , round(median(distance.probability.t$Mean.Probability),3) , round(median(distance.probability.t$Mean.Time),3) ) )
    row.names(summary.results) <- c("Distance","Probability","Time")
    write.csv(summary.results,file=paste0(connectivityExportDir,"summary95Percentile.csv"))
  
}

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------