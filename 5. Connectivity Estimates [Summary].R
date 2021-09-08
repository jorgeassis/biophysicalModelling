## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Config.R")
source("Dependences.R")
 
n.season <- "" # Spring; Summer; Autumn; Winter; "" for All

## ------------------------------------------------------------------------------------------------------------------

Connectivity <- read.big.matrix(paste0(results.folder,"/InternalProc/","connectivityEstimatesAveraged",n.season,".bm"))
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Time.max" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity

source.sink.xy <- read.big.matrix(paste0(results.folder,"/InternalProc/","source.sink.bm"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

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
lines( shortestPath(raster_tr_corrected, as.matrix(source.sink.xy[Pair == source.sink.xy$Pair[1],2:3]) , as.matrix(source.sink.xy[Pair == source.sink.xy$Pair[100],2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == source.sink.xy$Pair[1],2:3]) , as.matrix(source.sink.xy[Pair == source.sink.xy$Pair[100],2:3]) )

# ----------------------------------

n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores) ; registerDoParallel(cl.2)

marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  x.to <- x.to[x.to != 0]

  partial.distances <- sapply(x.to,function(x.i) { costDistance(raster_tr_corrected, as.matrix(source.sink.xy[ source.sink.xy$Pair == x , 2:3 ][1,]) , as.matrix( source.sink.xy[ source.sink.xy$Pair %in% x.i , 2:3 ][1,])) } )
  partial.distances <- data.frame(Pair.from=rep(x,length(x.to)),Pair.to=x.to,Distance=c(partial.distances)/1000)

  return( partial.distances )
  
}

stopCluster(cl.2) ; rm(cl.2)

head(marine.distances)

# ----------------------------------

distance.probability <- merge(Connectivity, marine.distances, by=c("Pair.from","Pair.to"))
distance.probability <- as.big.matrix(as.matrix(distance.probability))
write.big.matrix(distance.probability, paste0(results.folder,"/InternalProc/","connectivityEstimatesAveragedDistances",n.season,".bm"))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

distance.probability <- read.big.matrix(paste0(results.folder,"/InternalProc/","connectivityEstimatesAveragedDistances",n.season,".bm"))
distance.probability <- data.table(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")
head(distance.probability)

# ----------------------------------

n.repetitions <- 1:max(distance.probability$Time.max)

for( extract.simulation.days in n.repetitions) {
    
  n.days <- extract.simulation.days

  if( ! paste0("Sim",n.season,str_pad(n.days, 3, pad = "0"),"Days") %in% list.files(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/"))) ) { dir.create(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/","Sim",n.season,str_pad(n.days, 3, pad = "0"),"Days"))) }
  
  connectivityExportDir <- file.path(paste0(results.folder,"/connectivityExport/","Sim",n.season,str_pad(n.days, 3, pad = "0"),"Days/"))

  if( ! dir.exists(connectivityExportDir) ) { dir.create(file.path(connectivityExportDir), showWarnings = FALSE, recursive= TRUE) } 
  
  distance.probability.t <- distance.probability[Time.max <= extract.simulation.days,]
  
  x <- distance.probability.t$Distance
  y <- distance.probability.t$Probability
  
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

tMax <- max(distance.probability$Time.max)
resultsTime <- data.frame()

for( t in 1:tMax){
  
  distance.probability.t <- distance.probability[Time.max <= t,]
  
  x <- distance.probability.t$Distance
  y <- distance.probability.t$Probability
  
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
    
    distance.probability.t <- distance.probability[Time.max <= extract.simulation.days,]
    distance.probability.t[distance.probability.t >= 9e99999] <- NA
    
    summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance,na.rm=T),3) , round(max(distance.probability.t$Probability,na.rm=T),3) , round(max(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                                   Mean    = c( round(mean(distance.probability.t$Distance,na.rm=T),3) , round(mean(distance.probability.t$Probability,na.rm=T),3) , round(mean(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                                   SD      = c( round(sd(distance.probability.t$Distance,na.rm=T),3) , round(sd(distance.probability.t$Probability,na.rm=T),3) , round(sd(distance.probability.t$Mean.Time,na.rm=T),3) ) ,
                                   Median  = c( round(median(distance.probability.t$Distance,na.rm=T),3) , round(median(distance.probability.t$Probability,na.rm=T),3) , round(median(distance.probability.t$Mean.Time,na.rm=T),3) ) )
    row.names(summary.results) <- c("Distance","Probability","Time")
    
    write.csv(summary.results,file=paste0(connectivityExportDir,"summaryAverage.csv"))
    
    qt  <- quantile(distance.probability.t$Probability, probs = 0.95)
    distance.probability.t <- distance.probability.t[ Probability >= qt , ]
    
    summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance,na.rm=T),3) , round(max(distance.probability.t$Probability),3) , round(max(distance.probability.t$Mean.Time),3) ) ,
                                   Mean    = c( round(mean(distance.probability.t$Distance,na.rm=T),3) , round(mean(distance.probability.t$Probability),3) , round(mean(distance.probability.t$Mean.Time),3) ) ,
                                   SD      = c( round(sd(distance.probability.t$Distance,na.rm=T),3) , round(sd(distance.probability.t$Probability),3) , round(sd(distance.probability.t$Mean.Time),3) ) ,
                                   Median  = c( round(median(distance.probability.t$Distance,na.rm=T),3) , round(median(distance.probability.t$Probability),3) , round(median(distance.probability.t$Mean.Time),3) ) )
    row.names(summary.results) <- c("Distance","Probability","Time")
    write.csv(summary.results,file=paste0(connectivityExportDir,"summary95Percentile.csv"))
  
}

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
