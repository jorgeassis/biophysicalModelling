## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

number.cores <- 40

## ------------------------------------------------------------------------------------------------------------------

Connectivity <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.bm"))
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Time.max" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity

source.sink.xy <- read.big.matrix(paste0(project.folder,"/Results/source.sink.bm"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

# source.sink.xy <- source.sink.xy[Source == 1 & Pair %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## ------------------------------------------------------------------------------------------------------------------------------
## Marine Distances

cost.surface <- raster("Data/Rasters/Mask.tif")
# cost.surface <- aggregate(cost.surface, fact=2)

clipper <- as(extent(min(source.sink.xy[,2]) - 2,max(source.sink.xy[,2]) + 2,min(source.sink.xy[,3]) - 2,max(source.sink.xy[,3]) + 2), "SpatialPolygons")
plot(cost.surface)
lines(clipper)

cost.surface <- crop(cost.surface,clipper)
cost.surface[is.na(cost.surface)] <- 0
plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(source.sink.xy[Pair == 1167,2:3]) , as.matrix(source.sink.xy[Pair == 1,2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == 1167,2:3]) , as.matrix(source.sink.xy[Pair == 1,2:3]) )

# ----------------------------------

n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores) ; registerDoParallel(cl.2)

marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  x.to <- x.to[x.to != 0]
  partial.distances <- costDistance(raster_tr_corrected, as.matrix(source.sink.xy[ x , 2:3 ]) , as.matrix(source.sink.xy[ x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=rep(x,length(partial.distances)),Pair.to=x.to,Distance=c(partial.distances)/1000)
  
  zeros <- which(partial.distances$Distance == 0 & partial.distances$Pair.from != partial.distances$Pair.to)
  
  if( length(zeros) > 0 ) {
    
    for(z in 1:length(zeros)){
      
      partial.distances[zeros[z],3] <- spDistsN1( as.matrix(source.sink.xy[ partial.distances[zeros[z],1] , 2:3 ]), as.matrix(source.sink.xy[ partial.distances[zeros[z],2] , 2:3 ]), longlat=TRUE)
      
    }
    
  }
  
  
  return( partial.distances )
  
}

stopCluster(cl.2) ; rm(cl.2)

head(marine.distances)

# ----------------------------------

distance.probability <- merge(Connectivity, marine.distances, by=c("Pair.from","Pair.to"))
distance.probability <- as.big.matrix(as.matrix(distance.probability))
write.big.matrix(distance.probability, paste0(project.folder,"/Results/Connectivity.Distance.bm"))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

distance.probability <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.Distance.bm"))
distance.probability <- data.table(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

# ----------------------------------

extract.simulation.days <- 30

distance.probability.t <- distance.probability[Time.max <= extract.simulation.days,]
max(distance.probability.t$Time.max)

# Summary 10:8

ggplot(distance.probability.t , aes(x=Distance,y=Probability)) + 
  geom_point(alpha = 0.3) + 
  theme_bw(base_size = 14) + 
  labs(x = "Distance (km)" , y = "Mean probability of connectivity") +
  theme(panel.background = element_rect(colour = "black") )

# ----------------------------------

# Summary 1

summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance),3) , round(max(distance.probability.t$Probability),3) , round(max(distance.probability.t$Mean.Time),3) ) ,
                               Mean    = c( round(mean(distance.probability.t$Distance),3) , round(mean(distance.probability.t$Probability),3) , round(mean(distance.probability.t$Mean.Time),3) ) ,
                               SD      = c( round(sd(distance.probability.t$Distance),3) , round(sd(distance.probability.t$Probability),3) , round(sd(distance.probability.t$Mean.Time),3) ) ,
                               Median  = c( round(median(distance.probability.t$Distance),3) , round(median(distance.probability.t$Probability),3) , round(median(distance.probability.t$Mean.Time),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

qt  <- quantile(distance.probability.t$Probability, probs = 0.95)
distance.probability.t <- distance.probability.t[ Probability >= qt , ]

summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance),3) , round(max(distance.probability.t$Probability),3) , round(max(distance.probability.t$Mean.Time),3) ) ,
                               Mean    = c( round(mean(distance.probability.t$Distance),3) , round(mean(distance.probability.t$Probability),3) , round(mean(distance.probability.t$Mean.Time),3) ) ,
                               SD      = c( round(sd(distance.probability.t$Distance),3) , round(sd(distance.probability.t$Probability),3) , round(sd(distance.probability.t$Mean.Time),3) ) ,
                               Median  = c( round(median(distance.probability.t$Distance),3) , round(median(distance.probability.t$Probability),3) , round(median(distance.probability.t$Mean.Time),3) ) )
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
