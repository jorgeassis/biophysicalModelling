## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

source("0. Project Config.R")

## ------------------------------------------------------------------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

## -------------------

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## ------------------------------------------------------------------------------------------------------------------------------
## Prob. vs Distance Plot
##

clipper <- as(extent(min(source.sink.xy[,2] - 2),max(source.sink.xy[,2] + 2),min(source.sink.xy[,3] - 2),max(source.sink.xy[,3] + 2)), "SpatialPolygons")
crs(clipper) <- dt.projection


study.region <- crop(study.region, extent(new.extent.min.lon, new.extent.max.lon, new.extent.min.lat, new.extent.max.lat) + c(-resolution,+resolution,-resolution,+resolution) )
ocean.region <- crop(ocean.region, extent(new.extent.min.lon, new.extent.max.lon, new.extent.min.lat, new.extent.max.lat) + c(-resolution,+resolution,-resolution,+resolution) )

ocean.region[is.na(ocean.region)] <- 0
cost.surface <- ocean.region
cost.surface[cost.surface > 0] <- 1
plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(ocean.region,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(cells[cells[,1] == 1167,2:3]) , as.matrix(cells[cells[,1] == 8178,2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(cells[cells[,1] == 1167,2:3]), as.matrix(cells[cells[,1] == 8178,2:3]) )

# ----------------------------------

number.cores <- 2
n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores) ; registerDoParallel(cl.2)
marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  partial.distances <- costDistance(raster_tr_corrected, as.matrix(cells[ x , 2:3 ]) , as.matrix(cells[ x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=x,Pair.to=x.to,distance=c(partial.distances)/1000)
  return( partial.distances )
  
}
stopCluster(cl.2) ; rm(cl.2)
head(marine.distances)

# Save object

save(marine.distances,file=paste0(results.directory,"/marine.distances.RData"))
save(raster_tr_corrected,file=paste0(results.directory,"/cost.distance.raster.RData"))
