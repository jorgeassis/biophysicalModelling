## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=13) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)) ,
                   legend.title = element_blank() ,
                   legend.margin=margin(c(0.3,1,0.3,1), unit='lines') ,
                   legend.background = element_rect(fill="white", size=0.2, linetype="solid",  colour ="#979797"))

theme_map <- theme_minimal() + theme( text = element_text(family = "Helvetica", color = "#22211d"),
                                      axis.line = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.ticks = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      panel.grid.major = element_blank(), # element_line(color = "#979797", size = 0.05),
                                      panel.grid.minor = element_blank(), # element_line(color = "#ebebe5", size = 0.2),
                                      plot.background = element_rect(fill = "#FFFFFF", color = NA), 
                                      panel.background = element_rect(fill = "#FFFFFF", color = NA), 
                                      legend.background = element_rect(fill = "#FFFFFF", color = NA),
                                      panel.border = element_blank() )

## ----------------------------------------------------
## ----------------------------------------------------

project.name.c <- paste0(results.folder,"/Genetics/",mainfileName,"/",speciesName,"/",markerName,"/",toupper(index),"/Simulation",season,str_pad(pld.period, 3, pad = "0"),"Days/")
if( ! dir.exists(project.name.c) ) { dir.create(file.path(project.name.c),recursive = TRUE, showWarnings = FALSE) } 

## ----------------------------------------------------

# Open connectivity

source.sink.xy <- loadRData(paste0(results.folder,"/","sourceSinkSites.RData"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

# Open PD values to create function

PDObserved <- sort(read.csv("../Data/PDFunction.csv")[,1])

# Max observed PD
max.pd <- 110

PDObserved.Func <- density(PDObserved, kernel = c("gaussian"), adjust=15, from=0, to=max.pd)#max(PDObserved))
PDObserved.Func <- data.frame(x=PDObserved.Func$x,y=PDObserved.Func$y)
PDObserved.Func <- PDObserved.Func[PDObserved.Func$x > 0 & PDObserved.Func$x <= max.pd,]
plot(PDObserved.Func)

# Open connectivity

Connectivity.desc <- paste0(results.folder,"/","particlePairedConnectivityAveragedTable.desc")
load(file=paste0(results.folder,"/particlePairedConnectivityAveragedTableNames.RData"))
Connectivity <- attach.big.matrix(Connectivity.desc)
Connectivity <- as.data.table(Connectivity[])
colnames(Connectivity) <- particles.connectivity.names
Connectivity <- Connectivity[Connectivity$Max.Time <= max.pd,]

# Multiple by the function

for( i in 1:(max.pd+1)) {
  
  Connectivity[Connectivity$Max.Time >= i - 1 & Connectivity$Max.Time < i,"Mean.Probability"] <- Connectivity[Connectivity$Max.Time >= i - 1 & Connectivity$Max.Time < i,"Mean.Probability"] *  mean(PDObserved.Func[PDObserved.Func$x >= i - 1 & PDObserved.Func$x < i,"y"])
  
}

## ------------------------------------------------------------------------------------------------------------
## Read main sources

source.sink.xy.sp <- source.sink.xy
coordinates(source.sink.xy.sp) <- ~Lon+Lat
crs(source.sink.xy.sp) <- dt.projection

# Add missing connections

missingFrom <- source.sink.xy$Pair[which(! source.sink.xy$Pair %in% unique(Connectivity$Pair.from))]

if(length(missingFrom) > 0) { 
  
  Connectivity <- rbind(Connectivity,data.frame(Pair.from=missingFrom,Pair.to=missingFrom),fill=TRUE)
  Connectivity[is.na(Connectivity)] <- 0
  
}

missingTo <- source.sink.xy$Pair[which(! source.sink.xy$Pair %in% unique(Connectivity$Pair.to))]

if(length(missingTo) > 0) { 
  
  Connectivity <- rbind(Connectivity,data.frame(Pair.from=missingTo,Pair.to=missingTo),fill=TRUE)
  Connectivity[is.na(Connectivity)] <- 0
  
}

save(Connectivity,file=paste0("../Data/ConnectivityIntegrativePDFunction.RData"))

## --------------------------------------------------------------------
## --------------------------------------------------------------------
# Temporary Subset

subseter <- as.vector(extent(source.sink.xy.sp) + c(-0.5,0.5,-0.5,0.5))

## --------------------------------------------------------------------
## --------------------------------------------------------------------

robinson <- CRS("+proj=robin +over")

source.sink.xy <- source.sink.xy[source.sink.xy$Lon >= subseter[1] & source.sink.xy$Lon <= subseter[2] & source.sink.xy$Lat >= subseter[3] & source.sink.xy$Lat <= subseter[4], ]
source.sink.xy.sp <- crop(source.sink.xy.sp,extent(subseter))

if( c == 1 ) {
  
  hexagons.sourcesink.shp <- shapefile(paste0(results.folder,"/sourceSinkSites.shp"))
  
  worldMap <- ne_countries(scale = 10, returnclass = "sp")
  worldMap <- crop(worldMap,extent(subseter + c(-10,10,-10,10))  ) 
  
  bb <- sf::st_union(sf::st_make_grid(
    st_bbox(c(xmin = extent(worldMap)[1],
              xmax = extent(worldMap)[2],
              ymax = extent(worldMap)[3],
              ymin = extent(worldMap)[4]), crs = st_crs(4326)), n = 100))
  
  # bb_robinson <- st_transform(bb, as.character(robinson))
  # worldMap <- spTransform(worldMap, robinson)
  # hexagons.sourcesink.shp <- spTransform(hexagons.sourcesink.shp, robinson)
  
  mapRegion <- ggplot() +
    geom_polygon(data = worldMap , fill = "#CDCDCD", colour = "#CDCDCD" , size=0.25 ,  aes(long, lat, group = group))  + coord_map() +
    geom_polygon(data = hexagons.sourcesink.shp , fill = "gray", colour = "Black" , size=0.1 ,  aes(long, lat, group = group))  + coord_map() +
    coord_sf(xlim = c( extent(worldMap)[1], extent(worldMap)[2]), ylim = c(extent(worldMap)[3], extent(worldMap)[4]), expand = FALSE) +
    geom_sf(data=bb, colour='#686868', linetype='solid', fill = NA, size=0.1) +
    theme_map
  
}

# plot(source.sink.xy[,2:3])
# mapRegion

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Pairwise Connectivity estimates

popCoordinatesSS <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(popCoordinates),longlat = TRUE)
popCoordinatesSS <- apply(popCoordinatesSS,2,which.min)
popCoordinatesSS <- source.sink.xy[popCoordinatesSS,1]
popCoordinatesSS <- unlist(popCoordinatesSS)

# plot(popCoordinates)
# points(source.sink.xy[Pair %in% popCoordinatesSS,2:3],col="red")

unique.sites <- unique(popCoordinatesSS)
popDifferentiation.i <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))

for( u in 1:length(unique.sites)) {
  for( v in 1:length(unique.sites)) {
    
    pairs.u <- which( unique.sites[u] == popCoordinatesSS )
    pairs.v <- which( unique.sites[v] == popCoordinatesSS )
    popDifferentiation.i[u,v] <- mean(unlist(popDifferentiation[pairs.u,pairs.v]),na.rm=T)
    popDifferentiation.i[v,u] <- mean(unlist(popDifferentiation[pairs.u,pairs.v]),na.rm=T)
    
  } }

diag(popDifferentiation.i) <- NA
popDifferentiation <- popDifferentiation.i
popCoordinates <- data.frame(Lon=unlist(sapply(unique.sites,function(x) data.frame(source.sink.xy[Pair == x,2]) )),Lat=unlist(sapply(unique.sites,function(x) data.frame(source.sink.xy[Pair == x,3]) )))

if(TRUE %in% (duplicated(popCoordinatesSS))) { popCoordinates.n <- popCoordinates.n[-which(duplicated(popCoordinatesSS))] }
popCoordinatesSS <- unique.sites
colnames(popDifferentiation) <- popCoordinates.n

if(length(popCoordinatesSS) < 3) { next }

## ------------------------------------------
## ------------------------------------------

pdf(file=paste0(project.name.c,"RecAllSampling.pdf"), width=12)
print(
  mapRegion +
    geom_polygon(data = hexagons.sourcesink.shp[hexagons.sourcesink.shp$ID %in% popCoordinatesSS ,] , fill = "Red", colour = "Red" , size=1.5 ,  aes(long, lat, group = group))  + coord_map() +
    coord_sf(crs=robinson,xlim = c(min(popCoordinates$Lon) - 10,max(popCoordinates$Lon) + 10), ylim = c(min(popCoordinates$Lat) - 10,max(popCoordinates$Lat) + 10), expand = FALSE)
)
dev.off()

pdf(file=paste0(project.name.c,"RecSamplingLabels.pdf"), width=12)
print(
  mapRegion +
    geom_text(aes(x=popCoordinates$Lon,y=popCoordinates$Lat,label = popCoordinates.n), check_overlap = TRUE) +
    geom_label(aes(label=popCoordinates.n, x=popCoordinates$Lon,y=popCoordinates$Lat),check_overlap = T) +
    coord_sf(crs=robinson,xlim = c(min(popCoordinates$Lon) - 10,max(popCoordinates$Lon) + 10), ylim = c(min(popCoordinates$Lat) - 10,max(popCoordinates$Lat) + 10), expand = FALSE)
)
dev.off()

## ---------------------------------------------------
## ---------------------------------------------------
## All Pairs

comb <- Connectivity
comb <- as.data.frame( comb[ sort(comb$Mean.Probability , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Mean.Probability")] )
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
# E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
# E(graph.obj)$weight = comb[,3] 
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
if( forceUnidirectional ) { graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") } # min / mean / max
graph.obj <- simplify(graph.obj)

## ---------------------------------------------------

findNeightboursDist <- 16

cl.3 <- makeCluster(number.cores)
registerDoParallel(cl.3)

potential.connectivity <- foreach(from=popCoordinatesSS, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  res.connectivity.to <- numeric(0)
  
  for( to in popCoordinatesSS ) {
    
    from.i <- sort( spDists( as.matrix(source.sink.xy[,.(Lon,Lat)]) , as.matrix(source.sink.xy[Pair == from,.(Lon,Lat)]) , longlat = TRUE ) , index.return=TRUE)
    from.i <- source.sink.xy[from.i$ix[which(from.i$x <= findNeightboursDist )],Pair]
    to.i <- sort( spDists( as.matrix(source.sink.xy[,.(Lon,Lat)]) , as.matrix(source.sink.xy[Pair == to,.(Lon,Lat)]) , longlat = TRUE ) , index.return=TRUE)
    to.i <- source.sink.xy[to.i$ix[which(to.i$x <= findNeightboursDist )],Pair]
    
    trials <- expand.grid(from.i,to.i)
    
    for( tr in 1:nrow(trials)) {
      
      possible.paths.y <- get.shortest.paths(graph.obj,as.character( trials[tr,1] ) , as.character( trials[tr,2] ),mode="out")$vpath
      stones.t <- as.numeric(names(possible.paths.y[[1]]))
      stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
      
      path.values <- apply( stones.t.interm , 1 , function(z) { comb[ comb[,1] == z[1] & comb[,2] == z[2] , 3 ][1] }   )
      
      if( forceUnidirectional ) {
        
        path.values.i <- apply( stones.t.interm , 1 , function(z) { comb[ comb[,1] == z[1] & comb[,2] == z[2] , 3 ][1] }   )
        path.values.j <- apply( stones.t.interm , 1 , function(z) { comb[ comb[,1] == z[2] & comb[,2] == z[1] , 3 ][1] }   )
        path.values <- apply( cbind(path.values.i,path.values.j) , 1 , mean , na.rm=T)
        
      }
      
      if( length(path.values) > 0 & sum(is.na(path.values)) == 0 ) { break }
      
    }
    
    if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
    if( length(path.values) == 0) { path.values <- 0 }
    if( from == to ) { path.values <- 1 }
    
    res.connectivity.to <- c(res.connectivity.to,path.values)
    
  }
  
  ## ---------------------
  
  popDifferentiation.to.i <- sapply(popCoordinatesSS, function(x) { popDifferentiation[ which( popCoordinatesSS == from), which( popCoordinatesSS == x)] } )
  popDifferentiation.to.j <- sapply(popCoordinatesSS, function(x) { popDifferentiation[ which( popCoordinatesSS == x), which( popCoordinatesSS == from)] } )
  popDifferentiation.to <- apply(cbind(popDifferentiation.to.i,popDifferentiation.to.j), 1, mean, na.rm=T)
  popDifferentiation.to[is.na(popDifferentiation.to)] <- 0
  
  temp.res <- data.frame( pair.from = from,
                          pair.to = popCoordinatesSS , 
                          popDifferentiation = popDifferentiation.to, 
                          probability.ss = res.connectivity.to )
  
  return(temp.res)
  
  ## ---------------------
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

## ---------------------------------------------------

connectivity <- do.call(rbind,potential.connectivity)
connectivity <- connectivity[connectivity[,1] != connectivity[,2] ,]

if( sum(connectivity$probability.ss != 0) == 0 ) { connectivity$probability.ss <- 0.0001 }

connectivity <- connectivity[connectivity$probability.ss != 0 & ! is.na(connectivity$probability.ss),]
connectivity[,"probability.ss"] <- log(connectivity[,"probability.ss"])

## ---------------

norm <- expand.grid(unique(connectivity[,1]),unique(connectivity[,2]))
norm <- norm[norm[,1] != norm[,2] ,]
connectivity.final <- data.frame()

for( i in 1:nrow(norm)) {
  
  t.1 <- which( connectivity[,1] == norm[i,1] & connectivity[,2] == norm[i,2] )
  t.2 <- which( connectivity[,1] == norm[i,2] & connectivity[,2] == norm[i,1] )
  
  if( length(t.1) == 0 & length(t.2) == 0) { next }
  
  Differentiation <- mean(c( as.numeric(as.character(connectivity$popDifferentiation[t.1])), as.numeric(as.character(connectivity$popDifferentiation[t.2]))))
  Connectivity.max <- max(c( as.numeric(as.character(connectivity$probability.ss[t.1])), as.numeric(as.character(connectivity$probability.ss[t.2]))))
  Connectivity.mean <- mean(c( as.numeric(as.character(connectivity$probability.ss[t.1])), as.numeric(as.character(connectivity$probability.ss[t.2]))))
  Connectivity.min <- min(c( as.numeric(as.character(connectivity$probability.ss[t.1])), as.numeric(as.character(connectivity$probability.ss[t.2]))))
  
  connectivity.final <- rbind(connectivity.final,
                              data.frame(from=norm[i,1],
                                         to=norm[i,2],
                                         Differentiation = Differentiation,
                                         Connectivity.max = Connectivity.max ,
                                         Connectivity.mean = Connectivity.mean ,
                                         Connectivity.min = Connectivity.min
                                         , stringsAsFactors = FALSE ) )
}

write.csv(connectivity.final,paste0(project.name.c,"resultMatrix.csv"), row.names = FALSE)

## ------------------------------------
## ------------------------------------