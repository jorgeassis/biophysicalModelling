## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)

setwd("~/Desktop/Transport simulations/Git")

source("0. Project Config.R")
library(rgeos)
library(reshape2)
library(sf)

sql.project.name <- "MPA"
number.cores <- 8

mpa.shp.filename <- "/home/corallium/Desktop/Transport simulations/Data/Spatial/Europe/MPAs_Final.shp"
mpa.shp.filename.original <- "/home/corallium/Desktop/Transport simulations/Data/Spatial/Europe _ Original/European MPAs.shp"

sql.file <- "../Results/SQL/MPASimulationResults.sql"
bigmatrix.file <- "../InternalProc/particles.reference.desc"
sorce.sink.cells.file <- "../Results/source.sink.bm"

## ------------------------------------------------------------------------------------------------------------
## Read main sources

sql <- dbConnect(RSQLite::SQLite(), sql.file)
cell.to.process <- 1:nrow(dbReadTable(sql, "SourceSinkSites"))
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell[1]
n.new.particles.per.day <- dbReadTable(sql, "Parameters")$n.new.particles.per.day[1]
n.steps.per.day <- dbReadTable(sql, "Parameters")$n.hours.per.day[1]
dbDisconnect(sql)

source.sink.xy <- read.big.matrix("../Results/source.sink.bm")
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for all combinations of spawning months

# load("Connectivity.PERIOD.Rdata") OR run foreach section bellow

spawn.p <- c(1,2,3,4,5,6)

particles.reference.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.reference.desc"))

cl.2 <- makeCluster(number.cores , type="FORK")
registerDoParallel(cl.2)

Connectivity <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  particles.reference.bm.i <- attach.big.matrix(particles.reference.bm.desc)
  
  connectivity.temp.m <- particles.reference.bm.i[ mwhich(particles.reference.bm.i,2,list(cell.id.ref.f), list('eq')) , ]
  connectivity.temp.m <- data.frame(connectivity.temp.m)
  colnames(connectivity.temp.m) <- c("id","start.cell","start.year","start.month","start.day","pos.lon","pos.lat","pos.alt","state","t.start","t.finish","cell.rafted")
  
  connectivity.temp.m <- connectivity.temp.m[connectivity.temp.m$cell.rafted != 0 & connectivity.temp.m$state == 2,]
  connectivity.temp.m <- data.frame(connectivity.temp.m,travel.time= (1 + connectivity.temp.m$t.finish - connectivity.temp.m$t.start) / n.steps.per.day)
  connectivity.temp.m <- as.data.table(connectivity.temp.m)
  
  connectivity.temp.m <- connectivity.temp.m[start.month %in% spawn.p , ]
  
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

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity[is.na(Connectivity)] <- 0
Connectivity <- data.table(Connectivity[,])
Connectivity ; gc()

save(Connectivity,file="Connectivity.PERIOD.Rdata")

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Read main shapefiles

mpa.shp <- shapefile(mpa.shp.filename)
mpa.shp <-gBuffer(mpa.shp, byid=TRUE, width=0)
mpa.shp <- disaggregate(mpa.shp)

mpa.shp.original <- shapefile(mpa.shp.filename.original)
mpa.shp.original <-gBuffer(mpa.shp.original, byid=TRUE, width=0)

mpa.shp.notake <- which(mpa.shp.original$no_take == "All" | mpa.shp.original$no_take == "Part" )
mpa.shp.notake <- mpa.shp.original[mpa.shp.notake,]

## --------------------------------------------------------------------
## --------------------------------------------------------------------
# Subset temporario - Delete section

source.sink.xy <- source.sink.xy[source.sink.xy$Lon >= -20 & source.sink.xy$Lon <= 3 & source.sink.xy$Lat >= 20 & source.sink.xy$Lat <= 44.5, ]
plot(source.sink.xy[,2:3])

mpa.shp.notake <- crop(mpa.shp.notake,extent(-20,3,20,44.5))
plot(mpa.shp.notake)

## --------------------------------------------------------------------
## --------------------------------------------------------------------
## Identify MPAs in source sink sites

source.sink.xy <- cbind( source.sink.xy , data.frame(id.mpa=0,stringsAsFactors = FALSE) )

source.sink.xy.sp <- source.sink.xy[,2:3]
coordinates(source.sink.xy.sp) <- ~Lon+Lat
crs(source.sink.xy.sp) <- crs(mpa.shp)

# load("source.sink.xy.Rdata") OR run foreach section bellow

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

source.sink.xy.id <- foreach( source.sink.xy.i = 1:nrow(source.sink.xy) , .verbose=FALSE, .combine = rbind ,  .packages=c("rgeos","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  distances <- gDistance(source.sink.xy.sp[source.sink.xy.i], mpa.shp,byid=TRUE)
  return(which.min(as.vector(distances)))
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

source.sink.xy$id.mpa <- as.vector(unlist(source.sink.xy.id))
save(source.sink.xy,file="source.sink.xy.Rdata")

## ------------------------------------------------------------------------------------------

mpas.id <- unique(source.sink.xy$id.mpa)
mpa.shp <- mpa.shp[mpas.id,]
comb.pairs.mpas <- expand.grid(mpas.id,mpas.id)
nrow(comb.pairs.mpas)

plot(mpa.shp)

## ------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------

pld.period <- 30

## ------------------

# load("connectivity.mpas.Rdata") OR run foreach section bellow

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

connectivity <- foreach(pair.i=1:nrow(comb.pairs.mpas), .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN","bigmemory")) %dopar% { # 
  
  mpa.id.1 <- comb.pairs.mpas[pair.i,1]
  mpa.id.2 <- comb.pairs.mpas[pair.i,2]
  
  mpa.cells.1 <- as.vector(unlist(source.sink.xy[ id.mpa == mpa.id.1 , 1 ]))
  mpa.cells.2 <- as.vector(unlist(source.sink.xy[ id.mpa == mpa.id.2 , 1 ]))
  
  temp.result <- Connectivity[ Pair.from %in% mpa.cells.1 & Pair.to %in% mpa.cells.2 & Time.max <= pld.period , ]
  
  connectivity.pairs <- data.frame(  Pair.from = mpa.id.1,
                                     Pair.to = mpa.id.2,
                                     Number.events = mean(temp.result$Mean.events),
                                     Time.mean = mean(temp.result$Mean.Time),
                                     Time.max = mean(temp.result$Time.max),
                                     Time.sd = mean(temp.result$SD.Time),
                                     Probability = mean(temp.result$Probability) )
  
  return( connectivity.pairs )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

connectivity[ is.na(connectivity) ] <- 0
connectivity <- data.table(connectivity)
connectivity
save(connectivity,file="connectivity.mpas.Rdata")

## ------------------------------------------

# load("mpa.shp.notake.id.Rdata") OR run bellow

mpa.shp.centroids <- gCentroid(mpa.shp,byid=TRUE)

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

mpa.shp.notake.id <- foreach(i=1:length(mpa.shp.notake), .verbose=FALSE, .packages=c("sf","raster","rgeos","sp","bigmemory")) %dopar% { # 
  
  point.i <- gCentroid(mpa.shp.notake[i,],byid=TRUE)
  dist.i <- which.min(spDistsN1( as.matrix(as.data.frame(mpa.shp.centroids)) , as.matrix(as.data.frame(point.i)) , longlat =TRUE ))
  return( dist.i )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

mpa.shp.notake.id <- unlist(mpa.shp.notake.id)
mpa.shp.notake.id <- mpas.id[mpa.shp.notake.id]
save(mpa.shp.notake.id,file="mpa.shp.notake.id.Rdata")

## ------------------------------------------

plot(mpa.shp , col="gray",border="gray")
subset.mpa.shp <- mpa.shp[as.numeric(which(mpas.id %in% unique(mpa.shp.notake.id))),]
plot(subset.mpa.shp, col="black",add=TRUE)

## ------------------------------------------

## Choose first for all MPAs or second for FP MPAs

connectivity.matrix <- connectivity[,c(1,2,7)]
connectivity.matrix <- connectivity[ Pair.to %in% mpa.shp.notake.id & Pair.from %in% mpa.shp.notake.id , c(1,2,7) ]

connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )
diag(connectivity.matrix) <- 0

isolated.mpa <- colnames(connectivity.matrix)[which(apply(connectivity.matrix,1,sum,na.rm=T) == 0  & apply(connectivity.matrix,2,sum,na.rm=T) == 0)]

plot(mpa.shp , col="gray",border="gray")
subset.mpa.shp <- mpa.shp[as.numeric(which(mpas.id %in% isolated.mpa)),]
plot(subset.mpa.shp, col=sample(colours(), length(isolated.mpa))[1:length(isolated.mpa)],add=TRUE)

connect.index.1 <- data.frame(from=apply(connectivity.matrix,1,function(x){ sum(x != 0) } ) , to=apply(connectivity.matrix,2,function(x){ sum(x != 0) } ))

mean(unlist(connect.index.1))
sd(unlist(connect.index.1))
range(unlist(connect.index.1))

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

## Choose first for all MPAs or second for FP MPAs

comb <- connectivity[,c(1,2,7)]
comb <- connectivity[ Pair.to %in% mpa.shp.notake.id & Pair.from %in% mpa.shp.notake.id , c(1,2,7) ]

comb <- comb[ Pair.from != Pair.to ,]
comb <- as.data.frame( comb[ sort(comb[,Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )

# -------------

clustering.method <- "leading.eigenvector.community" # Uni: fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
graph.obj <- simplify(graph.obj)

clustering.method <- "clusters" # Uni: fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = comb[,3]
graph.obj=delete.edges(graph.obj, which(E(graph.obj)$weight ==0))

# -------------

membership.graph <- get(clustering.method)(graph.obj)$membership
modularity(graph.obj, membership.graph)
length(unique(membership.graph))

V(graph.obj)$color <- membership.graph

# 10 10 (maps 10 6)

plot(get(clustering.method)(graph.obj),graph.obj,vertex.label=NA)

subset.mpa.shp <- mpa.shp[ sapply( as.numeric(as_ids(V(graph.obj))) , function(x) { which(mpas.id %in% x ) }    ) ,]
subset.mpa.shp@data$COLOUR <- membership.graph

cols.to.use <- rainbow(length(unique(membership.graph)))[membership.graph]
plot(mpa.shp,col="gray",border="gray")
plot(subset.mpa.shp, col=cols.to.use,border=cols.to.use,add=T)

## ------------------------------------------

centrality <- eigen_centrality(graph.obj, directed = TRUE, scale = TRUE)$vector
centrality <- as.numeric(unlist(calculate_centralities(graph.obj, include = "Closeness Centrality (Freeman)")))

plot(mpa.shp, col="gray",border="gray")
subset.mpa.shp <- mpa.shp[as.numeric( which( centrality >=  as.numeric(quantile(centrality,0.95))  ) ) , ]
subset.mpa.shp@data$COLOUR <- 1
plot(subset.mpa.shp, col="red", border="red" , add=TRUE)

## ------------------------------------------

## Choose first for all MPAs or second for FP MPAs

comb <- connectivity[,c(1,2,7)]
comb <- connectivity[ Pair.to %in% mpa.shp.notake.id & Pair.from %in% mpa.shp.notake.id , c(1,2,7) ]

comb <- comb[ Pair.from != Pair.to ,]
comb <- as.data.frame( comb[ sort(comb[,Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )

# -------------

graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight == Inf))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
graph.obj <- simplify(graph.obj)

pairs.poly <- expand.grid(from=unique(mpa.shp.notake.id),to=unique(mpa.shp.notake.id))
pairs.poly <- pairs.poly[ pairs.poly$from != pairs.poly$to,]
pairs.poly$conn <- numeric(nrow(pairs.poly))

stones.list <- list()

for(p in 1:nrow(pairs.poly)) { 
  
  stones <- unlist(get.shortest.paths(graph.obj,as.character( pairs.poly[ p,1] ) , as.character( pairs.poly[ p,2] ),mode="out")$vpath)
  stones.list <- c(stones.list, list(names(stones)) )
  pairs.poly[ p,3] <- length(stones) - 1 
  
} 

sum(pairs.poly$conn == 0) / nrow(pairs.poly)
sum(pairs.poly$conn != 0) / nrow(pairs.poly)

stones <- unique(unlist(stones.list))
stones <- stones[ ! stones %in% unique(mpa.shp.notake.id) ] 
stones <- data.frame(stone=stones,times=sapply(stones , function(x) { sum(unlist(stones.list) == x  ) } ))
stones <- stones[sort(stones$times , index.return = T , decreasing=T)$ix,]

plot(mpa.shp, col="gray",border="gray")
subset.mpa.shp <- mpa.shp[ mpas.id %in% stones$stone, ]
plot(subset.mpa.shp, col="black", border="black",add=TRUE)
subset.mpa.shp <- mpa.shp[ mpas.id %in% mpa.shp.notake.id, ]
plot(subset.mpa.shp, col="red", border="red",add=TRUE)

plot(mpa.shp, col="gray",border="gray")
subset.mpa.shp <- mpa.shp[ mpas.id %in% stones$stone[stones$times >= quantile(stones$times,0.75)], ]
plot(subset.mpa.shp, col="black", border="black",add=TRUE)
subset.mpa.shp <- mpa.shp[ mpas.id %in% mpa.shp.notake.id, ]
plot(subset.mpa.shp, col="red", border="red",add=TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------------------------
