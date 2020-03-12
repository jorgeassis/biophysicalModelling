## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)
library(rnaturalearth)
library(geosphere)
source("0. Project Config.R")

sql.project.name <- "WAMPA"
number.cores <- 40

additional.landmass.shp <- "../Data/Spatial/surf_marine_2_new130520192.shp"
additional.landmass.shp <- gsub("Data/Spatial/","Final",additional.landmass.shp) 

sql.file <- "../Results/SQL/WAMPASimulationResults.sql"
bigmatrix.file <- "../InternalProc/particles.reference.desc"
sorce.sink.cells.file <- "../Results/source.sink.bm"

coordRef <- crs(shapefile(additional.landmass.shp))

## ------------------------------------------------------------------------------------------------------------
## Read main sources

sql <- dbConnect(RSQLite::SQLite(), sql.file)
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell[1]
n.new.particles.per.day <- dbReadTable(sql, "Parameters")$n.new.particles.per.day[1]
n.steps.per.day <- dbReadTable(sql, "Parameters")$n.hours.per.day[1]
dbDisconnect(sql)

source.sink.xy <- read.big.matrix("../Results/source.sink.bm")
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

source.sink.xy.sp <- source.sink.xy[,2:3]
coordinates(source.sink.xy.sp) <- ~Lon+Lat
crs(source.sink.xy.sp) <- coordRef

## ---------------------

polygon.shp <- shapefile(additional.landmass.shp)
polygon.shp <- gBuffer(polygon.shp, byid=TRUE, width=0)

worldMap <- ne_countries(scale = 50, returnclass = "sp")

plot(polygon.shp,col="Black",border="Black")
plot(worldMap,col="Gray",border="Gray",add=T)

plot(source.sink.xy.sp,col="Black",border="Black")

## --------------------------------------------------------------------
## --------------------------------------------------------------------
# Temporary Subset

subseter <- c(-23,40,20,47.5)

source.sink.xy <- source.sink.xy[source.sink.xy$Lon >= subseter[1] & source.sink.xy$Lon <= subseter[2] & source.sink.xy$Lat >= subseter[3] & source.sink.xy$Lat <= subseter[4], ]
plot(source.sink.xy[,2:3])

source.sink.xy.sp <- crop(source.sink.xy.sp,extent(subseter))
mpa.shp <- crop(mpa.shp,extent(source.sink.xy.sp))
polygon.shp <- crop(polygon.shp,extent(subseter))
plot(polygon.shp,col="Black",border="Black")
plot(source.sink.xy.sp,col="Black",border="Black")

## --------------------------------------------------------------------
## --------------------------------------------------------------------
## Produce main datasets for further queries

worldMap <- crop(worldMap,extent(polygon.shp) + c(-2.5,2.5,-3,3.75)) 
polygon.shp$ID <- 1:nrow(polygon.shp)

## ------------------------------------------------------------------------------------------------------------------------------
## Identify MPA id in source sink sites

# load OR run 1 foreach section bellow
# load("../Results/source.sink.xy.Rdata") OR run foreach section bellow

cl.2 <- makeCluster(20 , type="FORK")
registerDoParallel(cl.2)

source.sink.xy.polygon <- foreach( source.sink.xy.i = 1:nrow(source.sink.xy) , .verbose=FALSE, .combine = rbind ,  .packages=c("rgeos","raster","geosphere","FNN","bigmemory")) %dopar% { # 
  
  distances <- gDistance(source.sink.xy.sp[source.sink.xy.i], polygon.shp,byid=TRUE)
  allPolygoni <- polygon.shp[which.min(distances),]$ID
  
  return( data.frame( ID = allPolygoni ) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

source.sink.xy <- cbind(source.sink.xy[,.(Pair,Lon,Lat)],source.sink.xy.polygon)
head(source.sink.xy)

# save(source.sink.xy,file="../Results/source.sink.xy.Rdata")
# load("../Results/source.sink.xy.Rdata")

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

# List results

list.dirs(path = paste0("../Results"), recursive = FALSE)

season <- "YearRound" # c("YearRound","SeasonSummer","SeasonWinter")
pld.period <- 1:200 # c(10 , 30 , 90 , 120 , 200)
combinations <- expand.grid(season=season,pld.period=pld.period,stringsAsFactors = F)

## --------------------

polygonNames <- polygon.shp$ID
polygonNamesDuplicated <- which(duplicated(polygonNames))

for( i in 1:length(polygonNamesDuplicated) ) {
  
  t <- which(polygonNames %in% MPAnames[polygonNamesDuplicated[i]] )
  t <- t[-1]
  
  for( t.i in 1:length(t)) {
    
    polygonNames[t[t.i]] <- paste0(polygonNames[t[t.i]]," ", t.i,collapse = "")
    
  }
  
}

## --------------------

combResults <- data.frame()

isolatedResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(isolatedResults) <- polygonNames
colnames(isolatedResults) <- 1:nrow(combinations)

betweennessResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(betweennessResults) <- polygonNames
colnames(betweennessResults) <- 1:nrow(combinations)

higherBetweennessResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherBetweennessResults) <- polygonNames
colnames(higherBetweennessResults) <- 1:nrow(combinations)

eighenCentralityResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(eighenCentralityResults) <- polygonNames
colnames(eighenCentralityResults) <- 1:nrow(combinations)

highereighenCentralityResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(highereighenCentralityResults) <- polygonNames
colnames(highereighenCentralityResults) <- 1:nrow(combinations)

higherStoneResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherStoneResults) <- polygonNames
colnames(higherStoneResults) <- 1:nrow(combinations)

resistanceResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(resistanceResults) <- polygonNames
colnames(resistanceResults) <- 1:nrow(combinations)

higherResistanceResults <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherResistanceResults) <- polygonNames
colnames(higherResistanceResults) <- 1:nrow(combinations)

clusterAssignment <- data.frame(matrix(nrow=length(polygonNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(clusterAssignment) <- polygonNames
colnames(clusterAssignment) <- 1:nrow(combinations)

## ------------------------------------------------------------------------------

for( c in 1:nrow(combinations)){
  
  season <- combinations[c,1]
  pld.period <- combinations[c,2]
  
  if( season == "SeasonSummer" ) { spawn.p <- c(6,7,8,9) }
  if( season == "SeasonWinter" ) { spawn.p <- c(11,12,1,2) }
  if( season == "YearRound" ) { spawn.p <- 1:12 }
  
  project.name <- paste0(season,"_Pld",pld.period)
  
  ## ----------------------------------------------------
  
  if( ! dir.exists(paste0("../Results/",project.name)) ) { dir.create(file.path(paste0("../Results/",project.name)), showWarnings = FALSE) } 
  if( ! dir.exists(paste0("../Results/",project.name,"/Data")) ) { dir.create(file.path(paste0("../Results/",project.name,"/Data")), showWarnings = FALSE) } 
  if( ! dir.exists(paste0("../Results/",project.name,"/Maps")) ) { dir.create(file.path(paste0("../Results/",project.name,"/Maps")), showWarnings = FALSE) } 
  
  particles.reference.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.reference.desc"))
  cell.to.process <- unique(source.sink.xy$Pair)

  cl.2 <- makeCluster(40 , type="FORK")
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
  
  # save(connectivity.source.sink.xy,file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.xy.Rdata"))
  # project.name
  # load(file=paste0("../Results/YearRound_Pld144/Data/connectivity.source.sink.xy.Rdata"))
  
  ## ----------------------------------------------------
  
  polygonIDPairs <- expand.grid(from= polygon.shp$ID,to= polygon.shp$ID)
  polygonsCompute <- polygon.shp


  cl.2 <- makeCluster(40 , type="FORK")
  registerDoParallel(cl.2)
  
  connectivity <- foreach( pairs = 1:nrow(polygonIDPairs) , .verbose=FALSE, .combine = rbind ,  .packages=c("geosphere","rgeos","raster","data.table","FNN","bigmemory")) %dopar% { # 
    
    mpa.id.1 <- polygonIDPairs[pairs,1]
    mpa.id.2 <- polygonIDPairs[pairs,2]
    
    mpa.cells.1 <- as.vector(unlist(source.sink.xy[ ID == mpa.id.1 , 1 ]))
    mpa.cells.2 <- as.vector(unlist(source.sink.xy[ ID == mpa.id.2 , 1 ]))
    
    temp.result <- connectivity.source.sink.xy[ Pair.from %in% mpa.cells.1 & Pair.to %in% mpa.cells.2 , ]
    
    if(nrow(temp.result) == 0) {
      
      connectivity.pairs <- data.frame(  Pair.from = mpa.id.1,
                                         Pair.to = mpa.id.2,
                                         Number.events = 0,
                                         Time.mean = 0,
                                         Time.max = 0,
                                         Time.sd = 0,
                                         Probability = 0 )
      
    }
    
    if(nrow(temp.result) > 0) {
      
      connectivity.pairs <- data.frame(  Pair.from = mpa.id.1,
                                         Pair.to = mpa.id.2,
                                         Number.events = mean(temp.result$Mean.events),
                                         Time.mean = mean(temp.result$Mean.Time),
                                         Time.max = mean(temp.result$Max.Time),
                                         Time.sd = mean(temp.result$SD.Time),
                                         Probability = mean(temp.result$Mean.Probability) )
      
    }
    
    return( connectivity.pairs )
    
  }
  
  stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)
  
  # save(connectivity.notakeMPA,file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.notakeMPA.Rdata"))
  # project.name
  # load(file=paste0("../Results/YearRound_Pld144/Data/connectivity.source.sink.notakeMPA.Rdata"))
  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  
  connectivity.matrix <- connectivity[ ,c(1,2,7)]
  connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )
  diag(connectivity.matrix) <- 0
  
  isolated.polygon <- which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)
  isolated.polygon.id <- colnames(connectivity.matrix)[isolated.polygon]
  isolated.polygon.names <- polygonNames[which( colnames(connectivity.matrix) %in% isolated.polygon.id)]
  isolated.polygon.length <- length(isolated.polygon)
  
  isolatedResults[ which(rownames(isolatedResults) %in% isolated.polygon.names ),c] <- "TRUE"
  write.csv(isolatedResults,file="../Results/isolatedMPAs.csv")
  
  # Aggregation level (Proportion of non-isolated MPAs, at least one connection, in relation to the number of MPAs)
  aggregationAtLeastOne <- (nrow(connectivity.matrix)- isolated.polygon.length ) / nrow(connectivity.matrix)
  
  # Aggregation level (Based on overall connections)
  connectivity.matrix.binomial <- connectivity.matrix
  connectivity.matrix.binomial[connectivity.matrix.binomial != 0] <- 1
  aggregationAllConnections <- sum ( ( apply(connectivity.matrix.binomial,1,sum) + 1 ) / nrow(connectivity.matrix.binomial) ) / nrow(connectivity.matrix.binomial)
  
  # Resistance index
  resistance <- sapply( 1:nrow(connectivity.matrix.binomial) , function(res) { 1- length(unique(c(which(connectivity.matrix.binomial[res,] == 1 ) ,  which(connectivity.matrix.binomial[,res] == 1 )))) / nrow(connectivity.matrix.binomial) } )
  resistance[resistance == 1] <- NA
  resistanceResults[ ,c] <- resistance
  write.csv(resistanceResults,file="../Results/resistanceMPA.csv")
  
  # Those with higher / Percertil 95%
  temporaryRes <- polygonNames[ which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=TRUE))) ]
  higherResistanceResults[ as.numeric(unlist(sapply(  temporaryRes  , function(x)  which( rownames(higherResistanceResults) == x )))) ,c] <- "TRUE"
  write.csv(higherResistanceResults,file="../Results/higheResistanceMPA..csv")
  
  # Average connections between MPAs
  connect.index <- data.frame(MPA=colnames(connectivity.matrix),exportTo=apply(connectivity.matrix,1,function(x){ sum(x != 0) } ) , importFrom=apply(connectivity.matrix,2,function(x){ sum(x != 0) } ))
  averageConnections <- mean(unlist(connect.index[,-1]))
  maximumConnections <- max(apply(connect.index[,-1],1,max))
  
  # Export
  averageExport <- mean(connect.index[connect.index[,2] != 0,2])
  histData <- unlist(connect.index[,2])
  
  pdf(file=paste0("../Results/",project.name,"/histogramExportTo.pdf"), width=5, height=5)
  par(mar = c(4.5, 5.5, 4.5, 4.5))
  hist( histData[histData != 0], breaks=unique(round(seq(0,max(histData),length.out = 10))) , xlab="Number of MPAs",main=NULL)
  dev.off()
  
  # import
  averageImport <- mean(connect.index[connect.index[,3] != 0,3])
  histData <- unlist(connect.index[,3])
  
  pdf(file=paste0("../Results/",project.name,"/histogramImportFrom.pdf"), width=5, height=5)
  par(mar = c(4.5, 5.5, 4.5, 4.5))
  hist( histData[histData != 0], breaks=unique(round(seq(0,max(histData),length.out = 10))) , xlab="Number of MPAs",main=NULL)
  dev.off()
  
  ## --------------------------------------------------------------------------------
  ## --------------------------------------------------------------------------------
  
  comb <- connectivity
  comb <- comb[ comb$Pair.from != comb$Pair.to ,]
  comb <- as.data.frame( comb[ sort(comb[,"Probability"] , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Probability")] )
  
  # -------------
  
  graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
  graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
  graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
  graph.obj <- simplify(graph.obj)
  
  # -------------
  
  clustering.method <- "leading.eigenvector.community" # CLUSTERS ?? Uni: leading.eigenvector.community fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)
  clustering.graph <- get(clustering.method)(graph.obj,options=list(maxiter=1000000))
  membership.graph <- clustering.graph$membership
  # membership.graph <- membership.graph[sort(as.numeric( V(graph.obj) ),decreasing = FALSE,index.return=TRUE)$ix]
  
  modularity.val <- modularity(graph.obj, membership.graph) # Measuring goodness of fit of clustering (e.g, leading.eigenvector.community)
  
  clustersAll <- clusters(graph.obj)$membership
  membership.graph2 <- clustersAll
  clustersAll <- clustersAll[sort(as.numeric(names(clustersAll)),decreasing = FALSE,index.return=TRUE)$ix]
  
  clusterAssignment[ ,c] <- clustersAll
  write.csv(clusterAssignment,file="../Results/clusterAssignmentMPA.csv")
  
  # Number of clusters
  numberClustersLeading.eigenvector <- length(unique(membership.graph))
  numberClusters <- length(unique(clusters(graph.obj)$membership))
  
  # Aggregation factor based on clusters 
  # aggregationBasedOnLeading.eigenvector <- ( length(membership.graph) / numberClustersLeading.eigenvector ) / length(membership.graph)
  # aggregationBasedOnClusters <- ( length(membership.graph) / numberClusters ) / length(membership.graph)
  
  aggregationBasedOnLeading.eigenvector <- 1 - ( numberClustersLeading.eigenvector / length(membership.graph) )
  aggregationBasedOnClusters <- 1 - ( numberClusters / length(membership.graph) )
  
  V(graph.obj)$color <- "#000000"
  
  pdf(file=paste0("../Results/",project.name,"/Clustering.pdf"), width=11, height=11)
  plot(clustering.graph,graph.obj,vertex.label=NA,vertex.size=2,edge.curved = T  ) 
  title(paste0("Particle duration: ",c," day",ifelse(c==1,"","s")), line = 2,col="#797979" , cex=0.75)
  dev.off()
  
  polygon.shp.sorted <- polygon.shp[ sapply( as.numeric(as_ids(V(graph.obj))) , function(x) { which(polygon.shp$ID %in% x ) } ) ,]
  polygon.shp.sorted@data$COLOUR <- membership.graph2
  
  ## ------------------------------------------
  
  # Plot clusters with connections
  
  cols.to.use <- distinctColors(length(unique(membership.graph2)))[membership.graph2]
  cols.to.use[ which(names(membership.graph2) %in% isolated.polygon.id)] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapClusteringConnections.pdf"), width=11, height=8)
  plot(worldMap , col="#E8E8E8",border="#C9C9C9")
  text(-25.2, y = 24, labels = paste0("Particle duration: ",c," day",ifelse(c==1,"","s")  ) , col="#5E5E5E" , cex=0.9)
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(polygon.shp$ID == connected.pairs[i,1]),],
                                centroids[ which(polygon.shp$ID == connected.pairs[i,2]),],
                                n = 100, addStartEnd = TRUE, sp = TRUE)
    lines(  routes_sl , type="l" , col=colfunc(101)[round(strenght)])
  }
  
  subset.Centroids <- as.data.frame(gCentroid(polygon.shp.sorted,byid=TRUE),xy=T)
  points(subset.Centroids,col="#151515",bg=cols.to.use,pch=21, cex = 0.8)
  dev.off()
  
  ## ------------------------------------------
  
  # Plot eigenvector with connections
  
  cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
  cols.to.use[which(names(membership.graph2) %in% isolated.polygon.id)] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapClusteringEigenvectorConnections.pdf"), width=11, height=8)
  plot(worldMap , col="#E8E8E8",border="#C9C9C9")
  text(-25.2, y = 24, labels = paste0("Particle duration: ",c," day",ifelse(c==1,"","s")  ) , col="#5E5E5E" , cex=0.9)
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(polygon.shp$ID == connected.pairs[i,1]),],
                                centroids[ which(polygon.shp$ID == connected.pairs[i,2]),],
                                n = 100, addStartEnd = TRUE, sp = TRUE)
    lines(  routes_sl , type="l" , col=colfunc(101)[round(strenght)])
  }
  
  subset.Centroids <- as.data.frame(gCentroid(polygon.shp.sorted,byid=TRUE),xy=T)
  points(subset.Centroids,col="#151515",bg=cols.to.use,pch=21, cex = 0.8)
  dev.off()
  
  ## ------------------------------------------
  
  # Plot isolated with connections
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapIsolatedConnMPA.pdf"), width=11, height=8)
  plot(worldMap , col="#E8E8E8",border="#C9C9C9")
  text(-25.2, y = 24, labels = paste0("Particle duration: ",c," day",ifelse(c==1,"","s")  ) , col="#5E5E5E" , cex=0.9)
  non.isolated.mpaCentroids <- as.data.frame(gCentroid(polygon.shp[which(! colnames(connectivity.matrix) %in% isolated.polygon.id),],byid=TRUE),xy=T)
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  
  points(non.isolated.mpaCentroids,col="#151515",bg="white",pch=21, cex = 0.8)
  isolated.mpaCentroids <- as.data.frame(gCentroid(polygon.shp[which(colnames(connectivity.matrix) %in% isolated.polygon.id),],byid=TRUE),xy=T)
  points(isolated.mpaCentroids,col="#151515",bg="#9C2323",pch=21, cex = 0.8)
  dev.off()
  
  ## ------------------------------------------
  
  # Eigen centrality vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).
  # Cloness centrality measures how many steps is required to access every other vertex from a given vertex
  # Betweenness centrality is (roughly) defined by the number of geodesics (shortest paths) going through a vertex or an edge.
  
  betweennessIndex <- betweenness(graph.obj)
  betweennessIndex <- betweennessIndex[sort(as.numeric(names(betweennessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  averageBetweenness <- mean(betweennessIndex)
  
  eigen_centralityIndex <- eigen_centrality(graph.obj, directed = TRUE, scale = TRUE)$vector
  eigen_centralityIndex <- eigen_centralityIndex[sort(as.numeric(names(eigen_centralityIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  averageeigen_centrality <- mean(eigen_centralityIndex)
  
  closenessindex <- closeness(graph.obj)
  averageCloseness <- mean(closenessindex)
  
  betweennessResults[ ,c] <- betweennessIndex
  write.csv(betweennessResults,file="../Results/betweennessMPAs.csv")
  
  eighenCentralityResults[ ,c] <- eigen_centralityIndex
  write.csv(eighenCentralityResults,file="../Results/isolatedMPAs.csv")
  
  # Those with higher / Percertil 95%
  temporaryRes <- polygonNames[ as.numeric(names(which(betweennessIndex >=  as.numeric(quantile(betweennessIndex,0.95))))) ]
  higherBetweennessResults[ as.numeric(unlist(sapply(  temporaryRes  , function(x)  which( rownames(higherBetweennessResults) == x )))) ,c] <- "TRUE"
  write.csv(higherBetweennessResults,file="../Results/higherBetweennessMPA.csv")
  
  temporaryRes <- polygonNames[ as.numeric(names(which(eigen_centralityIndex >=  as.numeric(quantile(eigen_centralityIndex,0.95))))) ]
  highereighenCentralityResults[ as.numeric(unlist(sapply(  temporaryRes  , function(x)  which( rownames(highereighenCentralityResults) == x )))) ,c] <- "TRUE"
  write.csv(highereighenCentralityResults,file="../Results/higherBetweennessMPA.csv")
  
  # Plot Betweenness Connections
  
  polygon.shp.sorted <- polygon.shp[  as.numeric(names(which(betweennessIndex >=  as.numeric(quantile(betweennessIndex,0.95))))) , ]
  rbPal <- colorRampPalette(c('#BAE2FF','yellow','orange'))
  rbPal <- rbPal(100)[as.numeric(cut(betweennessIndex,breaks = 100))]
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapBetweennessConnections.pdf"), width=11, height=8)
  plot(worldMap , col="#E8E8E8",border="#C9C9C9")
  text(-25.2, y = 24, labels = paste0("Particle duration: ",c," day",ifelse(c==1,"","s")  ) , col="#5E5E5E" , cex=0.9)
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(polygon.shp$ID == connected.pairs[i,1]),],
                                centroids[ which(polygon.shp$ID == connected.pairs[i,2]),],
                                n = 100, addStartEnd = TRUE, sp = TRUE)
    lines(  routes_sl , type="l" , col=colfunc(101)[round(strenght)])
  }
  subset.shpCentroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  points(isolated.mpaCentroids,col="#151515",bg="white",pch=21, cex = 0.8)
  points(subset.shpCentroids,col="#151515",bg=rbPal,pch=21, cex = 0.8)
  subset.shpCentroids <- as.data.frame(gCentroid(polygon.shp.sorted,byid=TRUE),xy=T)
  points(subset.shpCentroids,col="#151515",bg="#9C2323",pch=21, cex = 0.8)
  dev.off()
  
  ## ------------------------------------------
  
  # Plot eighen Centrality Connections
  
  polygon.shp.sorted <- polygon.shp[  as.numeric(names(which(eigen_centralityIndex >=  as.numeric(quantile(eigen_centralityIndex,0.95))))) , ]
  rbPal <- colorRampPalette(c('#BAE2FF','yellow','orange'))
  rbPal <- rbPal(100)[as.numeric(cut(eigen_centralityIndex,breaks = 100))]
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapEighenCentralityConnections.pdf"), width=11, height=8)
  plot(worldMap , col="#E8E8E8",border="#C9C9C9")
  text(-25.2, y = 24, labels = paste0("Particle duration: ",c," day",ifelse(c==1,"","s")  ) , col="#5E5E5E" , cex=0.9)
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(polygon.shp$ID == connected.pairs[i,1]),],
                                centroids[ which(polygon.shp$ID == connected.pairs[i,2]),],
                                n = 100, addStartEnd = TRUE, sp = TRUE)
    lines(  routes_sl , type="l" , col=colfunc(101)[round(strenght)])
  }
  subset.shpCentroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  points(isolated.mpaCentroids,col="#151515",bg="white",pch=21, cex = 0.8)
  points(subset.shpCentroids,col="#151515",bg=rbPal,pch=21, cex = 0.8)
  subset.shpCentroids <- as.data.frame(gCentroid(polygon.shp.sorted,byid=TRUE),xy=T)
  points(subset.shpCentroids,col="#151515",bg="#9C2323",pch=21, cex = 0.8)
  dev.off()
  
  ## ------------------------------------------
  
  polygon.shp.sorted <- polygon.shp[ which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=T))) , ]
  rbPal <- colorRampPalette(c('#BAE2FF','yellow','orange'))
  rbPal <- rbPal(100)[as.numeric(cut(resistance,breaks = 100))]
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapResistance.pdf"), width=11, height=8)
  plot(worldMap , col="#E8E8E8",border="#C9C9C9")
  text(-25.2, y = 24, labels = paste0("Particle duration: ",c," day",ifelse(c==1,"","s")  ) , col="#5E5E5E" , cex=0.9)
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(polygon.shp$ID == connected.pairs[i,1]),],centroids[ which(polygon.shp$ID == connected.pairs[i,2]),],n = 100, addStartEnd = TRUE, sp = TRUE)
    lines(  routes_sl , type="l" , col=colfunc(101)[round(strenght)])
  }
  
  subset.shpCentroids <- as.data.frame(gCentroid(polygon.shp,byid=TRUE),xy=T)
  points(isolated.mpaCentroids,col="#151515",bg="white",pch=21, cex = 0.8)
  points(subset.shpCentroids,col="#151515",bg=rbPal,pch=21, cex = 0.8)
  subset.shpCentroids <- as.data.frame(gCentroid(polygon.shp.sorted,byid=TRUE),xy=T)
  points(subset.shpCentroids,col="#151515",bg="#9C2323",pch=21, cex = 0.8)
  dev.off()
  
  ## ------------------------------------------
  
  combResults <- rbind(combResults,
                       data.frame(pld=pld.period,
                                  n.isolated=isolated.polygon.length,
                                  aggregationAtLeastOne=aggregationAtLeastOne,
                                  aggregationAllConnections=aggregationAllConnections,
                                  averageConnections=averageConnections,
                                  maximumConnections=maximumConnections,
                                  averageResistance=mean(resistance,na.rm=T),
                                  averageExport=averageExport,
                                  averageImport=averageImport,
                                  numberClusters=numberClusters,
                                  numberClustersLeading.eigenvector=numberClustersLeading.eigenvector,
                                  modularity.val=modularity.val,
                                  aggregationBasedOnClusters=aggregationBasedOnClusters,
                                  aggregationBasedOnLeading.eigenvector=aggregationBasedOnLeading.eigenvector,
                                  averageeigen_centrality=averageeigen_centrality,
                                  averageBetweenness=averageBetweenness,
                                  averageCloseness=averageCloseness ) )
  
  write.csv(combResults,file="../Results/Results.csv")
  save(combResults,file=paste0("../Results/allPLDResults.Rdata"))
  
}

## ---------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------

x <- combResults$pld
x.lab <- "Propagule duration (day)"

colnames(combResults)

y <- combResults$n.isolated # 
y.lab <- "Isolated MPA (number)"
file.name <- "isolated"

y <- combResults$aggregationAtLeastOne # 
y.lab <- "Aggregation level (proportion of connection)"
file.name <- "aggregationAtLeastOne"

y <- combResults$aggregationAllConnections # 
y.lab <- "Aggregation level (proportion of connection)"
file.name <- "aggregationAllConnections"

y <- combResults$averageConnections # 
y.lab <- "Aggregation level (average number of connections)"
file.name <- "averageConnections"

y <- combResults$averageResistance # 
y.lab <- "Resistance level (average MPA resistance)"
file.name <- "averageResistance"

y <- combResults$numberClusters -  combResults$n.isolated # 
y.lab <- "Aggregation level (number of clusters)"
file.name <- "numberClusters"

y <- combResults$numberClustersLeading.eigenvector -  combResults$n.isolated# 
y.lab <- "Aggregation level (number of clusters based on leading eigenvector)"
file.name <- "numberClustersLeading.eigenvector"

y <- combResults$averageeigen_centrality # 
y.lab <- "Centrality (average MPA centrality)"
file.name <- "averageeigen_centrality"

y <- combResults$averageBetweenness # 
y.lab <- "Betweenness (average MPA betweenness)"
file.name <- "averageBetweenness"

y <- combResults$averageCloseness # 
y.lab <- "Closeness (average MPA closeness)"
file.name <- "averageCloseness"

# ---------

pdf(file=paste0("../Results/Plots/",file.name,".pdf"), width=12, height=5)
par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(x,y,pch=20,col="#A6A6A6", ylab="",xlab=x.lab,axes=FALSE)
title(ylab=y.lab, line=4)
lines(bezierCurve(x,y,100)$x,bezierCurve(x,y,100)$y,type="l", lwd=1, lty=2)
axis(2,las=2,col="White",col.ticks="Black", cex.axis=0.9)
axis(1,las=0,col="White",col.ticks="Black", cex.axis=0.9)
box()
dev.off()

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------

boxplotData <- data.frame( Area=area(polygon.shp) , Isolation=ifelse(isolatedResults[,ncol(isolatedResults)] == "TRUE","Isolated","Connected"))
boxplotData <- boxplotData[sort(boxplotData$Area,decreasing = TRUE,index.return=TRUE)$ix,]
boxplotData[,1] <- boxplotData[,1] / 1000000

mean(boxplotData[boxplotData$Isolation == "Isolated","Area"])
sd(boxplotData[boxplotData$Isolation == "Isolated","Area"])

mean(boxplotData[boxplotData$Isolation != "Isolated","Area"])
sd(boxplotData[boxplotData$Isolation != "Isolated","Area"])

kruskal.test( c(boxplotData[boxplotData$Isolation != "Isolated","Area"] , boxplotData[boxplotData$Isolation == "Isolated","Area"]) , c(rep(1,sum(boxplotData$Isolation != "Isolated")),rep(2,sum(boxplotData$Isolation == "Isolated"))) )
# Kruskal-Wallis chi-squared = 55.247, df = 1, p-value = 1.063e-13

p <- ggplot(boxplotData, aes(x=Isolation, y=Area)) + 
  geom_boxplot()
p

p <- ggplot(boxplotData, aes(x=Isolation, y=Area)) + 
  geom_boxplot() + ylim(0, 1100) + xlab("Isolation (level)") + ylab("MPA area (km2)")
p

pdf(file=paste0("../Results/Plots/AreaVsIsolationNoOut.pdf"), width=10, height=7)
ggplot(boxplotData, aes(x=Isolation, y=Area)) + 
geom_boxplot() + ylim(0, 1100) + xlab("Isolation (level)") + ylab("MPA area (km2)")
dev.off()

pdf(file=paste0("../Results/Plots/AreaVsIsolationAll.pdf"), width=10, height=7)
ggplot(boxplotData, aes(x=Isolation, y=Area)) + 
geom_boxplot() + xlab("Isolation (level)") + ylab("MPA area (km2)")
dev.off()
