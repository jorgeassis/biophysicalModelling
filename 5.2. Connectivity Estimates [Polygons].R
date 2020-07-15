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

sql.project.name <- "MPA"
number.cores <- 8

notakeMPA <- "../Data/notake_merged_Med_P.shp" 

sql.file <- "../Results/SQL/CentralityMPASimulationResults.sql"
bigmatrix.file <- "../InternalProc/particles.reference.desc"
sorce.sink.cells.file <- "../Results/source.sink.bm"

coordRef <- crs(shapefile(notakeMPA))

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

notakeMPA <- shapefile(notakeMPA)
notakeMPA <- gBuffer(notakeMPA, byid=TRUE, width=0)

worldMap <- ne_countries(scale = 10, returnclass = "sp")

## --------------------------------------------------------------------
## --------------------------------------------------------------------
# Temporary Subset

subseter <- c(-11.25,37.85,29.75,46.25)

## ------------

source.sink.xy <- source.sink.xy[source.sink.xy$Lon >= subseter[1] & source.sink.xy$Lon <= subseter[2] & source.sink.xy$Lat >= subseter[3] & source.sink.xy$Lat <= subseter[4], ]
plot(source.sink.xy[,2:3])

source.sink.xy.sp <- crop(source.sink.xy.sp,extent(subseter))
notakeMPA <- crop(notakeMPA,extent(subseter))
notakeMPA$ID <- 1:nrow(notakeMPA)

worldMap <- crop(worldMap,extent(subseter + c(-20,20,-20,20))  ) 

## -----------------

plot(notakeMPA , col="Black",border="Black")

## ------------------------------------------------------------------------------------------------------------------------------
## Identify MPA id in source sink sites

# load OR run 1 foreach section bellow
# load("../Results/source.sink.xy.Rdata") OR run foreach section bellow

cl.2 <- makeCluster(6 , type="FORK")
registerDoParallel(cl.2)

source.sink.xy.mpa <- foreach( source.sink.xy.i = 1:nrow(source.sink.xy) , .verbose=FALSE, .combine = rbind ,  .packages=c("rgeos","raster","geosphere","FNN","bigmemory")) %dopar% { # 
  
  distances <- gDistance(source.sink.xy.sp[source.sink.xy.i], notakeMPA,byid=TRUE)
  notakeMPAi <- notakeMPA[which.min(distances),]$ID
  
  return( data.frame( notakeMPAID = notakeMPAi  ) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

source.sink.xy <- cbind(source.sink.xy[,.(Pair,Lon,Lat)],source.sink.xy.mpa)
head(source.sink.xy)
save(source.sink.xy,file="../Results/source.sink.xy.Rdata")

load("../Results/source.sink.xy.Rdata")

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

# List results

list.dirs(path = paste0("../Results"), recursive = FALSE)

season <- "YearRound" # c("YearRound","SeasonSummer","SeasonWinter")
pld.period <- 1:200 # c(10 , 30 , 90 , 120 , 200)

combinations <- expand.grid(season=season,pld.period=pld.period,stringsAsFactors = F)

## --------------------

MPAnames <- notakeMPA$NAME

combResults <- data.frame()

isolatedResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(isolatedResults) <- MPAnames
colnames(isolatedResults) <- 1:nrow(combinations)

betweennessResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(betweennessResults) <- MPAnames
colnames(betweennessResults) <- 1:nrow(combinations)

higherBetweennessResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherBetweennessResults) <- MPAnames
colnames(higherBetweennessResults) <- 1:nrow(combinations)

eighenCentralityResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(eighenCentralityResults) <- MPAnames
colnames(eighenCentralityResults) <- 1:nrow(combinations)

highereighenCentralityResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(highereighenCentralityResults) <- MPAnames
colnames(highereighenCentralityResults) <- 1:nrow(combinations)

closenessResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(closenessResults) <- MPAnames
colnames(closenessResults) <- 1:nrow(combinations)

higherclosenessResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherclosenessResults) <- MPAnames
colnames(higherclosenessResults) <- 1:nrow(combinations)

clusterAssignment <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(clusterAssignment) <- MPAnames
colnames(clusterAssignment) <- 1:nrow(combinations)

resistanceResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(resistanceResults) <- MPAnames
colnames(resistanceResults) <- 1:nrow(combinations)

higherResistanceResults <- data.frame(matrix(nrow=length(MPAnames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherResistanceResults) <- MPAnames
colnames(higherResistanceResults) <- 1:nrow(combinations)

## ------------------------------------------------------------------------------

for( c in c(2,5,17,36) ){ # 1:nrow(combinations)
  
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
  
  cl.2 <- makeCluster(8 , type="FORK")
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
  save(connectivity.source.sink.xy,file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.xy.Rdata"))
  
  ## ----------------------------------------------------
  
  mpaIDPairs <- expand.grid(from=notakeMPA$ID,to=notakeMPA$ID)
  polygonsCompute <- mpaIDPairs
  
  cl.2 <- makeCluster(8 , type="FORK")
  registerDoParallel(cl.2)
  
  connectivity <- foreach( pairs = 1:nrow(mpaIDPairs) , .verbose=FALSE, .combine = rbind ,  .packages=c("geosphere","rgeos","raster","data.table","FNN","bigmemory")) %dopar% { # 
    
    mpa.id.1 <- mpaIDPairs[pairs,1]
    mpa.id.2 <- mpaIDPairs[pairs,2]
    
    mpa.cells.1 <- as.vector(unlist(source.sink.xy[ notakeMPAID == mpa.id.1 , 1 ]))
    mpa.cells.2 <- as.vector(unlist(source.sink.xy[ notakeMPAID == mpa.id.2 , 1 ]))
    
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
  
  save(connectivity,file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.notakeMPA.Rdata"))
  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  
  # project.name
  
  load(file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.notakeMPA.Rdata"))
  load(file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.xy.Rdata"))
  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  
  connectivity.matrix <- connectivity[ ,c(1,2,7)]
  connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )
  diag(connectivity.matrix) <- 0
  
  write.csv(connectivity.matrix,paste0("../Results/",project.name,"/connectivitymatrix.csv"))
  
  isolated.mpa <- which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)
  isolated.mpa.id <- colnames(connectivity.matrix)[isolated.mpa]
  isolated.mpa.names <- MPAnames[which( colnames(connectivity.matrix) %in% isolated.mpa.id)]
  isolated.mpa.length <- length(isolated.mpa)
  
  isolatedResults[ which(rownames(isolatedResults) %in% isolated.mpa.names ),c] <- "TRUE"
  write.csv(isolatedResults,file="../Results/isolatedMPAs.csv")
  
  # Aggregation level (Proportion of non-isolated MPAs, at least one connection, in relation to the number of MPAs)
  aggregationAtLeastOne <- (nrow(connectivity.matrix)-length(isolated.mpa)) / nrow(connectivity.matrix)
  
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
  temporaryRes <- MPAnames[ which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=TRUE))) ]
  higherResistanceResults[ as.numeric(unlist(sapply(  temporaryRes  , function(x)  which( rownames(higherResistanceResults) == x )))) ,c] <- "TRUE"
  write.csv(higherResistanceResults,file="../Results/higheResistanceMPA..csv")
  
  # Average connections between MPAs
  connect.index <- data.frame(MPA=colnames(connectivity.matrix),exportTo=apply(connectivity.matrix,1,function(x){ sum(x != 0) } ) , importFrom=apply(connectivity.matrix,2,function(x){ sum(x != 0) } ))
  
  averageConnections <- mean(unlist(connect.index[,-1]))
  sdConnections <- sd(unlist(connect.index[,-1]))
  maximumConnections <- max(apply(connect.index[,-1],1,max))
  
  ## --------------------------------------------------------------------------------
  ## --------------------------------------------------------------------------------
  
  comb <- connectivity
  comb <- comb[ comb$Pair.from != comb$Pair.to ,]
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
  
  mapSouthernEurope <- ggplot() +
    geom_polygon(data = worldMap , fill = "#C4C4C4", colour = "#ffffff" , size=0.25 ,  aes(long, lat, group = group)) +
    coord_map('lambert', lat0=29.75, lat1=45, xlim=c(-7.70, 35.75), ylim=c(29.75, 44.5)) + theme_map
  
  mapSouthernEuropeNet <- mapSouthernEurope
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(notakeMPA,byid=TRUE),xy=T)
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(notakeMPA$ID == connected.pairs[i,1]),],
                                centroids[ which(notakeMPA$ID == connected.pairs[i,2]),],
                                n = 100, addStartEnd = TRUE, sp = TRUE)
    SLDF = sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = NA), match.ID = F)
    mapSouthernEuropeNet <- mapSouthernEuropeNet + geom_path(data = SLDF, size=0.35 , aes(x = long, y = lat), col = "#797979") # colfunc(101)[round(strenght)]
  }
  mapSouthernEuropeNet
  
  ## --------------------------------------------------------------------------------
  
  graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  #E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
  E(graph.obj)$weight = comb[,3] # Hock, Karlo Mumby, Peter J 2015
  graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
  graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
  graph.obj <- simplify(graph.obj)
  
  # -------------
  
  membership.graph <- clusters(graph.obj)$membership
  membership.graph <- membership.graph[sort(as.numeric(names(membership.graph)),decreasing = FALSE,index.return=TRUE)$ix]
  
  clusterAssignment[ ,c] <- membership.graph
  write.csv(clusterAssignment,file="../Results/clusterAssignmentMPA.csv")
  
  # Number of clusters
  numberClusters <- length(unique(membership.graph)) - isolated.mpa.length
  
  aggregationBasedOnClusters <- 1 - ( numberClusters / length(membership.graph) )
  
  ## ------------------------------------------
  
  # Plot clusters with connections
  
  cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
  cols.to.use[which(names(membership.graph) %in% isolated.mpa.id)] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapClusteringConnections.pdf"), width=12)
  
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = 2.5, stroke = 0.35, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = 2.5, stroke = 0.35, alpha = 0.7)
  
  dev.off()
  
  ## ------------------------------------------
  
  # Plot isolated with connections
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapIsolatedConnMPA.pdf"), width=12)
  
  mapSouthernEuropeNet + 
    geom_point(data = centroids ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.7) +
    geom_point(data = centroids[as.numeric(isolated.mpa.id),] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "#9C2323", size = 2.5, stroke = 0.35, alpha = 0.9)
      
  dev.off()
  
  ## ------------------------------------------
  # 
  # betweennessIndex <- betweenness(graph.obj)
  # betweennessIndex <- betweennessIndex[sort(as.numeric(names(betweennessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  # averageBetweenness <- mean(betweennessIndex)
  # 
  # eigen_centralityIndex <- eigen_centrality(graph.obj, directed = FALSE, scale = FALSE)$vector
  # eigen_centralityIndex <- eigen_centralityIndex[sort(as.numeric(names(eigen_centralityIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  # averageeigen_centrality <- mean(eigen_centralityIndex)
  # 
  # closenessIndex <- closeness(graph.obj, mode="all")
  # closenessIndex <- closenessIndex[sort(as.numeric(names(closenessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  # averagecloseness <- mean(closenessIndex)

  eigen_centralityIndexQ95 <- rep(0,length(membership.graph))
  betweennessIndexQ95 <- rep(0,length(membership.graph))
  closenessIndexQ95 <- rep(0,length(membership.graph))
  
  for(ji in unique(membership.graph)) {
    
    subgraph.obj <- subgraph(graph.obj, names(which(membership.graph == ji)))
    
    if( length(V(subgraph.obj)) < 3 ) {
      
      eigen_centralityIndex <- rep(0,length(V(subgraph.obj)))
      betweennessIndex <- rep(0,length(V(subgraph.obj)))
      closenessIndex <- rep(0,length(V(subgraph.obj)))
      names(eigen_centralityIndex) <- names(V(subgraph.obj))
      names(betweennessIndex) <- names(V(subgraph.obj))
      names(closenessIndex) <- names(V(subgraph.obj))
      
      eigen_centralityIndex <- eigen_centralityIndex[sort(as.numeric(names(eigen_centralityIndex)),decreasing = FALSE,index.return=TRUE)$ix]
      betweennessIndex <- betweennessIndex[sort(as.numeric(names(betweennessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
      closenessIndex <- closenessIndex[sort(as.numeric(names(closenessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
      
    } 
    
    if( length(V(subgraph.obj)) > 2 ) {
      
      betweennessIndex <- betweenness(subgraph.obj)
      eigen_centralityIndex <- eigen_centrality(subgraph.obj, directed = FALSE, scale = FALSE)$vector
      closenessIndex <- closeness(subgraph.obj, mode="all")
    
      eigen_centralityIndex <- eigen_centralityIndex[sort(as.numeric(names(eigen_centralityIndex)),decreasing = FALSE,index.return=TRUE)$ix]
      betweennessIndex <- betweennessIndex[sort(as.numeric(names(betweennessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
      closenessIndex <- closenessIndex[sort(as.numeric(names(closenessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
      
      eigen_centralityIndexQ95[as.numeric(names(which(eigen_centralityIndex >= quantile(eigen_centralityIndex,probs=0.95))))] <- 1
      betweennessIndexQ95[as.numeric(names(which(betweennessIndex >= quantile(betweennessIndex,probs=0.95))))] <- 1
      closenessIndexQ95[as.numeric(names(which(closenessIndex >= quantile(closenessIndex,probs=0.95))))] <- 1
      
    }
    
    eighenCentralityResults[as.numeric(names(eigen_centralityIndex)) ,c] <- as.numeric(eigen_centralityIndex)
    betweennessResults[as.numeric(names(betweennessIndex)) ,c] <- as.numeric(betweennessIndex)
    closenessResults[as.numeric(names(closenessIndex)) ,c] <- as.numeric(closenessIndex)
    
  }
  
  write.csv(eighenCentralityResults,file="../Results/eighenCentralityMPAs.csv")
  write.csv(betweennessResults,file="../Results/betweennessMPAs.csv")
  write.csv(closenessResults,file="../Results/closenessMPAs.csv")
  
  averageEighenCentrality <- mean(as.numeric(eighenCentralityResults[,c]))
  averageBetweenness <- mean(as.numeric(betweennessResults[,c]))
  averageCloseness <- mean(as.numeric(closenessResults[,c]))
  
  sdEighenCentrality <- sd(as.numeric(eighenCentralityResults[,c]))
  sdBetweenness <- sd(as.numeric(betweennessResults[,c]))
  sdCloseness <- sd(as.numeric(closenessResults[,c]))
  
  highereighenCentralityResults[,c] <- eigen_centralityIndexQ95
  higherBetweennessResults[,c] <- betweennessIndexQ95
  higherclosenessResults[,c] <- closenessIndexQ95
    
  write.csv(highereighenCentralityResults,file="../Results/higherEighenCentrality.csv")
  write.csv(higherBetweennessResults,file="../Results/higherBetweenness.csv")
  write.csv(higherclosenessResults,file="../Results/highercloseness.csv")
  
  # Plot centrality indexes with Connections and clusters
  
  eighenCentralityIndexPlot <- as.numeric(eighenCentralityResults[,c])
  betweennessIndexPlot <- as.numeric(betweennessResults[,c])
  closenessIndexPlot <- as.numeric(closenessResults[,c])
  
  eighenCentralityIndexPlot <- (eighenCentralityIndexPlot / max(eighenCentralityIndexPlot) )
  eighenCentralityIndexPlot <- (eighenCentralityIndexPlot * 2) + 2.5
  
  betweennessIndexPlot <- (betweennessIndexPlot / max(betweennessIndexPlot) )
  betweennessIndexPlot <- (betweennessIndexPlot * 2) + 2.5
  
  closenessIndexPlot <- (closenessIndexPlot / max(closenessIndexPlot) )
  closenessIndexPlot <- (closenessIndexPlot * 2) + 2.5
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapEighenCentralityConnections.pdf"), width=12)
  
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = eighenCentralityIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = eighenCentralityIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[eigen_centralityIndexQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[eigen_centralityIndexQ95 == 1], size = eighenCentralityIndexPlot[eigen_centralityIndexQ95 == 1], stroke = 1.2, alpha = 0.7)
    
  dev.off()
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapBetweennessConnections.pdf"), width=12)
  
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = betweennessIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = betweennessIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[betweennessIndexQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[betweennessIndexQ95 == 1], size = betweennessIndexPlot[betweennessIndexQ95 == 1], stroke = 1.2, alpha = 0.7)
  
  dev.off()
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapClosenessConnections.pdf"), width=12)
  
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = closenessIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = closenessIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[closenessIndexQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[closenessIndexQ95 == 1], size = closenessIndexPlot[closenessIndexQ95 == 1], stroke = 1.2, alpha = 0.7)
  
  dev.off()
    
  # ----------------------------------
  # eighenvector vs area
  
  par(mar=c(5,5,3,3),bg = 'white')  
  
  data <- data.frame(Type=ifelse(eigen_centralityIndexQ95 == 1 , "Hub" , "non-Hub"),data=area(notakeMPA)/1000000)
  data <- data[ data$data != max(data$data) ,]
  
  pdf(file=paste0("../Results/",project.name,"/areaEffectEigencentrality.pdf"),  width=10)
  
  boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
          names = c("Hub", "non-Hub"),
          las = 2,
          bg="#9A9A9A",
          col = c("#9A9A9A","#D4D4D4"),
          border = "#505050",
          horizontal = TRUE ,
          medcol = c("#9A9A9A","#D4D4D4"),
          xlab = "Marine reserve area (x 10^6 km2)" )
  
  points(mean(data[data$Type == "Hub",2]),1,pch=19)
  points(mean(data[data$Type != "Hub",2]),2,pch=19)
  text(max(data$data)-3,2.3,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(data$Type))$p.value, digits=4)))
  
  dev.off()
  
  data <- data.frame(Type=ifelse(betweennessIndexQ95 == 1 , "Hub" , "non-Hub"),data=area(notakeMPA)/1000000)
  data <- data[ data$data != max(data$data) ,]
  
  pdf(file=paste0("../Results/",project.name,"/areaEffectBetweenness.pdf"), width=10)
  
  boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
          names = c("Hub", "non-Hub"),
          las = 2,
          bg="#9A9A9A",
          col = c("#9A9A9A","#D4D4D4"),
          border = "#505050",
          horizontal = TRUE ,
          medcol = c("#9A9A9A","#D4D4D4"),
          xlab = "Marine reserve area (x 10^6 km2)" )
  
  points(mean(data[data$Type == "Hub",2]),1,pch=19)
  points(mean(data[data$Type != "Hub",2]),2,pch=19)
  text(max(data$data)-3,2.3,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(data$Type))$p.value, digits=4)))
  
  dev.off()
  
  data <- data.frame(Type=ifelse(closenessIndexQ95 == 1 , "Hub" , "non-Hub"),data=area(notakeMPA)/1000000)
  data <- data[ data$data != max(data$data) ,]
  
  pdf(file=paste0("../Results/",project.name,"/areaEffectCloseness.pdf"), width=10)
  
  boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
          names = c("Hub", "non-Hub"),
          las = 2,
          bg="#9A9A9A",
          col = c("#9A9A9A","#D4D4D4"),
          border = "#505050",
          horizontal = TRUE ,
          medcol = c("#9A9A9A","#D4D4D4"),
          xlab = "Marine reserve area (x 10^6 km2)" )
  
  points(mean(data[data$Type == "Hub",2]),1,pch=19)
  points(mean(data[data$Type != "Hub",2]),2,pch=19)
  text(max(data$data)-3,2.3,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(data$Type))$p.value, digits=4)))
  
  dev.off()
  
  # ----------------------------------
  
  combResults <- rbind(combResults,
                       data.frame(pld=pld.period,
                                  n.isolated.mpa=isolated.mpa.length,
                                  aggregationAtLeastOne=aggregationAtLeastOne,
                                  aggregationAllConnections=aggregationAllConnections,
                                  averageConnections=averageConnections,
                                  sdConnections=sdConnections,
                                  maximumConnections=maximumConnections,
                                  numberClusters=numberClusters,
                                  aggregationBasedOnClusters=aggregationBasedOnClusters,
                                  averageEighenCentrality=averageEighenCentrality,
                                  sdEighenCentrality=sdEighenCentrality,
                                  averageBetweenness=averageBetweenness,
                                  sdBetweenness=sdBetweenness,
                                  averageCloseness=averageCloseness,
                                  sdCloseness=sdCloseness) )
  
  write.csv(combResults,file="../Results/Results.csv")
  save(combResults,file=paste0("../Results/allPLDResults.Rdata"))
  dev.off()
  
}

## ---------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------

x <- combResults$pld
x.lab <- "Propagule duration (day)"

y <- combResults$averageeigen_centrality # colnames(combResults)
y.lab <- "Eigen centrality (average)"

pdf(file=paste0("../Results/Plots/averageEigen_centrality.pdf"), width=12, height=5)

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