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
library(rgeos)

source("0. Project Config.R")

regionsOfInterest <- "../Data/mainRegions.shp" 

bigmatrix.file <- paste0(project.folder,"/Results/",project.name,"/InternalProc/","particles.reference.desc")
sorce.sink.cells.file <- paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm")
coordRef <- crs(shapefile(regionsOfInterest))

## ------------------------------------------------------------------------------------------------------------
## Read main sources

load(paste0(project.folder,"/Results/",project.name,"/InternalProc/","Parameters.RData"))

sim.extent <-unique(as.numeric(unlist(strsplit(global.simulation.parameters$extent, split=","))))
months <- unique(as.numeric(unlist(strsplit(global.simulation.parameters$sim.months , split=","))))
n.hours.per.day <- global.simulation.parameters$n.hours.per.day
n.particles.per.cell <- global.simulation.parameters$n.particles.per.cell
n.new.particles.per.day <- global.simulation.parameters$n.new.particles.per.day
n.steps.per.day <- global.simulation.parameters$n.hours.per.day

source.sink.xy <- read.big.matrix(sorce.sink.cells.file)
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

source.sink.xy.sp <- source.sink.xy[,2:3]
coordinates(source.sink.xy.sp) <- ~Lon+Lat
crs(source.sink.xy.sp) <- coordRef

## ---------------------

regionsOfInterest <- shapefile(regionsOfInterest)
regionsOfInterest <- gBuffer(regionsOfInterest, byid=TRUE, width=0)
worldMap <- ne_countries(scale = 10, returnclass = "sp")

regionsOfInterest[1,"Region"] <- "Acores"
regionsOfInterest[8,"Region"] <- "Mediterranean"

# Aggregate based on names

regionsOfInterest$Region[which(grepl("Spain",regionsOfInterest$Region))] <- regionsOfInterest$Region[which(grepl("Spain",regionsOfInterest$Region))[1]]

regionsOfInterest <- aggregate(regionsOfInterest, by = list(RegionFinal = regionsOfInterest$Region), FUN=mean, dissolve = TRUE, areaWeighted = FALSE)
regionsOfInterest$Area <- round(raster::area(regionsOfInterest)/ 1000000,digits=3)

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
## Identify MPA id in source sink sites

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
## Revert regions to landmasses

RegionNames <- regionsOfInterest$RegionFinal
RegionID <- regionsOfInterest$ID
regionsOfInterest <- gIntersection(regionsOfInterest, gUnaryUnion(worldMap), byid = TRUE, id = sapply(regionsOfInterest@polygons, slot, name = "ID"))
regionsOfInterest$Names <- RegionNames
regionsOfInterest$ID <- RegionID
regionsOfInterest
plot(regionsOfInterest)

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

# List results

list.dirs(path = paste0("../Results"), recursive = FALSE)

season <- "YearRound" # c("YearRound","SeasonSummer","SeasonWinter")
pld.period <- 30 # 1:120 c(10 , 30 , 90 , 120 , 200)
combinations <- expand.grid(season=season,pld.period=pld.period,stringsAsFactors = F)

## --------------------

combResults <- data.frame()

isolatedResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(isolatedResults) <- RegionNames
colnames(isolatedResults) <- 1:nrow(combinations)

betweennessResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(betweennessResults) <- RegionNames
colnames(betweennessResults) <- 1:nrow(combinations)

higherBetweennessResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherBetweennessResults) <- RegionNames
colnames(higherBetweennessResults) <- 1:nrow(combinations)

eighenCentralityResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(eighenCentralityResults) <- RegionNames
colnames(eighenCentralityResults) <- 1:nrow(combinations)

highereighenCentralityResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(highereighenCentralityResults) <- RegionNames
colnames(highereighenCentralityResults) <- 1:nrow(combinations)

closenessResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(closenessResults) <- RegionNames
colnames(closenessResults) <- 1:nrow(combinations)

higherclosenessResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherclosenessResults) <- RegionNames
colnames(higherclosenessResults) <- 1:nrow(combinations)

clusterAssignment <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(clusterAssignment) <- RegionNames
colnames(clusterAssignment) <- 1:nrow(combinations)

resistanceResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(resistanceResults) <- RegionNames
colnames(resistanceResults) <- 1:nrow(combinations)

higherResistanceResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherResistanceResults) <- RegionNames
colnames(higherResistanceResults) <- 1:nrow(combinations)

outDegreeResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(outDegreeResults) <- RegionNames
colnames(outDegreeResults) <- 1:nrow(combinations)

higherOutDegreeResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherOutDegreeResults) <- RegionNames
colnames(higherOutDegreeResults) <- 1:nrow(combinations)

selfRecruitmentResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(selfRecruitmentResults) <- RegionNames
colnames(selfRecruitmentResults) <- 1:nrow(combinations)

higherSelfRecruitmentResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
rownames(higherSelfRecruitmentResults) <- RegionNames
colnames(higherSelfRecruitmentResults) <- 1:nrow(combinations)

## ------------------------------------------------------------------------------

dev.off()
doParallelCalculations <- TRUE

for( c in 1:nrow(combinations) ){ #  
  
  cat(c,"\n")
  gc(reset=TRUE)
  
  season <- combinations[c,1]
  pld.period <- combinations[c,2]
  
  if( season == "SeasonSummer" ) { spawn.p <- c(6,7,8,9) }
  if( season == "SeasonWinter" ) { spawn.p <- c(11,12,1,2) }
  if( season == "YearRound" ) { spawn.p <- 1:12 }
  
  project.name.c <- paste0(project.name,"/",season,"_Pld",pld.period)
  
  ## ----------------------------------------------------
  
  if( ! dir.exists(paste0("../Results/",project.name.c)) ) { dir.create(file.path(paste0("../Results/",project.name.c)), showWarnings = FALSE) } 
  if( ! dir.exists(paste0("../Results/",project.name.c,"/Data")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Data")), showWarnings = FALSE) } 
  if( ! dir.exists(paste0("../Results/",project.name.c,"/Maps")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Maps")), showWarnings = FALSE) } 
  if( ! dir.exists(paste0("../Results/",project.name.c,"/Networks")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Networks")), showWarnings = FALSE) } 
  
  if(doParallelCalculations) {
    
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
    
    save(connectivity,file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.notakeMPA.Rdata"))
    
  }
  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  
  load(file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.notakeMPA.Rdata"))
  load(file=paste0("../Results/",project.name.c,"/Data/connectivity.source.sink.xy.Rdata"))
  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  
  connectivity.matrix <- connectivity[ ,c(1,2,7)]
  connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )
  
  ## ----------------------------------
  
  colnames(connectivity.matrix) <- RegionNames
  rownames(connectivity.matrix) <- RegionNames
  
  ## ----------------------------------
  
  retention <- diag(connectivity.matrix)
  sumRows <- apply(connectivity.matrix,1,sum)
  SR <- retention / sumRows
  SR[is.na(SR)] <- 0
  
  selfRecruitmentResults[ ,c] <- SR
  write.csv(selfRecruitmentResults,file="../Results/selfRecruitment.csv")
  
  higherSelfRecruitment <- which(SR >=  as.numeric(quantile(SR,0.95,na.rm=TRUE))) 
  higherSelfRecruitmentResults[ higherSelfRecruitment ,c] <- "TRUE"
  write.csv(higherSelfRecruitmentResults,file="../Results/higherSelfRecruitment.csv.csv")
  
  diag(connectivity.matrix) <- 0
  
  write.csv(connectivity.matrix,paste0("../Results/",project.name.c,"/connectivitymatrix.csv"))
  
  isolated.mpa <- which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)
  isolatedResults[ isolated.mpa,c] <- "TRUE"
  write.csv(isolatedResults,file="../Results/isolatedMPAs.csv")
  
  # Aggregation level (Proportion of non-isolated MPAs, at least one connection, in relation to the number of MPAs)
  aggregationAtLeastOne <- (nrow(connectivity.matrix)-length(isolated.mpa)) / nrow(connectivity.matrix)
  
  # Aggregation level (Based on overall connections)
  connectivity.matrix.binomial <- connectivity.matrix
  connectivity.matrix.binomial[connectivity.matrix.binomial != 0] <- 1
  aggregationAllConnections <- sum ( ( apply(connectivity.matrix.binomial,1,sum) + 1 ) / nrow(connectivity.matrix.binomial) ) / nrow(connectivity.matrix.binomial)
  
  # Resistance index
  resistance <- sapply( 1:nrow(connectivity.matrix.binomial) , function(res) { 1- length(unique(c(which(connectivity.matrix.binomial[res,] == 1 ) ,  which(connectivity.matrix.binomial[,res] == 1 )))) / nrow(connectivity.matrix.binomial) } )
  resistanceResults[ ,c] <- resistance
  write.csv(resistanceResults,file="../Results/resistanceMPA.csv")
  
  # Those with higher / Percertil 95%
  higherResistance <- which(resistance >=  as.numeric(quantile(resistance[resistance !=1 ],0.95,na.rm=TRUE))) 
  higherResistanceResults[ higherResistance ,c] <- "TRUE"
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
  
  # extent(worldMap) + c(20,-20,20,-20)
  
  mapSouthernEurope <- ggplot() +
    geom_polygon(data = worldMap , fill = "#C4C4C4", colour = "#ffffff" , size=0.25 ,  aes(long, lat, group = group)) +
    coord_map('lambert', lat0=29.75, lat1=45, xlim=c(-30, 10), ylim=c(10, 52)) + theme_map
  
  mapSouthernEuropeNet <- mapSouthernEurope
  connected.pairs <- comb[comb$Probability > 0,]
  centroids <- as.data.frame(gCentroid(regionsOfInterest,byid=TRUE),xy=T)
  colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))
  for( i in nrow(connected.pairs):1 ){
    strenght <- (connected.pairs[i,3] * 100) + 1 
    routes_sl <- gcIntermediate(centroids[ which(regionsOfInterest$ID == connected.pairs[i,1]),],
                                centroids[ which(regionsOfInterest$ID == connected.pairs[i,2]),],
                                n = 100, addStartEnd = TRUE, sp = TRUE)
    SLDF = sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = NA), match.ID = F)
    mapSouthernEuropeNet <- mapSouthernEuropeNet + geom_path(data = SLDF, size=0.35 , aes(x = long, y = lat), col = "#797979") # colfunc(101)[round(strenght)]
  }
  
  ## --------------------------------------------------------------------------------
  
  graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  #E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
  E(graph.obj)$weight = comb[,3] # Hock, Karlo Mumby, Peter J 2015
  graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
  graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
  graph.obj <- simplify(graph.obj)
  
  # -------------
  
  graph.obj.ss <- graph.obj
  E(graph.obj.ss)$weight <- 1 - E(graph.obj.ss)$weight
  connectivity.matrix.ss <- connectivity.matrix
  for(i in 1:nrow(connectivity.matrix.ss)) {
    for(j in 1:nrow(connectivity.matrix.ss)) {
      connectivity.matrix.ss[i,j] <- distances(graph.obj.ss, as.character(i), as.character(j))
    }
  }
  connectivity.matrix.ss[connectivity.matrix.ss == Inf] <- 1
  unLinked <- which(apply(connectivity.matrix.ss,1,sum) == nrow(connectivity.matrix.ss) - 1 & apply(connectivity.matrix.ss,2,sum) == nrow(connectivity.matrix.ss) - 1)
  connectivity.matrix.ss <- connectivity.matrix.ss[-unLinked,-unLinked]
  
  # Euclidean distance
  
  dist <- as.dist(connectivity.matrix.ss)
  dist[7] <- 1.1
  # Hierarchical Clustering with hclust
  hc <- hclust(dist, method = "ward.D")
  
  # Plot the result
  plot(hc)
  
  # -------------
  
  membership.graph <- clusters(graph.obj)$membership
  membership.graph <- membership.graph[sort(as.numeric(names(membership.graph)),decreasing = FALSE,index.return=TRUE)$ix]
  
  clusterAssignment[ ,c] <- membership.graph
  write.csv(clusterAssignment,file="../Results/clusterAssignmentMPA.csv")
  
  # Number of clusters
  numberClusters <- length(unique(membership.graph)) - length(isolated.mpa)
  
  aggregationBasedOnClusters <- 1 - ( numberClusters / length(membership.graph) )
  
  ## ------------------------------------------
  
  # Plot clusters
  
  cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
  cols.to.use[isolated.mpa] <- "white"
  
  l <- layout.fruchterman.reingold(graph.obj)
  reducedNames <- sapply(names(V(graph.obj)),function(x) { regionsOfInterest$Names[regionsOfInterest$ID == x][1] })
  
  pdf(file=paste0("../Results/",project.name,"/Networks/MapClusteringConnections.pdf"), width=12)
  plot(graph.obj,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
  dev.off()
  
  ## ------------------------------------------
  
  # Plot clusters with connections
  
  cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
  cols.to.use[isolated.mpa] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapClusteringConnections.pdf"), width=12)
  print(
    mapSouthernEuropeNet + 
      geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = 2.5, stroke = 0.35, alpha = 0.9) +
      geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = 2.5, stroke = 0.35, alpha = 0.7)
  )
  dev.off()
  
  ## ------------------------------------------
  
  # Plot isolated with connections
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapIsolatedConnMPA.pdf"), width=12)
  print(
    mapSouthernEuropeNet + 
      geom_point(data = centroids ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.7) +
      geom_point(data = centroids[isolated.mpa,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "#9C2323", size = 2.5, stroke = 0.35, alpha = 0.9)
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
  cols.to.use[isolated.mpa] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapBetweennessConnections.pdf"), width=12)
  print(
    mapSouthernEuropeNet + 
      geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = betweennessIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
      geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = betweennessIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
      geom_point(data = centroids[betweennessQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "#9C2323", size = betweennessIndexPlot[betweennessQ95 == 1], stroke = 1.2, alpha = 0.7)
    + labs(caption = paste0("95th: ", round(quantile(betweennessIndex,probs=0.95),2)))
  )
  dev.off()
  
  cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
  cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(outDegreeIndexPlot),breaks = 20))]
  cols.to.use[isolated.mpa] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapStrengthOut.pdf"), width=12)
  print(
    mapSouthernEuropeNet + 
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
  cols.to.use[isolated.mpa] <- "white"
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapSRConnections.pdf"), width=12)
  print(
    mapSouthernEuropeNet + 
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
                                  n.isolated.mpa=length(isolated.mpa),
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
  
}

## ---------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------

library(ggplot2)
names(combResults)

x <- combResults$pld
x.lab <- "Propagule duration (day)"
y <- combResults$n.isolated.mpa
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
ggplot(combResults, aes(x = pld, y = n.isolated.mpa)) +
  geom_line() + 
  ylab("Isolation degree (number of reserves)") +
  xlab("Propagule duration (day)")
dev.off()

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------