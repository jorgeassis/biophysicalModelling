## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("../Project Config 0.R")
source("Dependences.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

if( ! "connectivityExport" %in% file.path(paste0(project.folder,"/Results/",project.name)) ) { dir.create(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport"))) }
connectivityExportDir <- file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/"))

## -------------------------

additional.source.sink <- shapefile(additional.source.sink.shp)
additional.source.sink <- additional.source.sink[,"ID"]
writeOGR(obj=additional.source.sink, dsn=paste0(project.folder,"/Results/",project.name,"/"), layer="sourceSinkPolyg", driver="ESRI Shapefile") 

## -------------------------

distance.probability <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveraged.bm") )
distance.probability <- as.data.frame(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events")
distance.probability

source.sink.xy <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm") )
source.sink.xy <- as.data.frame(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy

## -------------------------

## Sort connectivity estimates by polygon order

sorted.ids <- data.frame(polygons.id=additional.source.sink$ID,pair=NA)

for( i in 1:nrow(additional.source.sink)) {
  
  polygon.i <- additional.source.sink[i,]
  polygon.i.centroid <- coordinates(polygon.i)
  closest.match <- which.min(spDistsN1(as.matrix(source.sink.xy[,c("Lon","Lat")]),polygon.i.centroid))
  sorted.ids[i,2] <- source.sink.xy[closest.match,"Pair"]
  
}

head(sorted.ids)

## -------------------------
## Force all vectors

missing.Pair.from <- sorted.ids$pair[which( ! sorted.ids$pair %in% distance.probability$Pair.from)]
if( length(missing.Pair.from) > 0 ) {
  for( i in 1:length(missing.Pair.from)) {
    distance.probability <- rbind(distance.probability,distance.probability[1,])
    distance.probability[nrow(distance.probability),"Pair.from"] <- missing.Pair.from[i]
    distance.probability[nrow(distance.probability),3:11] <- 0
  }
}

missing.Pair.to <- sorted.ids$pair[which( ! sorted.ids$pair %in% distance.probability$Pair.to)]
if( length(missing.Pair.to) > 0 ) {
  for( i in 1:length(missing.Pair.to)) {
    distance.probability <- rbind(distance.probability,distance.probability[1,])
    distance.probability[nrow(distance.probability),"Pair.to"] <- missing.Pair.to[i]
    distance.probability[nrow(distance.probability),3:11] <- 0
  }
}

## --------------------------------------------------------------------
## --------------------------------------------------------------------

## Assign IDs [e.g., region of interest] to sink.source sites

source.sink.xy.ids <- assignIDs(source.sink.xy[,c("Lon","Lat")])
sum(is.na(source.sink.xy.ids[,1]))
save(source.sink.xy.ids, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.ids.RData"))

source.sink.xy <- cbind(source.sink.xy,source.sink.xy.ids)

unique(source.sink.xy.ids$Name)
coords.centroid.ids <- data.frame(Name=unique(source.sink.xy.ids$Name),
                                  Lon=c(-31.203339, -31.105835, -28.701468, -28.329306, -28.058767, -28.012076, -25.477667, -25.102516, -27.222877),
                                  Lat=c(39.440657, 39.696841, 38.573669, 38.467299, 38.654148, 39.054118, 37.772072, 36.975784, 38.723658))

save(coords.centroid.ids, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","coords.centroid.ids.RData"))

## --------------------------------------------------------------------
## --------------------------------------------------------------------

## Attribute pairwise connectivity estimates for each polygon 

n.days <- 30

new.extent <- c(min(source.sink.xy[,2]),max(source.sink.xy[,2]),min(source.sink.xy[,3]),max(source.sink.xy[,3]))
network <- produce.network("Prob",distance.probability,n.days,FALSE,10,source.sink.xy,new.extent)

## -------------------------
## Direct Connectivity

library(reshape2)
matrixConnectivityP.square <- acast(distance.probability[,c("Pair.from","Pair.to","Probability")], formula=Pair.from~Pair.to, fun.aggregate= mean , value.var="Probability",drop = FALSE)
matrixConnectivityT.square <- acast(distance.probability[,c("Pair.from","Pair.to","Mean.Time")], formula=Pair.from~Pair.to,  fun.aggregate=mean, value.var="Mean.Time",drop = FALSE)

dim(matrixConnectivityP)
dim(matrixConnectivityT)

matrixConnectivityP <- melt(matrixConnectivityP.square)
matrixConnectivityT <- melt(matrixConnectivityT.square)

## -------------------------
## Stepping-stone Connectivity

cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)

matrixConnectivitySSP <- foreach(from=sorted.ids$pair, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  network.x <- network[[2]]
  connectivity.x <- network[[1]]
  res.connectivity.to <- numeric(0)
  
  for( to in sorted.ids$pair ) {
    
    possible.paths.y <- get.shortest.paths(network.x,as.character( from ) , as.character( to ),mode="out")$vpath
    stones.t <- as.numeric(names(possible.paths.y[[1]]))
    stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
    path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
    
    if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
    if( length(path.values) == 0) { path.values <- 0 }
    if( from == to ) { path.values <- connectivity.x[ connectivity.x[,1] == from & connectivity.x[,2] == to , 3 ][1]  }
    
    res.connectivity.to <- c(res.connectivity.to,path.values)
    
  }
  
  temp.res <- data.frame( pair.from = from,
                          pair.to = sorted.ids$pair , 
                          connectivity = res.connectivity.to)
  return(temp.res)
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

matrixConnectivitySSP <- do.call(rbind,matrixConnectivitySSP)
matrixConnectivitySSP.square <- acast(matrixConnectivitySSP, formula=pair.from~pair.to, fun.aggregate= mean , value.var="connectivity",drop = FALSE)

## -----------------------------------------
## Export Connectivity

# Square Matrix

dim(matrixConnectivityP.square)
dim(matrixConnectivityT.square)
dim(matrixConnectivitySSP.square)

matrixConnectivityP.square[is.na(matrixConnectivityP.square)] <- 0
matrixConnectivityT.square[is.na(matrixConnectivityT.square)] <- 0
matrixConnectivitySSP.square[is.na(matrixConnectivitySSP.square)] <- 0

write.table(file=paste0(connectivityExportDir,"matrixConnectivityP.square.csv"),x=matrixConnectivityP.square,row.names=FALSE,col.names=FALSE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityT.square.csv"),x=matrixConnectivityT.square,row.names=FALSE,col.names=FALSE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivitySSP.square.csv"),x=matrixConnectivitySSP.square,row.names=FALSE,col.names=FALSE,sep=";",dec=".")

# Pairs

matrixConnectivityP[is.na(matrixConnectivityP)] <- 0
matrixConnectivityT[is.na(matrixConnectivityT)] <- 0
matrixConnectivitySSP[is.na(matrixConnectivitySSP)] <- 0

colnames(matrixConnectivityP) <- c("From","To","Probability")
colnames(matrixConnectivityT) <- c("From","To","Time")
colnames(matrixConnectivitySSP) <- c("From","To","Probability")

write.table(file=paste0(connectivityExportDir,"matrixConnectivityP.csv"),x=matrixConnectivityP,row.names=FALSE,col.names=FALSE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityT.csv"),x=matrixConnectivityT,row.names=FALSE,col.names=FALSE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivitySSP.csv"),x=matrixConnectivitySSP,row.names=FALSE,col.names=FALSE,sep=";",dec=".")

## -----------------------------------------
## Plot connectivity [networks]

graph.i <- network[[2]]

## -----------------------------------------
## Aggregate and plot network assignments

matrixConnectivityP.agg <- data.frame(expand.grid(from=unique(source.sink.xy$Name),to=unique(source.sink.xy$Name)),val=NA)
matrixConnectivityPSS.agg <- data.frame(expand.grid(from=unique(source.sink.xy$Name),to=unique(source.sink.xy$Name)),val=NA)
matrixConnectivityT.agg <- data.frame(expand.grid(from=unique(source.sink.xy$Name),to=unique(source.sink.xy$Name)),val=NA)

for ( i in 1:nrow(matrixConnectivityP.agg)) {
  
  from <- matrixConnectivityP.agg[i,1]
  from <- source.sink.xy[source.sink.xy$Name == from,"Pair"]
  to <- matrixConnectivityP.agg[i,2]
  to <- source.sink.xy[source.sink.xy$Name == to,"Pair"]
  
  matrixConnectivityP.agg[i,3] <- mean(matrixConnectivityP[ matrixConnectivityP$From %in% from & matrixConnectivityP$To %in% to , "Probability" ])
  matrixConnectivityPSS.agg[i,3] <- mean(matrixConnectivitySSP[ matrixConnectivitySSP$From %in% from & matrixConnectivitySSP$To %in% to , "Probability" ])
  
  timeVal <- matrixConnectivityT[ matrixConnectivityT$From %in% from & matrixConnectivityT$To %in% to , "Time" ]
  timeVal <- mean(timeVal[timeVal != 0])
  matrixConnectivityT.agg[i,3] <- timeVal

}

matrixConnectivityP.agg <- matrixConnectivityP.agg[ sort( matrixConnectivityP.agg$val , decreasing = TRUE, index.return =TRUE)$ix , ]
matrixConnectivityPSS.agg <- matrixConnectivityPSS.agg[ sort( matrixConnectivityPSS.agg$val , decreasing = TRUE, index.return =TRUE)$ix , ]
matrixConnectivityT.agg <- matrixConnectivityT.agg[ sort( matrixConnectivityT.agg$val , decreasing = TRUE, index.return =TRUE)$ix , ]

# ----------

comb <- matrixConnectivityP.agg
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight <- comb[,3]
E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,1/-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015

graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "min") # min / mean / max
graph.obj <- simplify(graph.obj)

clustering.graph <- cluster_leading_eigen(graph.obj,options=list(maxiter=1000000))
membership.graph <- clustering.graph$membership

cols.to.use <- distinctColors(length(unique(membership.graph)))
cols.to.use <- cols.to.use[membership.graph]
V(graph.obj)$color <- cols.to.use
l <- layout.fruchterman.reingold(graph.obj)

lineThinkness <- E(graph.obj)$weight
lineThinkness <- lineThinkness / max(lineThinkness)

lineThinkness[ lineThinkness >= 0.5 ] <- 1
lineThinkness[ lineThinkness < 0.5 & lineThinkness >= 0.1 ] <- 0.5
lineThinkness[ lineThinkness < 0.1 & lineThinkness >= 0.01 ] <- 0.25
lineThinkness[ lineThinkness < 0.01 & lineThinkness >= 0.001 ] <- 0.1
lineThinkness[ lineThinkness < 0.001 & lineThinkness >= 0.0001 ] <- 0.05
lineThinkness[ lineThinkness < 0.0001 ] <- 0.01

pdf( file=paste0(connectivityExportDir,"idsNetworkClustering.pdf") , width = 10 )
plot(graph.obj,edge.width=lineThinkness,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

## -----------------------------------------
## Map network assignments

coords.centroid.ids

get.vertex.attribute(graph.obj)

## -----------------------------------------
## Map connectivity [clusters]





## -------------------------
## -------------------------

Get no data polygons and get closest one - polygons without hexagon assigned!
  
  
  
rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)

library(rnaturalearth)
library(geosphere)
library(rgeos)

source("../Project Config 0.R")

sql.project.name <- "MPA"
number.cores <- 8

notakeMPA <- "../Data/notake_merged_Med_P.shp" 

sql.file <- "../Results/SQL/CentralityMPASimulationResults.sql"
bigmatrix.file <- "../InternalProc/particles.reference.desc"
sorce.sink.cells.file <- "../Results/source.sink.bm"
coordRef <- crs(shapefile(notakeMPA))

## ------------------------------------------------------------------------------------------------------------
## Read main sources


## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

# List results

list.dirs(path = paste0("../Results"), recursive = FALSE)

season <- "YearRound" # c("YearRound","SeasonSummer","SeasonWinter")
pld.period <- 1:120 # c(10 , 30 , 90 , 120 , 200)
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

dev.off()
doParallelCalculations <- TRUE

for( c in 77:nrow(combinations) ){ #  
  
  cat(c,"\n")
  gc(reset=TRUE)
  
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
  
  if(doParallelCalculations) {
      
        particles.reference.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.reference.desc"))
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
        save(connectivity.source.sink.xy,file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.xy.Rdata"))
        
        ## ----------------------------------------------------
        
        mpaIDPairs <- expand.grid(from=notakeMPA$ID,to=notakeMPA$ID)
        polygonsCompute <- mpaIDPairs
        
        cl.2 <- makeCluster(number.cores , type="FORK")
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
        
  }
  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  
  # load(file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.notakeMPA.Rdata"))
  # load(file=paste0("../Results/",project.name,"/Data/connectivity.source.sink.xy.Rdata"))
  
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
  resistanceResults[ ,c] <- resistance
  write.csv(resistanceResults,file="../Results/resistanceMPA.csv")
  
  # Those with higher / Percertil 95%
  higherResistance <- which(resistance >=  as.numeric(quantile(resistance[resistance !=1 ],0.95,na.rm=TRUE))) 
  temporaryRes <- MPAnames[ higherResistance ]
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
    geom_point(data = centroids[as.numeric(isolated.mpa.id),] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = "#9C2323", size = 2.5, stroke = 0.35, alpha = 0.9)
  )
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
  print(
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = eighenCentralityIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = eighenCentralityIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[eigen_centralityIndexQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[eigen_centralityIndexQ95 == 1], size = eighenCentralityIndexPlot[eigen_centralityIndexQ95 == 1], stroke = 1.2, alpha = 0.7)
  )
  dev.off()
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapBetweennessConnections.pdf"), width=12)
  print(
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = betweennessIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = betweennessIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[betweennessIndexQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[betweennessIndexQ95 == 1], size = betweennessIndexPlot[betweennessIndexQ95 == 1], stroke = 1.2, alpha = 0.7)
  )
  dev.off()
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapClosenessConnections.pdf"), width=12)
  print(
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = closenessIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = closenessIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[closenessIndexQ95 == 1,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[closenessIndexQ95 == 1], size = closenessIndexPlot[closenessIndexQ95 == 1], stroke = 1.2, alpha = 0.7)
  )
  dev.off()
    
  # ----------------------------------
  # Plot resistance
  
  resistanceIndexPlot <- (resistance / max(resistance) )
  resistanceIndexPlot <- (resistanceIndexPlot * 2) + 2.5
  
  pdf(file=paste0("../Results/",project.name,"/Maps/MapResistanceConnections.pdf"), width=12)
  print(
  mapSouthernEuropeNet + 
    geom_point(data = centroids[cols.to.use == "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use == "white"], size = resistanceIndexPlot[cols.to.use == "white"], stroke = 0.25, alpha = 0.9) +
    geom_point(data = centroids[cols.to.use != "white",] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[cols.to.use != "white"], size = resistanceIndexPlot[cols.to.use != "white"], stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroids[higherResistance,] ,  aes(x = x, y = y) , shape = 21, colour = "black", fill = cols.to.use[higherResistance], size = resistanceIndexPlot[higherResistance], stroke = 1.2, alpha = 0.7)
  )
  dev.off()
  
  # ----------------------------------
  # centrality vs area vs shape (corrected perimeter-area ratio; https://doi.org/10.1016/j.gecco.2018.e00504)
  
  shapeIndex <- perimeter(notakeMPA) / sqrt( 4*pi*area(notakeMPA))
  
  par(mar=c(5,5,3,3),bg = 'white')  
  
  data <- data.frame(Type=ifelse(eigen_centralityIndexQ95 == 1 , "Hub" , "non-Hub"),data=shapeIndex)
  data <- data[ data$data != max(data$data) ,]
  
  pdf(file=paste0("../Results/",project.name,"/shapeEffectEigencentrality.pdf"),  width=10)
  
  boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
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
  text(max(data$data)-0.3,2.3,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(data$Type))$p.value, digits=4)))
  shapeEffectEigencentrality <- kruskal.test(data[,2],as.numeric(data$Type))$p.value < 0.05
  
  dev.off()
  
  data <- data.frame(Type=ifelse(betweennessIndexQ95 == 1 , "Hub" , "non-Hub"),data=shapeIndex)
  data <- data[ data$data != max(data$data) ,]
  
  pdf(file=paste0("../Results/",project.name,"/shapeEffectBetweenness.pdf"), width=10)
  
  boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
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
  text(max(data$data)-0.3,2.3,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(data$Type))$p.value, digits=4)))
  shapeEffectBetweenness <- kruskal.test(data[,2],as.numeric(data$Type))$p.value < 0.05
  
  dev.off()
  
  data <- data.frame(Type=ifelse(closenessIndexQ95 == 1 , "Hub" , "non-Hub"),data=shapeIndex)
  data <- data[ data$data != max(data$data) ,]
  
  pdf(file=paste0("../Results/",project.name,"/shapeEffectCloseness.pdf"), width=10)
  
  boxplot(data[data$Type == "Hub",2], data[data$Type == "non-Hub",2],
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
  text(max(data$data)-0.3,2.3,paste0("p-value: ",round(kruskal.test(data[,2],as.numeric(data$Type))$p.value, digits=4)))
  shapeEffectCloseness <- kruskal.test(data[,2],as.numeric(data$Type))$p.value < 0.05
  
  dev.off()
  
  # --------------
  
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
  areaEffectEigencentrality <- kruskal.test(data[,2],as.numeric(data$Type))$p.value < 0.05
  
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
  areaEffectBetweenness <- kruskal.test(data[,2],as.numeric(data$Type))$p.value < 0.05
  
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
  areaEffectCloseness <- kruskal.test(data[,2],as.numeric(data$Type))$p.value < 0.05
  
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
                                  averageResistance=mean(resistance,na.rm=T),
                                  averageEighenCentrality=averageEighenCentrality,
                                  sdEighenCentrality=sdEighenCentrality,
                                  averageBetweenness=averageBetweenness,
                                  sdBetweenness=sdBetweenness,
                                  averageCloseness=averageCloseness,
                                  sdCloseness=sdCloseness,
                                  shapeEffectEigencentrality=shapeEffectEigencentrality,
                                  shapeEffectBetweenness=shapeEffectBetweenness,
                                  shapeEffectCloseness=shapeEffectCloseness,
                                  areaEffectEigencentrality=areaEffectEigencentrality,
                                  areaEffectBetweenness=areaEffectBetweenness,
                                  areaEffectCloseness=areaEffectCloseness ))
  
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
y <- combResults$areaEffectBetweenness
y.lab <- "areaEffectEigencentrality"

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(x,y,pch=20,col="#A6A6A6", ylab="",xlab=x.lab,axes=FALSE)
title(ylab=y.lab, line=4)
lines(bezierCurve(x,y,100)$x,bezierCurve(x,y,100)$y,type="l", lwd=1, lty=2)
axis(2,las=2,col="White",col.ticks="Black", cex.axis=0.9)
axis(1,las=0,col="White",col.ticks="Black", cex.axis=0.9)
box()

pdf(file=paste0("../Results/aggregation degree.pdf"), width=12)

ggplot(combResults, aes(x = pld, y = averageConnections)) +
  geom_line() +
  geom_ribbon(aes(ymin = averageConnections - sdConnections/2,
                  ymax = averageConnections + sdConnections/2), alpha = 0.2) + xlab("Propagule duration (day)") + ylab("Network connections (average ± standard deviation)") + xlim(0,120)

dev.off()

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------