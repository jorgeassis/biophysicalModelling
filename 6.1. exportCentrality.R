
rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)

setwd("/Volumes/StingRay/Dropbox/Manuscripts/Oceanographic connectivity explains the distribution of genetic diversity of marine forests/Git")

source("0. Config.R")
source("Dependences.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

pipeLiner <- TRUE
doParallelCalculations <- TRUE # repeat all parallel computations
type <- "points" # points polygons

global.simulation.parameters <- loadRData(paste0(results.folder,"/","modelParameters.RData"))
n.days.max <- global.simulation.parameters$particle.max.duration
n.days.max

## ----------------------------------------------------------------

pld.period <- c(0) # 1:n.days.max # 
n.seasons <- "" # c("","Spring","Summer","Autumn","Winter")
combinations <- expand.grid(season=n.seasons,pld.period=pld.period,stringsAsFactors = F)
season <- combinations[1,1]

## ------------

source.sink.xy <- loadRData(paste0(results.folder,"/","sourceSinkSites.RData"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

# Open connectivity

Connectivity.desc <- paste0(results.folder,"/","particlePairedConnectivityAveragedTable.desc")
load(file=paste0(results.folder,"/particlePairedConnectivityAveragedTableNames.RData"))
ConnectivityAll <- attach.big.matrix(Connectivity.desc)
ConnectivityAll <- as.data.table(ConnectivityAll[])
colnames(ConnectivityAll) <- particles.connectivity.names

number.cores <- 3

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

marine.distances <- foreach(c=1:nrow(combinations), .combine='rbind', .verbose=FALSE, .packages=c("igraph","data.table","stringr")) %dopar% {
  
  pld.period <- combinations[c,2]
  
  # Read from external Function if PD is zero
  if( pld.period == 0 ) { Connectivity <- loadRData(paste0("../Data/ConnectivityIntegrativePDFunction.RData")) }
  
  if( pld.period != 0 ) { Connectivity <- ConnectivityAll[ConnectivityAll$Max.Time <= pld.period,] }
  
  ## ------------------------------------------------------------------------------------------------------------
  
  project.name.c <- paste0(results.folder,"/connectivityExport/","Sim",season,str_pad(pld.period, 3, pad = "0"),"Days/")
  
  ## ----------------------------------------------------
  
  if( ! dir.exists(project.name.c) ) { dir.create(file.path(project.name.c), showWarnings = FALSE) } 
  if( ! dir.exists(paste0(project.name.c,"/Data")) ) { dir.create(file.path(paste0(project.name.c,"/Data")), showWarnings = FALSE) } 
  if( ! dir.exists(paste0(project.name.c,"/Maps")) ) { dir.create(file.path(paste0(project.name.c,"/Maps")), showWarnings = FALSE) } 
  if( ! dir.exists(paste0(project.name.c,"/Networks")) ) { dir.create(file.path(paste0(project.name.c,"/Networks")), showWarnings = FALSE) } 
  
  write.csv(source.sink.xy,paste0(project.name.c,"/sourceSinkXY.csv"), row.names = FALSE)
  
  comb <- as.data.frame( Connectivity )
  comb <- comb[ which(comb$Pair.from != comb$Pair.to) ,]
  comb <- comb[ sort(comb$Mean.Probability , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Mean.Probability")] 
  
  graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  #E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
  E(graph.obj)$weight = comb[,3] # Hock, Karlo Mumby, Peter J 2015
  graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
  graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
  graph.obj <- simplify(graph.obj)
  
  closenessIndex <- closeness(graph.obj)
  closenessIndex <- closenessIndex[sort(as.numeric(names(closenessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  closenessIndex <- (closenessIndex - min(closenessIndex)) / (max(closenessIndex) - min(closenessIndex)) 
  write.csv(data.frame(Pair=names(closenessIndex),Closeness=closenessIndex),paste0(project.name.c,"/closenessIndex.csv"), row.names = FALSE)
  
  harmonicIndex <- harmonic_centrality(graph.obj)
  harmonicIndex <- harmonicIndex[sort(as.numeric(names(harmonicIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  harmonicIndex <- (harmonicIndex - min(harmonicIndex)) / (max(harmonicIndex) - min(harmonicIndex)) 
  write.csv(data.frame(Pair=names(harmonicIndex),Closeness=harmonicIndex),paste0(project.name.c,"/harmonicIndex.csv"), row.names = FALSE)
  
  betweennessIndex <- betweenness(graph.obj)
  betweennessIndex <- betweennessIndex[sort(as.numeric(names(betweennessIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  betweennessIndex <- (betweennessIndex - min(betweennessIndex)) / (max(betweennessIndex) - min(betweennessIndex)) 
  write.csv(data.frame(Pair=names(betweennessIndex),Betweenness=betweennessIndex),paste0(project.name.c,"/betweennessIndex.csv"), row.names = FALSE)
  
  degreeIndex <- degree(graph.obj)
  degreeIndex <- degreeIndex[sort(as.numeric(names(degreeIndex)),decreasing = FALSE,index.return=TRUE)$ix]
  write.csv(data.frame(Pair=names(degreeIndex),Degree=degreeIndex),paste0(project.name.c,"/degreeIndex.csv"), row.names = FALSE)
  
  return( NULL )
  
}

stopCluster(cl.2) ; rm(cl.2)
