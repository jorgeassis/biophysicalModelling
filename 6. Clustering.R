## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

source("0. Project Config.R")

## ---------------

produce.network <- function(network.type,connectivity,extract.simulation.days,crop.network,buffer,cells) {
  
  Connectivity <- as.data.table(connectivity)
  
  if(crop.network) {  final.cells <- which( cells[,2] >= (min( cells[,2] ) - buffer) & 
                                            cells[,2] <= (max( cells[,2] ) + buffer) & 
                                            cells[,3] >= (min( cells[,3] ) - buffer) & 
                                            cells[,3] <= (max( cells[,3] ) + buffer) )   
  
  final.cells <- cells[final.cells,1]
  
  }
  
  if( ! crop.network ) {  final.cells <- cells[,1]  }
  
  if( network.type == "Prob" ) {
    
    comb <- Connectivity[Time.max <= extract.simulation.days,.(Pair.from,Pair.to,Probability)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )
     
    # norm <- t(combn(position.matrix, 2))
    # 
    # for( i in 1:nrow(norm)) {
    #   
    #   t.1 <- which( comb[,1] == norm[i,1] & comb[,2] == norm[i,2] )
    #   
    #   if( length(t.1) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,1] , Pair.to = norm[i,2] ,  Mean.Probability = 0)) }
    #   
    #   t.2 <- which( comb[,1] == norm[i,2] & comb[,2] == norm[i,1] )
    #   
    #   if( length(t.2) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,2] , Pair.to = norm[i,1] ,  Mean.Probability = 0)) }
    #   
    # }
    
    net.function <<- prod
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    # E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
    E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  if( network.type == "Time" ) {
    
    comb <- Connectivity[Time.max <= extract.simulation.days,.(Pair.from,Pair.to,Time.mean)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Time.mean] , decreasing = TRUE, index.return =TRUE)$ix , ] )
    
    # norm <- t(combn(position.matrix, 2))
    # 
    # for( i in 1:nrow(norm)) {
    #   
    #   t.1 <- which( comb[,1] == norm[i,1] & comb[,2] == norm[i,2] )
    #   
    #   if( length(t.1) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,1] , Pair.to = norm[i,2] ,  Mean.Time = 9e9999)) }
    #   
    #   t.2 <- which( comb[,1] == norm[i,2] & comb[,2] == norm[i,1] )
    #   
    #   if( length(t.2) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,2] , Pair.to = norm[i,1] ,  Mean.Time = 9e9999)) }
    #   
    #   
    # }
    
    net.function <<- sum
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    E(graph.obj)$weight = comb[,3]
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  return(list(comb,graph.obj))
  
}

## ------------------------------------------------------------------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

## -------------------

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## -------------------

max(Connectivity$Time.max)

network <- produce.network("Prob",Connectivity,30,TRUE,2,source.sink.xy)
g2 <- network[[2]]

gs <- g2
gs <- as.undirected(g2, mode = "collapse", edge.attr.comb = "max") # For Probabilities

clustering.method <- "edge.betweenness.community" # fastgreedy.community** leading.eigenvector.community walktrap.community fastgreedy.community clusters edge.betweenness.community(slow)

length.of.tests <- 100
tests.modularity <- round(seq(from=1,to=ecount(gs),length.out=length.of.tests))
e.weight <- edge.attributes(gs)$weight
e.weight <- sort(e.weight, decreasing = FALSE, index.return =TRUE)$ix

cl <- makeCluster(10) ; registerDoParallel(cl)
mods <- foreach( i = tests.modularity, .verbose=FALSE, .packages=c("igraph") ) %dopar% {
  graph <- delete.edges( gs, e.weight[seq(length=i)] )
  try( membership.graph <- get(clustering.method)(graph)$membership , silent=TRUE )
  if(exists("membership.graph")) { return( c( modularity(graph, membership.graph) , length(unique(membership.graph)) ) ) } else { return( c(NA,NA) )  }
}
stopCluster(cl) ; rm(cl)

ModTab <- cbind(tests.modularity,do.call("rbind",mods))
head(ModTab)

par(mar=c(5,4,4,5))
plot(ModTab[,1],ModTab[,2],type="l",col="black", lty=1,ylab="Modularity",xlab="Removed edges", axes = FALSE, ylim=c(0,1))
axis( 1, lwd = 1) ; axis( 2, lwd = 1)
par(new=TRUE)
plot(ModTab[,1],ModTab[,3],type="l",col="black", lty=2, axes = FALSE,ylim=c(0,100),xlab="",ylab="")
axis(4)
mtext("Number of clusters",side=4,line=3)
legend("topleft",col=c("black","black"),lty=c(1,2),legend=c("modularity","clusters"),bty ="n")

## --------------------------------

ModTab <- ModTab[which(ModTab[,3] <= 30),]
best.edges <- 2533

graph <- delete.edges(gs, e.weight[seq(length=best.edges)] )
membership.graph <- get(clustering.method)(graph)$membership
n.clusters <- length(unique(membership.graph))
n.clusters
modularity(graph,get(clustering.method)(graph)$membership) 

library(randomcoloR)
#landmass <- crop(shapefile(landmass.shp),extent(-62,-55,-60.25,-50))
plot(landmass, col="grey" , border="grey")
points(source.sink.xy[as.numeric(vertex.attributes(graph)$name),2:3],pch=20,cex=0.4,col=distinctColorPalette(length(unique(membership.graph)))[membership.graph])

raster.memberships <- rasterFromXYZ(cbind(source.sink.xy[as.numeric(vertex.attributes(graph)$name),2:3],membership.graph))
projection(raster.memberships) <- CRS("+proj=longlat +datum=WGS84")
plot(raster.memberships,col=sample(rainbow(n.clusters),replace=F),box=FALSE,legend=FALSE)
writeRaster(raster.memberships,filename=paste0(results.directory,"clustering7"),format="GTiff",overwrite=T) 

## --------------------------------

n.cells <- length(get(clustering.method)(graph)$membership)
unique.members <- sort(unique(get(clustering.method)(graph)$membership))

permutat <- 4999
mods <- sapply(1:permutat, function(i){
  modularity( graph , sample( unique.members , n.cells , replace = TRUE) )
})

signif(sum( mods > modularity(graph,get(clustering.method)(graph)$membership) ) / permutat,digits=4)

# ----------

V(graph)$size <- 16 * (evcent(graph)$vector) # Centrality score # edge.betweenness(graph) degree(graph)    
V(graph)$label <- 1:npops
V(graph)$color <- membership.graph

labels <- rep(NA,npops)
labels <- 1:npops

comm <- get(clustering.method)(graph)
plot(comm,graph, vertex.label=1:24,vertex.size=1, vertex.color=membership.graph)

centrality <- alpha_centrality(graph, nodes = V(graph), alpha = 1, loops = FALSE, exo = 1, weights = NULL, tol = 1e-07, sparse = TRUE)

## -----------------------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------
