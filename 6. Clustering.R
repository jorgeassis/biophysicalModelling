## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

source("0. Project Config.R")

sql.project.name <- "Atlantic"

number.cores <- 5

## ------------------------------------------------------------------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",sql.project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity[is.na(Connectivity)] <- 0

Connectivity ; gc()

## -------------------

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## -------------------

max(Connectivity$Max.Time)
n.days <- 30

comb <- Connectivity
comb[ which(comb[,"Max.Time"] > n.days) , "Mean.Probability" ] <- 0
comb <- comb[,c("Pair.from","Pair.to","Mean.Probability")]
comb <- comb[comb$Pair.from != comb$Pair.to,]
comb <- as.data.frame( comb[ sort( as.vector(unlist(comb[,"Mean.Probability"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )

net.function <<- prod
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )

# E(graph.obj)$weight = ifelse(log(comb[,3]) == -Inf,0, comb[,3]) # Hock, Karlo Mumby, Peter J 2015

# E(graph.obj)$weight = comb[,3] # The wheight has a negative impact on finding the closest path
 E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,1/-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015

graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "min") # min / mean / max
graph.obj <- simplify(graph.obj, remove.multiple = TRUE, remove.loops = TRUE)

clustering.method <- "walktrap.community" # Uni: fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)

## --------------------------------
# 
# length.of.tests <- 100
# tests.modularity <- round(seq(from=1,to=ecount(gs),length.out=length.of.tests))
# e.weight <- edge.attributes(gs)$weight
# e.weight <- sort(e.weight, decreasing = FALSE, index.return =TRUE)$ix
# 
# cl <- makeCluster(10) ; registerDoParallel(cl)
# mods <- foreach( i = tests.modularity, .verbose=FALSE, .packages=c("igraph") ) %dopar% {
#   graph <- delete.edges( gs, e.weight[seq(length=i)] )
#   try( membership.graph <- get(clustering.method)(graph)$membership , silent=TRUE )
#   if(exists("membership.graph")) { return( c( modularity(graph, membership.graph) , length(unique(membership.graph)) ) ) } else { return( c(NA,NA) )  }
# }
# stopCluster(cl) ; rm(cl)
# 
# ModTab <- cbind(tests.modularity,do.call("rbind",mods))
# head(ModTab)
# 
# par(mar=c(5,4,4,5))
# plot(ModTab[,1],ModTab[,2],type="l",col="black", lty=1,ylab="Modularity",xlab="Removed edges", axes = FALSE, ylim=c(0,1))
# axis( 1, lwd = 1) ; axis( 2, lwd = 1)
# par(new=TRUE)
# plot(ModTab[,1],ModTab[,3],type="l",col="black", lty=2, axes = FALSE,ylim=c(0,100),xlab="",ylab="")
# axis(4)
# mtext("Number of clusters",side=4,line=3)
# legend("topleft",col=c("black","black"),lty=c(1,2),legend=c("modularity","clusters"),bty ="n")## --------------------------------
# 
# head(ModTab)
# 
# ModTab <- ModTab[ModTab[,3] == min(ModTab[,3]),]
# ModTab <- ModTab[ModTab[,2] == max(ModTab[,2]),]
# ModTab
# 
# best.edges <- 1
# graph <- gs
# graph <- delete.edges(gs, e.weight[seq(length=best.edges)] )
# graph <- simplify(graph,remove.multiple = TRUE)

membership.graph <- get(clustering.method)(graph.obj)$membership
membership.graph.cells <- data.frame(Cell=source.sink.xy[as.numeric(vertex.attributes(graph.obj)$name),1],source.sink.xy[as.numeric(vertex.attributes(graph.obj)$name),2:3],Membership=membership.graph)

# Remove isolated cells

# membership.graph.t <- aggregate(membership.graph.cells$Membership, by=list(Category=membership.graph.cells$Membership), FUN=sum)
# membership.graph.t <- membership.graph.t[sort(membership.graph.t$x,decreasing = FALSE,index.return = T)$ix,]
# hist(membership.graph.cells$Membership)
# isolated.cells <- membership.graph.t[which(membership.graph.t$x == 1),1]
# membership.graph.cells <- membership.graph.cells[ ! membership.graph.cells$Cell %in% isolated.cells, ]
# 
# ## --------------------------------
# 
# landmass <- shapefile(landmass.shp)
# crs(landmass) <- dt.projection
# landmass <- gBuffer(landmass, byid=TRUE, width=0)
# clipper <- as(extent(min(source.sink.xy[,2] - 2),max(source.sink.xy[,2] + 2),min(source.sink.xy[,3] - 2),max(source.sink.xy[,3] + 2)), "SpatialPolygons")
# crs(clipper) <- dt.projection
# landmass <- gIntersection(landmass, clipper, byid=TRUE)

## --------------------------------

# 2000px | 12 12 Clustering_1Day

plot(landmass, col="grey" , border="grey")
#points(membership.graph.cells[,2:3],pch=20,cex=0.4,col= sample(rainbow(max(membership.graph.cells$Membership)) , length(rainbow(max(membership.graph.cells$Membership))) , replace=FALSE)[membership.graph.cells$Membership])
points(membership.graph.cells[,2:3],pch=20,cex=0.4,col=distinctColorPalette(length(unique(membership.graph.cells$Membership)))[membership.graph.cells$Membership])

modularity(graph.obj,get(clustering.method)(graph.obj)$membership) 
n.clusters <- length(unique(membership.graph.cells$Membership))
n.clusters

## --------------------------------
## Tweek

position.matrix <- spDists(as.matrix(membership.graph.cells[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
membership.graph.cells[position.matrix,]

# membership.graph.cells[membership.graph.cells$y >= & membership.graph.cells$y >= , "Membership"] 
# membership.graph.cells[membership.graph.cells$y >= & membership.graph.cells$y >= , "Membership"] <- 

## --------------------------------

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

## ----------------------------------------

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
