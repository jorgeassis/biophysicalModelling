## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")
source("Dependences.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

if( ! exists("pipeLiner") ) { 
  
  n.days <- 30
  
}

## -------------------------

if( ! paste0("connectivityExport","_",n.days,"Days") %in% file.path(paste0(project.folder,"/Results/",project.name)) ) { dir.create(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport","_",n.days,"Days/"))) }
connectivityExportDir <- file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport","_",n.days,"Days/"))

## -------------------------

sampling.sites <- read.csv("../Data/coords.csv",sep=";",dec=",",header=T)
colnames(sampling.sites) <- c("Code","Lon","Lat")

## -------------------------

Connectivity <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveragedPolys.bm") )
Connectivity <- as.data.frame(distance.probability[,])
colnames(Connectivity) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events")
Connectivity[ which(Connectivity[,"Time.max"] > n.days) , "Probability" ] <- 0

source.sink.xy <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm") )
source.sink.xy <- as.data.frame(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

## -------------------------

sampling.sites.address <- geo_to_h3(sampling.sites[,c("Lat","Lon")], sim.resolution )
sampling.sites.pair <- sapply(1:nrow(sampling.sites), function(x) { source.sink.xy$Pair[which.min(spDistsN1(as.matrix(source.sink.xy[,c("Lon","Lat")]),as.matrix(sampling.sites[x,c("Lon","Lat")],ncol=2),longlat = TRUE))] } )
sampling.sites <- data.frame(h3Address=sampling.sites.address,pair=sampling.sites.pair,sampling.sites)

write.table(file=paste0(connectivityExportDir,"samplingSites.csv"),x=sampling.sites,row.names=FALSE,col.names=TRUE,sep=";",dec=".")

## -------------------------

new.extent <- c(min(source.sink.xy[,2]),max(source.sink.xy[,2]),min(source.sink.xy[,3]),max(source.sink.xy[,3]))
network <- produce.network("Prob",distance.probability,n.days,FALSE,10,source.sink.xy,new.extent)

## -------------------------
## Direct Connectivity

matrixConnectivityP <- matrix(NA,ncol=nrow(sampling.sites),nrow=nrow(sampling.sites))
matrixConnectivityT <- matrix(NA,ncol=nrow(sampling.sites),nrow=nrow(sampling.sites))

for( i in 1:nrow(sampling.sites) ) {
  
  for( j in 1:nrow(sampling.sites) ) {
    
    i.pair <- sampling.sites$pair[i]
    j.pair <- sampling.sites$pair[j]
    
    if(length(which(distance.probability$Pair.from == i.pair & distance.probability$Pair.to == j.pair)) > 0) {
        matrixConnectivityP[i,j] <- distance.probability[distance.probability$Pair.from == i.pair & distance.probability$Pair.to == j.pair, "Probability"]
        matrixConnectivityT[i,j] <- distance.probability[distance.probability$Pair.from == i.pair & distance.probability$Pair.to == j.pair, "Mean.Time"]
    }
  }
}

colnames(matrixConnectivityP) <- sampling.sites$pair
rownames(matrixConnectivityP) <- sampling.sites$pair
colnames(matrixConnectivityT) <- sampling.sites$pair
rownames(matrixConnectivityT) <- sampling.sites$pair

matrixConnectivityP <- melt(matrixConnectivityP)
matrixConnectivityT <- melt(matrixConnectivityT)

## -------------------------
## Stepping-stone Connectivity

cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)

potential.connectivity <- foreach(from=sampling.sites$pair, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  network.x <- network[[2]]
  connectivity.x <- network[[1]]
  res.connectivity.to <- numeric(0)

  for( to in sampling.sites$pair ) {
    
    possible.paths.y <- get.shortest.paths(network.x,as.character( from ) , as.character( to ),mode="out")$vpath
    stones.t <- as.numeric(names(possible.paths.y[[1]]))
    stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
    path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
    
    if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
    if( length(path.values) == 0) { path.values <- 0 }
    if( from == to ) { path.values <- 1 }
    
    res.connectivity.to <- c(res.connectivity.to,path.values)

  }
  
  temp.res <- data.frame( pair.from = from,
                          pair.to = sampling.sites$pair , 
                          connectivity = res.connectivity.to)
  return(temp.res)

}

stopCluster(cl.3) ; rm(cl.3) ; gc()

connectivity <- do.call(rbind,potential.connectivity)

## -------------------------
## Export Connectivity

## Matrix

write.table(file=paste0(connectivityExportDir,"matrixDirectConnectivityP.csv"),x=matrixConnectivityP,row.names=FALSE,col.names=FALSE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixDirectConnectivityT.csv"),x=matrixConnectivityT,row.names=FALSE,col.names=FALSE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixSSConnectivityP.csv"),x=connectivity,row.names=FALSE,col.names=FALSE,sep=";",dec=".")

## Graph

graph.obj <- graph.edgelist( cbind( as.character( connectivity$pair.from) , as.character(connectivity$pair.to) ) , directed = TRUE )

# E(graph.obj)$weight = connectivity$connectivity
# E(graph.obj)$weight = ifelse(-log(connectivity$connectivity) == Inf,0,-log(connectivity$connectivity)) # Hock, Karlo Mumby, Peter J 2015
E(graph.obj)$weight = ifelse(-log(connectivity$connectivity) == Inf,0,1/-log(connectivity$connectivity)) # Hock, Karlo Mumby, Peter J 2015

graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "min") # min / mean / max
graph.obj <- simplify(graph.obj,remove.multiple = TRUE, remove.loops = TRUE)

# clustering.graph <- clusters(graph.obj)
# clustering.graph <- walktrap.community(graph.obj)
# clustering.graph <- cluster_fast_greedy(graph.obj)
clustering.graph <- cluster_leading_eigen(graph.obj,options=list(maxiter=1000000))

membership.graph <- clustering.graph$membership
reducedNames <- sapply(names(V(graph.obj)),function(x) { sampling.sites$Code[sampling.sites$pair == x][1] })

cols.to.use <- distinctColors(length(unique(membership.graph)))
cols.to.use <- cols.to.use[membership.graph]
V(graph.obj)$color <- cols.to.use
l <- layout.fruchterman.reingold(graph.obj)

pdf( file=paste0(connectivityExportDir,"sitesNetworkClustering.pdf") , width = 10 )
plot(graph.obj,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

# To implete Yet !

weights <- (((-log(comb[,3]) - (max(-log(comb[,3])))) * (-1)) / 100)
weights <- comb[,3]
weights.reclass <- weights / max(weights)
weights.reclass[weights.reclass >= 0.75] <- 8
weights.reclass[weights.reclass >= 0.5 & weights.reclass < 0.75] <- 4
weights.reclass[weights.reclass >= 0.25 & weights.reclass < 0.5] <- 1
weights.reclass[weights.reclass < 0.25 ] <- 0.25
plot(get(clustering.method)(graph.obj),graph.obj,vertex.label.color="Black",vertex.label.family="Helvetica",edge.width=weights.reclass,edge.color="Black") 


# Simplify

graph.obj.simp <- graph.obj

repeat {
  graph.obj.simp.temp <- delete_edges(graph.obj.simp, which.min(E(graph.obj.simp)$weight))
  if(clusters(graph.obj.simp.temp)$no == 1) { graph.obj.simp <- graph.obj.simp.temp }
  if(clusters(graph.obj.simp.temp)$no != 1) { break }
}

l <- layout.fruchterman.reingold(graph.obj.simp)

pdf( file=paste0(connectivityExportDir,"sitesNetworkClusteringReduced.pdf") , width = 10 )
plot(graph.obj.simp,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , color=cols.to.use )
dev.off()

## Clustering [Map]

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
    panel.grid.major = element_line(color = "#979797", size = 0.05),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank()
  )

boxes <- data.frame(maxlat = max.lat+2,minlat = min.lat-2,maxlong = max.lon+2,minlong = min.lon-2, id="1")
boxes <- transform(boxes, laby=(min.lat + min.lat )/2, labx=(max.lon+min.lon )/2)

boxesSim <- data.frame(maxlat = max.lat,minlat = min.lat,maxlong = max.lon,minlong = min.lon, id="1")
boxesSim <- transform(boxesSim, laby=(min.lat + min.lat )/2, labx=(max.lon+min.lon )/2)

if( is.null(landmass.shp)) { land.polygon <- getMap(resolution = "high") }
if( ! is.null(landmass.shp)) { land.polygon <- shapefile(landmass.shp) }
crs(land.polygon) <- dt.projection 
land.polygon <- crop(land.polygon,extent(c(min.lon, max.lon,  min.lat,max.lat)))
land.polygon@bbox <- as.matrix(extent(c(min.lon-1, max.lon+1,  min.lat+1,max.lat+1)))

mainTitle <- "Laminaria pallida :: SW Africa"

map <- ggplot() +
  geom_polygon(data = land.polygon , fill = "#C4C4C4", colour = "#ffffff" , size=0.15 ,  aes(long, lat, group = group)) +
  #geom_rect(data=boxes, aes(xmin=min.lon , xmax=max.lon, ymin=min.lat, ymax=max.lat ), color="transparent", fill="transparent") +
  coord_equal() + theme_map # coord_map(projection = "mercator")
plot(map)

map.links <- map
connected.pairs <- as_long_data_frame(graph.obj.simp)
connected.pairs <- connected.pairs[sort(connected.pairs$weight, index.return=T , decreasing = FALSE )$ix,]

for( i in 1:nrow(connected.pairs) ){
  strenght <- (connected.pairs[i,3] * 100) + 1 
  sampling.sites.from <- sampling.sites[which(sampling.sites$pair == connected.pairs[i,"from_name"])[1],c("Lon","Lat")]
  coordinates(sampling.sites.from) <- ~Lon+Lat
  sampling.sites.to <- sampling.sites[which(sampling.sites$pair == connected.pairs[i,"to_name"])[1],c("Lon","Lat")]
  coordinates(sampling.sites.to) <- ~Lon+Lat
  routes_sl <- gcIntermediate(sampling.sites.from,sampling.sites.to,n = 100, addStartEnd = TRUE, sp = TRUE)
  SLDF = sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = NA), match.ID = F)
  map.links <- map.links + geom_path(data = SLDF, size=0.2 , aes(x = long, y = lat), col = "#797979") # colfunc(101)[round(strenght)]
}

centroids <- sampling.sites[ sapply(as.numeric(names(V(graph.obj.simp))), function(x) which(sampling.sites$pair == x )[1] ) ,  c("Lon","Lat") ]

pdf( file=paste0(connectivityExportDir,"mapSitesNetworkClusteringLink.pdf") , width = 10 )
print(
  map.links + geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = cols.to.use, size = 2.5, stroke = 0.35, alpha = 0.9)
  )
dev.off()

pdf( file=paste0(connectivityExportDir,"mapSitesNetworkClustering.pdf") , width = 10 )
print(
  map + geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = cols.to.use, size = 3, stroke = 0.35, alpha = 0.9)
  )
dev.off()
