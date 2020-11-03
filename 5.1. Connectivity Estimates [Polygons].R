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

# n.days <- 30

## -------------------------

if( ! paste0("connectivityExport") %in% list.files(file.path(paste0(project.folder,"/Results/",project.name))) ) { dir.create(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/"))) }
if( ! paste0("Sim",str_pad(n.days, 3, pad = "0"),"Days") %in% list.files(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/"))) ) { dir.create(file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/","Sim",str_pad(n.days, 3, pad = "0"),"Days"))) }
connectivityExportDir <- file.path(paste0(project.folder,"/Results/",project.name,"/connectivityExport/","Sim",str_pad(n.days, 3, pad = "0"),"Days/"))

## -------------------------

coords.centroid.ids <- data.frame(Name=c("Flores","Corvo","Faial","Pico","Sao Jorge","Graciosa","Sao Miguel","Sta Maria", "Terceira"),
                                  Lon=c(-31.203339, -31.105835, -28.701468, -28.329306, -28.058767, -28.012076, -25.477667, -25.102516, -27.222877),
                                  Lat=c(39.440657, 39.696841, 38.573669, 38.467299, 38.654148, 39.054118, 37.772072, 36.975784, 38.723658))

save(coords.centroid.ids, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","coords.centroid.ids.RData"))

## -------------------------

additional.source.sink <- shapefile(additional.source.sink.shp)
additional.source.sink <- additional.source.sink[,"ID"]
writeOGR(obj=additional.source.sink, dsn=paste0(project.folder,"/Results/",project.name,"/"), layer="sourceSinkPolys", driver="ESRI Shapefile") 

## -------------------------

Connectivity <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveragedPolys.bm") )
Connectivity <- as.data.frame(Connectivity[,])
colnames(Connectivity) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events")
Connectivity[ which(Connectivity[,"Time.max"] > n.days) , "Probability" ] <- 0
head(Connectivity)

source.sink.xy <- read.big.matrix( paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.Polys.bm") )
source.sink.xy <- as.data.frame(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
head(source.sink.xy)

## -------------------------

if( is.null(landmass.shp)) { land.polygon <- getMap(resolution = "high") }
if( ! is.null(landmass.shp)) { land.polygon <- shapefile(landmass.shp) }
land.polygon <- disaggregate(land.polygon)
crs(land.polygon) <- dt.projection 
land.polygon <- crop(land.polygon,extent(min.lon,max.lon,min.lat,max.lat))

theme_map <- 
  theme_minimal() +
  theme(text = element_text(family = "Helvetica", color = "#22211d"),
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
        panel.border = element_blank())

mapBuffer <- 0.5
boxes <- data.frame(maxlat = extent(land.polygon)[4]+mapBuffer,minlat = extent(land.polygon)[3]-mapBuffer,maxlong = extent(land.polygon)[2]+mapBuffer,minlong = extent(land.polygon)[1]-mapBuffer, id="1")

map <- ggplot() +
  geom_polygon(data = land.polygon , fill = "#A6A6A6", colour = "#ffffff" , size=0.15 ,  aes(long, lat, group = group)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") +
  coord_equal() + theme_map # coord_map(projection = "mercator")
map

## -----------------------------------------------------------------
## -----------------------------------------------------------------

## Assign IDs [e.g., region of interest] to sink.source sites

head(coords.centroid.ids)
coords.centroid.ids.pts <- coords.centroid.ids
coordinates(coords.centroid.ids.pts) <- ~Lon+Lat
crs(coords.centroid.ids.pts) <- crs(land.polygon)

land.polygon <- land.polygon[ sapply( 1:length(coords.centroid.ids.pts) , function(x) { which(!is.na(over(land.polygon,coords.centroid.ids.pts[x,]))) } ) , ]
land.polygon$ROI <- coords.centroid.ids$Name
writeOGR(obj=land.polygon, dsn=paste0(project.folder,"/Results/",project.name,"/"), layer="sourceSinkPolygons", driver="ESRI Shapefile") 

source.sink.xy.pts <- source.sink.xy
coordinates(source.sink.xy.pts) <- ~Lon+Lat
crs(source.sink.xy.pts) <- crs(land.polygon)

source.sink.xy.ids <- sapply( 1:length(source.sink.xy.pts) , function(x) { coords.centroid.ids[as.numeric(apply(gDistance(source.sink.xy.pts[x,], land.polygon,byid=TRUE),2,which.min)),"Name"] } )
save(source.sink.xy.ids, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.ids.RData"))
source.sink.xy <- data.frame(source.sink.xy,Name=source.sink.xy.ids)
head(source.sink.xy)

## --------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------
## Force all Pairs in Matrices

missing.Pair.from <- source.sink.xy$Pair[which( ! source.sink.xy$Pair %in% Connectivity$Pair.from)]
if( length(missing.Pair.from) > 0 ) {
  for( i in 1:length(missing.Pair.from)) {
    Connectivity <- rbind(Connectivity,Connectivity[1,])
    Connectivity[nrow(Connectivity),"Pair.from"] <- missing.Pair.from[i]
    Connectivity[nrow(Connectivity),3:11] <- 0
  }
}

missing.Pair.to <- source.sink.xy$Pair[which( ! source.sink.xy$Pair %in% Connectivity$Pair.to)]
if( length(missing.Pair.to) > 0 ) {
  for( i in 1:length(missing.Pair.to)) {
    Connectivity <- rbind(Connectivity,Connectivity[1,])
    Connectivity[nrow(Connectivity),"Pair.to"] <- missing.Pair.to[i]
    Connectivity[nrow(Connectivity),3:11] <- 0
  }
}

## --------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------
## Pairwise connectivity estimates for polygons
  
## -------------------------
## Direct Connectivity

matrixConnectivityP.square <- acast(Connectivity[,c("Pair.from","Pair.to","Probability")], formula=Pair.from~Pair.to, fun.aggregate= mean , value.var="Probability",drop = FALSE)
matrixConnectivityT.square <- acast(Connectivity[,c("Pair.from","Pair.to","Mean.Time")], formula=Pair.from~Pair.to,  fun.aggregate=mean, value.var="Mean.Time",drop = FALSE)

dim(matrixConnectivityP.square)
dim(matrixConnectivityT.square)

matrixConnectivityP <- melt(matrixConnectivityP.square)
colnames(matrixConnectivityP) <- c("Pair.from","Pair.to","Probability")
matrixConnectivityT <- melt(matrixConnectivityT.square)
colnames(matrixConnectivityT) <- c("Pair.from","Pair.to","Mean.Time")
  
## -------------------------

## Stepping-stone Connectivity

## Produce network to estimate least cost distances

networkData <- matrixConnectivityP
networkData <- networkData[!is.na(networkData$Probability),]
networkData <- as.data.frame( networkData[ sort( as.vector(unlist(networkData[,"Probability"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )

network <- graph.edgelist( cbind( as.character( networkData[,1]) , as.character(networkData[,2]) ) , directed = TRUE )
E(network)$weight = ifelse(-log(networkData[,3]) == Inf,0,-log(networkData[,3])) # Hock, Karlo Mumby, Peter J 2015
network <- delete.edges(network, which(E(network)$weight ==0))
network <- simplify(network)

cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)

matrixConnectivitySSP <- foreach(from=source.sink.xy$Pair, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  possible.paths <- get.shortest.paths(network,as.character( from ) ,mode="out")
  possible.paths <- possible.paths[[1]]

  res.connectivity.to <- data.frame( pair.From=from,pair.To=source.sink.xy$Pair, connectivity=NA)
                                     
  for( i in 1:length(possible.paths) ) {
    
    possible.paths.i <- as.numeric(names(possible.paths[i][[1]]))
    to <- possible.paths.i[length(possible.paths.i)]
    to.i <- which(res.connectivity.to$pair.To == to)
    
    if( from == to ) { path.values <- networkData[ networkData[,1] == from & networkData[,2] == to , 3 ][1]  }
    
    if( from != to ) {   
          
          stones.t.interm <- cbind(possible.paths.i[-length(possible.paths.i)],possible.paths.i[-1])
          path.values <- apply( stones.t.interm , 1 , function(z) { networkData[ networkData[,1] == z[1] & networkData[,2] == z[2] , 3 ][1] }   )
          
          if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
          if( length(path.values) == 0) { path.values <- 0 }
          
    }
    
    res.connectivity.to[to.i,2] <- to
    res.connectivity.to[to.i,3] <- path.values
    
  }
  
  return( res.connectivity.to )
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

matrixConnectivitySSP <- do.call(rbind,matrixConnectivitySSP)
colnames(matrixConnectivitySSP) <- c("From","To","Probability")
matrixConnectivitySSP.square <- acast(matrixConnectivitySSP, formula=From~To, fun.aggregate= mean , value.var="Probability",drop = FALSE)

## -----------------------------------------
## Export Connectivity

# Square Matrix

dim(matrixConnectivityP.square)
dim(matrixConnectivityT.square)
dim(matrixConnectivitySSP.square)

all.equal(colnames(matrixConnectivityP.square),colnames(matrixConnectivitySSP.square))

matrixConnectivityP.square[is.na(matrixConnectivityP.square)] <- 0
matrixConnectivitySSP.square[is.na(matrixConnectivitySSP.square)] <- 0

write.table(file=paste0(connectivityExportDir,"matrixConnectivityP.square.csv"),x=matrixConnectivityP.square,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityT.square.csv"),x=matrixConnectivityT.square,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivitySSP.square.csv"),x=matrixConnectivitySSP.square,row.names=TRUE,col.names=TRUE,sep=";",dec=".")

matrixConnectivityP.square <- read.csv(paste0(connectivityExportDir,"matrixConnectivityP.square.csv"),sep=";" , header=T,check.names=FALSE)
matrixConnectivityT.square <- read.csv(paste0(connectivityExportDir,"matrixConnectivityT.square.csv"),sep=";" , header=T,check.names=FALSE)
matrixConnectivitySSP.square <- read.csv(paste0(connectivityExportDir,"matrixConnectivitySSP.square.csv"),sep=";" , header=T,check.names=FALSE)

# Pairs

matrixConnectivityP[is.na(matrixConnectivityP)] <- 0
matrixConnectivitySSP[is.na(matrixConnectivitySSP)] <- 0

write.table(file=paste0(connectivityExportDir,"matrixConnectivityP.csv"),x=matrixConnectivityP,row.names=FALSE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityT.csv"),x=matrixConnectivityT,row.names=FALSE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivitySSP.csv"),x=matrixConnectivitySSP,row.names=FALSE,col.names=TRUE,sep=";",dec=".")

matrixConnectivityP <- read.csv(paste0(connectivityExportDir,"matrixConnectivityP.csv"),sep=";" , header=F,check.names=FALSE)
matrixConnectivityT <- read.csv(paste0(connectivityExportDir,"matrixConnectivityT.csv"),sep=";" , header=F,check.names=FALSE)
matrixConnectivitySSP <- read.csv(paste0(connectivityExportDir,"matrixConnectivitySSP.csv"),sep=";" , header=F,check.names=FALSE)

## --------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------
## Plot connectivity [networks]

## -------------------------
## Produce network

networkData <- matrixConnectivityP
networkData <- as.data.frame( networkData[ sort( as.vector(unlist(networkData[,"Probability"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )

network <- graph.edgelist( cbind( as.character( networkData[,1]) , as.character(networkData[,2]) ) , directed = TRUE )
E(network)$weight <- networkData[,3]
network <- delete.edges(network, which(E(network)$weight ==0))
network <- as.undirected(network, mode = "collapse", edge.attr.comb = "min") # min / mean / max
network <- simplify(network)

network.clustering <- cluster_walktrap(network) # network.clustering <- cluster_leading_eigen(network,options=list(maxiter=1000000))
network.membership <- network.clustering$membership
cols.to.use <- distinctColors(length(unique(network.membership)))
cols.to.use <- cols.to.use[network.membership]
V(network)$color <- cols.to.use

network.modularity <- modularity(network,network.membership)

lineThinkness <- reclassVals(E(network)$weight,min(E(network)$weight),max(E(network)$weight))

l <- source.sink.xy[sapply(get.vertex.attribute(network)$name,function(x) { which( source.sink.xy$Pair == x) } ),c("Lon","Lat")]
l <- as.matrix(l)

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
plot(network,edge.width=lineThinkness,vertex.label=NA,vertex.size=2,edge.curved = F , color=cols.to.use , layout=l,pch=20 )
dev.off()

additional.source.sink$Membership <- sapply(additional.source.sink$ID,function(x) { network.membership[which(get.vertex.attribute(network)$name == x)] })
additional.source.sink.sf <- st_as_sf(additional.source.sink)

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(Membership),color = as.factor(Membership)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(network.membership)), values = unique(cols.to.use))
dev.off()

## -------------------------------------------------------
## -------------------------------------------------------
## Export major results

connectivityMatrix <- matrixConnectivityP.square
diag(connectivityMatrix) <- 0

connectivityMatrixreclass <- connectivityMatrix
connectivityMatrixreclass[ connectivityMatrixreclass != 0] <- 1

isolated <- which(apply(connectivityMatrix,1,sum,na.rm=T) == 0 & apply(connectivityMatrix,2,sum,na.rm=T) == 0)
isolated.id <- colnames(connectivityMatrix)[isolated]
isolated.length <- length(isolated)

linkage.export <- apply(connectivityMatrixreclass,1,sum,na.rm=T)
linkage.import <- apply(connectivityMatrixreclass,2,sum,na.rm=T)

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=12) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) )

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
ggplot(data=data.frame(links=unique(linkage.export),sum=sapply(unique(linkage.export),function(x) { sum( linkage.export == x)  })), aes(x=links, y=sum)) +
  geom_bar(stat="identity", fill="#236292") + ylab("Number of sites") + xlab("Number of connections between sites [export]") + mainTheme
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
ggplot(data=data.frame(links=unique(linkage.import),sum=sapply(unique(linkage.import),function(x) { sum( linkage.import == x)  })), aes(x=links, y=sum)) +
  geom_bar(stat="identity", fill="#236292") + ylab("Number of sites") + xlab("Number of connections between sites [import]") + mainTheme
dev.off()

linkage <- apply(rbind( linkage.export , linkage.import),2,mean)
higherLinkage <- which(linkage >=  as.numeric(quantile(linkage,0.95,na.rm=TRUE))) 
linkage.on.average <- mean(linkage)

mean.probability <- mean( apply(rbind( apply(connectivityMatrix,1,mean,na.rm=T) , apply(connectivityMatrix,2,mean,na.rm=T)),2,mean) )

numberClusters <- length(unique(network.membership)) - length(isolated)
aggregationBasedOnClusters <- 1 - ( numberClusters / length(network.membership) )

resistance <- sapply( 1:nrow(connectivityMatrixreclass) , function(res) { 1- length(unique(c(which(connectivityMatrixreclass[res,] == 1 ) ,  which(connectivityMatrixreclass[,res] == 1 )))) / nrow(connectivityMatrixreclass) } )
higherResistance <- which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=TRUE)))

betweennessIndex <- betweenness(network)
higherBetweenness <- which(betweennessIndex >=  as.numeric(quantile(betweennessIndex,0.95,na.rm=TRUE)))
average.betweenness <- mean(betweennessIndex)

eigencentralityIndex <- eigen_centrality(network, directed = FALSE, scale = FALSE)$vector
higherEigencentrality <- which(eigencentralityIndex >=  as.numeric(quantile(eigencentralityIndex,0.95,na.rm=TRUE)))
average.eigencentrality <- mean(eigencentralityIndex)

closenessIndex <- closeness(network, mode="all")
higherCloseness <- which(closenessIndex >=  as.numeric(quantile(closenessIndex,0.95,na.rm=TRUE)))
average.closeness <- mean(closenessIndex)

# ---------------------

if( exists("pipeLiner") ) {
  
  c <- which(n.repetitions == n.days)
  
  if( ! exists("isolatedAll") ) { 
    summaryAll <- data.frame()
    isolatedAll <- data.frame( matrix(0,nrow=nrow(connectivityMatrix),ncol=length(n.repetitions)), row.names=colnames(connectivityMatrix))
    colnames(isolatedAll) <- n.repetitions
    higherresistanceAll <- resistanceAll <- clusterAssignmentAll <- higherclosenessAll <- closenessAll <- highereighenCentralityAll <- eighenCentralityAll <- higherBetweennessAll <- betweennessAll <- higherLinkageLevelAll <- linkageLevelAll <- isolatedAll
  }
  
  isolatedAll[ isolated ,c] <- 1
  write.csv(isolatedAll,file=paste0(connectivityExportDir,"isolatedAll.csv"))
  linkageLevelAll[,c] <- linkage
  write.csv(linkageLevelAll,file=paste0(connectivityExportDir,"linkageLevelAll.csv"))
  higherLinkageLevelAll[higherLinkage,c] <- 1
  write.csv(higherLinkageLevelAll,file=paste0(connectivityExportDir,"higherLinkageLevelAll.csv"))
  clusterAssignmentAll[,c] <- network.membership
  write.csv(clusterAssignmentAll,file=paste0(connectivityExportDir,"clusterAssignmentAll.csv"))
  betweennessAll[,c] <- betweennessIndex
  write.csv(betweennessAll,file=paste0(connectivityExportDir,"betweennessAll.csv"))
  higherBetweennessAll[higherBetweenness,c] <- 1
  write.csv(higherBetweennessAll,file=paste0(connectivityExportDir,"higherBetweennessAll.csv"))
  eighenCentralityAll[,c] <- eigencentralityIndex
  write.csv(eighenCentralityAll,file=paste0(connectivityExportDir,"eighenCentralityAll.csv"))
  highereighenCentralityAll[higherEigencentrality,c] <- 1
  write.csv(highereighenCentralityAll,file=paste0(connectivityExportDir,"highereighenCentralityAll.csv"))
  closenessAll[,c] <- closenessIndex
  write.csv(closenessAll,file=paste0(connectivityExportDir,"closenessAll.csv"))
  higherclosenessAll[higherCloseness,c] <- 1
  write.csv(higherclosenessAll,file=paste0(connectivityExportDir,"higherclosenessAll.csv"))
  resistanceAll[,c] <- resistance
  write.csv(resistanceAll,file=paste0(connectivityExportDir,"resistanceAll.csv"))
  higherresistanceAll[higherResistance,c] <- 1
  write.csv(higherresistanceAll,file=paste0(connectivityExportDir,"higherresistanceAll.csv"))
  
}

additional.source.sink$Isolated <- isolatedAll[,c]
additional.source.sink$Links <- linkage
additional.source.sink$HigherLinks <- higherLinkageLevelAll[,c]
additional.source.sink$Betweenness <- betweennessAll[,c]
additional.source.sink$HigherBetweenness <- higherBetweennessAll[,c]
additional.source.sink$EighenCentrality <- eighenCentralityAll[,c]
additional.source.sink$HighereighenCentrality <- higherLinkageLevelAll[,c]
additional.source.sink$Closeness <- closenessAll[,c]
additional.source.sink$HigherCloseness <- higherclosenessAll[,c]
additional.source.sink$Resistance <- resistanceAll[,c]
additional.source.sink$HigherResistance <- higherresistanceAll[,c]
writeOGR(obj=additional.source.sink, dsn=connectivityExportDir, layer="sourceSinkPolygonsData", driver="ESRI Shapefile") 

summaryAll <- rbind(summaryAll,data.frame( day.sim=n.days,
                                           numberClusters=numberClusters,
                                           aggregation=aggregationBasedOnClusters,
                                           n.isolated=isolated.length,
                                           links.on.average=linkage.on.average,
                                           probability.on.average=mean.probability,
                                           betweenness.on.average=average.betweenness,
                                           eigencentrality.on.average=average.eigencentrality,
                                           closeness.on.average=average.closeness,
                                           modularity=network.modularity ))


write.csv(summaryAll,paste0(project.folder,"/Results/",project.name,"/connectivityExport/summaryAll.csv"))

# ---------------------

additional.source.sink.sf <- st_as_sf(additional.source.sink)

cols.to.use <- distinctColors(length(unique(additional.source.sink.sf$Isolated)))
cols.to.use <- c("white","red")

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(Isolated),color = as.factor(Isolated)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(additional.source.sink.sf$Isolated)), values = unique(cols.to.use))
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(HigherLinks),color = as.factor(HigherLinks)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(additional.source.sink.sf$HigherLinks)), values = unique(cols.to.use))
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(HigherBetweenness),color = as.factor(HigherBetweenness)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(additional.source.sink.sf$HigherBetweenness)), values = unique(cols.to.use))
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(HighereighenCentrality),color = as.factor(HighereighenCentrality)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(additional.source.sink.sf$HighereighenCentrality)), values = unique(cols.to.use))
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(higherclosenessAll),color = as.factor(higherclosenessAll)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(additional.source.sink.sf$higherclosenessAll)), values = unique(cols.to.use))
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
map + geom_sf(data=additional.source.sink.sf,aes(fill = as.factor(HigherResistance),color = as.factor(HigherResistance)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(additional.source.sink.sf$HigherResistance)), values = unique(cols.to.use))
dev.off()

## --------------------------------------------------------------------
## --------------------------------------------------------------------
## Aggregate and plot network assignments

matrixConnectivityP.agg <- data.frame(expand.grid(From=unique(source.sink.xy$Name),To=unique(source.sink.xy$Name)),Val=NA)
matrixConnectivityPSS.agg <- data.frame(expand.grid(From=unique(source.sink.xy$Name),To=unique(source.sink.xy$Name)),Val=NA)
matrixConnectivityT.agg <- data.frame(expand.grid(From=unique(source.sink.xy$Name),To=unique(source.sink.xy$Name)),Val=NA)

for ( i in 1:nrow(matrixConnectivityP.agg)) {
  
  from <- as.character(matrixConnectivityP.agg[i,1])
  from <- source.sink.xy[which(source.sink.xy$Name == from),"Pair"]
  to <- as.character(matrixConnectivityP.agg[i,2])
  to <- source.sink.xy[which(source.sink.xy$Name == to),"Pair"]
  
  matrixConnectivityP.agg[i,3] <- mean(matrixConnectivityP[ matrixConnectivityP[,1] %in% from & matrixConnectivityP[,2] %in% to , "Probability" ])
  matrixConnectivityPSS.agg[i,3] <- mean(matrixConnectivitySSP[ matrixConnectivitySSP[,1] %in% from & matrixConnectivitySSP[,2] %in% to , "Probability" ])
  
  timeVal <- matrixConnectivityT[ matrixConnectivityT[,1] %in% from & matrixConnectivityT[,2] %in% to , "Mean.Time" ]
  timeVal <- mean(timeVal[timeVal != 0],na.rm=T)
  matrixConnectivityT.agg[i,3] <- timeVal

}

write.table(file=paste0(connectivityExportDir,"matrixConnectivityP.agg.csv"),x=matrixConnectivityP.agg,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityPSS.agg.csv"),x=matrixConnectivityPSS.agg,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityT.agg.csv"),x=matrixConnectivityT.agg,row.names=TRUE,col.names=TRUE,sep=";",dec=".")

matrixConnectivityP.agg.square <- acast(matrixConnectivityP.agg, formula=From~To, fun.aggregate= mean , value.var="Val",drop = FALSE)
matrixConnectivityPSS.agg.square <- acast(matrixConnectivityPSS.agg, formula=From~To, fun.aggregate= mean , value.var="Val",drop = FALSE)
matrixConnectivityT.agg.square <- acast(matrixConnectivityT.agg, formula=From~To, fun.aggregate= mean , value.var="Val",drop = FALSE)

write.table(file=paste0(connectivityExportDir,"matrixConnectivityP.agg.square.csv"),x=matrixConnectivityP.agg.square,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityPSS.agg.square.csv"),x=matrixConnectivityPSS.agg.square,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(connectivityExportDir,"matrixConnectivityT.agg.square.csv"),x=matrixConnectivityT.agg.square,row.names=TRUE,col.names=TRUE,sep=";",dec=".")

# ----------

matrixConnectivityP.agg <- matrixConnectivityP.agg[ sort( matrixConnectivityP.agg$Val , decreasing = TRUE, index.return =TRUE)$ix , ]
networkData <- matrixConnectivityP.agg

network <- graph.edgelist( cbind( as.character( networkData[,1]) , as.character(networkData[,2]) ) , directed = TRUE )
E(network)$weight <- networkData[,3]
network <- delete.edges(network, which(E(network)$weight ==0))
network <- as.undirected(network, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
network <- simplify(network)
network.clustering <- cluster_walktrap(network) # network.clustering <- cluster_leading_eigen(network)
network.membership <- network.clustering$membership

cols.to.use <- distinctColors(length(unique(network.membership)))
cols.to.use <- cols.to.use[network.membership]
V(network)$color <- cols.to.use
l <- layout.fruchterman.reingold(network)
l <- layout_nicely(network)

lineThinkness <- reclassVals(E(network)$weight,min(E(network)$weight),max(E(network)$weight))
  
pdf( file=paste0(connectivityExportDir,"idsNetworkClustering.pdf") , width = 10 )
plot(network,edge.width=lineThinkness,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

# ----------

l <- coords.centroid.ids[sapply(get.vertex.attribute(network)$name,function(x) { which(coords.centroid.ids == x) } ),c("Lon","Lat")]
l <- as.matrix(l)

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
plot(network,edge.width=lineThinkness,vertex.label.dist=2,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

## -------------------------------------------------------------
## Map network assignments

networkData <- matrixConnectivityP.agg
networkData <- networkData[networkData$Val != 0,]
mapNet <- ggplot()

for( i in nrow(comb):1 ){
  lineThinkness <- reclassVals(comb[i,3],min(comb[,3]),max(comb[,3]))
  routes_sl <- gcIntermediate(coords.centroid.ids[coords.centroid.ids$Name == comb[i,1],c("Lon","Lat")],
                              coords.centroid.ids[coords.centroid.ids$Name == comb[i,2],c("Lon","Lat")],
                              n = 100, addStartEnd = TRUE, sp = TRUE)
  SLDF = sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = NA), match.ID = F)
  mapNet <- mapNet + geom_path(data = SLDF, size=lineThinkness , aes(x = long, y = lat), col = "#797979") # colfunc(101)[round(strenght)]
}

coords.centroid.ids.pts <- coords.centroid.ids
coordinates(coords.centroid.ids.pts) <- ~Lon+Lat
crs(coords.centroid.ids.pts) <- crs(land.polygon)
land.polygon.ids.c <- get.vertex.attribute(network)$color

mapToPlot <- mapNet

for( i in 1:nrow(coords.centroid.ids)) {
  
  ids.i <- get.vertex.attribute(network)$name[i]
  land.polygon.ids <- land.polygon[ which( land.polygon$ROI == ids.i),]
  
  mapToPlot <- mapToPlot +
    geom_polygon(data = land.polygon.ids , fill = land.polygon.ids.c[i], colour = "black" , size=0.1 ,  aes(long, lat, group = group)) +
    geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent")
}

mapToPlot <- mapToPlot + coord_equal() + theme_map # coord_map(projection = "mercator")

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoordsMapped.pdf") , width = 10 )
mapToPlot
dev.off()

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

## Export major results

connectivityMatrix <- matrixConnectivityP.agg.square
diag(connectivityMatrix) <- 0

connectivityMatrixreclass <- connectivityMatrix
connectivityMatrixreclass[ connectivityMatrixreclass != 0] <- 1

isolated <- which(apply(connectivityMatrix,1,sum,na.rm=T) == 0 & apply(connectivityMatrix,2,sum,na.rm=T) == 0)
isolated.id <- colnames(connectivityMatrix)[isolated]
isolated.length <- length(isolated)

linkage.export <- apply(connectivityMatrixreclass,1,sum,na.rm=T)
linkage.import <- apply(connectivityMatrixreclass,2,sum,na.rm=T)

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=12) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) )

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
ggplot(data=data.frame(links=unique(linkage.export),sum=sapply(unique(linkage.export),function(x) { sum( linkage.export == x)  })), aes(x=links, y=sum)) +
  geom_bar(stat="identity", fill="#236292") + ylab("Number of sites") + xlab("Number of connections between sites [export]") + mainTheme
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
ggplot(data=data.frame(links=unique(linkage.import),sum=sapply(unique(linkage.import),function(x) { sum( linkage.import == x)  })), aes(x=links, y=sum)) +
  geom_bar(stat="identity", fill="#236292") + ylab("Number of sites") + xlab("Number of connections between sites [import]") + mainTheme
dev.off()

linkage <- apply(rbind( linkage.export , linkage.import),2,mean)
higherLinkage <- which(linkage >=  as.numeric(quantile(linkage,0.95,na.rm=TRUE))) 
linkage.on.average <- mean(linkage)

mean.probability <- mean( apply(rbind( apply(connectivityMatrix,1,mean,na.rm=T) , apply(connectivityMatrix,2,mean,na.rm=T)),2,mean) )

numberClusters <- length(unique(network.membership)) - length(isolated)
aggregationBasedOnClusters <- 1 - ( numberClusters / length(network.membership) )

resistance <- sapply( 1:nrow(connectivityMatrixreclass) , function(res) { 1- length(unique(c(which(connectivityMatrixreclass[res,] == 1 ) ,  which(connectivityMatrixreclass[,res] == 1 )))) / nrow(connectivityMatrixreclass) } )
higherResistance <- which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=TRUE)))

betweennessIndex <- betweenness(network)
higherBetweenness <- which(betweennessIndex >=  as.numeric(quantile(betweennessIndex,0.95,na.rm=TRUE)))
average.betweenness <- mean(betweennessIndex)

eigencentralityIndex <- eigen_centrality(network, directed = FALSE, scale = FALSE)$vector
higherEigencentrality <- which(eigencentralityIndex >=  as.numeric(quantile(eigencentralityIndex,0.95,na.rm=TRUE)))
average.eigencentrality <- mean(eigencentralityIndex)

closenessIndex <- closeness(network, mode="all")
higherCloseness <- which(closenessIndex >=  as.numeric(quantile(closenessIndex,0.95,na.rm=TRUE)))
average.closeness <- mean(closenessIndex)

# ---------------------

if( exists("pipeLiner") ) {
  
  c <- which(n.repetitions == n.days)
  
  if( ! exists("isolatedAll") ) { 
    summaryAgg <- data.frame()
    isolatedAgg<- data.frame( matrix(0,nrow=nrow(connectivityMatrix),ncol=length(n.repetitions)), row.names=colnames(connectivityMatrix))
    colnames(isolatedAgg) <- n.repetitions
    higherresistanceAgg <- resistanceAgg <- clusterAssignmentAgg <- higherclosenessAgg <- closenessAgg <- highereighenCentralityAgg <- eighenCentralityAgg <- higherBetweennessAgg <- betweennessAgg <- higherLinkageLevelAgg <- linkageLevelAgg <- isolatedAgg
  }
  
  isolatedAgg[ isolated ,c] <- 1
  write.csv(isolatedAgg,file=paste0(connectivityExportDir,"isolatedAgg.csv"))
  linkageLevelAgg[,c] <- linkage
  write.csv(linkageLevelAgg,file=paste0(connectivityExportDir,"linkageLevelAgg.csv"))
  higherLinkageLevelAgg[higherLinkage,c] <- 1
  write.csv(higherLinkageLevelAgg,file=paste0(connectivityExportDir,"higherLinkageLevelAgg.csv"))
  clusterAssignmentAgg[,c] <- network.membership
  write.csv(clusterAssignmentAgg,file=paste0(connectivityExportDir,"clusterAssignmentAgg.csv"))
  betweennessAgg[,c] <- betweennessIndex
  write.csv(betweennessAgg,file=paste0(connectivityExportDir,"betweennessAgg.csv"))
  higherBetweennessAgg[higherBetweenness,c] <- 1
  write.csv(higherBetweennessAgg,file=paste0(connectivityExportDir,"higherBetweennessAgg.csv"))
  eighenCentralityAgg[,c] <- eigencentralityIndex
  write.csv(eighenCentralityAgg,file=paste0(connectivityExportDir,"eighenCentralityAgg.csv"))
  highereighenCentralityAgg[higherEigencentrality,c] <- 1
  write.csv(highereighenCentralityAgg,file=paste0(connectivityExportDir,"highereighenCentralityAgg.csv"))
  closenessAgg[,c] <- closenessIndex
  write.csv(closenessAgg,file=paste0(connectivityExportDir,"closenessAgg.csv"))
  higherclosenessAgg[higherCloseness,c] <- 1
  write.csv(higherclosenessAgg,file=paste0(connectivityExportDir,"higherclosenessAgg.csv"))
  resistanceAgg[,c] <- resistance
  write.csv(resistanceAgg,file=paste0(connectivityExportDir,"resistanceAgg.csv"))
  higherresistanceAgg[higherResistance,c] <- 1
  write.csv(higherresistanceAgg,file=paste0(connectivityExportDir,"higherresistanceAgg.csv"))
  
}

if( ! all.equal(rownames(isolatedAgg),land.polygon$ROI)) { stop("Error :: 0003") }

land.polygon$Isolated <- isolatedAgg[,c]
land.polygon$Links <- linkage
land.polygon$HigherLinks <- higherLinkageLevelAgg[,c]
land.polygon$Betweenness <- betweennessAgg[,c]
land.polygon$HigherBetweenness <- higherBetweennessAgg[,c]
land.polygon$EighenCentrality <- eighenCentralityAgg[,c]
land.polygon$HighereighenCentrality <- higherLinkageLevelAgg[,c]
land.polygon$Closeness <- closenessAgg[,c]
land.polygon$HigherCloseness <- higherclosenessAgg[,c]
land.polygon$Resistance <- resistanceAgg[,c]
land.polygon$HigherResistance <- higherresistanceAgg[,c]
writeOGR(obj=land.polygon, dsn=connectivityExportDir, layer="sourceSinkPolygonsData", driver="ESRI Shapefile") 

summaryAgg <- rbind(summaryAgg,data.frame( day.sim=n.days,
                                           numberClusters=numberClusters,
                                           aggregation=aggregationBasedOnClusters,
                                           n.isolated=isolated.length,
                                           links.on.average=linkage.on.average,
                                           probability.on.average=mean.probability,
                                           betweenness.on.average=average.betweenness,
                                           eigencentrality.on.average=average.eigencentrality,
                                           closeness.on.average=average.closeness,
                                           modularity=network.modularity ))


write.csv(summaryAll,paste0(project.folder,"/Results/",project.name,"/connectivityExport/summaryAll.csv"))

# ---------------------

land.polygon.sf <- st_as_sf(land.polygon)

cols.to.use <- distinctColors(length(unique(additional.source.sink.sf$Isolated)))
cols.to.use <- c("white","red")



pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
mapNet +
  geom_sf(data=land.polygon.sf,aes(fill = as.factor(Isolated),color = as.factor(Isolated)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(land.polygon.sf$Isolated)), values = unique(cols.to.use)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") + theme_map + theme(legend.position = "none")
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
mapNet +
  geom_sf(data=land.polygon.sf,aes(fill = as.factor(HigherLinks),color = as.factor(HigherLinks)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(land.polygon.sf$HigherLinks)), values = unique(cols.to.use)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") + theme_map + theme(legend.position = "none")
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
mapNet +
  geom_sf(data=land.polygon.sf,aes(fill = as.factor(HigherBetweenness),color = as.factor(HigherBetweenness)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(land.polygon.sf$HigherBetweenness)), values = unique(cols.to.use)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") + theme_map + theme(legend.position = "none")
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
mapNet +
  geom_sf(data=land.polygon.sf,aes(fill = as.factor(HighereighenCentrality),color = as.factor(HighereighenCentrality)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(land.polygon.sf$HighereighenCentrality)), values = unique(cols.to.use)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") + theme_map + theme(legend.position = "none")
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
mapNet +
  geom_sf(data=land.polygon.sf,aes(fill = as.factor(HigherCloseness),color = as.factor(HigherCloseness)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(land.polygon.sf$HigherCloseness)), values = unique(cols.to.use)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") + theme_map + theme(legend.position = "none")
dev.off()

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
mapNet +
  geom_sf(data=land.polygon.sf,aes(fill = as.factor(HigherResistance),color = as.factor(HigherResistance)),size=0.5) + theme(legend.position = "none") + scale_fill_manual(limits = as.factor(unique(land.polygon.sf$HigherResistance)), values = unique(cols.to.use)) +
  geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent") + theme_map + theme(legend.position = "none")
dev.off()

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------