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
## Produce network to estimate least cost distances

networkData <- matrixConnectivityP
networkData <- as.data.frame( networkData[ sort( as.vector(unlist(networkData[,"Probability"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )

network <- graph.edgelist( cbind( as.character( networkData[,1]) , as.character(networkData) ) , directed = TRUE )
E(network)$weight <- networkData[,3]
network <- delete.edges(network, which(E(network)$weight ==0))
network <- as.undirected(network, mode = "collapse", edge.attr.comb = "min") # min / mean / max
network <- simplify(network)
network.clustering <- cluster_leading_eigen(network,options=list(maxiter=1000000))
network.membership <- network.clustering$membership

cols.to.use <- distinctColors(length(unique(membership.graph)))
cols.to.use <- cols.to.use[membership.graph]
V(graph.obj)$color <- cols.to.use

lineThinkness <- reclassVals(E(graph.obj)$weight,min(E(graph.obj)$weight),max(E(graph.obj)$weight))

sorted.ids
source.sink.xy.ids

l <- source.sink.xy[sapply(get.vertex.attribute(graph.obj)$name,function(x) { which( as.character(sorted.ids$polygons.id) == x) } ),c("Lon","Lat")]
l <- as.matrix(l)

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
plot(graph.obj,edge.width=lineThinkness,vertex.label.dist=2,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

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

write.table(file=paste0(matrixConnectivityP.agg,"matrixConnectivityP.agg.csv"),x=matrixConnectivityP,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(matrixConnectivityPSS.agg,"matrixConnectivityPSS.agg.csv"),x=matrixConnectivityP,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(matrixConnectivityT.agg,"matrixConnectivityT.agg.csv"),x=matrixConnectivityP,row.names=TRUE,col.names=TRUE,sep=";",dec=".")

matrixConnectivityP.agg.square <- acast(matrixConnectivityP.agg[,c("from","to","val")], formula=from~to, fun.aggregate= mean , value.var="val",drop = FALSE)
matrixConnectivityPSS.agg.square <- acast(matrixConnectivityPSS.agg[,c("from","to","val")], formula=from~to, fun.aggregate= mean , value.var="val",drop = FALSE)
matrixConnectivityT.agg.square <- acast(matrixConnectivityT.agg[,c("from","to","val")], formula=from~to, fun.aggregate= mean , value.var="val",drop = FALSE)

write.table(file=paste0(matrixConnectivityP.agg.square,"matrixConnectivityP.agg.square.csv"),x=matrixConnectivityP,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(matrixConnectivityPSS.agg.square,"matrixConnectivityPSS.agg.square.csv"),x=matrixConnectivityP,row.names=TRUE,col.names=TRUE,sep=";",dec=".")
write.table(file=paste0(matrixConnectivityT.agg.square,"matrixConnectivityT.agg.square.csv"),x=matrixConnectivityP,row.names=TRUE,col.names=TRUE,sep=";",dec=".")

# ----------

comb <- matrixConnectivityP.agg

graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight <- comb[,3]
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "min") # min / mean / max
graph.obj <- simplify(graph.obj)
clustering.graph <- cluster_leading_eigen(graph.obj,options=list(maxiter=1000000))
membership.graph <- clustering.graph$membership

cols.to.use <- distinctColors(length(unique(membership.graph)))
cols.to.use <- cols.to.use[membership.graph]
V(graph.obj)$color <- cols.to.use
l <- layout.fruchterman.reingold(graph.obj)
l <- layout_nicely(graph.obj)

lineThinkness <- reclassVals(E(graph.obj)$weight,min(E(graph.obj)$weight),max(E(graph.obj)$weight))
  
pdf( file=paste0(connectivityExportDir,"idsNetworkClustering.pdf") , width = 10 )
plot(graph.obj,edge.width=lineThinkness,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

# ----------

l <- coords.centroid.ids[sapply(get.vertex.attribute(graph.obj)$name,function(x) { which(coords.centroid.ids == x) } ),c("Lon","Lat")]
l <- as.matrix(l)

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoords.pdf") , width = 10 )
plot(graph.obj,edge.width=lineThinkness,vertex.label.dist=2,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

## -------------------------------------------------------------
## Map network assignments

comb <- comb[comb$val != 0,]
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
land.polygon.ids.c <- get.vertex.attribute(graph.obj)$color

mapToPlot <- mapNet

for( i in 1:nrow(coords.centroid.ids)) {
  
  land.polygon.ids <- land.polygon[ which(!is.na(unlist(over(land.polygon,coords.centroid.ids.pts[i,])))),]
  
  mapToPlot <- mapToPlot +
    geom_polygon(data = land.polygon.ids , fill = land.polygon.ids.c[i], colour = land.polygon.ids.c[i] , size=0.15 ,  aes(long, lat, group = group)) +
    geom_rect(data=boxes, mapping=aes(xmin=minlong, xmax=maxlong, ymin=minlat, ymax=maxlat), color="transparent", fill="transparent")
}

mapToPlot <- mapToPlot + coord_equal() + theme_map # coord_map(projection = "mercator")

pdf( file=paste0(connectivityExportDir,"idsNetworkClusteringCoordsMapped.pdf") , width = 10 )
mapToPlot
dev.off()

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

connectivity.matrix <- connectivity[ ,c(1,2,7)]
connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )
diag(connectivity.matrix) <- 0

write.csv(connectivity.matrix,paste0("../Results/",project.name,"/connectivitymatrix.csv"))

isolated.mpa <- which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)
isolated.mpa.id <- colnames(connectivity.matrix)[isolated.mpa]
isolated.mpa.names <- MPAnames[which( colnames(connectivity.matrix) %in% isolated.mpa.id)]
isolated.mpa.length <- length(isolated.mpa)

isolatedResults[ which(rownames(isolatedResults) %in% isolated.mpa.names ),c] <- "TRUE"

if( exists("pipeLiner")) {
  
  isolatedAll
  write.csv(isolatedAll,file="../Results/isolatedMPAs.csv")
  
  isolatedAll <- data.frame()
  betweennessAll <- data.frame()
  higherBetweennessAll <- data.frame()
  eighenCentralityAll <- data.frame()
  highereighenCentralityAll <- data.frame()
  closenessAll <- data.frame()
  higherclosenessAll <- data.frame()
  clusterAssignmentAll <- data.frame()
  resistanceAll <- data.frame()
  higherresistanceAll <- data.frame()
  
  
  
  
}


  
  ## ------------------------------------------------------------------------------
  ## ------------------------------------------------------------------------------
  

  
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