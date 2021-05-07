## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

if( exists("pipeLiner") ) {
  
  rm( list=(ls()[ls()!="v"]) )
  gc(reset=TRUE)

  source("0. Project Config.R")
  source("Dependences.R")
  
  list.dirs(path = paste0("../Results"), recursive = FALSE)
  season <- "YearRound" # c("YearRound","SeasonSummer","SeasonWinter")
  spawn.p <- 1:12  # spawn.p <- c(6,7,8,9)
  pld.period <- 30
  c <- 1
  
}

worldMap <- ne_countries(scale = 10, returnclass = "sp")

bigmatrix.file <- paste0(project.folder,"/Results/",project.name,"/InternalProc/","particles.reference.desc")
sorce.sink.cells.file <- paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm")

# Open connectivity

n.season <- "" # Spring; Summer; Autumn; Winter; "" for All

Connectivity <- read.big.matrix(paste0(project.folder,"/Results/",project.name,"/InternalProc/","connectivityEstimatesAveraged",n.season,".bm"))
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Time.max" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity

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
crs(source.sink.xy.sp) <- dt.projection
source.sink.xy.sp$Pair <- source.sink.xy$Pair
  
# Add missing connections

missingFrom <- source.sink.xy$Pair[which(! source.sink.xy$Pair %in% unique(Connectivity$Pair.from))]

if(length(missingFrom) > 0) { 
  
  Connectivity <- rbind(Connectivity,data.frame(Pair.from=missingFrom,Pair.to=missingFrom),fill=TRUE)
  Connectivity[is.na(Connectivity)] <- 0
  
}

missingTo <- source.sink.xy$Pair[which(! source.sink.xy$Pair %in% unique(Connectivity$Pair.to))]

if(length(missingTo) > 0) { 
  
  Connectivity <- rbind(Connectivity,data.frame(Pair.from=missingTo,Pair.to=missingTo),fill=TRUE)
  Connectivity[is.na(Connectivity)] <- 0
  
}

## --------------------------------------------------------------------
## --------------------------------------------------------------------
# Temporary Subset

subseter <- as.vector(extent(source.sink.xy.sp) + c(-0.1,0.1,-0.1,0.1)) # c(-11.25,37.85,29.75,46.25)

## --------------------------------------------------------------------
## --------------------------------------------------------------------

source.sink.xy <- source.sink.xy[source.sink.xy$Lon >= subseter[1] & source.sink.xy$Lon <= subseter[2] & source.sink.xy$Lat >= subseter[3] & source.sink.xy$Lat <= subseter[4], ]
plot(source.sink.xy[,2:3])

source.sink.xy.sp <- crop(source.sink.xy.sp,extent(subseter))
worldMap <- crop(worldMap,extent(subseter + c(-10,10,-10,10))  ) 

## -----------------

plot(worldMap , col="Black",border="Black")
plot(source.sink.xy.sp , col="Red",border="Black", add=T)

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

combResults <- data.frame()
RegionNames <- as.character(source.sink.xy$Pair)

if( ! exists("isolatedResults")) {
  
  isolatedResults <- data.frame(matrix(nrow=length(RegionNames),ncol=nrow(combinations),""),stringsAsFactors = FALSE)
  rownames(isolatedResults) <- RegionNames
  colnames(isolatedResults) <- 1:nrow(combinations)
  betweennessResults <- higherBetweennessResults <- eighenCentralityResults <- highereighenCentralityResults <- closenessResults <- higherclosenessResults <- clusterAssignment <- resistanceResults <- higherResistanceResults <- outDegreeResults <- higherOutDegreeResults <- selfRecruitmentResults <- higherSelfRecruitmentResults <- isolatedResults

}

## ------------------------------------------------------------------------------

dev.off()
cat(c,"\n")
gc(reset=TRUE)

project.name.c <- paste0(project.name,"/",season,"_Pld",pld.period)

## ----------------------------------------------------

if( ! dir.exists(paste0("../Results/",project.name.c)) ) { dir.create(file.path(paste0("../Results/",project.name.c)), showWarnings = FALSE) } 
if( ! dir.exists(paste0("../Results/",project.name.c,"/Data")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Data")), showWarnings = FALSE) } 
if( ! dir.exists(paste0("../Results/",project.name.c,"/Maps")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Maps")), showWarnings = FALSE) } 
if( ! dir.exists(paste0("../Results/",project.name.c,"/Networks")) ) { dir.create(file.path(paste0("../Results/",project.name.c,"/Networks")), showWarnings = FALSE) } 

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

connectivity.matrix <- Connectivity[ ,c(1,2,3)]
connectivity.matrix[is.na(connectivity.matrix)] <- 0
connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )

write.csv(connectivity.matrix,paste0("../Results/",project.name.c,"/connectivitymatrix.csv"))

## ----------------------------------


retention <- diag(connectivity.matrix)
sumRows <- apply(connectivity.matrix,1,sum,na.rm=T)
retention[ is.na(retention) ] <- 0
selfRecruitment <- retention / sumRows
selfRecruitment[ is.na(selfRecruitment) ] <- 0

selfRecruitmentResults[ sapply(names(selfRecruitment),function(x) which( rownames(selfRecruitmentResults) == x)) ,c] <- selfRecruitment
write.csv(selfRecruitmentResults,file="../Results/selfRecruitment.csv")

higherSelfRecruitment.calc <- which(selfRecruitment >=  as.numeric(quantile(selfRecruitment,0.95,na.rm=TRUE)))
higherSelfRecruitmentResults[which(rownames(higherSelfRecruitmentResults) %in% names(higherSelfRecruitment.calc)),c] <- 1
write.csv(higherSelfRecruitmentResults,file="../Results/higherSelfRecruitment.csv")

diag(connectivity.matrix) <- 0

isolated.sourceSink <- which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)
isolatedResults[ sapply(names(isolated.sourceSink),function(x) which( rownames(isolatedResults) == x)),c] <- 1
write.csv(isolatedResults,file="../Results/isolated.csv")

# Aggregation level (Proportion of non-isolated MPAs, at least one connection, in relation to the number of MPAs)
aggregationAtLeastOne <- (nrow(connectivity.matrix)-length(isolated.sourceSink)) / nrow(connectivity.matrix)

# Aggregation level (Based on overall connections)
connectivity.matrix.binomial <- connectivity.matrix
connectivity.matrix.binomial[connectivity.matrix.binomial != 0] <- 1
aggregationAllConnections <- sum ( ( apply(connectivity.matrix.binomial,1,sum,na.rm=T) + 1 ) / nrow(connectivity.matrix.binomial) ) / nrow(connectivity.matrix.binomial)

# Resistance index
resistance <- sapply( 1:nrow(connectivity.matrix.binomial) , function(res) { 1- length(unique(c(which(connectivity.matrix.binomial[res,] == 1 ) ,  which(connectivity.matrix.binomial[,res] == 1 )))) / nrow(connectivity.matrix.binomial) } )
names(resistance) <- colnames(connectivity.matrix.binomial)

resistanceResults[ sapply(names(resistance),function(x) which( rownames(resistanceResults) == x)) ,c] <- resistance
write.csv(resistanceResults,file="../Results/resistance.csv")

# Those with higher / Percertil 95%
higherResistance.calc <- which(resistance >=  as.numeric(quantile(resistance,0.95,na.rm=TRUE)))
higherResistanceResults[which(rownames(higherResistanceResults) %in% names(higherResistance.calc)),c] <- 1
write.csv(higherResistanceResults,file="../Results/higherSelfRecruitment.csv")

# Average connections between MPAs
connect.index <- data.frame(SSSites=colnames(connectivity.matrix),exportTo=apply(connectivity.matrix,1,function(x){ sum(x != 0,na.rm=T) } ) , importFrom=apply(connectivity.matrix,2,function(x){ sum(x != 0,na.rm=T) } ))

averageConnections <- mean(unlist(connect.index[,-1]),na.rm=T)
sdConnections <- sd(unlist(connect.index[,-1]),na.rm=T)
maximumConnections <- max(apply(connect.index[,-1],1,max))

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

comb <- Connectivity
comb <- comb[ comb$Pair.from != comb$Pair.to ,]
comb <- as.data.frame( comb[ sort(comb$Probability , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Probability")] )

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

mapRegion <- ggplot() +
  geom_polygon(data = worldMap , fill = "#C4C4C4", colour = "#ffffff" , size=0.25 ,  aes(long, lat, group = group))  + coord_map() +
  coord_sf(xlim = c( extent(worldMap)[1], extent(worldMap)[2]), ylim = c(extent(worldMap)[3], extent(worldMap)[4]), expand = FALSE)
#  + coord_map('lambert', lat0=(extent(worldMap) + c(10,-10,10,-10))[3] + (extent(worldMap) + c(10,-10,10,-10))[4] / 2, lat1=45, xlim=c((extent(worldMap) + c(10,-10,10,-10))[1], (extent(worldMap) + c(10,-10,10,-10))[2]), ylim=c((extent(worldMap) + c(10,-10,10,-10))[3], (extent(worldMap) + c(10,-10,10,-10))[4])) + theme_map

mapRegionNet <- mapRegion
connected.pairs <- comb[comb$Probability > 0,]
colfunc <- colorRampPalette(c("#6C6C6C", "#CC6633","#C40F0F"))

for( i in nrow(connected.pairs):1 ){
  strenght <- (connected.pairs[i,3] * 100) + 1 
  routes_sl.1 <- which(source.sink.xy.sp$Pair == connected.pairs[i,1])
  routes_sl.2 <- which(source.sink.xy.sp$Pair == connected.pairs[i,2])
  routes_sl <- gcIntermediate(source.sink.xy.sp[routes_sl.1,],source.sink.xy.sp[routes_sl.2,],n = 100, addStartEnd = TRUE, sp = TRUE)
  SLDF = sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = NA), match.ID = F)
  mapRegionNet <- mapRegionNet + geom_path(data = SLDF, size=0.35 , aes(x = long, y = lat), col = colfunc(101)[round(strenght)] ) # "#797979"
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
for(i in 1:length(V(graph.obj.ss)$name)) {
  for(j in 1:length(V(graph.obj.ss)$name)) {
    connectivity.matrix.ss[i,j] <- distances(graph.obj.ss, V(graph.obj.ss)$name[i], V(graph.obj.ss)$name[j])
  }
}

connectivity.matrix.ss[connectivity.matrix.ss == Inf] <- max(connectivity.matrix.ss[connectivity.matrix.ss!=Inf],na.rm=T)
connectivity.matrix.ss[ is.na(connectivity.matrix.ss)] <- max(connectivity.matrix.ss[connectivity.matrix.ss!=Inf],na.rm=T)

#unLinked <- which(apply(connectivity.matrix.ss,1,sum) == nrow(connectivity.matrix.ss) - 1 & apply(connectivity.matrix.ss,2,sum) == nrow(connectivity.matrix.ss) - 1)
#connectivity.matrix.ss <- connectivity.matrix.ss[-unLinked,-unLinked]

# Euclidean distance

dist <- as.dist(connectivity.matrix.ss) # dist[7] <- 1.1 # Canarias vs. Selvagens 

# Hierarchical Clustering with hclust
hc <- hclust(dist, method = "ward.D")

# Plot the result
plot(hc)

# -------------

membership.graph <- clusters(graph.obj)$membership
clusterAssignment[ sapply(names(membership.graph),function(x) which(row.names(clusterAssignment) == x) ),c] <- membership.graph
write.csv(clusterAssignment,file="../Results/clusterAssignment.csv")

# Number of clusters
numberClusters <- length(unique(membership.graph)) - length(isolated.sourceSink)
aggregationBasedOnClusters <- 1 - ( numberClusters / length(membership.graph) )

## ------------------------------------------

# Plot clusters

cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]
isolated.sourceSink.col <- sapply(names(isolated.sourceSink),function(x) names( rownames(membership.graph) == x))
cols.to.use[ifelse(is.list(isolated.sourceSink.col),numeric(0),isolated.sourceSink.col)] <- "white"

l <- layout.fruchterman.reingold(graph.obj)
reducedNames <- RegionNames

pdf(file=paste0("../Results/",project.name.c,"/Networks/MapClusteringConnections.pdf"), width=12)
plot(graph.obj,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , color=cols.to.use , layout=l )
dev.off()

## ------------------------------------------

# Plot clusters with connections

cols.to.use <- distinctColors(length(unique(membership.graph)))[membership.graph]

centroids <- source.sink.xy.sp[sapply(names(membership.graph),function(x) which(source.sink.xy.sp$Pair == x)),]
centroids <- data.frame(centroids)

if(length(isolated.sourceSink) > 0) {
  centroidsIsolated <- source.sink.xy.sp[sapply(names(isolated.sourceSink),function(x) which(source.sink.xy.sp$Pair == x)),]
  centroidsIsolated <- data.frame(centroidsIsolated)
} else { centroidsIsolated <- centroids[1,]; centroidsIsolated$Lon <- -181; centroidsIsolated$Lat <- 91 }

pdf(file=paste0("../Results/",project.name.c,"/Maps/clusteringConnections.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroidsIsolated ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.9) +
    geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = cols.to.use, size = 2.5, stroke = 0.35, alpha = 0.7)
)
dev.off()

## ------------------------------------------

# Plot isolated with connections

pdf(file=paste0("../Results/",project.name.c,"/Maps/isolatedConn.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.7) +
    geom_point(data = centroidsIsolated ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "#9C2323", size = 2.5, stroke = 0.35, alpha = 0.9)
)
dev.off()

## ------------------------------------------

betweennessIndex <- betweenness(graph.obj)
betweennessIndex <- (betweennessIndex - min(betweennessIndex)) / (max(betweennessIndex) - min(betweennessIndex)) 
averageBetweenness <- mean(betweennessIndex)
sdBetweenness <- sd(betweennessIndex)

betweennessQ95 <- betweennessIndex
betweennessQ95[!is.na(betweennessQ95)] <- 0
betweennessQ95[which(betweennessIndex >= quantile(betweennessIndex,probs=0.95))] <- 1

betweennessResults[ sapply(names(betweennessIndex),function(x) which( rownames(betweennessResults) == x)) ,c] <- betweennessIndex
higherBetweennessResults[ sapply(names(betweennessIndex),function(x) which( rownames(betweennessResults) == x)),c] <- betweennessQ95
write.csv(betweennessResults,file="../Results/betweenness.csv")
write.csv(higherBetweennessResults,file="../Results/higherBetweenness.csv")

outDegree <- strength(graph.obj, mode = c("out"))
outDegree <- (outDegree - min(outDegree)) / (max(outDegree) - min(outDegree)) 
averageOutDegree <- mean(outDegree)
sdOutDegree <- sd(outDegree)

outDegreeQ95 <- outDegree
outDegreeQ95[!is.na(outDegreeQ95)] <- 0
outDegreeQ95[which(outDegree >= quantile(outDegree,probs=0.95))] <- 1

outDegreeResults[ sapply(names(outDegree),function(x) which( rownames(outDegreeResults) == x)),c] <- outDegree
higherOutDegreeResults[ sapply(names(outDegreeQ95),function(x) which( rownames(higherOutDegreeResults) == x)) ,c] <- outDegreeQ95

write.csv(outDegreeResults,file="../Results/outDegree.csv")
write.csv(higherOutDegreeResults,file="../Results/higherOutDegree.csv")

# Plot centrality indexes with Connections and clusters

betweennessIndexPlot <- betweennessIndex
betweennessIndexPlot <- (betweennessIndexPlot * 1.25) + 2

cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(betweennessIndexPlot),breaks = 20))]

centroids <- source.sink.xy.sp[sapply(names(betweennessIndexPlot),function(x) which(source.sink.xy.sp$Pair == x)),]
centroids <- data.frame(centroids)
centroidsHiger <- source.sink.xy.sp[sapply(names(betweennessQ95[betweennessQ95 == 1 ]),function(x) which(source.sink.xy.sp$Pair == x)),]
centroidsHiger <- data.frame(centroidsHiger)

pdf(file=paste0("../Results/",project.name.c,"/Maps/betweennessConnections.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroidsIsolated ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.9) +
    geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = cols.to.use, size = betweennessIndexPlot, stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroidsHiger ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "#9C2323", size = max(betweennessIndexPlot), stroke = 1.2, alpha = 0.7)
)
dev.off()
write.csv(data.frame(Q95=paste0("95th of betweenness centrality: ", round(quantile(betweennessIndexPlot,probs=0.95),2))),file=paste0("../Results/",project.name.c,"/Maps/betweennessConnections.csv"))

# Plot centrality indexes with Connections and clusters

outDegreeIndexPlot <- outDegree
outDegreeIndexPlot <- (outDegreeIndexPlot * 1.25) + 2

cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(outDegreeIndexPlot),breaks = 20))]

centroids <- source.sink.xy.sp[sapply(names(outDegreeIndexPlot),function(x) which(source.sink.xy.sp$Pair == x)),]
centroids <- data.frame(centroids)
centroidsHiger <- source.sink.xy.sp[sapply(names(outDegreeQ95[outDegreeQ95 == 1 ]),function(x) which(source.sink.xy.sp$Pair == x)),]
centroidsHiger <- data.frame(centroidsHiger)

pdf(file=paste0("../Results/",project.name.c,"/Maps/outDegreeConnections.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroidsIsolated ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.9) +
    geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = cols.to.use, size = outDegreeIndexPlot, stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroidsHiger ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "#9C2323", size = max(outDegreeIndexPlot), stroke = 1.2, alpha = 0.7)
)
dev.off()
write.csv(data.frame(Q95=paste0("95th of out degree: ", round(quantile(outDegreeIndexPlot,probs=0.95),2))),file=paste0("../Results/",project.name.c,"/Maps/outDegreeConnections.csv"))

## ------------------------------------------

# Plot centrality indexes with Connections and clusters

selfRecruitmentPlot <- selfRecruitment
selfRecruitmentPlot <- (selfRecruitmentPlot * 1.25) + 2

cols.to.use <- colorRampPalette(c('#BAE2FF','yellow','orange','#9C2323'))
cols.to.use <- cols.to.use(20)[as.numeric(cut(as.numeric(selfRecruitmentPlot),breaks = 20))]

centroids <- source.sink.xy.sp[sapply(names(selfRecruitmentPlot),function(x) which(source.sink.xy.sp$Pair == x)),]
centroids <- data.frame(centroids)
centroidsHiger <- source.sink.xy.sp[sapply(names(higherSelfRecruitment[higherSelfRecruitment == 1 ]),function(x) which(source.sink.xy.sp$Pair == x)),]
centroidsHiger <- data.frame(centroidsHiger)

pdf(file=paste0("../Results/",project.name.c,"/Maps/selfRecruitment.pdf"), width=12)
print(
  mapRegionNet + 
    geom_point(data = centroidsIsolated ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "white", size = 2.5, stroke = 0.35, alpha = 0.9) +
    geom_point(data = centroids ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = cols.to.use, size = selfRecruitmentPlot, stroke = 0.25, alpha = 0.7) +
    geom_point(data = centroidsHiger ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "#9C2323", size = max(selfRecruitmentPlot), stroke = 1.2, alpha = 0.7)
)
dev.off()
write.csv(data.frame(Q95=paste0("95th of self recruitment: ", round(quantile(selfRecruitmentPlot,probs=0.95),2))),file=paste0("../Results/",project.name.c,"/Maps/selfRecruitment.csv"))

# ----------------------------------
# ----------------------------------

combResults <- rbind(combResults,
                     data.frame(pld=pld.period,
                                n.isolated.sourceSink=length(isolated.sourceSink),
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
                                sdOutDegree=sdOutDegree
                            ))

write.csv(combResults,file="../Results/Results.csv")
save(combResults,file=paste0("../Results/allPLDResults.Rdata"))

# list.memory()
rm(connectivity.source.sink.xy )

## ---------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------

if( exists("pipeLiner") ) {
  
  names(combResults)
  
  x <- combResults$pld
  x.lab <- "Propagule duration (day)"
  y <- combResults$n.isolated.sourceSink
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
  ggplot(combResults, aes(x = pld, y = n.isolated.sourceSink)) +
    geom_line() + 
    ylab("Isolation degree (number of reserves)") +
    xlab("Propagule duration (day)")
  dev.off()

}

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------
