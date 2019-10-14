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
source("Dependences.R")

sql.file <- "../Results/SQL/MPASimulationResults.sql"
bigmatrix.file <- "../InternalProc/particles.reference.desc"
sorce.sink.cells.file <- "../Results/source.sink.bm"
mpa.shp.notake.filename <- "../Data/Shapefiles/noTake.shp" 

sql.project.name <- "noTakeMPA"
number.cores <- 8

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

## --------------------------------------------------------------------
## --------------------------------------------------------------------

load("../Results/allPLDIndex.Rdata")

x <- results$pld
y <- results$n.clusters # colnames(results)
y.lab <- "Number of clusters"
x.lab <- "Propagule duration (day)"

par(mar = c(4.5, 5.5, 4.5, 4.5))
plot(x,y,pch=20,col="#A6A6A6", ylab="",xlab=x.lab,axes=FALSE)
title(ylab=y.lab, line=4)
lines(bezierCurve(x,y,100)$x,bezierCurve(x,y,100)$y,type="l", lwd=1,axes=FALSE, lty=2)
axis(2,las=2,col="White",col.ticks="Black", cex.axis=0.9)
axis(1,las=0,col="White",col.ticks="Black", cex.axis=0.9)
box()

## --------------------------------------------------------------------
## --------------------------------------------------------------------

load("../Results/notakeMPA.Rdata") 

worldMap <- ne_countries(scale = 10, returnclass = "sp")
worldMap <- crop(worldMap,extent(c(min(source.sink.xy[,2]) - 2.5 , max(source.sink.xy[,2]) + 2.5 , min(source.sink.xy[,3]) - 2.5 ,max(source.sink.xy[,3]) + 2.5 )))

allMPA <- shapefile(mpa.shp.notake.filename)

## --------------------------------------------------------------------
## --------------------------------------------------------------------

load("../Results/source.sink.xy.Rdata")

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Produce connectivity for different spawning months and pld periods

# List results

list.dirs(path = paste0("../Results"), recursive = FALSE)

season <- "YearRound" # SeasonSummer SeasonWinter YearRound
pld.period <- 200 # 10 30 90 120 200
type <- "notakeMPA" # allMPA notakeMPA notakeMPAConnectness

## --------------------

if( season == "SeasonSummer" ) { spawn.p <- c(6,7,8,9) }
if( season == "SeasonWinter" ) { spawn.p <- c(11,12,1,2) }
if( season == "YearRound" ) { spawn.p <- 1:12 }

project.name <- paste0(season,"_Pld",pld.period)

## --------------------

load(paste0("../Results/",project.name,"/connectivity.source.sink.xy.Rdata"))
load(paste0("../Results/",project.name,"/connectivity.source.sink.notakeMPA.Rdata")) 

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------

connectivity <- get(paste0("connectivity.",type))
connectivity.matrix <- connectivity[ ,c(1,2,7)] 

connectivity.matrix <- acast(connectivity.matrix, Pair.from ~ Pair.to )
diag(connectivity.matrix) <- 0

isolated.mpa <- colnames(connectivity.matrix)[which(apply(connectivity.matrix,1,sum,na.rm=T) == 0 & apply(connectivity.matrix,2,sum,na.rm=T) == 0)]
isolated.mpa

# Isolated MPAs
length(isolated.mpa)
get(type)$name[which( colnames(connectivity.matrix) %in% isolated.mpa)]

plot(worldMap , col="#E8E8E8",border="#C9C9C9")
plot(allMPA , col="gray",border="gray",add=T)
plot(get(type)[which(colnames(connectivity.matrix) %in% isolated.mpa),] , col="black",border="black",add=T)

# Aggregation level (Proportion of non-isolated MPAs, at least one connection, in relation to the number of MPAs)
(nrow(connectivity.matrix)-length(isolated.mpa)) / nrow(connectivity.matrix)

# Aggregation level (Based on overaal connections)
connectivity.matrix.binomial <- connectivity.matrix
connectivity.matrix.binomial[connectivity.matrix.binomial != 0] <- 1
sum ( ( apply(connectivity.matrix.binomial,1,sum) + 1 ) / nrow(connectivity.matrix.binomial) ) / nrow(connectivity.matrix.binomial)

# Average connections between MPAs
connect.index <- data.frame(MPA=colnames(connectivity.matrix),exportTo=apply(connectivity.matrix,1,function(x){ sum(x != 0) } ) , importFrom=apply(connectivity.matrix,2,function(x){ sum(x != 0) } ))

mean(unlist(connect.index[,-1]))
sd(unlist(connect.index[,-1]))
min(apply(connect.index[,-1],1,sum))
max(apply(connect.index[,-1],1,max))

# Export

mean(unlist(connect.index[,2]))
sd(unlist(connect.index[,2]))
min(unlist(connect.index[,2]))
max(unlist(connect.index[,2]))

histData <- unlist(connect.index[,2])
hist( histData[histData != 0], breaks=10)

# import

mean(unlist(connect.index[,3]))
sd(unlist(connect.index[,3]))
min(unlist(connect.index[,3]))
max(unlist(connect.index[,3]))

histData <- unlist(connect.index[,3])
hist( histData[histData != 0], breaks=10)

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

comb <- get(paste0("connectivity.",type))
comb <- comb[ comb$Pair.from != comb$Pair.to ,]
comb <- as.data.frame( comb[ sort(comb[,"Probability"] , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Probability")] )

# -------------

clustering.method <- "leading.eigenvector.community" # CLUSTERS ?? Uni: leading.eigenvector.community fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
graph.obj <- simplify(graph.obj)

# -------------

clustering.graph <- get(clustering.method)(graph.obj,options=list(maxiter=1000000))
membership.graph <- clustering.graph$membership
modularity(graph.obj, membership.graph) # Measuring goodness of fit of clustering (e.g, leading.eigenvector.community)

# Number of clusters
length(unique(membership.graph))

# Aggregation factor based on clusters 

( length(membership.graph) / length(unique(membership.graph)) ) / length(membership.graph)

V(graph.obj)$color <- "#000000"

plot(clustering.graph,graph.obj,vertex.size=1) 
plot(clustering.graph,graph.obj,vertex.label=NA,vertex.size=2) 

subset.mpa.shp <- get(type)[ sapply( as.numeric(as_ids(V(graph.obj))) , function(x) { which(get(type)$ID %in% x ) }    ) ,]
subset.mpa.shp@data$COLOUR <- membership.graph
cols.to.use <- rainbow(length(unique(membership.graph)))[membership.graph]

plot(worldMap , col="#E8E8E8",border="#C9C9C9")
plot(allMPA , col="gray",border="gray",add=T)
plot(subset.mpa.shp, col=cols.to.use,border=cols.to.use,add=T)

# Plot connections

connected.pairs <- comb[comb$Probability > 0,]
connected.pairs[,3] <- log(connected.pairs[,3] )
connected.pairs[,3] <- (1/connected.pairs[,3]) *(-1) * 100
connected.pairs[,3] <- connected.pairs[,3] - min(connected.pairs[,3])
connected.pairs[,3] <- connected.pairs[,3] / max(connected.pairs[,3])

# connected.pairs <- connected.pairs[ connected.pairs[,3] != 0 ,1:3]

centroids <- as.data.frame(gCentroid(get(type),byid=TRUE),xy=T)

# 1000px figure

plot(worldMap , col="#E8E8E8",border="#C9C9C9")
plot(get(type) , col="black",border="black",add=T)

colfunc <- colorRampPalette(c("Gray", "#CC6633","#C40F0F"))

for( i in nrow(connected.pairs):1 ){
  
  strenght <- (connected.pairs[i,3] * 100) + 1 
  routes_sl <- gcIntermediate(centroids[ which(get(type)$ID == connected.pairs[i,1]),],
                              centroids[ which(get(type)$ID == connected.pairs[i,2]),],
                              n = 100, addStartEnd = TRUE, sp = TRUE)
  
  lines(  routes_sl , type="l" , col=colfunc(101)[strenght])
  
}

## ------------------------------------------

# Eigen centrality vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).
# Cloness centrality measures how many steps is required to access every other vertex from a given vertex
# Betweenness centrality is (roughly) defined by the number of geodesics (shortest paths) going through a vertex or an edge.

vertexIndex <- closeness(graph.obj)
vertexIndex <- eigen_centrality(graph.obj, directed = TRUE, scale = TRUE)$vector
vertexIndex <- betweenness(graph.obj)

# Those with higher / Percertil 95%
temporaryRes <- get(type)$name[which( vertexIndex >=  as.numeric(quantile(vertexIndex,0.95)))]
temporaryRes

write.csv(temporaryRes,file="../Results/temporaryRes.csv")
write.csv(temporaryRes,file="../Results/temporaryRes.csv")

plot(worldMap , col="#E8E8E8",border="#C9C9C9")
plot(allMPA, col="gray",border="gray",add=T)
subset.mpa.shp <- get(type)[as.numeric( which( vertexIndex >=  as.numeric(quantile(vertexIndex,0.95))  ) ) , ]
subset.mpa.shp@data$COLOUR <- 1
plot(subset.mpa.shp, col="black", border="black" , add=TRUE)

## ------------------------------------------

pairs.poly <- expand.grid(MPA.from=get(type)$ID,MPA.to=get(type)$ID)
pairs.poly <- pairs.poly[ pairs.poly$MPA.from != pairs.poly$MPA.to,]
pairs.poly$Number.steps <- numeric(nrow(pairs.poly))

stones.list <- list()

for(p in 1:nrow(pairs.poly)) { 
  
  stones <- get.shortest.paths(graph.obj,as.character( pairs.poly[ p,1] ) , as.character( pairs.poly[ p,2] ),mode="out")$vpath
  stones <- names(unlist(stones))
  
  if(length(stones) == 1) {
    pairs.poly[ p,3] <- 0
  }
  if(length(stones) > 1) {
    pairs.poly[ p,3] <- 0
    stones.list <- c(stones.list, list(stones[-c(1,length(stones))]) )
    pairs.poly[ p,3] <- length(stones) - 1 
  }
}

head(pairs.poly)

# # Proportion of MPA pairs disconnected
# sum(pairs.poly$Number.steps == 0) / nrow(pairs.poly)
# 
# # Proportion of MPA pairs connected
# sum(pairs.poly$Number.steps != 0) / nrow(pairs.poly)

# Most important MPAs allowing connectivity
stones <- unique(unlist(stones.list))
stones <- data.frame(stone=stones,times=sapply(stones , function(x) { sum(unlist(stones.list) == x  ) } ))
stones <- stones[sort(stones$times , index.return = T , decreasing=T)$ix,]

head(stones)

plot(worldMap , col="#E8E8E8",border="#C9C9C9")
plot(allMPA, col="gray",border="gray",add=T)

# Connected (again)

subset.mpa.shp <- get(type)[ get(type)$ID %in% stones$stone, ]
plot(subset.mpa.shp, col="black", border="black",add=TRUE)

quantileThreshold <- 0.95

# Those with higher importance 
get(type)$name[get(type)$ID %in% stones$stone[stones$times >= quantile(stones$times,quantileThreshold)]]

plot(worldMap , col="#E8E8E8",border="#C9C9C9")
plot(allMPA, col="gray",border="gray",add=T)
subset.mpa.shp <- get(type)[ get(type)$ID %in% stones$stone[stones$times >= quantile(stones$times,quantileThreshold)], ]
plot(subset.mpa.shp, col="black", border="black",add=TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------------------------