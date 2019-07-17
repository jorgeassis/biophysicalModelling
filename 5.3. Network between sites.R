## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

number.cores <- 6

distance.probability <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.Distance.bm"))
distance.probability <- data.table(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

source.sink.xy <- read.big.matrix(paste0(project.folder,"/Results/source.sink.bm"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

clipper <- as(extent(min(source.sink.xy[,2]) - 2,max(source.sink.xy[,2]) + 2,min(source.sink.xy[,3]) - 2,max(source.sink.xy[,3]) + 2), "SpatialPolygons")

cost.surface <- raster("Data/Rasters/Mask.tif")
cost.surface <- crop(cost.surface,clipper)
cost.surface[is.na(cost.surface)] <- 0
plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

## ------------------------------------------------------------------------------------------------------------------------------
## Stepping stone network between sampling sites

file.sampling.sites<- "/Volumes/Jellyfish/Dropbox/Manuscripts/_ Under Revision/Rejection of the Abundant Centre Hypothesis in marine mussels/Data/Sampling_sites.txt"
sampling.sites <- read.table(file.sampling.sites,sep=",",header=T)
sampling.sites.names <- as.character(sampling.sites$Site)
sampling.sites <- sapply(1:nrow(sampling.sites), function(x) { which.min(spDistsN1(as.matrix(source.sink.xy[,.(Lon,Lat)]),as.matrix(sampling.sites[x,2:3],ncol=2),longlat = TRUE)) } )

new.extent <- c(min(source.sink.xy[,2]),max(source.sink.xy[,2]),min(source.sink.xy[,3]),max(source.sink.xy[,3]))
network <- produce.network("Prob",as.data.frame(distance.probability),30,TRUE,5,as.data.frame(source.sink.xy),new.extent)
network.x <- network[[2]]
connectivity.x <- network[[1]]

network.sites <- expand.grid(from=sampling.sites,to=sampling.sites)
network.sites$probability <- 0

network.sites.direct <- expand.grid(from=sampling.sites,to=sampling.sites)
network.sites.direct$probability <- 0

for(i in 1:length(sampling.sites)) {
  for(j in 1:length(sampling.sites)) {
    
    possible.paths <- get.shortest.paths(network.x,as.character( sampling.sites[i] ) , as.character( sampling.sites[j] ),mode="out")$vpath
    stones.t <- as.numeric(names(possible.paths[[1]]))
    stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
    path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
    
    if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
    if( length(path.values) == 0) { path.values <- 0 }
    if( i == j ) { path.values <- 1 }
    
    network.sites[network.sites[,1] == sampling.sites[i] & network.sites[,2] == sampling.sites[j],3] <- path.values
    network.sites.direct[network.sites[,1] == sampling.sites[i] & network.sites[,2] == sampling.sites[j],3] <- connectivity.x[ connectivity.x[,1] == sampling.sites[i] & connectivity.x[,2] == sampling.sites[j] , 3 ][1]
    
  }
}

comb <- network.sites[-which(network.sites[,1] == network.sites[,2]),]

comb.matrix <- acast(network.sites, from ~ to , value.var = "probability",fill=0)
colnames(comb.matrix) <- c(paste0("P",1:9),"M3")
rownames(comb.matrix) <- c(paste0("P",1:9),"M3")
comb.matrix <- comb.matrix[c(1,2,3,4,10,5,6,7,8,9),c(1,2,3,4,10,5,6,7,8,9)]

comb.matrix.direct <- acast(network.sites.direct, from ~ to , value.var = "probability",fill=0)
colnames(comb.matrix.direct) <- c(paste0("P",1:9),"M3")
rownames(comb.matrix.direct) <- c(paste0("P",1:9),"M3")
comb.matrix.direct <- comb.matrix.direct[c(1,2,3,4,10,5,6,7,8,9),c(1,2,3,4,10,5,6,7,8,9)]

# -------------

clustering.method <- "cluster_optimal" # Uni: cluster_spinglass cluster_optimal cluster_louvain fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )

E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
E(graph.obj)$weight = 1 / -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
 
#graph.obj <- set.vertex.attribute(graph.obj, "name", value=sampling.sites.names)
graph.obj <- set.vertex.attribute(graph.obj, "name", value=c(paste0("P",1:9),"M3"))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "min") # min / mean / max
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight == 0))
graph.obj <- simplify(graph.obj)

membership.graph <- get(clustering.method)(graph.obj)$membership
modularity(graph.obj, membership.graph)

# --------------------------------------

weights <- (((-log(comb[,3]) - (max(-log(comb[,3])))) * (-1)) / 100)
weights <- comb[,3]
weights.reclass <- weights / max(weights)
weights.reclass[weights >= quantile(weights)[1] & weights < quantile(weights)[2]] <- 0.5
weights.reclass[weights >= quantile(weights)[2] & weights < quantile(weights)[3]] <- 0.5
weights.reclass[weights >= quantile(weights)[3] & weights < quantile(weights)[4]] <- 0.5
weights.reclass[weights >= quantile(weights, prob=0.75)] <- 2

plot(get(clustering.method)(graph.obj),graph.obj,vertex.label.color="Black",vertex.label.family="Helvetica",edge.width=weights.reclass,edge.color="Black") 
