## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

if( ! exists("pipeLiner") ) {
  
  rm( list=(ls()[ls()!="v"]) )
  gc(reset=TRUE)
  
  source("0. Project Config.R")
  source("Dependences.R")
  
  list.dirs(path = paste0("../Results"), recursive = FALSE)
  pld.period <- 30
  
  popCoordinates <- read.csv("../Data/DiffRecords.csv", sep=";")
  popDifferentiation <- read.csv("../Data/DiffFST.csv", sep=";", header=FALSE)
  
}

## --------------

transformFST <- TRUE

worldMap <- ne_countries(scale = 10, returnclass = "sp")
worldMap <- crop(worldMap,extent(min(popCoordinates[,2])-10,max(popCoordinates[,2])+10,min(popCoordinates[,3])-10,max(popCoordinates[,3])+10))

Connectivity.file <- paste0(project.folder,"/Results/",project.name,"/InternalProc","/connectivityEstimatesAveragedDistances.bm")
Connectivity <- read.big.matrix(Connectivity.file)
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

sorce.sink.cells.file <- paste0(project.folder,"/Results/",project.name,"/InternalProc/","source.sink.bm")
source.sink.xy <- read.big.matrix(sorce.sink.cells.file)
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

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

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Pairwise Connectivity estimates

popCoordinates.n <- popCoordinates[,1]
popCoordinates <- popCoordinates[,c(2,3)]

if( nrow(popDifferentiation) != length(popCoordinates.n) ) { stop("Files do not match in extent") }

## ---------------

popCoordinatesSS <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(popCoordinates),longlat = TRUE)
popCoordinatesSS <- apply(popCoordinatesSS,2,which.min)
popCoordinatesSS <- source.sink.xy[popCoordinatesSS,1]
popCoordinatesSS <- unlist(popCoordinatesSS)

plot(popCoordinates)
points(source.sink.xy[Pair %in% popCoordinatesSS,2:3],col="red")

if( length(popCoordinatesSS) != length(unique(popCoordinatesSS)) ) { 
  
  unique.sites <- unique(popCoordinatesSS)
  popDifferentiation.i <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  
  for( u in 1:length(unique.sites)) {
    for( v in 1:length(unique.sites)) {
      
      pairs.u <- which( unique.sites[u] == popCoordinatesSS )
      pairs.v <- which( unique.sites[v] == popCoordinatesSS )
      popDifferentiation.i[u,v] <- mean(unlist(popDifferentiation[pairs.u,pairs.v]),na.rm=T)
      
    }  }
  
  diag(popDifferentiation.i) <- NA
  
  popDifferentiation <- popDifferentiation.i
  popCoordinates <- data.frame(Lon=unlist(sapply(unique.sites,function(x) data.frame(source.sink.xy[Pair == x,2]) )),Lat=unlist(sapply(unique.sites,function(x) data.frame(source.sink.xy[Pair == x,3]) )))
  
  popCoordinates.n <- popCoordinates.n[-which(duplicated(popCoordinatesSS))]
  popCoordinatesSS <- unique.sites
  
}

ggplot() +
  geom_polygon(data = worldMap , fill = "#C4C4C4", colour = "#ffffff" , size=0.25 ,  aes(long, lat, group = group))  + coord_map() +
  geom_point(data = popCoordinates  ,  aes(x = Lon, y = Lat) , shape = 21, colour = "black", fill = "Red", size = 2.5, stroke = 0.35, alpha = 0.9) +
  geom_label(aes(label=popCoordinates.n, x=popCoordinates$Lon,y=popCoordinates$Lat),check_overlap = T) 
  # geom_text(aes(label=popCoordinates.n, x=popCoordinatesMap$Lon,y=popCoordinatesMap$Lat),check_overlap = T) 

## ---------------------------------------------------
## ---------------------------------------------------

clipper <- as(extent(min(source.sink.xy[,2]) - 5,max(source.sink.xy[,2]) + 5,min(source.sink.xy[,3]) - 5,max(source.sink.xy[,3]) + 5), "SpatialPolygons")

costSurface <- raster("Data/Rasters/Mask.tif")
costSurface <- crop(costSurface,clipper)
costSurface[is.na(costSurface)] <- 0
plot(costSurface,box=FALSE,legend=FALSE,col=c("black","white"))
raster_tr <- transition(costSurface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

## ---------------------------------------------------
## ---------------------------------------------------
## All Pairs

comb <- Connectivity
comb[ which(comb[,"Time.max"] > pld.period) , "Probability" ] <- 0
comb <- as.data.frame( comb[ sort(comb$Probability , decreasing = TRUE, index.return =TRUE)$ix , c("Pair.from","Pair.to","Probability")] )
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
# E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
# E(graph.obj)$weight = comb[,3] 
E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015

#graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "mean") # min / mean / max
graph.obj <- simplify(graph.obj)

## ---------------------------------------------------
## Sampled Pairs

graph.obj.sampled <- delete_vertices(graph.obj, names(V(graph.obj))[! names(V(graph.obj)) %in% as.character(popCoordinatesSS)])
membership.graph <- clusters(graph.obj.sampled)$membership

cols.to.use <- distinctColors(length(membership.graph))
cols.to.use <- cols.to.use[membership.graph]

l <- layout.fruchterman.reingold(graph.obj.sampled)
reducedNames <- unlist(sapply( names(V(graph.obj.sampled)),function(x) { popCoordinates.n[which( popCoordinatesSS == x)] } ))
plot(graph.obj.sampled,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,layout = l,edge.curved = F , vertex.color=cols.to.use )

pdf( file=paste0(project.folder,"Results/Clustering Sampled Sites.pdf") , width = 10, height = 8 )
plot(graph.obj.sampled,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,layout = l,edge.curved = F , vertex.color=cols.to.use )
dev.off()

## ---------------------------------------------------

cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)

potential.connectivity <- foreach(from=popCoordinatesSS, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  res.connectivity.to <- numeric(0)
  res.distance <- numeric(0)
  
  for( to in popCoordinatesSS ) {

    possible.paths.y <- get.shortest.paths(graph.obj,as.character( from ) , as.character( to ),mode="out")$vpath
    stones.t <- as.numeric(names(possible.paths.y[[1]]))
    stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
    path.values <- apply( stones.t.interm , 1 , function(z) { comb[ comb[,1] == z[1] & comb[,2] == z[2] , 3 ][1] }   )
    
    if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
    if( length(path.values) == 0) { path.values <- 0 }
    if( from == to ) { path.values <- 1 }
    
    res.connectivity.to <- c(res.connectivity.to,path.values)
    res.distance.t <- costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == from,2:3]) , as.matrix(source.sink.xy[Pair == to,2:3]) )
    
    if(res.distance.t == 0 | res.distance.t == Inf) { res.distance.t <- spDistsN1( as.matrix(source.sink.xy[ Pair == from , 2:3 ]), as.matrix(source.sink.xy[ Pair == to , 2:3 ]), longlat=TRUE) }
    
    res.distance <- c(res.distance, res.distance.t)
    
  }
    
  ## ---------------------
  
  popDifferentiation.to <- unlist(popDifferentiation[ which( popCoordinatesSS == from), ])
  popDifferentiation.to[is.na(popDifferentiation.to)] <- 0
  
  temp.res <- data.frame( pair.from = from,
                          pair.to = popCoordinatesSS , 
                          popDifferentiation = popDifferentiation.to, 
                          distance = res.distance, 
                          probability.ss = res.connectivity.to )
                          
  return(temp.res)
          
  ## ---------------------
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

## ---------------------------------------------------
  



# REVIEW FROM THIS POINT ONWARDS



connectivity <- do.call(rbind,potential.connectivity)
connectivity <- connectivity[connectivity[,1] != connectivity[,2] ,]
connectivity[connectivity[,5] == 0,5] <- 1e-299

## ---------------

connectivity[,5] <- log(connectivity[,5])
connectivity[ connectivity < -9e10 ] <- NA
connectivity[ connectivity > 9e10 ] <- NA

norm <- t(combn(popCoordinatesSS, 2))
norm <- norm[norm[,1] != norm[,2] ,]

connectivity.final <- data.frame()

for( i in 1:nrow(norm)) {
  
  t.1 <- which( connectivity[,1] == norm[i,1] & connectivity[,2] == norm[i,2] )
  t.2 <- which( connectivity[,1] == norm[i,2] & connectivity[,2] == norm[i,1] )
  
  connectivity.final <- rbind(connectivity.final,
                              data.frame(from=norm[i,1],
                                         to=norm[i,2],
                                         Differantiation = as.numeric(as.character(connectivity$popDifferentiation[t.1])),
                                         Distance = connectivity$distance[t.1],
                                         Connectivity.max = min(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2])) ,
                                         Connectivity.mean = mean(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2])) ,
                                         Connectivity.min = max(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2]))
                              , stringsAsFactors = FALSE ) )
}

connectivity.final <- connectivity.final[connectivity.final[,1] != connectivity.final[,2] ,]

if(transform.fst) { connectivity.final[,3] <- connectivity.final$Differantiation/(1-connectivity.final$Differantiation) }

connectivity.final <- connectivity.final[which(complete.cases(connectivity.final)),]

cor.ibd <- cor(connectivity.final$Distance,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.ibd <- lm(Differantiation ~ Distance, data=connectivity.final , na.action = na.omit)

cor.min <- cor(connectivity.final$Connectivity.min,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.min <- lm(Differantiation ~ Connectivity.min, data=connectivity.final)

cor.max <- cor(connectivity.final$Connectivity.max,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.max <- lm(Differantiation ~ Connectivity.max, data=connectivity.final)

cor.mean <- cor(connectivity.final$Connectivity.mean,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.mean <- lm(Differantiation ~ Connectivity.mean, data=connectivity.final)

r2.ibd = summary(fit.ibd)$adj.r.squared
p.ibd=summary(fit.ibd)$coefficients
p.ibd= ifelse( nrow(p.ibd) == 2 , p.ibd[2,4] , NA)

r2.min = summary(fit.min)$adj.r.squared
p.min=summary(fit.min)$coefficients
p.min= ifelse( nrow(p.min) == 2 , p.min[2,4] , NA)

r2.mean = summary(fit.mean)$adj.r.squared
p.mean=summary(fit.mean)$coefficients
p.mean= ifelse( nrow(p.mean) == 2 , p.mean[2,4] , NA)

r2.max = summary(fit.max)$adj.r.squared
p.max=summary(fit.max)$coefficients
p.max= ifelse( nrow(p.max) == 2 , p.max[2,4] , NA)

data.frame(day = n.days,
          aic.ibd= AIC(fit.ibd),
          r2.ibd = r2.ibd,
          cor.ibd = cor.ibd,
          p.ibd=p.ibd,
          aic.min= AIC(fit.min),
          r2.min = r2.min,
          cor.min = cor.min,
          p.min=p.min,
          aic.mean= AIC(fit.mean),
          r2.mean = r2.mean,
          cor.mean = cor.mean,
          p.mean=p.mean,
          aic.max= AIC(fit.max),
          r2.max = r2.max,
          cor.max = cor.max,
          p.max=p.max ) 
      
## ------------------------------------
## ------------------------------------

pdf( file=paste0(project.folder,"Results/",project.name,".lm.fit.pdf") , width = 9, height = 9 )

par(mfrow=c(2,2))

plot(connectivity.final$Distance,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="Marine distance (km)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic popDifferentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Distance),to=max(connectivity.final$Distance),length.out=100),predict(fit.ibd, data.frame(Distance=seq(from=min(connectivity.final$Distance),to=max(connectivity.final$Distance),length.out=100) )),lty=2,col="#902828")

plot(connectivity.final$Connectivity.min,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="-log(Min. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic popDifferentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Connectivity.min,na.rm=T),to=max(connectivity.final$Connectivity.min,na.rm=T),length.out=100),predict(fit.min, data.frame(Connectivity.min=seq(from=min(connectivity.final$Connectivity.min,na.rm=T),to=max(connectivity.final$Connectivity.min,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(connectivity.final$Connectivity.mean,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="-log(Mean. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic popDifferentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Connectivity.mean,na.rm=T),to=max(connectivity.final$Connectivity.mean,na.rm=T),length.out=100),predict(fit.mean, data.frame(Connectivity.mean=seq(from=min(connectivity.final$Connectivity.mean,na.rm=T),to=max(connectivity.final$Connectivity.mean,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(connectivity.final$Connectivity.max,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="-log(Max. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic popDifferentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Connectivity.max,na.rm=T),to=max(connectivity.final$Connectivity.max,na.rm=T),length.out=100),predict(fit.max, data.frame(Connectivity.max=seq(from=min(connectivity.final$Connectivity.max,na.rm=T),to=max(connectivity.final$Connectivity.max,na.rm=T),length.out=100) )),lty=2,col="#902828")

dev.off()
      
## ------------------------------------------------------------------------------------------------------------------------------

# Plot IBD vs Best Ocean Model
# 8:8

fit.mix <- lm(Differantiation ~ Connectivity.mean+Distance, data=connectivity.final)
summary(fit.ibd)
summary(fit.mean)
summary(fit.mix)
AIC(fit.mix) ; AIC(fit.ibd) ; AIC(fit.mean)
anova(fit.ibd,fit.mean)
anova(fit.mix,fit.ibd)
anova(fit.mix,fit.mean)
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))

plot(connectivity.final$Differantiation,predict(fit.mix),lty=1,col="#5E5E5E",ylab="",xlab="Observed genetic popDifferentiation",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Predicted genetic popDifferentiation",mgp=c(4,1,0)) 
fit.mix.line <- lm(Pred ~ Differantiation, data=data.frame(Differantiation=connectivity.final$Differantiation,Pred=predict(fit.mix)), na.action = na.omit)
lines(seq(min(connectivity.final$Differantiation),max(connectivity.final$Differantiation),length.out = 100),predict(fit.mix.line, data.frame(Differantiation=seq(min(connectivity.final$Differantiation),max(connectivity.final$Differantiation),length.out = 100))) ,lty=2,col="#902828")

## ------------------------------------------------------------------------------------------------------------------------------

comb <- do.call(rbind,potential.connectivity)[,c(1,2,5)]
comb <- as.data.frame( comb[ sort( as.vector(unlist(comb[,"probability.ss"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )

net.function <<- prod
graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
E(graph.obj)$weight = comb[,3] # Hock, Karlo Mumby, Peter J 2015
graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))

network <- simplify(graph.obj,remove.multiple = TRUE, remove.loops = TRUE)

clustering.graph <- cluster_leading_eigen(network,options=list(maxiter=1000000))
# clustering.graph <- cluster_edge_betweenness(network)
clustering.graph <- cluster_walktrap(network)
clustering.graph <- clusters(network)
 
membership.graph <- clustering.graph$membership

distinctColors <- function(n) {
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(sample(col_vector, n))
}

cols.to.use <- distinctColors(length(unique(membership.graph)))
cols.to.use <- cols.to.use[membership.graph]
V(network)$color <- cols.to.use
# l <- layout.fruchterman.reingold(network)
reducedNames <- unlist(sapply( names(V(network)),function(x) { popCoordinates.n[which( popCoordinatesSS == x)] } ))
plot(network,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , color=cols.to.use , layout=layout.fruchterman.reingold(network) ) # layout = l,

pdf( file=paste0(project.folder,"Results/Clustering SS PD 60.pdf") , width = 10, height = 8 )
plot(network,vertex.label.dist=1.5,vertex.label.family="Helvetica",vertex.label.color="Black",vertex.label.cex=0.75,vertex.label=reducedNames,vertex.size=10,edge.curved = F , color=cols.to.use , layout=layout.fruchterman.reingold(network) ) # layout = l,
dev.off()
