## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

sql.project.name <- "Atlantic"
number.cores <- 40

## ------------------------------------------------------------------------------------------------------------------

Connectivity <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.bm"))
Connectivity <- data.table(Connectivity[,])
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Time.max" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity

source.sink.xy <- read.big.matrix(paste0(project.folder,"/Results/source.sink.bm"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

# source.sink.xy <- source.sink.xy[Source == 1 & Pair %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## ------------------------------------------------------------------------------------------------------------------------------
## Marine Distances

cost.surface <- raster("Data/Rasters/Mask.tif")
# cost.surface <- aggregate(cost.surface, fact=2)

clipper <- as(extent(min(source.sink.xy[,2]) - 2,max(source.sink.xy[,2]) + 2,min(source.sink.xy[,3]) - 2,max(source.sink.xy[,3]) + 2), "SpatialPolygons")
plot(cost.surface)
lines(clipper)

cost.surface <- crop(cost.surface,clipper)
cost.surface[is.na(cost.surface)] <- 0
plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(source.sink.xy[Pair == 1167,2:3]) , as.matrix(source.sink.xy[Pair == 1,2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == 1167,2:3]) , as.matrix(source.sink.xy[Pair == 1,2:3]) )

# ----------------------------------

n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores) ; registerDoParallel(cl.2)

marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  x.to <- x.to[x.to != 0]
  partial.distances <- costDistance(raster_tr_corrected, as.matrix(source.sink.xy[ x , 2:3 ]) , as.matrix(source.sink.xy[ x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=rep(x,length(partial.distances)),Pair.to=x.to,Distance=c(partial.distances)/1000)
  return( partial.distances )
  
}

stopCluster(cl.2) ; rm(cl.2)

head(marine.distances)

# ----------------------------------

distance.probability <- merge(Connectivity, marine.distances, by=c("Pair.from","Pair.to"))
distance.probability <- as.big.matrix(as.matrix(distance.probability))
write.big.matrix(distance.probability, paste0(project.folder,"/Results/Connectivity.Distance.bm"))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

distance.probability <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.Distance.bm"))
distance.probability <- data.table(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

# ----------------------------------

extract.simulation.days <- 5

distance.probability.t <- distance.probability[Time.max <= extract.simulation.days,]
max(distance.probability.t$Time.max)

# Summary 10:8

ggplot(distance.probability.t , aes(x=Distance,y=Probability)) + 
  geom_point(alpha = 0.3) + 
  theme_bw(base_size = 14) + 
  labs(x = "Distance (km)" , y = "Mean probability of connectivity") +
  theme(panel.background = element_rect(colour = "black") )

# ----------------------------------

# Summary 1

summary.results <- data.frame( Max     = c( round(max(distance.probability$Distance),3) , round(max(distance.probability$Probability),3) , round(max(distance.probability$Time.mean),3) ) ,
                               Mean    = c( round(mean(distance.probability$Distance),3) , round(mean(distance.probability$Probability),3) , round(mean(distance.probability$Time.mean),3) ) ,
                               SD      = c( round(sd(distance.probability$Distance),3) , round(sd(distance.probability$Probability),3) , round(sd(distance.probability$Time.mean),3) ) ,
                               Median  = c( round(median(distance.probability$Distance),3) , round(median(distance.probability$Probability),3) , round(median(distance.probability$Time.mean),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

qt  <- quantile(distance.probability$Probability, probs = 0.95)
distance.probability.t <- distance.probability[ Probability >= qt , ]

summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance),3) , round(max(distance.probability.t$Probability),3) , round(max(distance.probability.t$Time.mean),3) ) ,
                               Mean    = c( round(mean(distance.probability.t$Distance),3) , round(mean(distance.probability.t$Probability),3) , round(mean(distance.probability.t$Time.mean),3) ) ,
                               SD      = c( round(sd(distance.probability.t$Distance),3) , round(sd(distance.probability.t$Probability),3) , round(sd(distance.probability.t$Time.mean),3) ) ,
                               Median  = c( round(median(distance.probability.t$Distance),3) , round(median(distance.probability.t$Probability),3) , round(median(distance.probability.t$Time.mean),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

# Identify which have a high threshold

land.surface <- raster("Data/Rasters/Mask.tif")

cells.i <- distance.probability[ Time.max <= 5 & Distance >= 200 & Distance < 1000000000, Pair.from   ]
cells.j <- distance.probability[ Time.max <= 5 & Distance >= 200 & Distance < 1000000000, Pair.to  ]

clipper <- as(extent(min(source.sink.xy[unique(c(cells.i,cells.j)),2]) - 2 , max(source.sink.xy[unique(c(cells.i,cells.j)),2]) + 2,min(source.sink.xy[unique(c(cells.i,cells.j)),3]) - 2,max(source.sink.xy[unique(c(cells.i,cells.j)),3]) + 2), "SpatialPolygons")
land.surface <- crop(land.surface,clipper)

plot(land.surface,box=FALSE,legend=FALSE,col=c("black","Gray"))
points(source.sink.xy[cells.i,2:3],col="red")
points(source.sink.xy[cells.j,2:3],col="green")

distance.probability[ Pair.from %in% cells.i , ]
reference.table[ cell %in% cells.i & cell.rafted %in% cells.j , ]

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Pairwise Connectivity estimates

file.sampling.sites <- paste0(project.folder,"/Data/Differentiation/ID#0_Coords.txt")
file.differentiation <- paste0(project.folder,"/Data/Differentiation/ID#0_FST.txt")

transform.fst <- TRUE

if( transform.fst ) { project.name <- paste0(project.name,".FstT") }

sampling.sites <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,2:3] 
sampling.sites.n <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,1] 
differentiation <- read.table(file.differentiation,header = T,sep=";",stringsAsFactors=F)[,-1]
differentiation[lower.tri(differentiation)] <- t(differentiation)[lower.tri(differentiation)]

differentiation <- read.table(file.differentiation,header = T,sep=";",stringsAsFactors=F)[,-1]

if( nrow(differentiation) != length(sampling.sites.n) | sum(is.na(differentiation)) != length(sampling.sites.n) ) { next }

## ---------------

subseter <- which(sampling.sites$Lon >= min(source.sink.xy[,2]) & sampling.sites$Lon <= max(source.sink.xy[,2]),
                  sampling.sites$Lat >= min(source.sink.xy[,3]) & sampling.sites$Lat <= max(source.sink.xy[,3]))

subseter <- which(sampling.sites$Lat  < 50)

sampling.sites <- sampling.sites[subseter,]
sampling.sites.n <- sampling.sites.n[subseter]
differentiation <- differentiation[subseter,subseter]

## ---------------

position.matrix <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- source.sink.xy[position.matrix,1]

if( length(position.matrix) != length(unique(position.matrix)) ) { 
  
  unique.sites <- unique(position.matrix)
  differentiation.i <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  
  for( u in 1:length(unique.sites)) {
    for( v in 1:length(unique.sites)) {
      
      pairs.u <- which( unique.sites[u] == position.matrix )
      pairs.v <- which( unique.sites[v] == position.matrix )
      differentiation.i[u,v] <- mean(unlist(differentiation[pairs.u,pairs.v]),na.rm=T)
      
    }  }
  
  diag(differentiation.i) <- NA
  
  differentiation <- differentiation.i
  sampling.sites <- unique(sampling.sites)
  position.matrix <- unique(position.matrix)
  
}

position.matrix <- as.vector(unlist(position.matrix))

par(mfrow=c(1,1),mar = c(1,1,1,1))
plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
points(sampling.sites)
setkey(source.sink.xy,Pair)
points(source.sink.xy[Pair %in% position.matrix,2:3],col="red")

## ---------------

n.days <- 30

## ---------------------------------------------------

new.extent <- c(min(source.sink.xy[position.matrix,2]),max(source.sink.xy[position.matrix,2]),min(source.sink.xy[position.matrix,3]),max(source.sink.xy[position.matrix,3]))
network <- produce.network("Prob",distance.probability,n.days,TRUE,5,source.sink.xy,new.extent)

## ---------------------------------------------------

cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)

potential.connectivity <- foreach(from=position.matrix, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  network.x <- network[[2]]
  connectivity.x <- network[[1]]
  res.connectivity.to <- numeric(0)
  res.distance <- numeric(0)
  
  for( to in position.matrix ) {

    possible.paths.y <- get.shortest.paths(network.x,as.character( from ) , as.character( to ),mode="out")$vpath
    stones.t <- as.numeric(names(possible.paths.y[[1]]))
    stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
    path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
    
    if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
    if( length(path.values) == 0) { path.values <- 0 }
    if( from == to ) { path.values <- 1 }
    
    res.connectivity.to <- c(res.connectivity.to,path.values)
    res.distance <- c(res.distance, costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == from,2:3]) , as.matrix(source.sink.xy[Pair == to,2:3]) ))
    
  }
    
  ## ---------------------
  
  differentiation.to <- unlist(differentiation[ which( position.matrix == from), ])
  differentiation.to[is.na(differentiation.to)] <- 0
  
  temp.res <- data.frame( pair.from = from,
                          pair.to = position.matrix , 
                          differentiation = differentiation.to, 
                          distance = res.distance, 
                          probability.ss = res.connectivity.to )
                          
  return(temp.res)
          
  ## ---------------------
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

## ---------------------------------------------------
  
connectivity <- do.call(rbind,potential.connectivity)
connectivity <- connectivity[connectivity[,1] != connectivity[,2] ,]
connectivity[connectivity[,5] == 0,5] <- 1e-299

## ---------------

connectivity[,5] <- log(connectivity[,5])
connectivity[ connectivity < -9e10 ] <- NA
connectivity[ connectivity > 9e10 ] <- NA

norm <- t(combn(position.matrix, 2))
norm <- norm[norm[,1] != norm[,2] ,]

connectivity.final <- data.frame()

for( i in 1:nrow(norm)) {
  
  t.1 <- which( connectivity[,1] == norm[i,1] & connectivity[,2] == norm[i,2] )
  t.2 <- which( connectivity[,1] == norm[i,2] & connectivity[,2] == norm[i,1] )
  
  connectivity.final <- rbind(connectivity.final,
                              data.frame(from=norm[i,1],
                                         to=norm[i,2],
                                         Differantiation = as.numeric(as.character(connectivity$differentiation[t.1])),
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
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Distance),to=max(connectivity.final$Distance),length.out=100),predict(fit.ibd, data.frame(Distance=seq(from=min(connectivity.final$Distance),to=max(connectivity.final$Distance),length.out=100) )),lty=2,col="#902828")

plot(connectivity.final$Connectivity.min,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="-log(Min. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Connectivity.min,na.rm=T),to=max(connectivity.final$Connectivity.min,na.rm=T),length.out=100),predict(fit.min, data.frame(Connectivity.min=seq(from=min(connectivity.final$Connectivity.min,na.rm=T),to=max(connectivity.final$Connectivity.min,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(connectivity.final$Connectivity.mean,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="-log(Mean. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(connectivity.final$Connectivity.mean,na.rm=T),to=max(connectivity.final$Connectivity.mean,na.rm=T),length.out=100),predict(fit.mean, data.frame(Connectivity.mean=seq(from=min(connectivity.final$Connectivity.mean,na.rm=T),to=max(connectivity.final$Connectivity.mean,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(connectivity.final$Connectivity.max,connectivity.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="-log(Max. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
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

plot(connectivity.final$Differantiation,predict(fit.mix),lty=1,col="#5E5E5E",ylab="",xlab="Observed genetic differentiation",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Predicted genetic differentiation",mgp=c(4,1,0)) 
fit.mix.line <- lm(Pred ~ Differantiation, data=data.frame(Differantiation=connectivity.final$Differantiation,Pred=predict(fit.mix)), na.action = na.omit)
lines(seq(min(connectivity.final$Differantiation),max(connectivity.final$Differantiation),length.out = 100),predict(fit.mix.line, data.frame(Differantiation=seq(min(connectivity.final$Differantiation),max(connectivity.final$Differantiation),length.out = 100))) ,lty=2,col="#902828")

