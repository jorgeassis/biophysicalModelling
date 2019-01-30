## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

source("0. Project Config.R")

## ------------------------------------------------------------------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

## -------------------

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## ------------------------------------------------------------------------------------------------------------------------------
## Prob. vs Distance Plot
##

cost.surface <- raster("Data/Rasters/Mask.tif")
clipper <- as(extent(min(source.sink.xy[,2] - 2),max(source.sink.xy[,2] + 2),min(source.sink.xy[,3] - 2),max(source.sink.xy[,3] + 2)), "SpatialPolygons")
cost.surface <- crop(cost.surface,clipper)

cost.surface[is.na(cost.surface)] <- 0
plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(source.sink.xy[source.sink.xy$cells.id == 1167,2:3]) , as.matrix(source.sink.xy[source.sink.xy$cells.id == 1,2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(source.sink.xy[source.sink.xy$cells.id == 1167,2:3]), as.matrix(source.sink.xy[source.sink.xy$cells.id == 1,2:3]) )

# ----------------------------------

# REVISE!

number.cores.t <- 2
n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores.t) ; registerDoParallel(cl.2)
marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  partial.distances <- costDistance(raster_tr_corrected, as.matrix(cells[ x , 2:3 ]) , as.matrix(cells[ x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=x,Pair.to=x.to,distance=c(partial.distances)/1000)
  return( partial.distances )
  
}
stopCluster(cl.2) ; rm(cl.2)
head(marine.distances)

# Save object

save(marine.distances,file=paste0(project.folder,"/Data/marine.distances.RData"))
save(raster_tr_corrected,file=paste0(project.folder,"/Data/cost.distance.raster.RData"))

extract.simulation.days <- 10

distance.probability <- merge( marine.distances.DT, Connectivity.DT, by.x = c("Pair.from","Pair.to"), by.y = c("Pair.from","Pair.to")) 
colnames(distance.probability) <- c("Pair.from" , "Pair.to" ,"Distance" , "Mean.Probability" , "Max.Probability" , "Mean.Time" , "Max.Time" , "N.events" )

distance.probability <- distance.probability[Max.Time <= extract.simulation.days & Distance < 1500,]

max(distance.probability$Max.Time)

# Summary 0

ggplot(distance.probability , aes(x=Distance,y=Mean.Probability)) + 
  geom_point(alpha = 0.3) + 
  theme_bw(base_size = 14) + 
  labs(x = "Distance (km)" , y = "Mean probability of connectivity") +
  theme(panel.background = element_rect(colour = "black") )

# ----------------------------------

# Summary 1

summary.results <- data.frame( Max     = c( round(max(distance.probability$Distance),3) , round(max(distance.probability$Mean.Probability),3) , round(max(distance.probability$Mean.Time),3) ) ,
                               Mean    = c( round(mean(distance.probability$Distance),3) , round(mean(distance.probability$Mean.Probability),3) , round(mean(distance.probability$Mean.Time),3) ) ,
                               SD      = c( round(sd(distance.probability$Distance),3) , round(sd(distance.probability$Mean.Probability),3) , round(sd(distance.probability$Mean.Time),3) ) ,
                               Median  = c( round(median(distance.probability$Distance),3) , round(median(distance.probability$Mean.Probability),3) , round(median(distance.probability$Mean.Time),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

qt  <- quantile(distance.probability$Mean.Probability, probs = 0.95)
distance.probability.t <- distance.probability[ Mean.Probability >= qt , ]

summary.results <- data.frame( Max     = c( round(max(distance.probability.t$Distance),3) , round(max(distance.probability.t$Mean.Probability),3) , round(max(distance.probability.t$Mean.Time),3) ) ,
                               Mean    = c( round(mean(distance.probability.t$Distance),3) , round(mean(distance.probability.t$Mean.Probability),3) , round(mean(distance.probability.t$Mean.Time),3) ) ,
                               SD      = c( round(sd(distance.probability.t$Distance),3) , round(sd(distance.probability.t$Mean.Probability),3) , round(sd(distance.probability.t$Mean.Time),3) ) ,
                               Median  = c( round(median(distance.probability.t$Distance),3) , round(median(distance.probability.t$Mean.Probability),3) , round(median(distance.probability.t$Mean.Time),3) ) )
row.names(summary.results) <- c("Distance","Probability","Time")
summary.results

# Identify which have a high threshold

cells.i <- distance.probability[ Distance >= 1300, Pair.from  ]
cells.j <- distance.probability[ Distance >= 1300, Pair.to  ]

plot(ocean.region,box=FALSE,legend=FALSE,col=c("black"))
points(cells[cells.i,2:3],col="red")
points(cells[cells.j,2:3],col="green")

distance.probability[ Pair.from %in% cells.i & Pair.to %in% cells.j , ]
reference.table[ cell %in% cells.i & cell.rafted %in% cells.j , ]

## ------------------------------------------------------------------------------------------------------------------------------
## Pairwise Connectivity estimates

file.sampling.sites <- paste0(project.folder,"/Data/2_Coords.txt")
file.differentiation <- paste0(project.folder,"/Data/2_FST.txt")

sampling.sites <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,c(3,2)]
sampling.sites.n <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,1] 
differentiation <- read.table(file.differentiation,header = T,sep=";",stringsAsFactors=F)[,-1]
differentiation[lower.tri(differentiation)] <- t(differentiation)[lower.tri(differentiation)]
data.frame(colnames(differentiation.r) , sampling.sites.n)

## ---------------

position.matrix <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- source.sink.xy[position.matrix,1]

plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
points(sampling.sites)
points(source.sink.xy[position.matrix,2:3],col="red")

## ---------------

n.days <- 30
new.extent <- c(min(source.sink.xy[position.matrix,2]),max(source.sink.xy[position.matrix,2]),min(source.sink.xy[position.matrix,3]),max(source.sink.xy[position.matrix,3]))

network <- produce.network("Prob",Connectivity,n.days,TRUE,5,source.sink.xy,new.extent)

## ---------------------------------------------------

cl.3 <- makeCluster(2) ; registerDoParallel(cl.3)

connectivity <- foreach(x=1:length(position.matrix), .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  options(warn=-1)
  
  g2 <- network[[2]]
  comb <- network[[1]]
  
  possible.paths <- get.shortest.paths(g2,as.character( position.matrix[x] ),as.character( position.matrix ),mode="out")$vpath
  temp.res.ss <- numeric(0)
  
  for( p in 1:length(possible.paths) ) {
    
    stones.t <- as.numeric(names(possible.paths[[p]]))
    stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
    path.values <- apply( stones.t.interm , 1 , function(z) { comb[ comb[,1] == z[1] & comb[,2] == z[2] , 3 ][1] }   )
    temp.res.ss <- c(temp.res.ss,apply( t(path.values) , 1 , prod ))
    
  }
  
  temp.res.distance <- as.numeric( costDistance(raster_tr_corrected, as.matrix(cells[ position.matrix[x] , 2:3 ]) , as.matrix(cells[ position.matrix , 2:3 ]) )  ) / 1000
  
  temp.res <- data.frame( pair.from=position.matrix[x],
                          pair.to=position.matrix , 
                          differentiation= unlist(differentiation[ x, ]), 
                          distance= temp.res.distance, 
                          probability.ss = temp.res.ss
  )
  
  options(warn=0)
  
  return( temp.res ) }

stopCluster(cl.3) ; rm(cl.3)

connectivity <- do.call(rbind,connectivity)
connectivity <- connectivity[connectivity[,1] != connectivity[,2] ,]

## ---------------

connectivity[,5] <- log(connectivity[,5])
connectivity[ connectivity < -9e10 ] <- NA
connectivity[ connectivity > 9e10 ] <- NA

norm <- t(combn(position.matrix, 2))
connectivity.final <- data.frame()

for( i in 1:nrow(norm)) {
  
  t.1 <- which( connectivity[,1] == norm[i,1] & connectivity[,2] == norm[i,2] )
  t.2 <- which( connectivity[,1] == norm[i,2] & connectivity[,2] == norm[i,1] )
  
  connectivity.final <- rbind(connectivity.final,
                              data.frame(from=norm[i,1],
                                         to=norm[i,2],
                                         Differantiation = connectivity$differentiation[t.1],
                                         Distance = connectivity$distance[t.1],
                                         Connectivity.max = min(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2])) ,
                                         Connectivity.mean = mean(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2])) ,
                                         Connectivity.min = max(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2]))
                              )
  )
  
}

connectivity.final[,3] <- connectivity.final$Differantiation/(1-connectivity.final$Differantiation)
connectivity.final <- connectivity.final[complete.cases(connectivity.final),]

cor.ibd <- cor(connectivity.final$Distance,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.ibd <- lm(Differantiation ~ Distance, data=connectivity.final , na.action = na.omit)

cor.min <- cor(connectivity.final$Connectivity.min,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.min <- lm(Differantiation ~ Connectivity.min, data=connectivity.final)

cor.max <- cor(connectivity.final$Connectivity.max,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.max <- lm(Differantiation ~ Connectivity.max, data=connectivity.final)

cor.mean <- cor(connectivity.final$Connectivity.mean,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
fit.mean <- lm(Differantiation ~ Connectivity.mean, data=connectivity.final)

data.frame( day = n.days,
            aic.ibd= AIC(fit.ibd),
            r2.ibd = summary(fit.ibd)$adj.r.squared,
            cor.ibd = cor.ibd,
            aic.min= AIC(fit.min),
            r2.min = summary(fit.min)$adj.r.squared,
            cor.min = cor.min,
            aic.mean= AIC(fit.mean),
            r2.mean = summary(fit.mean)$adj.r.squared,
            cor.mean = cor.mean,
            aic.max= AIC(fit.max),
            r2.max = summary(fit.max)$adj.r.squared,
            cor.max = cor.max) 

plot(connectivity.final$Connectivity.max,connectivity.final$Differantiation,xlab="Connectivity.max" , ylab="Differantiation")
plot(connectivity.final$Distance,connectivity.final$Differantiation,xlab="Connectivity.max" , ylab="Differantiation")

## ------------------------------------------------------------------------------------------------------------------------------