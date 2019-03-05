## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

sql.project.name <- "SouthAfrica"

## ------------------------------------------------------------------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",sql.project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Probability" , "Max.Probability" , "Time.mean" , "Time.max" , "Number.events" )
Connectivity

## -------------------

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## ------------------------------------------------------------------------------------------------------------------------------
## Prob. vs Distance Plot
##

cost.surface <- raster("Data/Rasters/Mask.tif")

# cost.surface <- disaggregate(cost.surface, fact=4)

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

number.cores.t <- 4
n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores.t) ; registerDoParallel(cl.2)
marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  x.to <- x.to[x.to != 0]
  partial.distances <- costDistance(raster_tr_corrected, as.matrix(source.sink.xy[ x , 2:3 ]) , as.matrix(source.sink.xy[ x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=rep(x,length(partial.distances)),Pair.to=x.to,distance=c(partial.distances)/1000)
  return( partial.distances )
  
}
stopCluster(cl.2) ; rm(cl.2)
head(marine.distances)

# ----------------------------------

distance.probability <- cbind( Connectivity, marine.distances$distance) 
colnames(distance.probability) <- c("Pair.from" , "Pair.to" , "Probability", "Max.Probability", "Time.mean", "Time.max", "Number.events","Distance")

distance.probability <- distance.probability[distance.probability$Distance != Inf,]

save(marine.distances,file=paste0(project.folder,"/Data/marine.distances.RData"))

# ----------------------------------

extract.simulation.days <- 30

distance.probability <- distance.probability[Time.max <= extract.simulation.days,]
max(distance.probability$Time.max)

# Summary 0

ggplot(distance.probability , aes(x=Distance,y=Probability)) + 
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

cells.i <- distance.probability[ Distance >= 20 & Probability > 0.065, Pair.from  ]
cells.j <- distance.probability[ Distance >= 20 & Probability > 0.065, Pair.to  ]

land.surface <- raster("Data/Rasters/Mask.tif")
plot(land.surface,box=FALSE,legend=FALSE,col=c("black"))
points(source.sink.xy[cells.i,2:3],col="red")
points(source.sink.xy[cells.j,2:3],col="green")

distance.probability[ Pair.from %in% cells.i & Pair.to %in% cells.j , ]
reference.table[ cell %in% cells.i & cell.rafted %in% cells.j , ]

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Pairwise Connectivity estimates

# for(sp in 3:31) {

file.sampling.sites <- paste0(project.folder,"/Connectivity of Laminaria Pallida/Data/Coords.csv")
file.differentiation <- paste0(project.folder,"/Connectivity of Laminaria Pallida/Data/FST.csv")

transform.fst <- TRUE
if(transform.fst) { project.name <- paste0(project.name,".FstT")}

sampling.sites <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,2:3] 
sampling.sites.n <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,1] 
differentiation <- read.table(file.differentiation,header = T,sep=";",stringsAsFactors=F)[,-1]
differentiation[upper.tri(differentiation)] <- t(differentiation)[upper.tri(differentiation)]

if( nrow(differentiation) != length(sampling.sites.n) | sum(is.na(differentiation)) != length(sampling.sites.n) ) { next }

## ---------------

subseter <- which(sampling.sites$Lon >= min(source.sink.xy[,2]) & sampling.sites$Lon <= max(source.sink.xy[,2]),
                  sampling.sites$Lat >= min(source.sink.xy[,3]) & sampling.sites$Lat <= max(source.sink.xy[,3]))

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

par(mfrow=c(1,1),mar = c(1,1,1,1))
plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
points(sampling.sites)
points(source.sink.xy[position.matrix,2:3],col="red")

## ---------------

max.days.sim <- 30
connectivity.per.days <- data.frame()
connectivity.per.days.matrices <- list()

# n.days <- 30

for(n.days in 1:max.days.sim) {
  
    sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",sql.project.name,"SimulationResults.sql"))
    Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
    source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
    dbDisconnect(sql)
    
    ## ---------------------------------------------------
    
    Connectivity.av <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
    colnames(Connectivity.av) <- c("Pair.from" , "Pair.to" , "Probability" , "Max.Probability" , "Time.mean" , "Time.max" , "Number.events" )
    Connectivity.av
    
    ## ---------------------------------------------------
    
    time.i <- Sys.time()
    if(n.days == 1) { time.f <- time.i}
    progress.percent <- round((n.days / max.days.sim) * 100)
    time.take.step.min <- round(as.numeric(difftime(time.i, time.f, units = "mins")))
    
    cat('\014')
    cat('\n')
    cat('\n Running step #',n.days,'| Time taken',time.take.step.min,'mins.')
    cat('\n')
    cat('\n',paste0(rep("-",100),collapse = ""))
    cat('\n',paste0(rep("-",progress.percent),collapse = ""),"||",progress.percent,"%")
    cat('\n',paste0(rep("-",100),collapse = ""))
    
    ## ---------------------------------------------------
      
    new.extent <- c(min(source.sink.xy[position.matrix,2]),max(source.sink.xy[position.matrix,2]),min(source.sink.xy[position.matrix,3]),max(source.sink.xy[position.matrix,3]))
    network <- produce.network("Prob",Connectivity.av,n.days,FALSE,5,source.sink.xy,new.extent)
    
    ## ---------------------------------------------------
    
    cl.3 <- makeCluster(2) ; registerDoParallel(cl.3)
    
    connectivity <- foreach(from=position.matrix, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
      
      neightbors.n <- 4
      
      options(warn=-1)
      
      network.x <- network[[2]]
      connectivity.x <- network[[1]]
      
      res.connectivity.to <- numeric(0)
      
      for( to in position.matrix ) {
        
        # HAVe A BETTEr APPROCH TO INCREASE PROBS!
        # Try with those that are zero
        
        neightbors <- spDistsN1(as.matrix(source.sink.xy[,2:3]),as.matrix(source.sink.xy[from,2:3]),longlat = TRUE)
        neightbors <- sort(neightbors,decreasing = FALSE,index.return=TRUE)$ix[2:(neightbors.n+1)] 
        
        res.connectivity.to.t <- sapply(1:neightbors.n,function(x) {
                  
                possible.paths.y <- get.shortest.paths(network.x,as.character( neightbors[x] ) , as.character( to ),mode="out")$vpath
                stones.t <- as.numeric(names(possible.paths.y[[1]]))
                stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
                path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
                
                if( length(path.values) > 0 ) { return(apply( t(path.values) , 1 , prod )) }
                if( length(path.values) == 0) { return(0) }
        })
        
        res.connectivity.to <- c(res.connectivity.to,max(res.connectivity.to.t))
        
      }

      temp.res.distance <- as.numeric( costDistance(raster_tr_corrected, as.matrix(source.sink.xy[ from , 2:3 ]) , as.matrix(source.sink.xy[ position.matrix , 2:3 ]) )  ) / 1000
      
      temp.res <- data.frame( pair.from=from,
                              pair.to=position.matrix , 
                              differentiation= unlist(differentiation[ which( position.matrix == from), ]), 
                              distance= temp.res.distance, 
                              probability.ss = res.connectivity.to
      )
      
      options(warn=0)
      
      return( temp.res ) }
    
    stopCluster(cl.3) ; rm(cl.3)
    
    connectivity <- do.call(rbind,connectivity)
    connectivity <- connectivity[connectivity[,1] != connectivity[,2] ,]
    connectivity[connectivity[,5] == 0,5] <- 1e-299

    if( sum(!is.na(connectivity$probability.ss) & connectivity$probability.ss != 0) > 2 ) {
            
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
          
          connectivity.per.days.matrices <- c(connectivity.per.days.matrices,list(connectivity.final))
          
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

          connectivity.per.days <- rbind(connectivity.per.days, data.frame( day = n.days,
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
                                                                            p.max=p.max ) )
          
    }
    
    if( sum(!is.na(connectivity$probability.ss & connectivity$probability.ss != 0)) <= 2 ) {
      
      cor.ibd <- cor.min <- cor.max <- cor.mean <- NA
      fit.ibd <- fit.min <- fit.max <- fit.mean <- NA
      
      connectivity.per.days.matrices <- c(connectivity.per.days.matrices,list(NA))

      connectivity.per.days <- rbind(connectivity.per.days, data.frame( day = n.days,
                                                                        aic.ibd= NA,
                                                                        r2.ibd = NA,
                                                                        cor.ibd = NA,
                                                                        p.ibd=NA,
                                                                        aic.min= NA,
                                                                        r2.min = NA,
                                                                        cor.min = NA,
                                                                        p.min=NA,
                                                                        aic.mean= NA,
                                                                        r2.mean = NA,
                                                                        cor.mean = NA,
                                                                        p.mean=NA,
                                                                        aic.max= NA,
                                                                        r2.max = NA,
                                                                        cor.max = NA,
                                                                        p.max=NA) )
    }
    
}

## ---------------

connectivity.per.days <- connectivity.per.days[complete.cases(connectivity.per.days),]
limits.plot <- c(min(connectivity.per.days[,c("aic.ibd","aic.min","aic.mean","aic.max")]),max(connectivity.per.days[,c("aic.ibd","aic.min","aic.mean","aic.max")]))
threshold <- which(connectivity.per.days[,c("aic.min","aic.mean","aic.max")] == min(connectivity.per.days[,c("aic.min","aic.mean","aic.max")]) , arr.ind = TRUE)

connectivity.per.days[threshold[1,1],]

write.table( connectivity.per.days[threshold[1,1],] , file =paste0(project.folder,"Results/",project.name,".result.txt") , col.names = TRUE , sep=";" , dec="." ,  quote = FALSE)

par(mfrow=c(1,1),mar = c(1, 1, 1, 1))

pdf( file=paste0(project.folder,"Results/",project.name,".test.days.pdf") , width = 9, height = 6 )

par(mar = c(5, 5.5, 3, 3))
plot(connectivity.per.days$day,connectivity.per.days$aic.min,ylim=limits.plot,lty=1,col="#5E5E5E",type="l",ylab="",xlab="Dispersal period (day)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Model performance (AIC)",mgp=c(4,1,0)) 

# lines(connectivity.per.days$day,connectivity.per.days$aic.min,ylim=limits.plot,lty=2,col="#5E5E5E")
lines(connectivity.per.days$day,connectivity.per.days$aic.mean,ylim=limits.plot,lty=3,col="#5E5E5E")
lines(connectivity.per.days$day,connectivity.per.days$aic.max,ylim=limits.plot,lty=4,col="#5E5E5E")
abline(h = min(connectivity.per.days[,2]), lty=3, col="#902828")
points(connectivity.per.days$day[threshold[1,1]],connectivity.per.days[threshold[1,1],c(2,6,10,14)[threshold[1,2]]],pch=21,bg="#902828",col="black")
legend(49, limits.plot[1] + (limits.plot[2] - limits.plot[1])/5, legend=c("Distance", "Ocean Min.", "Ocean Mean.", "Ocean Max."),border = "gray",col=c("#902828","#5E5E5E", "#5E5E5E","#5E5E5E","#5E5E5E"), lty=c(3,1,2,3,4), cex=0.8)

dev.off()

# Plot IBD vs Best Ocean Model

connectivity.final <- connectivity.per.days.matrices[[connectivity.per.days[threshold[1,1],1]]]

fit.ibd <- lm(Differantiation ~ Distance, data=connectivity.final , na.action = na.omit)
fit.min <- lm(Differantiation ~ Connectivity.min, data=connectivity.final)
fit.mean <- lm(Differantiation ~ Connectivity.mean, data=connectivity.final)
fit.max <- lm(Differantiation ~ Connectivity.max, data=connectivity.final)

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
      
# }

## ------------------------------------------------------------------------------------------------------------------------------

