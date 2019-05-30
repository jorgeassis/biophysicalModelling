## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

number.cores <- 40

distance.probability <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.Distance.bm"))
distance.probability <- data.table(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

source.sink.xy <- read.big.matrix(paste0(project.folder,"/Results/source.sink.bm"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

clipper <- as(extent(min(source.sink.xy[,2]) - 2,max(source.sink.xy[,2]) + 2,min(source.sink.xy[,3]) - 2,max(source.sink.xy[,3]) + 2), "SpatialPolygons")

# Distance 

cost.surface <- raster("Data/Rasters/Mask.tif")
cost.surface <- crop(cost.surface,clipper)
cost.surface[is.na(cost.surface)] <- 0

plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

# Temperature 

cost.surface.temp <- raster("../Data/Temperature.tif")
cost.surface.temp <- crop(cost.surface.temp,clipper)
cost.surface.temp <- 1 - (cost.surface.temp / max(getValues(cost.surface.temp),na.rm=T))
cost.surface.temp[is.na(cost.surface.temp)] <- 0

plot(cost.surface.temp,box=FALSE,legend=FALSE)

raster_tr_temp <- transition(cost.surface.temp, mean, directions=8)
raster_tr_temp_corrected <- geoCorrection(raster_tr_temp, type="c", multpl=FALSE)

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## Pairwise Connectivity estimates

rocky.bottoms <- raster("../Data/Rocky.tif")

subsetter <- which(extract(rocky.bottoms,source.sink.xy[,.(Lon,Lat)]) == 1)
source.sink.xy <- source.sink.xy[subsetter,]
distance.probability <- distance.probability[ Pair.from %in% source.sink.xy[,Pair] &  Pair.to %in% source.sink.xy[,Pair],]

## -----------------------------

file.sampling.sites <- paste0(project.folder,"/Data/Differentiation/ID#0_Coords.txt")
file.differentiation <- paste0(project.folder,"/Data/Differentiation/ID#0_FST.txt")

transform.fst <- TRUE

if( transform.fst ) { project.name <- paste0(project.name,".FstT") }

sampling.sites <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,2:3] 
sampling.sites.n <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,1] 
differentiation <- read.table(file.differentiation,header = T,sep=";",stringsAsFactors=F)[,-1]

# differentiation[lower.tri(differentiation)] <- t(differentiation)[lower.tri(differentiation)]

differentiation[upper.tri(differentiation)] <- t(differentiation)[upper.tri(differentiation)]

if( nrow(differentiation) != length(sampling.sites.n) | sum(is.na(differentiation)) != length(sampling.sites.n) ) { next }

## ---------------

# subseter <- which(sampling.sites$Lon >= min(source.sink.xy[,2]) & sampling.sites$Lon <= max(source.sink.xy[,2]),
#                   sampling.sites$Lat >= min(source.sink.xy[,3]) & sampling.sites$Lat <= max(source.sink.xy[,3]))
# 
# subseter <- which(sampling.sites$Lat > 50 & sampling.sites$Lon < 1.5)
# 
# sampling.sites <- sampling.sites[-subseter,]
# sampling.sites.n <- sampling.sites.n[-subseter]
# differentiation <- differentiation[-subseter,-subseter]

## ---------------

position.matrix <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- source.sink.xy[position.matrix,1]
position.matrix <- unlist(position.matrix)

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

## ---------------------------------------------

max.days.sim <- 60
connectivity.per.days <- data.frame()
connectivity.per.days.matrices <- list()

for(n.days in 1:max.days.sim) {
  
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
    network <- produce.network("Prob",distance.probability,n.days,FALSE,5,source.sink.xy,new.extent)
    
    ## ---------------------------------------------------
    
    cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)
    
    potential.connectivity <- foreach(from=position.matrix, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
      
      network.x <- network[[2]]
      connectivity.x <- network[[1]]
      res.connectivity.to <- numeric(0)
      res.distance <- numeric(0)
      res.temp <- numeric(0)
      
      for( to in position.matrix ) {

        if( to %in% V(network.x)$name & from %in% V(network.x)$name ) {
          
          possible.paths.y <- get.shortest.paths(network.x,as.character( from ) , as.character( to ),mode="out")$vpath
          stones.t <- as.numeric(names(possible.paths.y[[1]]))
          stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
          path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
          
          if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) }
          if( length(path.values) == 0) { path.values <- 0 }
          if( from == to ) { path.values <- 1 }
          
          res.connectivity.to <- c(res.connectivity.to,path.values)
          
        } else {
          
          res.connectivity.to <- c(res.connectivity.to,0)
          
        }
        
        res.distance <- c(res.distance, costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == from,2:3]) , as.matrix(source.sink.xy[Pair == to,2:3]) ))
        res.temp <- c(res.temp, costDistance(raster_tr_temp_corrected, as.matrix(source.sink.xy[Pair == from,2:3]) , as.matrix(source.sink.xy[Pair == to,2:3]) ))
        
      }
        
      ## ---------------------
    
      zeros <- which(res.distance == 0 & position.matrix != from )
    
      if( length(zeros) > 0 ) {
    
        for(z in 1:length(zeros)){
    
          res.distance[zeros[z]] <- spDistsN1( as.matrix(source.sink.xy[ Pair == from , 2:3 ]), as.matrix(source.sink.xy[ Pair == position.matrix[z] , 2:3 ]), longlat=TRUE)
    
        }
    
      }
      
      ## ---------------------
      
      differentiation.to <- unlist(differentiation[ which( position.matrix == from), ])
      differentiation.to[is.na(differentiation.to)] <- 0
      
      temp.res <- data.frame( pair.from = from,
                              pair.to = position.matrix , 
                              differentiation = differentiation.to, 
                              distance = res.distance, 
                              temperature = res.temp, 
                              probability.ss = res.connectivity.to )
                              
      return(temp.res)
              
      ## ---------------------
      
    }
    
    stopCluster(cl.3) ; rm(cl.3) ; gc()
    
    ## ---------------------------------------------------
      
    connectivity <- do.call(rbind,potential.connectivity)
    connectivity <- connectivity[connectivity[,1] != connectivity[,2] ,]
    
    connectivity[connectivity[,"probability.ss"] == 0,5] <- NA
    connectivity <- connectivity[complete.cases(connectivity),]
    
    #connectivity[connectivity[,5] == 0,5] <- 1e-299
    
    ## ---------------
    
    connectivity[,"probability.ss"] <- log(connectivity[,"probability.ss"])
    #connectivity[ connectivity < -9e10 ] <- NA
    #connectivity[ connectivity > 9e10 ] <- NA
    
    norm <- connectivity[,1:2]
    
    connectivity.final <- data.frame()
    
    for( i in 1:nrow(norm)) {
      
      t.1 <- which( connectivity[,1] == norm[i,1] & connectivity[,2] == norm[i,2] )
      t.2 <- which( connectivity[,1] == norm[i,2] & connectivity[,2] == norm[i,1] )
      
      connectivity.final <- rbind(connectivity.final,
                                  data.frame(from=norm[i,1],
                                             to=norm[i,2],
                                             Differantiation = as.numeric(as.character(connectivity$differentiation[t.1])),
                                             Distance = connectivity$distance[t.1],
                                             Temperature = connectivity$temperature[t.1],
                                             Connectivity.max = min(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2])) ,
                                             Connectivity.mean = mean(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2])) ,
                                             Connectivity.min = max(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2]))
                                  , stringsAsFactors = FALSE ) )
    }
    
    connectivity.final <- connectivity.final[connectivity.final[,1] != connectivity.final[,2] ,]
    
    if(transform.fst) { connectivity.final[,"Differantiation"] <- connectivity.final$Differantiation/(1-connectivity.final$Differantiation) }
    
    connectivity.final <- connectivity.final[which(complete.cases(connectivity.final)),]
    
    connectivity.per.days.matrices <- c(connectivity.per.days.matrices,list(connectivity.final))
    
    cor.ibd <- cor(connectivity.final$Distance,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
    fit.ibd <- lm(Differantiation ~ Distance, data=connectivity.final , na.action = na.omit)
    
    cor.temp <- cor(connectivity.final$Temperature,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
    fit.temp <- lm(Differantiation ~ Temperature, data=connectivity.final , na.action = na.omit)
    
    cor.min <- cor(connectivity.final$Connectivity.min,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
    fit.min <- lm(Differantiation ~ Connectivity.min, data=connectivity.final)
    
    cor.max <- cor(connectivity.final$Connectivity.max,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
    fit.max <- lm(Differantiation ~ Connectivity.max, data=connectivity.final)
    
    cor.mean <- cor(connectivity.final$Connectivity.mean,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
    fit.mean <- lm(Differantiation ~ Connectivity.mean, data=connectivity.final)
    
    r2.ibd = summary(fit.ibd)$adj.r.squared
    p.ibd=summary(fit.ibd)$coefficients
    p.ibd= ifelse( nrow(p.ibd) == 2 , p.ibd[2,4] , NA)
    
    r2.temp = summary(fit.temp)$adj.r.squared
    p.temp=summary(fit.temp)$coefficients
    p.temp= ifelse( nrow(p.temp) == 2 , p.temp[2,4] , NA)
    
    r2.min = summary(fit.min)$adj.r.squared
    p.min=summary(fit.min)$coefficients
    p.min= ifelse( nrow(p.min) == 2 , p.min[2,4] , NA)
    
    r2.mean = summary(fit.mean)$adj.r.squared
    p.mean=summary(fit.mean)$coefficients
    p.mean= ifelse( nrow(p.mean) == 2 , p.mean[2,4] , NA)
    
    r2.max = summary(fit.max)$adj.r.squared
    p.max=summary(fit.max)$coefficients
    p.max= ifelse( nrow(p.max) == 2 , p.max[2,4] , NA)
    
    connectivity.per.days <- rbind(connectivity.per.days,     data.frame(day = n.days,
                                                                         aic.ibd= AIC(fit.ibd),
                                                                         r2.ibd = r2.ibd,
                                                                         cor.ibd = cor.ibd,
                                                                         p.ibd=p.ibd,
                                                                         
                                                                         aic.temp= AIC(fit.temp),
                                                                         r2.temp = r2.temp,
                                                                         cor.temp = cor.temp,
                                                                         p.temp=p.temp,
                                                                         
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
                                                                         p.max=p.max )  )
    
      
}

## ---------------------------------------------

save( connectivity.per.days.matrices,connectivity.per.days ,file= paste0(project.folder,"Results/",project.name,".Results.RData") )

## ---------------------------------------------
## ---------------------------------------------

complete.cases(connectivity.per.days)

connectivity.per.days <- connectivity.per.days[,]
limits.plot <- c(min(connectivity.per.days[,c("aic.ibd","aic.min","aic.mean","aic.max")]),max(connectivity.per.days[,c("aic.ibd","aic.min","aic.mean","aic.max")]))
threshold <- which(connectivity.per.days[,c("aic.min","aic.mean","aic.max")] == min(connectivity.per.days[,c("aic.min","aic.mean","aic.max")]) , arr.ind = TRUE)

connectivity.per.days[threshold[1,1],]

## ------------------------------------

pdf( file=paste0(project.folder,"Results/",project.name,".test.days.pdf") , width = 9, height = 6 )

par(mfrow=c(1,1),mar = c(1, 1, 1, 1))

par(mar = c(5, 5.5, 3, 3))
plot(connectivity.per.days$day,connectivity.per.days$aic.min,ylim=limits.plot,lty=1,col="#5E5E5E",type="l",ylab="",xlab="Dispersal period (day)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Model performance (AIC)",mgp=c(4,1,0)) 

lines(connectivity.per.days$day,connectivity.per.days$aic.mean,ylim=limits.plot,lty=3,col="#5E5E5E")
lines(connectivity.per.days$day,connectivity.per.days$aic.max,ylim=limits.plot,lty=4,col="#5E5E5E")
abline(h = min(connectivity.per.days[,2]), lty=3, col="#902828")
points(connectivity.per.days$day[threshold[1,1]],connectivity.per.days[threshold[1,1],c("aic.min","aic.men","aic.max")[threshold[1,2]]],pch=21,bg="#902828",col="black")
legend(49, limits.plot[1] + (limits.plot[2] - limits.plot[1]), legend=c("Distance", "Ocean Min.", "Ocean Mean.", "Ocean Max."),border = "gray",col=c("#902828","#5E5E5E", "#5E5E5E","#5E5E5E","#5E5E5E"), lty=c(3,1,2,3,4), cex=0.8)

dev.off()

## ------------------------------------
## ------------------------------------

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
      
## ------------------------------------------------------------------------------------------------------------------------------

# Plot IBD vs Best Ocean Model
# 8:8

fit.mix1 <- lm(Differantiation ~ Connectivity.min+Distance, data=connectivity.final)
fit.mix2 <- lm(Differantiation ~ Connectivity.min+Temperature, data=connectivity.final)
fit.mix <- lm(Differantiation ~ Connectivity.min+Temperature+Distance, data=connectivity.final)
summary(fit.mix1)
summary(fit.mix2)
summary(fit.mix)
AIC(fit.mix1) ; AIC(fit.mix2) ; AIC(fit.mix)

summary(fit.ibd)
summary(fit.temp)
summary(fit.min)
summary(fit.max)
summary(fit.mix1)
summary(fit.mix2)
AIC(fit.mix) ; AIC(fit.ibd) ; AIC(fit.min)
anova(fit.temp,fit.ibd)
anova(fit.max,fit.min)
anova(fit.temp,fit.max)
anova(fit.temp,fit.mean)

par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
plot(connectivity.final$Differantiation,predict(fit.mix),lty=1,col="#5E5E5E",ylab="",xlab="Observed genetic differentiation",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Predicted genetic differentiation",mgp=c(4,1,0)) 
fit.mix.line <- lm(Pred ~ Differantiation, data=data.frame(Differantiation=connectivity.final$Differantiation,Pred=predict(fit.mix)), na.action = na.omit)
lines(seq(min(connectivity.final$Differantiation),max(connectivity.final$Differantiation),length.out = 100),predict(fit.mix.line, data.frame(Differantiation=seq(min(connectivity.final$Differantiation),max(connectivity.final$Differantiation),length.out = 100))) ,lty=2,col="#902828")

## ------------------------------------------------------------------------------------------------------------------------------

