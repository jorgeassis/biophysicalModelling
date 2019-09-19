rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

number.cores <- 30

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

# Include Cystoseiras and reprocess everything for final publication

for( spi in 0:34 ) {
  
              results.name <- paste0("SP#",spi)
              
              file.sampling.sites <- paste0(project.folder,"/Data/Differentiation/ID#",spi,"_Coords.txt")
              file.differentiation <- paste0(project.folder,"/Data/Differentiation/ID#",spi,"_FST.txt")
              
              transform.fst <- TRUE
              
              sampling.sites <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,2:3] 
              sampling.sites.n <- read.table(file.sampling.sites,header = T,sep=";",stringsAsFactors=F)[,1] 
              differentiation <- read.table(file.differentiation,header = T,sep=";",stringsAsFactors=F)[,-1]
              
              if( is.na(differentiation[nrow(differentiation),1]) ) { differentiation[lower.tri(differentiation)] <- t(differentiation)[lower.tri(differentiation)] }
              if( is.na(differentiation[1,ncol(differentiation)]) ) { differentiation[upper.tri(differentiation)] <- t(differentiation)[upper.tri(differentiation)] }
              
              if( nrow(differentiation) != length(sampling.sites.n) ) { next }
              
              differentiation
              
              ## ---------------
              
              # subseter <- which(sampling.sites$Lon >= min(source.sink.xy[,2]) & sampling.sites$Lon <= max(source.sink.xy[,2]),
              #                   sampling.sites$Lat >= min(source.sink.xy[,3]) & sampling.sites$Lat <= max(source.sink.xy[,3]))
              # 
              # subseter <- which(sampling.sites$Lat  < 50)
              # 
              # sampling.sites <- sampling.sites[subseter,]
              # sampling.sites.n <- sampling.sites.n[subseter]
              # differentiation <- differentiation[subseter,subseter]
              
              ## ---------------
              
              # rocky.bottoms <- raster("../Data/Rocky.tif")
              
              # subsetter <- which(extract(rocky.bottoms,source.sink.xy[,.(Lon,Lat)]) == 1)
              # source.sink.xy <- source.sink.xy[subsetter,]
              # distance.probability <- distance.probability[ Pair.from %in% source.sink.xy[,Pair] &  Pair.to %in% source.sink.xy[,Pair],]
              
              ## ---------------
              
              position.matrix <- spDists(as.matrix(source.sink.xy[source.sink.xy$Source == 1,2:3]),as.matrix(sampling.sites),longlat = TRUE)
              
              tofar.to.remove <- apply(position.matrix,2,min)
              
              position.matrix <- apply(position.matrix,2,which.min)
              position.matrix <- source.sink.xy[position.matrix,1]
              position.matrix <- unlist(position.matrix)
              
              if( sum(tofar.to.remove > 100) > 0 ) {
                
                to.remove <- which(tofar.to.remove > 100)
                
                position.matrix <- position.matrix[-to.remove]
                differentiation <- differentiation[-to.remove,-to.remove]
                sampling.sites <- sampling.sites[-to.remove,]
                sampling.sites.n <- sampling.sites.n[-to.remove]
                
              }
              
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
              
              pdf( file=paste0(project.folder,"Results/",results.name," map.pdf") , width = 9, height = 6 )
              
              par(mfrow=c(1,1),mar = c(1,1,1,1))
              plot(cost.surface,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
              points(sampling.sites)
              setkey(source.sink.xy,Pair)
              points(source.sink.xy[Pair %in% position.matrix,2:3],col="red")
              
              dev.off()
              
              ## ---------------------------------------------------
              ## ---------------------------------------------------
              
              max.days.sim <- 30
              no.zeros.conn <- numeric(30)
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
                network <- produce.network("Prob",distance.probability,n.days,FALSE,10,source.sink.xy,new.extent)
                
                ## ---------------------------------------------------
                
                cl.3 <- makeCluster(number.cores) ; registerDoParallel(cl.3)
                
                potential.connectivity <- foreach(from=position.matrix, .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
                  
                  network.x <- network[[2]]
                  connectivity.x <- network[[1]]
                  res.connectivity.to <- numeric(0)
                  res.distance <- numeric(0)
                  
                  for( to in position.matrix ) {
                  
                    sampling.sites.i <- sampling.sites[which( position.matrix == from),]
                    position.matrix.i <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(sampling.sites.i),longlat = TRUE)
                    position.matrix.i <- sort(c(unlist(position.matrix.i)),decreasing = FALSE, index.return=TRUE)$ix[1:3]
                    position.matrix.i <- source.sink.xy[position.matrix.i,1]
                    position.matrix.i <- unlist(position.matrix.i)
                    
                    sampling.sites.j <- sampling.sites[which( position.matrix == to),]
                    position.matrix.j <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(sampling.sites.j),longlat = TRUE)
                    position.matrix.j <- sort(c(unlist(position.matrix.j)),decreasing = FALSE, index.return=TRUE)$ix[1:3]
                    position.matrix.j <- source.sink.xy[position.matrix.j,1]
                    position.matrix.j <- unlist(position.matrix.j)
                    
                    positions <- expand.grid(position.matrix.i,position.matrix.j)
                    
                    for( p.i in 1:nrow(positions) ){
                      
                      if( from == to ) { path.values <- 1 ; break }
                      
                      possible.paths.y <- tryCatch({
                        get.shortest.paths(network.x,as.character( positions[p.i,1] ) , as.character( positions[p.i,2]  ),mode="out")$vpath
                        
                      }, error = function(e) {
                        error <- TRUE
                      } )

                      if( ! class(possible.paths.y) == "list" ) { possible.paths.y <- 0 }
                      
                      stones.t <- as.numeric(names(possible.paths.y[[1]]))
                      stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
                      path.values <- apply( stones.t.interm , 1 , function(z) { connectivity.x[ connectivity.x[,1] == z[1] & connectivity.x[,2] == z[2] , 3 ][1] }   )
                      
                      if( length(path.values) > 0 ) { path.values <- apply( t(path.values) , 1 , prod ) ; break }
                      if( length(path.values) == 0) { path.values <- 0 }
                      
                    }
                   
                    res.connectivity.to <- c(res.connectivity.to,path.values)
                    res.distance.t <- costDistance(raster_tr_corrected, as.matrix(source.sink.xy[Pair == from,2:3]) , as.matrix(source.sink.xy[Pair == to,2:3]) )
                    
                    if(res.distance.t == 0 | res.distance.t == Inf) { res.distance.t <- spDistsN1( as.matrix(source.sink.xy[ Pair == from , 2:3 ]), as.matrix(source.sink.xy[ Pair == to , 2:3 ]), longlat=TRUE) }
                    
                    res.distance <- c(res.distance, res.distance.t)
                    
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
                
                ## ---------------
                
                connectivity[,5] <- log(connectivity[,5])
                connectivity[ connectivity < -9e10 ] <- NA
                connectivity[ connectivity > 9e10 ] <- NA
                
                norm <- t(combn( unique(c(connectivity$pair.from,connectivity$pair.to)) , 2))
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
                                                         Connectivity.min = min(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2]),na.rm=T) ,
                                                         Connectivity.mean = mean(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2]),na.rm=T) ,
                                                         Connectivity.max = max(c(connectivity$probability.ss[t.1],connectivity$probability.ss[t.2]),na.rm=T) 
                                                         , stringsAsFactors = FALSE ) )
                }
                
                connectivity.final <- data.frame(connectivity.final,Connectivity.bi = ( connectivity.final$Connectivity.min - connectivity.final$Connectivity.max ) / connectivity.final$Connectivity.min)
                connectivity.final <- connectivity.final[connectivity.final[,1] != connectivity.final[,2] ,]
                nrow.i <- nrow(connectivity.final)
                
                connectivity.final <- connectivity.final[complete.cases(connectivity.final),]
                nrow.j <- nrow(connectivity.final)
                
                if( nrow.i == nrow.j ) { no.zeros.conn[n.days] <- 1 }

                if( nrow.j == 0 ) {  next }
                
                # connectivity[connectivity[,5] == 0,5] <- 1e-299

                if(transform.fst) { connectivity.final[,3] <- connectivity.final$Differantiation/(1-connectivity.final$Differantiation) }
                
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
                
                connectivity.per.days <- rbind(connectivity.per.days,     data.frame(day = n.days,
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
                                                                                     p.max=p.max )  )
                
              }
              
              ## ------------
              
              save( no.zeros.conn, connectivity.per.days.matrices,connectivity.per.days ,file= paste0(project.folder,"Results/Raw/",results.name,".MatricesResults.RData") )
              
              ## ------------------------------------
              ## ------------------------------------
              
              # ADD DIRECT
              
              connectivity.per.days <- connectivity.per.days[,]
              
              all.connected <- min(which(no.zeros.conn == 1))
              
              connectivity.per.days <- connectivity.per.days[ -(1:(all.connected-1)),]
              
              limits.plot <- c(min(connectivity.per.days[,c("aic.ibd","aic.min","aic.mean","aic.max")]),max(connectivity.per.days[,c("aic.ibd","aic.min","aic.mean","aic.max")]))
              threshold <- which(connectivity.per.days[,c("aic.min","aic.mean","aic.max")] == min(connectivity.per.days[,c("aic.min","aic.mean","aic.max")],na.rm=T) , arr.ind = TRUE)
              
              best.combination <- data.frame(day.conn=all.connected,connectivity.per.days[threshold[1,1],])
              best.combination
              
              write.csv(best.combination,file=paste0(project.folder,"Results/",results.name," Summary Model.csv"))
              
              ## ------------------------------------
              
              pdf( file=paste0(project.folder,"Results/",results.name," test days.pdf") , width = 9, height = 6 )
              
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
              points(connectivity.per.days$day[threshold[1,1]],connectivity.per.days[threshold[1,1],c("aic.min","aic.mean","aic.max")[threshold[1,2]]],pch=21,bg="#902828",col="black")
              legend(49, limits.plot[1] + (limits.plot[2] - limits.plot[1]), legend=c("Distance", "Ocean Min.", "Ocean Mean.", "Ocean Max."),border = "gray",col=c("#902828","#5E5E5E", "#5E5E5E","#5E5E5E","#5E5E5E"), lty=c(3,1,2,3,4), cex=0.8)
              
              dev.off()
              
              ## ------------------------------------
              
              connectivity.final <- connectivity.per.days.matrices[[ connectivity.per.days[threshold[1,1],1] - sum(! 1:30 %in% 1:length(connectivity.per.days.matrices)) ]]
              
              ## ------------------------------------
              
              fit.ibd <- lm(Differantiation ~ Distance, data=connectivity.final , na.action = na.omit)
              fit.min <- lm(Differantiation ~ Connectivity.min, data=connectivity.final)
              fit.mean <- lm(Differantiation ~ Connectivity.mean, data=connectivity.final)
              fit.max <- lm(Differantiation ~ Connectivity.max, data=connectivity.final)
              
              fit.min.bi <- lm(Differantiation ~ Connectivity.min + Connectivity.bi, data=connectivity.final)
              fit.mean.bi <- lm(Differantiation ~ Connectivity.mean + Connectivity.bi, data=connectivity.final)
              fit.max.bi <- lm(Differantiation ~ Connectivity.max + Connectivity.bi, data=connectivity.final)
              
              best.combination.2 <- data.frame(min=AIC(fit.min),min.bi=AIC(fit.min.bi),min.bi.r2=summary(fit.min.bi)$adj.r.squared,mean=AIC(fit.mean),mean.bi=AIC(fit.mean.bi),mean.bi.r2=summary(fit.mean.bi)$adj.r.squared,max=AIC(fit.max),max.bi=AIC(fit.max.bi),max.bi.r2=summary(fit.max.bi)$adj.r.squared)
              best.combination.2
              
              write.csv(best.combination.2,file=paste0(project.folder,"Results/",results.name," Summary Model Dirrect.csv"))
              
              pdf( file=paste0(project.folder,"Results/",results.name," lm best fit.pdf") , width = 9, height = 9 )
              
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

}

## ------------------------------------------------------------------------------------------------------------------------------