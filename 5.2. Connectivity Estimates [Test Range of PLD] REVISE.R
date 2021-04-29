
result <- data.frame()

for( n.days in 1:60 ) {
  
new.extent <- c(min(source.sink.xy[position.matrix,2]),max(source.sink.xy[position.matrix,2]),min(source.sink.xy[position.matrix,3]),max(source.sink.xy[position.matrix,3]))
network <- produce.network("Prob",distance.probability,n.days,FALSE,10,source.sink.xy,new.extent)

## ------------

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

result <- rbind(result,
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
           p.max=p.max ) )

}
## ------------------------------------
## ------------------------------------

result

