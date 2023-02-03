
clipper <- as(extent(min(mainfile$Lon) - 10,max(mainfile$Lon) + 10,min(mainfile$Lat) - 10,max(mainfile$Lat) + 10), "SpatialPolygons")

costSurface <- raster("Data/Rasters/Mask.tif")
costSurface <- crop(costSurface,clipper)
costSurface[is.na(costSurface)] <- 0

costSurface.i <- crop(costSurface,clipper)
raster_tr <- transition(costSurface.i, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

cl.3 <- makeCluster(number.cores)
registerDoParallel(cl.3)

res.distance.i <- foreach(from=1:nrow(mainfile), .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
  
  res.distance <- numeric(0)
  
  for( to in 1:nrow(mainfile) ) {
    
    res.distance.t <- costDistance(raster_tr_corrected, as.matrix(mainfile[from,c("Lon","Lat")]) , as.matrix(mainfile[to,c("Lon","Lat")]) ) / 1000
    res.distance <- c(res.distance, res.distance.t)
    
  }
  
  ## ---------------------
  
  return(res.distance)
  
  ## ---------------------
  
}

stopCluster(cl.3) ; rm(cl.3) ; gc()

# --------

popIBD <- expand.grid(from=unique.sites,to=unique.sites,distance=NA,differentiation=NA)

for( u in 1:length(unique.sites)) {
  for( v in 1:length(unique.sites)) {
    
    pairs.u <- which( unique.sites[u] == popCoordinatesSS )
    pairs.v <- which( unique.sites[v] == popCoordinatesSS )
    popIBD[which(popIBD$from == unique.sites[u] & popIBD$to == unique.sites[v]),"distance"] <- res.distance.i[[u]][v]
    popIBD[which(popIBD$from == unique.sites[u]& popIBD$to == unique.sites[v]),"differentiation"] <- mean(unlist(popDifferentiation[pairs.u,pairs.v]),na.rm=T)

  } }

popIBD <- popIBD[complete.cases(popIBD),]

plot1 <- ggplot(popIBD, aes(x=distance, y=differentiation, group = 1)) +
  geom_point(size=1,color="#5B5B5B", alpha = 0.5) +
  geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=TRUE,size=0.5) +
  xlab(paste0("Minimum marine distance [km]")) + ylab("Genetic differentiation") + mainTheme +
  theme_minimal(base_size = 14) + mainTheme

project.name.c <- paste0(results.folder,"/Genetics/",speciesName,"/",markerName,"/IBD/")
if( ! dir.exists(project.name.c) ) { dir.create(file.path(project.name.c),recursive = TRUE, showWarnings = FALSE) }

pdf( file=paste0(project.name.c,"IBDFit.pdf") , width = 9, height = 9 )
print(plot1)
dev.off()

cor.ibd <- cor(popIBD$distance,popIBD$differentiation , use = "complete.obs",method="pearson")
fit.ibd <- lm(differentiation ~ distance, data=popIBD)
r2.ibd = summary(fit.ibd)$adj.r.squared
p.ibd = summary(fit.ibd)$coefficients
p.ibd = ifelse( nrow(p.ibd) == 2 , p.ibd[2,4] , NA)
aic.ibd = AIC(fit.ibd)

resultsIBD <- data.frame(study = mainfileName,
                                speciesName = speciesName,
                                markerName = markerName,
                                nPop = nPop,
                                nPairsRegression = nrow(connectivity.final),
                                day = pld.period,
                                aic.ibd= aic.ibd,
                                r2.ibd = r2.ibd,
                                cor.ibd = cor.ibd,
                                p.ibd=p.ibd)

write.csv(resultsIBD,file=paste0(project.name.c,"regressionIBD.csv") , sep=";")