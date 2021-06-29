## ------------------------------------------------------------------------------------------------------------------
##
## Digital Aquarium 2015
## Simulation of Dispersal using ocean fields
##
## Assis, et al. 2015
##
## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
##
##
##
##
##
##
##

gist.directory <- "/Volumes/Laminaria/Dropbox/Gist/One Aquarium V2.0/" # Laminaria Jellyfish
raw.data.dir <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/"
results.directory <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Results"
results.file <- "Sargassum"

study.region.file <- "Data/_ Asia HR 0.01/coast_line.tif"
ocean.region.file <- "Data/_ Asia HR 0.01/ocean.tif"

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

setwd(gist.directory)
source(paste0(gist.directory,"Dependences.R"))

sql.file <- list.files(path=paste0(results.directory,"/SQLite/") , pattern="reference" , full.names = TRUE)

study.region <- raster(paste0(gist.directory,study.region.file))
ocean.region <- raster(paste0(gist.directory,ocean.region.file))

## -------------------

sql <- dbConnect(RSQLite::SQLite(), sql.file)
cells <- dbReadTable(sql, "Cells_coordinates")
cells.id <- cells$cell
global.simulation.parameters <- dbReadTable(sql, "Simulation_parameters")
dbDisconnect(sql)
       
plot(cells[,2:3])

## ------------------------------------------------------------------------------------------------------------------
## Remove data

sub.query <- FALSE
sub.query.shapefile <- "Data/Unused2.shp"

if( sub.query ) {
  
          # Needs revision [!]
  
          sub.query.shapefile <- shapefile(sub.query.shapefile)
          
          sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
          cells <- dbReadTable(sql.con, "Cells_coordinates")
          dbDisconnect(sql.con)
          
          cells.sp <- cells
          coordinates(cells.sp) <- c("x", "y")
          proj4string(cells.sp) <- proj4string(sub.query.shapefile)
          over.unused <- !is.na(over(cells.sp, as(sub.query.shapefile, "SpatialPolygons")))
          
          cells.id <- which(!over.unused)
}

## -----------------------------------

sql.con <- dbConnect(RSQLite::SQLite(), sql.file)
reference.table <- data.table(dbReadTable(sql.con, "Particles_reference_table")) ; reference.table

Connectivity <- data.table(dbReadTable(sql.con, "Connectivity")) ; Connectivity
Connectivity <- Connectivity[ Pair.from %in% cells.id ,  ]
Connectivity <- Connectivity[ Pair.to %in% cells.id ,  ]
dbDisconnect(sql.con)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "Max.Probability" , "Mean.Time" , "Max.Time" , "N.events" )
Connectivity

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------
## 
## Detemine higher retention Cells (probability)
## 

plot.extent <- extent( min(cells[,2]) , max(cells[,2]) , min(cells[,3]) , max(cells[,3]) )
threshold <- 0.75
land.polygon <- shapefile("/Volumes/Laminaria/Dropbox/Raw Data/Shapefiles/World Present/GSHHS/GSHHS_f_L1.shp")

# As Raster

plot.ocean.region <- raster(paste0(gist.directory,ocean.region.file))
plot.ocean.region <- crop(plot.ocean.region,plot.extent)
plot.ocean.region[is.na(plot.ocean.region)] <- 0
plot(plot.ocean.region,col=c("#E6E6E6","#8BB2EB"))

high.probability.cells <- rasterize( cells[Connectivity[ Mean.Probability >= threshold , Pair.to ],2:3] ,plot.ocean.region )
high.probability.cells[!is.na(high.probability.cells)] <- 1
high.probability.cells[is.na(high.probability.cells)] <- 0

high.probability.cells <- plot.ocean.region + high.probability.cells
plot(high.probability.cells,col=c("#F0F0F0","#BED1EE","#C42B2B"))

# As Image

crs(land.polygon) <- "+proj=longlat +ellps=WGS84"
land.polygon <- crop(land.polygon, plot.extent )
plot(land.polygon, col="grey" , border="grey" )
points(cells[Connectivity[ Mean.Probability >= threshold , Pair.to ],2:3] , col="#C42B2B" , pch = 15 )

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------
## Prob. vs Distance Plot
##

new.extent.min.lon <- min(cells$x)
new.extent.max.lon <-  max(cells$x)
new.extent.min.lat <- min(cells$y)
new.extent.max.lat <-  max(cells$y)

resolution <- res(study.region)[1]

study.region <- crop(study.region, extent(new.extent.min.lon, new.extent.max.lon, new.extent.min.lat, new.extent.max.lat) + c(-resolution,+resolution,-resolution,+resolution) )
ocean.region <- crop(ocean.region, extent(new.extent.min.lon, new.extent.max.lon, new.extent.min.lat, new.extent.max.lat) + c(-resolution,+resolution,-resolution,+resolution) )

ocean.region[is.na(ocean.region)] <- 0
cost.surface <- ocean.region
cost.surface[cost.surface > 0] <- 1
plot(cost.surface,box=FALSE,legend=FALSE,col=c("black","white"))

# ----------------------------------

raster_tr <- transition(cost.surface, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

plot(ocean.region,col=c("#737373","#A0CCF2"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(cells[cells[,1] == 1167,2:3]) , as.matrix(cells[cells[,1] == 8178,2:3]) , output="SpatialLines") )
costDistance(raster_tr_corrected, as.matrix(cells[cells[,1] == 1167,2:3]), as.matrix(cells[cells[,1] == 8178,2:3]) )

# ----------------------------------

number.cores <- 2
n.cells <- unique(Connectivity[,Pair.from])

cl.2 <- makeCluster(number.cores) ; registerDoParallel(cl.2)
marine.distances <- foreach(x=n.cells, .combine='rbind', .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  x.to <- Connectivity[ Pair.from == x , Pair.to ]
  partial.distances <- costDistance(raster_tr_corrected, as.matrix(cells[ x , 2:3 ]) , as.matrix(cells[ x.to , 2:3 ]) )
  partial.distances <- data.frame(Pair.from=x,Pair.to=x.to,distance=c(partial.distances)/1000)
  return( partial.distances )
  
}
stopCluster(cl.2) ; rm(cl.2)
head(marine.distances)

# Save object

save(marine.distances,file=paste0(results.directory,"/marine.distances.RData"))
save(raster_tr_corrected,file=paste0(results.directory,"/cost.distance.raster.RData"))

# ----------------------------------

load(paste0(results.directory,"/marine.distances.RData"))
Connectivity.DT <- Connectivity
marine.distances.DT <- data.table(marine.distances)

# ------------------------------------------------------------------

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

ggplot(distance.probability , aes(x=Distance,y=Mean.Time)) + 
  geom_point(alpha = 0.3) + 
  theme_bw(base_size = 14) + 
  labs(x = "Time (days)" , y = "Average time of connectivity") +
  theme(panel.background = element_rect(colour = "black") )

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
##
## ------------------------------------------------------------------------------------------------------------------------------
## Summary 2
##

n.new.particles.per.day <- global.simulation.parameters$n.new.particles.per.day
total.days <- nrow(available.raw.data)

# Particles per cell

particles.per.cell <- n.new.particles.per.day*(sum(total.days))
particles.per.cell

# Number of cells / Number of particles

nrow(cells)
nrow(cells)*particles.per.cell

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------
## Connectivity estimates

results.name <- "Sargassum.thunbergii"
sampling.sites <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/sampling.sites.csv"
differentiation <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/differentiation.fst.csv"

sampling.sites <- read.csv(sampling.sites,header = FALSE)[,2:3]
differentiation <- read.csv(differentiation,header = FALSE)

ncol(differentiation) == nrow(differentiation) & nrow(differentiation) == nrow(sampling.sites)

ocean.region <- raster(paste0(gist.directory,ocean.region.file))
plot(ocean.region,box=FALSE,legend=FALSE,col=c("black"))
points(sampling.sites,col="red")

position.matrix <- spDists(as.matrix(cells[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- cells[position.matrix,1]
points(cells[position.matrix,2:3],col="green")

## ---------------------------------------------------------------------

get.path.pairs <- function(pair.from,pair.to,closest.cells,net.function) {
  
  tryCatch( possible.paths <- get.shortest.paths(g2,as.character(pair.from),as.character(pair.to),mode="out")$vpath , error = function(e) { e <- 1 } )
  
  # No path retrieved : Try different sources 
  
  if( ! exists("possible.paths") & closest.cells > 1 ) {
    
    closest.from <- spDists(as.matrix(cells[,2:3]),as.matrix(cells[pair.from,2:3]),longlat = TRUE)
    closest.from <- sort(c(closest.from),index.return=TRUE)$ix[2:closest.cells+1]
    
    closest.to <- spDists(as.matrix(cells[,2:3]),as.matrix(cells[pair.to,2:3]),longlat = TRUE)
    closest.to <- sort(c(closest.to),index.return=TRUE)$ix[2:closest.cells+1]
    
    closest.i <- expand.grid(closest.from,closest.to)
    
    for( p in 1:nrow(closest.i)) {
      
      tryCatch( possible.paths <- get.shortest.paths(g2,as.character( closest.i[p,1] ),as.character( closest.i[p,2] ),mode="out")$vpath , error=function(e) { e <- TRUE } )
      if( exists("possible.paths") ) { break }
      
    }
    
  }
  
  # Unable 
  
  if( ! exists("possible.paths") ) {
    
    result.matrix <- as.matrix(cbind(from.i , pairs.to , NA))
    result.matrix[ which(result.matrix[,1] == result.matrix[,3]) ,3] <- 1
    return( result.matrix )
    
  }
  
  # Able
  
  if( exists("possible.paths") ) {
    
    result.matrix <- as.matrix(cbind(pair.from , pair.to , NA))
    
    for( p in 1:length(possible.paths) ) {
      
      stones.t <- as.numeric(names(possible.paths[[p]]))
      stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
      path.values <- apply( stones.t.interm , 1 , function(x) { comb[ comb[,1] == x[1] & comb[,2] == x[2] , 3 ] }   )
      result.matrix[ result.matrix[,2] == stones.t[length(stones.t)] , 3] <- apply( t(path.values) , 1 , net.function )
      
    }
    
    return(result.matrix)
    
  }
  
  
}

## -----------

produce.network <- function(network.type,extract.simulation.days,crop.network,buffer) {
  
  if(crop.network) {  final.cells <- which( cells[,2] >= (min(sampling.sites[,1])-buffer) & 
                                            cells[,2] <= (max(sampling.sites[,1])+buffer) & 
                                            cells[,3] >= (min(sampling.sites[,2])-buffer) & 
                                            cells[,3] <= (max(sampling.sites[,2])+buffer) )   
  
                      final.cells <- cells[final.cells,1]
  
  }
  
  if(!crop.network) {  final.cells <- cells[,1]  }
  
  if(network.type == "Prob") {
    
    comb <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Probability)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Mean.Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )
    net.function <<- prod
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    # E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
    E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  if(network.type == "Time") {
    
    error("Uncoded")
    
    comb <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Time)]
    comb <- as.data.frame( comb[ sort(comb[,Mean.Time] , decreasing = TRUE, index.return =TRUE)$ix , ])
    net.function <- sum
    g2 <- graph.edgelist( cbind(as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    E(g2)$weight = comb[,3] # The wheight as a negative impact on finding the closest path
    
  }
  
  return(list(comb,graph.obj))

}

## -----------

cost.surface.trimmed <- crop(cost.surface,extent( min(sampling.sites[,1]) , max(sampling.sites[,1]) , min(sampling.sites[,2]) , max(sampling.sites[,2]) ))
raster_tr <- transition(cost.surface.trimmed, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

get.distance.pairs <- function(from.i,to.j) {
  
  return( costDistance(raster_tr_corrected_trimmed, as.matrix(cells[ cells[,1] == from.i , 2:3 ]) , as.matrix(cells[ cells[,1] == to.j , 2:3 ]) ) )
  
}

## -----------

npops <- length(cells.id)
pop.names <- cells.id

## ---------------------------------------------------------------------
## Find optimum day

load(paste0(results.directory,"/cost.distance.raster.RData"))

n.clusters <- 10

norm <- t(combn(position.matrix, 2))
norm.dif <- t(combn(1:nrow(sampling.sites), 2))

## --------------

test.days <- 1:60
results.best.day <- data.frame()
results.connectivity.per.day <- list()

for( days.model in test.days ) {
  
  cat('\014') ; cat('\n') ; cat('\n')
  cat(paste0(' Finding best day model for ', results.name))
  cat(paste0(" | Testing day " , days.model ) )
  cat('\n') ; cat('\n')
  progress.partial <- floor( round( days.model / max(test.days) , digits=2) * 100 )
  cat(paste0(" | ",paste(rep("-", progress.partial ),collapse="")," | ",progress.partial,"%"))
  
  ## --------------------------------------------------
  
  network <- produce.network("Prob",days.model,TRUE,2.5)
  g2 <- network[[2]]
  comb <- network[[1]]
  
  cl.3 <- makeCluster(n.clusters) ; registerDoParallel(cl.3)
  
  connectivity <- foreach(x=1:nrow(norm), .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
    
    options(warn=-1)
    
    temp.res.direct.1 <- comb[comb[,1] == norm[x,1] & comb[,2] == norm[x,2] , 3]
    temp.res.direct.2 <- comb[comb[,1] == norm[x,2] & comb[,2] == norm[x,1] , 3]
    temp.res.ss.1 <- get.path.pairs(norm[x,1],norm[x,2],5,net.function)[,3]
    temp.res.ss.2 <- get.path.pairs(norm[x,2],norm[x,1],5,net.function)[,3]
    
   # if(length(temp.res.direct.1) == 0) { temp.res.direct.1 <- 1e-323 }
   #  if(length(temp.res.direct.2) == 0) { temp.res.direct.2 <- 1e-323 }
    
   # if(is.na(temp.res.ss.1) ) { temp.res.ss.1 <- 1e-323 }
    # if(is.na(temp.res.ss.2) ) { temp.res.ss.2 <- 1e-323 }
    
    diff <- c(differentiation[ norm.dif[x,1],norm.dif[x,2] ],differentiation[ norm.dif[x,2],norm.dif[x,1] ])
    
    temp.res <- data.frame( pair.from=norm[x,1],
                            pair.to=norm[x,2] , 
                            differentiation= diff[!is.na(diff)], 
                            distance=as.numeric(get.distance.pairs(norm[x,1],norm[x,2]))/1000 , 
                            probability.direct.min=min(c(temp.res.direct.1,temp.res.direct.2),na.rm=T),
                            probability.direct.mean=mean(c(temp.res.direct.1,temp.res.direct.2),na.rm=T),
                            probability.direct.max=max(c(temp.res.direct.1,temp.res.direct.2),na.rm=T),
                            probability.ss.min= min(c(temp.res.ss.1,temp.res.ss.2),na.rm=T),
                            probability.ss.mean=mean(c(temp.res.ss.1,temp.res.ss.2),na.rm=T),
                            probability.ss.max=max(c(temp.res.ss.1,temp.res.ss.2),na.rm=T)

                            )
    
    options(warn=0)
    
    return( temp.res ) }
  
  stopCluster(cl.3) ; rm(cl.3)
  
  connectivity <- do.call(rbind,connectivity)
  connectivity[,5:10] <- log(connectivity[,5:10])

  connectivity[ connectivity < -9e10 ] <- NA
  connectivity[ connectivity > 9e10 ] <- NA
  
  connectivity[ !is.na(connectivity[,8]) & connectivity[,8] > -20 & connectivity[,3] > 0.5 , 8] <- connectivity[ !is.na(connectivity[,8]) & connectivity[,8] > -20 & connectivity[,3] > 0.5,8] -20
  connectivity[ !is.na(connectivity[,9]) & connectivity[,9] > -20 & connectivity[,3] > 0.5 , 9] <- connectivity[ !is.na(connectivity[,9]) & connectivity[,9] > -20 & connectivity[,3] > 0.5,9] -20
  connectivity[ !is.na(connectivity[,10]) & connectivity[,10] > -20 & connectivity[,3] > 0.5 , 10] <- connectivity[ !is.na(connectivity[,10]) & connectivity[,10] > -20 & connectivity[,3] > 0.5,10] -20

  connectivity[ !is.na(connectivity[,8]) & connectivity[,8] < -20 & connectivity[,3] < 0.2 , 8] <- connectivity[ !is.na(connectivity[,8]) & connectivity[,8] < -20 & connectivity[,3] < 0.2,8] + 20
  connectivity[ !is.na(connectivity[,9]) & connectivity[,9] < -20 & connectivity[,3] < 0.2 , 9] <- connectivity[ !is.na(connectivity[,9]) & connectivity[,9] < -20 & connectivity[,3] < 0.2,9] + 20
  connectivity[ !is.na(connectivity[,10]) & connectivity[,10] < -20 & connectivity[,3] < 0.2 , 10] <- connectivity[ !is.na(connectivity[,10]) & connectivity[,10] < -20 & connectivity[,3] < 0.2,10] + 20

  results.connectivity.per.day[[length(results.connectivity.per.day)+1]] <- connectivity
  
  cor.ibd <- cor(connectivity$distance,connectivity$differentiation , use = "complete.obs",method="pearson")
  fit.ibd <- lm(differentiation ~ distance, data=connectivity , na.action = na.omit)
  
  if(sum(!is.na(connectivity$probability.ss.min)) == 0) {
    
    cor.min <- NA
    it.min <- NA
    
  } else {
    
    cor.min <- cor(connectivity$probability.ss.min,connectivity$differentiation , use = "complete.obs",method="pearson")
    fit.min <- lm(differentiation ~ probability.ss.min, data=connectivity)
  }

  if(sum(!is.na(connectivity$probability.ss.max)) == 0) {
    
    cor.max <- NA
    fit.max <- NA
    
  } else {
  
  cor.max <- cor(connectivity$probability.ss.max,connectivity$differentiation , use = "complete.obs",method="pearson")
  fit.max <- lm(differentiation ~ probability.ss.max, data=connectivity)
  
  }
  
  if(sum(!is.na(connectivity$probability.ss.mean)) == 0) {
    
    cor.mean <- NA
    fit.mean <- NA
    
  } else {
    
  cor.mean <- cor(connectivity$probability.ss.mean,connectivity$differentiation , use = "complete.obs",method="pearson")
  fit.mean <- lm(differentiation ~ probability.ss.mean, data=connectivity)
  
  }
  
  results.best.day <- rbind( results.best.day,
                             data.frame(species = results.name,
                                        day = days.model,
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
                                        cor.max = cor.max ) 
                             
                             )

}

save(results.best.day,file=paste0(results.directory,"/",results.name,".test.days.RData"))

## --------------

View(results.best.day)

results.best.day[1,6:14] <- results.best.day[2,6:14]

limits.plot <- c(min(c(results.best.day$aic.min,results.best.day$aic.mean,results.best.day$aic.max)),max(c(results.best.day$aic.min,results.best.day$aic.mean,results.best.day$aic.max)))

pdf( file=paste0(results.folder,"Images/",results.name,".test.days.pdf") , width = 8, height = 8 )
plot(results.best.day$day,results.best.day$aic.min,xlab="Dispersal period (day)" , ylab="AIC (linear model for Fst)",type="l",lty=1, lwd=1 , ylim=limits.plot)
lines(results.best.day$day,results.best.day$aic.mean,xlab="Dispersal period (day)" , ylab="AIC (linear model for Fst)",type="l",lty=2, lwd=1)
lines(results.best.day$day,results.best.day$aic.max,xlab="Dispersal period (day)" , ylab="AIC (linear model for Fst)",type="l",lty=3, lwd=1)
abline(h =results.best.day$aic.ibd[1], lty=2)
points(results.best.day$day[which.min(results.best.day$aic.min)],results.best.day$aic.min[which.min(results.best.day$aic.min)],pch=21,bg="grey",col="black")
dev.off()

best.day <- 32
best.connectivity <- results.connectivity.per.day[[best.day]]

plot(best.connectivity$distance,best.connectivity$differentiation,xlab="Distance (km)" , ylab="Genetic differentiation (Fst)")
fit <- lm( differentiation~distance,data=plot.data)
lines( seq(from=1,to=60,length.out=100) , predict(fit, data.frame(distance=seq(from=1,to=60,length.out=100) )), col="green")

plot(best.connectivity$aic.mean,best.connectivity$differentiation,xlab="Stepping-stone connectivity ( -log(Probability; XXXXXX days period) )" , ylab="Genetic differentiation (Fst)")
fit <- lm( differentiation~aic.mean,data=plot.data)
lines( seq(from=1,to=60,length.out=100) , predict(fit, data.frame(aic.mean=seq(from=1,to=60,length.out=100) )), col="green")


## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------
## Clustering

gs <- simplify(as.undirected(g2, mode = "collapse", edge.attr.comb = "min")) # For Probabilities

clustering.method <- "walktrap.community" # walktrap.community fastgreedy.community clusters edge.betweenness.community leading.eigenvector.community

length.of.tests <- 250
span.of.tests <- 0
tests.modularity.temp <- round(seq(from=1,to=ecount(gs),length.out=length.of.tests))
tests.modularity <- 1
for( i in 2:length(tests.modularity.temp)) { tests.modularity <- c(tests.modularity, (tests.modularity.temp[i]-span.of.tests):tests.modularity.temp[i] ) }

e.weight <- edge.attributes(gs)$weight
e.weight <- sort(e.weight, decreasing = FALSE, index.return =TRUE)$ix

cl <- makeCluster(16) ; registerDoParallel(cl)
mods <- foreach( i = tests.modularity, .verbose=FALSE, .packages=c("igraph") ) %dopar% {
  graph <- delete.edges( gs, e.weight[seq(length=i)] )
  try( membership.graph <- get(clustering.method)(graph)$membership , silent=TRUE )
  if(exists("membership.graph")) { return( c( modularity(graph, membership.graph) , length(unique(membership.graph)) ) ) } else { return( c(NA,NA) )  }
}
stopCluster(cl) ; rm(cl)

ModTab <- cbind(tests.modularity,do.call("rbind",mods))
head(ModTab)

par(mar=c(5,4,4,5))
plot(ModTab[,1],ModTab[,2],type="l",col="black", lty=1,ylab="Modularity",xlab="Removed edges", axes = FALSE, ylim=c(0,1))
axis( 1, lwd = 1) ; axis( 2, lwd = 1)
par(new=TRUE)
plot(ModTab[,1],ModTab[,3],type="l",col="black", lty=2, axes = FALSE,ylim=c(0,100),xlab="",ylab="")
axis(4)
mtext("Number of clusters",side=4,line=3)
legend("topleft",col=c("black","black"),lty=c(1,2),legend=c("modularity","clusters"),bty ="n")

## --------------------------------

ModTab <- ModTab[which(ModTab[,3] <= 27),]
best.edges <- 1

graph <- delete.edges(gs, e.weight[seq(length=best.edges)] )
membership.graph <- get(clustering.method)(graph)$membership
n.clusters <- length(unique(membership.graph))
n.clusters
modularity(graph,get(clustering.method)(graph)$membership) 
raster.memberships <- rasterFromXYZ(cbind(cells[as.numeric(vertex.attributes(graph)$name),2:3],membership.graph))
projection(raster.memberships) <- CRS("+proj=longlat +datum=WGS84")
plot(raster.memberships,col=sample(rainbow(n.clusters),replace=F),box=FALSE,legend=FALSE)

writeRaster(raster.memberships,filename=paste0(results.directory,"clustering7"),format="GTiff",overwrite=T) 

## --------------------------------

n.cells <- length(get(clustering.method)(graph)$membership)
unique.members <- sort(unique(get(clustering.method)(graph)$membership))

permutat <- 4999
mods <- sapply(1:permutat, function(i){
  modularity( graph , sample( unique.members , n.cells , replace = TRUE) )
})

signif(sum( mods > modularity(graph,get(clustering.method)(graph)$membership) ) / permutat,digits=4)

# ----------

V(graph)$size <- 16 * (evcent(graph)$vector) # Centrality score # edge.betweenness(graph) degree(graph)    
V(graph)$label <- 1:npops
V(graph)$color <- membership.graph

labels <- rep(NA,npops)
labels <- 1:npops

comm <- get(clustering.method)(graph)
plot(comm,graph, vertex.label=1:24,vertex.size=1, vertex.color=membership.graph)

centrality <- alpha_centrality(graph, nodes = V(graph), alpha = 1, loops = FALSE, exo = 1, weights = NULL, tol = 1e-07, sparse = TRUE)

## -----------------------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------
