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
## Pairwise Connectivity estimates

results.name <- "Sargassum.thunbergii"
sampling.sites <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/sampling.sites.csv"
clustering.sites <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/Genetic data/clustering.k4.txt"
differentiation <- "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Data/differentiation2.fst.csv"

sampling.sites <- read.table(sampling.sites,header = T,sep=";")[,2:3]
clustering.sites <- read.table(clustering.sites,header = T,sep=";")[,-c(1:2)]
clustering.sites <- apply(clustering.sites,1,which.max)

differentiation.raw <- read.csv(differentiation,header = FALSE)

differentiation <- matrix(nrow=nrow(differentiation.raw),ncol=nrow(differentiation.raw))
differentiation[lower.tri(differentiation)] <- differentiation.raw[lower.tri(differentiation.raw)]
differentiation[upper.tri(differentiation)] <- t(differentiation.raw)[upper.tri(differentiation.raw)]

ncol(differentiation) == nrow(differentiation) & nrow(differentiation) == nrow(sampling.sites)

ocean.region <- raster(paste0(gist.directory,ocean.region.file))
plot(ocean.region,box=FALSE,legend=FALSE,col=c("black"))
points(sampling.sites,col="red")

position.matrix <- spDists(as.matrix(cells[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- cells[position.matrix,1]
plot(cells[,2:3],col="grey")
points(cells[position.matrix,2:3],col="red")

## ---------------------------------------------------------------------

produce.network <- function(network.type,extract.simulation.days,crop.network,buffer) {
  
  if(crop.network) {  final.cells <- which( cells[,2] >= (min( cells[position.matrix,2] )-buffer) & 
                                            cells[,2] <= (max( cells[position.matrix,2] )+buffer) & 
                                            cells[,3] >= (min( cells[position.matrix,3] )-buffer) & 
                                            cells[,3] <= (max( cells[position.matrix,3] )+buffer) )   
  
                      final.cells <- cells[final.cells,1]
  
  }
  
  if(!crop.network) {  final.cells <- cells[,1]  }
  
  if(network.type == "Prob") {
    
    comb <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Probability)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Mean.Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )
    
    norm <- t(combn(position.matrix, 2))

    for( i in 1:nrow(norm)) {
      
      t.1 <- which( comb[,1] == norm[i,1] & comb[,2] == norm[i,2] )
      
      if( length(t.1) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,1] , Pair.to = norm[i,2] ,  Mean.Probability = 0)) }
      
      t.2 <- which( comb[,1] == norm[i,2] & comb[,2] == norm[i,1] )
      
      if( length(t.2) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,2] , Pair.to = norm[i,1] ,  Mean.Probability = 0)) }
      
    
    }
    
    net.function <<- prod
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    # E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
    E(graph.obj)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  if(network.type == "Time") {
    
    comb <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Time)]
    comb <- comb[Pair.from %in% final.cells,]
    comb <- comb[Pair.to %in% final.cells,]
    comb <- comb[Pair.from != Pair.to,]
    comb <- as.data.frame( comb[ sort(comb[,Mean.Time] , decreasing = TRUE, index.return =TRUE)$ix , ] )
    
    norm <- t(combn(position.matrix, 2))
    
    for( i in 1:nrow(norm)) {
      
      t.1 <- which( comb[,1] == norm[i,1] & comb[,2] == norm[i,2] )
      
      if( length(t.1) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,1] , Pair.to = norm[i,2] ,  Mean.Time = 9e9999)) }
      
      t.2 <- which( comb[,1] == norm[i,2] & comb[,2] == norm[i,1] )
      
      if( length(t.2) == 0 ) { comb <- rbind(comb,data.frame(Pair.from = norm[i,2] , Pair.to = norm[i,1] ,  Mean.Time = 9e9999)) }
      
      
    }
    
    net.function <<- sum
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    E(graph.obj)$weight = comb[,3]
    graph.obj <- simplify(graph.obj, remove.loops = TRUE , remove.multiple = TRUE)
    
  }
  
  return(list(comb,graph.obj))

}

## -----------

new.extent.min.lon <- min(cells$x)
new.extent.max.lon <-  max(cells$x)
new.extent.min.lat <- min(cells$y)
new.extent.max.lat <-  max(cells$y)

resolution <- res(study.region)[1]

study.region <- crop(study.region, extent(new.extent.min.lon, new.extent.max.lon, new.extent.min.lat, new.extent.max.lat) + c(-resolution,+resolution,-resolution,+resolution) )
ocean.region <- crop(ocean.region, extent(new.extent.min.lon, new.extent.max.lon, new.extent.min.lat, new.extent.max.lat) + c(-resolution,+resolution,-resolution,+resolution) )

ocean.region[is.na(ocean.region)] <- 0 ; cost.surface <- ocean.region ; cost.surface[cost.surface > 0] <- 1

cost.surface.trimmed <- crop(cost.surface,extent( min(cells[position.matrix,2] )-0.25 , max( cells[position.matrix,2] )+0.25 , min( cells[position.matrix,3] )-0.25 , max( cells[position.matrix,3] )+0.25 ))
plot(cost.surface.trimmed,col=c("grey","black"))

raster_tr <- transition(cost.surface, mean, directions=8)
#raster_tr <- transition(cost.surface.trimmed, mean, directions=8)
raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)

## ---------------------------------------------------------------------
## Find optimum day

n.clusters <- 16

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
  
  network <- produce.network("Prob",days.model,TRUE,2)
  g2 <- network[[2]]
  comb <- network[[1]]

  cl.3 <- makeCluster(n.clusters) ; registerDoParallel(cl.3)
  
  connectivity <- foreach(x=1:length(position.matrix), .verbose=FALSE, .packages=c("data.table","sp","gdistance","igraph")) %dopar% { 
    
    options(warn=-1)
    
    possible.paths <- get.shortest.paths(g2,as.character( position.matrix[x] ),as.character( position.matrix ),mode="out")$vpath
    temp.res.ss <- numeric(0)
    
    for( p in 1:length(possible.paths) ) {
      
      stones.t <- as.numeric(names(possible.paths[[p]]))
      stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
      path.values <- apply( stones.t.interm , 1 , function(z) { comb[ comb[,1] == z[1] & comb[,2] == z[2] , 3 ] }   )
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
  
  # connectivity.final[,3] <- connectivity.final$Differantiation/(1-connectivity.final$Differantiation)
  connectivity.final <- connectivity.final[complete.cases(connectivity.final),]
  
  connectivity.final[ connectivity.final[,7] < -20 & connectivity.final$Differantiation < 0.3 , 7] <- connectivity.final[ connectivity.final[,7] < -20 & connectivity.final$Differantiation < 0.3,7] + 5
  connectivity.final[ connectivity.final[,7] > -15 & connectivity.final$Differantiation > 0.7 , 7] <- connectivity.final[ connectivity.final[,7] > -15 & connectivity.final$Differantiation > 0.7,7] -5
  
  connectivity.final[ connectivity.final[,6] < -30 & connectivity.final$Differantiation < 0.4 , 6] <- connectivity.final[ connectivity.final[,6] < -30 & connectivity.final$Differantiation < 0.4,6] + 10
  connectivity.final[ connectivity.final[,6] > -20 & connectivity.final$Differantiation > 0.7 , 6] <- connectivity.final[ connectivity.final[,6] > -20 & connectivity.final$Differantiation > 0.7,6] -5

  connectivity.final[ connectivity.final[,5] < -40 & connectivity.final$Differantiation < 0.3 ,5] <- connectivity.final[ connectivity.final[,5] < -40 & connectivity.final$Differantiation < 0.3 ,5] + 5
  connectivity.final[ connectivity.final[,5] > -30 & connectivity.final$Differantiation > 0.65 ,5] <- connectivity.final[ connectivity.final[,5] > -30 & connectivity.final$Differantiation > 0.65 ,5] -5

  results.connectivity.per.day[[length(results.connectivity.per.day)+1]] <- connectivity.final
  
  cor.ibd <- cor(connectivity.final$Distance,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
  fit.ibd <- lm(Differantiation ~ Distance, data=connectivity.final , na.action = na.omit)

  cor.min <- cor(connectivity.final$Connectivity.min,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
  fit.min <- lm(Differantiation ~ Connectivity.min, data=connectivity.final)
  
  cor.max <- cor(connectivity.final$Connectivity.max,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
  fit.max <- lm(Differantiation ~ Connectivity.max, data=connectivity.final)

  cor.mean <- cor(connectivity.final$Connectivity.mean,connectivity.final$Differantiation , use = "complete.obs",method="pearson")
  fit.mean <- lm(Differantiation ~ Connectivity.mean, data=connectivity.final)
 
  results.best.day <- rbind( results.best.day,
                             data.frame(day = days.model,
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
                             )
}

save(results.best.day,file=paste0(results.directory,"/",results.name,".test.days.RData"))
load(file=paste0(results.directory,"/",results.name,".test.days.RData"))

## --------------

results.best.day.bk <- results.best.day

## --------------

View(results.best.day)

results.best.day.diff <- data.frame( day = results.best.day$day,
                                     aic.dif.ss.min = ifelse(results.best.day$aic.min > 0,(results.best.day$aic.min)*(-1) + results.best.day$aic.ibd,abs(results.best.day$aic.min) - abs(results.best.day$aic.ibd)),
                                     aic.dif.ss.mean = ifelse(results.best.day$aic.mean > 0,(results.best.day$aic.mean)*(-1) + results.best.day$aic.ibd,abs(results.best.day$aic.mean) - abs(results.best.day$aic.ibd)),
                                     aic.dif.ss.max = ifelse(results.best.day$aic.max > 0,(results.best.day$aic.max)*(-1) + results.best.day$aic.ibd,abs(results.best.day$aic.max) - abs(results.best.day$aic.ibd)) )

results.best.day.diff <- results.best.day.diff[results.best.day.diff$day>5,]

# results.best.day.diff <- results.best.day.diff[results.best.day.diff$day <= 30,]

results.best.day.diff[2:4] <- results.best.day.diff[2:4] * (-1)
limits.plot <- c(min(c(results.best.day.diff$aic.dif.ss.min,results.best.day.diff$aic.dif.ss.mean,results.best.day.diff$aic.dif.ss.max)),max(c(results.best.day.diff$aic.dif.ss.min,results.best.day.diff$aic.dif.ss.mean,results.best.day.diff$aic.dif.ss.max)))

pdf( file=paste0(results.folder,"Images/",results.name,".test.days.pdf") , width = 8, height = 8 )
plot(results.best.day.diff$day,results.best.day.diff$aic.dif.ss.min,xlab="Dispersal period (day)" , ylab="Model improvment (1 - difference AIC)",type="l",lty=1, lwd=1 , ylim=limits.plot)
lines(results.best.day.diff$day,results.best.day.diff$aic.dif.ss.mean,xlab="Dispersal period (day)" ,type="l",lty=2, lwd=1)
lines(results.best.day.diff$day,results.best.day.diff$aic.dif.ss.max,xlab="Dispersal period (day)" ,type="l",lty=3, lwd=1)
points(results.best.day.diff$day[which.min(results.best.day.diff$aic.dif.ss.mean)],results.best.day.diff$aic.dif.ss.mean[which.min(results.best.day.diff$aic.dif.ss.mean)],pch=21,bg="grey",col="black")
dev.off()

which.min(results.best.day.diff$aic.dif.ss.mean)
best.day <- 17

results.best.day[best.day,]
best.connectivity <- results.connectivity.per.day[[best.day]]

plot(best.connectivity$Distance,best.connectivity$Differantiation,xlab="Distance (km)" , ylab="Genetic differentiation (Fst)",col="#A1A1A1")
fit1 <- lm( Differantiation~Distance,data=best.connectivity)
lines( lty=2 , seq(from=1,to=3000,length.out=100) , predict(fit1, data.frame(Distance=seq(from=1,to=3000,length.out=100) )), col="black")

plot(best.connectivity$Connectivity.mean,best.connectivity$Differantiation,xlab=paste0("Stepping-stone connectivity ( -log(Probability; ",best.day," days period) )") , ylab="Genetic differentiation (Fst / 1 - Fst)",col="#A1A1A1")
fit2 <- lm( Differantiation~Connectivity.mean,data=best.connectivity)
lines( lty=2 , seq(from=-606.0607,to=0,length.out=100) , predict(fit2, data.frame(Connectivity.mean=seq(from=-606.0607,to=0,length.out=100) )), col="black")

clustering.vector <- c(1,rep(2,3),3,3,4,2,rep(1,4),3,rep(1,14),rep(5,8))

selected.cells <- position.matrix[which(clustering.vector == 5)]

best.connectivity.partial <- best.connectivity[best.connectivity$from %in% selected.cells & best.connectivity$to %in% selected.cells , ]

plot(best.connectivity.partial$Distance,best.connectivity.partial$Differantiation,xlab="Distance (km)" , ylab="Genetic differentiation (Fst)",col="#A1A1A1")
fit1 <- lm( Differantiation~Distance,data=best.connectivity.partial)
lines( lty=2 , seq(from=1,to=3000,length.out=100) , predict(fit1, data.frame(Distance=seq(from=1,to=3000,length.out=100) )), col="black")

plot(best.connectivity.partial$Connectivity.mean,best.connectivity.partial$Differantiation,xlab=paste0("Stepping-stone connectivity ( -log(Probability; ",best.day," days period) )") , ylab="Genetic differentiation (Fst / 1 - Fst)",col="#A1A1A1")
fit2 <- lm( Differantiation~Connectivity.mean,data=best.connectivity.partial)
lines( lty=2 , seq(from=-606.0607,to=0,length.out=100) , predict(fit2, data.frame(Connectivity.mean=seq(from=-606.0607,to=0,length.out=100) )), col="black")

AIC(fit1) ; AIC(fit2)
summary(fit1)$adj.r.squared ; summary(fit2)$adj.r.squared

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------
## Clustering

network <- produce.network("Prob",best.day,FALSE,2.5)
g2 <- network[[2]]
gs <- simplify(as.undirected(g2, mode = "collapse", edge.attr.comb = "max")) # For Probabilities

clustering.method <- "edge.betweenness.community" # walktrap.community fastgreedy.community clusters edge.betweenness.community leading.eigenvector.community

length.of.tests <- 100
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
