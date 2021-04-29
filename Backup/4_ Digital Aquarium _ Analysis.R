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
       
## -------------------

plot(cells[,2:3])
points(cells[8994,2:3],col="Red")
points(cells[8359,2:3],col="Green")

plot.set.points <- c(8936 , 8939 , 8940 , 8942 , 8963 , 8994 , 9012 , 9013 , 9017 , 9021 , 9032 , 9077 , 9125 , 9329)
plot.set.points <- unique(c(cells.i,cells.j))

for(i in plot.set.points ) {
  
  plot(cells[,2:3])
  points(cells[i,2:3],col="Red")
  points(cells[Connectivity[Pair.from == i & Pair.to %in% plot.set.points & Year == 2012,Pair.to],2:3],col="Green")
  readline(prompt="Press [enter] to continue")
  
}

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
## Detemine higher retention Cells (probability)
## 
## ------------------------------------------------------------------------------------------------------------------------------

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
## Prob. vs Distance Plot
## 
## ------------------------------------------------------------------------------------------------------------------------------

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

# ------------------------------------------------------------------

# Summary 2

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}

# -------------

raw.data.currents.files <- list.files(paste0(raw.data.dir,"/","Currents"),full.names = TRUE,pattern="nc")
available.raw.data <- data.frame()

for( file in 1:length(raw.data.currents.files)) {
  nc <- nc_open( raw.data.currents.files[file] , verbose=FALSE )
  nc.date <- as.Date(ncvar_get( nc, "Time"), origin = "1970-01-01")
  nc_close( nc )
  available.raw.data <- rbind( available.raw.data, data.frame( simulation=file, year=substr(nc.date, 1, 4) , month=substr(nc.date, 6, 7) , day=substr(nc.date, 9, 10)) )
}

# -------------

months <- min(as.numeric(as.character(available.raw.data$month))):max(as.numeric(as.character(available.raw.data$month)))
years <- min(as.numeric(as.character(available.raw.data$year))):max(as.numeric(as.character(available.raw.data$year)))
years <- 2003:2011
n.new.particles.per.day <- global.simulation.parameters$n.new.particles.per.day

# Days of simulation

total.days <- apply( expand.grid(years,months) , 1 , function(x){ numberOfDays(as.Date(  paste(x[1],"-",x[2],"-","01",sep="")  , "%Y-%m-%d") ) } ) 
sum(total.days)

# Particles per cell

particles.per.cell <- n.new.particles.per.day*(sum(total.days)-1)
particles.per.cell

# Number of cells / Number of particles

nrow(cells)
nrow(cells)*particles.per.cell

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

# Baseline for analyses

npops <- length(cells.id)
pop.names <- cells.id

extract.simulation.days <- 30

network.type <- "Prob" # Time

if(network.type == "Prob") {
  
  comb.direct <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Probability)]
  
  comb <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Probability)]
  comb <- as.data.frame( comb[ sort(comb[,Mean.Probability] , decreasing = TRUE, index.return =TRUE)$ix , ] )
  net.function <- prod
  g2 <- graph.edgelist( cbind(as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  # E(g2)$weight = 1 - comb[,3] # The wheight as a negative impact on finding the closest path
  E(g2)$weight = -log(comb[,3]) # Hock, Karlo Mumby, Peter J 2015
}

if(network.type == "Time") {
  comb <- Connectivity[Max.Time <= extract.simulation.days,.(Pair.from,Pair.to,Mean.Time)]
  comb <- as.data.frame( comb[ sort(comb[,Mean.Time] , decreasing = TRUE, index.return =TRUE)$ix , ])
  net.function <- sum
  g2 <- graph.edgelist( cbind(as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
  E(g2)$weight = comb[,3] # The wheight as a negative impact on finding the closest path
}

# -------------------------------------------------------
# Function of stepping-stones

get.path.pairs <- function(from.i,pairs.to,closest.cells) {
  
      if( exists("possible.paths") ) { rm(possible.paths) }

      tryCatch( possible.paths <- get.all.shortest.paths(g2,as.character(from.i),as.character(pairs.to),mode="out")$res , error = function(e) { e <- 1 } )
  
      if( ! exists("possible.paths") ) { possible.paths <- 1 }
  
      # No path retrieved : Try different sources 
    
      if( length(possible.paths) == 1 ) {
            
            closest.i <- spDists(as.matrix(cells[,2:3]),as.matrix(cells[from.i,2:3]),longlat = TRUE)
            closest.i <- sort(c(closest.i),index.return=TRUE)$ix[2:closest.cells+1]
            
            for( p in 1:length(closest.i)) {
              
              tryCatch( possible.paths <- get.all.shortest.paths(g2,as.character(closest.i[p]),as.character(pairs.to),mode="out")$res , error=function(e) { e <- TRUE } )
              
              if( ! exists("possible.paths") ) { possible.paths <- 1 }
              
              if( length(possible.paths) > 1 ) { break }
              
            }
          
      }
      
      # Unable 
      
      if( length(possible.paths) <= 1 ) {
        
        result.matrix <- as.matrix(cbind(from.i , pairs.to , NA))
        result.matrix[ which(result.matrix[,1] == result.matrix[,3]) ,3] <- 1
        return( result.matrix )
    
      }
      
      # Able
    
      if( length(possible.paths) > 1 ) {
        
        result.matrix <- as.matrix(cbind(from.i , pairs.to , NA))
        
        for( p in 1:length(possible.paths) ) {
          
          stones.t <- as.numeric(names(possible.paths[[p]]))
          
          if( ! is.na( result.matrix[ result.matrix[,2] == stones.t[length(stones.t)] , 3 ] ) ) { next }
          
          stones.t.interm <- cbind(stones.t[-length(stones.t)],stones.t[-1])
          path.values <- apply( stones.t.interm , 1 , function(x) { comb[ comb[,1] == x[1] & comb[,2] == x[2] , 3 ] }   )
          result.matrix[ result.matrix[,2] == stones.t[length(stones.t)] , 3] <- apply( t(path.values) , 1 , net.function )
          
        }
        
        return(result.matrix)
        
      }
      

}

# -------------------------------------------------------
# Files and names

# Saccorhiza.polyschides 2.xlsx
# Fucus.vesiculosus 29.xlsx

results.name <- "Fucus.vesiculosus"
sampling.sites.directory <- "/Volumes/Laminaria/Dropbox/Manuscripts/Phylogeographic patterns in the North Atlantic and Adjacent Seas/Dispersal simulations/Data/Species"
sampling.sites.file <- "29.xlsx"

sampling.sites <- as.data.frame(read.xls(paste0(sampling.sites.directory,"/",sampling.sites.file)))[,2:3]
sampling.sites.names <- as.data.frame(read.xls(paste0(sampling.sites.directory,"/",sampling.sites.file)))[,4]
plot(ocean.region,box=FALSE,legend=FALSE,col=c("black"))
points(sampling.sites,col="red")

# -------------------------------------------------------
# Connectivity Matrix between Pairs

position.matrix <- spDists(as.matrix(cells[,2:3]),as.matrix(sampling.sites),longlat = TRUE)
position.matrix <- apply(position.matrix,2,which.min)
position.matrix <- cells[position.matrix,1]
points(cells[position.matrix,2:3],col="green")

norm <- expand.grid(from=position.matrix ,to=position.matrix ) 
norm.names <- expand.grid(from=sampling.sites.names,to=sampling.sites.names)

# --------------------------------------

# Direct connectivity

cl.2 <- makeCluster(16) ; registerDoParallel(cl.2)
connectivity.direct <- foreach(x=1:nrow(norm), .verbose=FALSE, .packages=c("data.table","sp")) %dopar% { 
  
      prob <- comb.direct[ Pair.from == norm[x,1] & Pair.to == norm[x,2] , Mean.Probability ]
      if( length(prob) == 0 ) { prob = 0 } else { prob = 1 - prob }
      return( prob ) }

stopCluster(cl.2) ; rm(cl.2)

connectivity.direct.matrix <- data.frame(from=norm.names$from , to=norm.names$to , probability=unlist(connectivity.direct))
connectivity.direct.matrix <- acast(connectivity.direct.matrix, from~to, value.var="probability")

matcher <- order(match(colnames(connectivity.direct.matrix),as.character(sampling.sites.names)))

write.table(connectivity.direct.matrix[matcher,matcher], file = paste0(results.directory,"/conn.direct.",results.name,".",extract.simulation.days,".days.txt"), sep = ";", row.names = FALSE, col.names = FALSE, na = "NA", dec = ".")

# --------------------------------------

# Stepping Stone connectivity

cl.2 <- makeCluster(16) ; registerDoParallel(cl.2)
connectivity.stepping.stone <- foreach(x=position.matrix, .combine=rbind , .verbose=FALSE, .packages=c("igraph","sp")) %dopar% { get.path.pairs( x , position.matrix , 10 ) } #
stopCluster(cl.2) ; rm(cl.2)

colnames(connectivity.stepping.stone) <- c("from","to","prob")
connectivity.stepping.stone <- connectivity.stepping.stone[,c("to","from","prob")]

connectivity.stepping.stone <- data.frame(norm.names,probability=unlist(connectivity.stepping.stone))
connectivity.stepping.stone.matrix <- acast(connectivity.stepping.stone, from~to, value.var="probability")
write.table(connectivity.stepping.stone.matrix[matcher,matcher], file = paste0(results.directory,"/conn.ss.",results.name,".",extract.simulation.days,".days.txt"), sep = ";", row.names = FALSE, col.names = FALSE, na = "NA", dec = ".")

# -----------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------
# Clustering

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
