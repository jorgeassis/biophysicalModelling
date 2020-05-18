## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")
source("Dependences.R")

## ------------------------------------------------------------------------------------------------------------------------------
## 
## # Video with Particle Flow
## # Video with assignments
## 
## ------------------------------------------------------------------------------------------------------------------------------

# Video with Particle Flow

simulation.name <- project.name

particles.video.location.x.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.video.location.x.desc"))
particles.video.location.y.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.video.location.y.desc"))

particles.lon <- attach.big.matrix(particles.video.location.x.bm.desc)
particles.lat <- attach.big.matrix(particles.video.location.y.bm.desc)

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

if(class(movie.sites.xy) == "character") { movie.sites.xy <- as.data.frame(shapefile(paste0(project.folder,movie.sites.xy)))[,2:3] }

movie.sites.xy <- movie.sites.xy[complete.cases(movie.sites.xy),]
movie.sites.id <- sort( as.vector(get.knnx( source.sink.xy[ source.sink.xy[,4] == 1,c("x","y") ] , movie.sites.xy , k = 1 + movie.sites.buffer , algorithm="kd_tree" )$nn.index) )
movie.sites.id <-  source.sink.xy[ source.sink.xy[,4] == 1, "cells.id" ][movie.sites.id]

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
sim.extent <- unique(as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$extent, split=","))))
movie.year <- unique(dbReadTable(sql, "Parameters")$movie.year)
months <- unique(as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$sim.months , split=","))))
#particles.to.sql.id <- unique(as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$particles.to.sql.id , split=","))) )
movie.sites.id <- unique(as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$movie.sites.id , split=","))))
n.hours.per.day <- unique(dbReadTable(sql, "Parameters")$n.hours.per.day)
n.particles.per.cell <- unique(dbReadTable(sql, "Parameters")$n.particles.per.cell)
source.sink.xy <- unique(dbReadTable(sql, "SourceSinkSites"))
dbDisconnect(sql)

# --------------------------------------

show.polygon.region.interest <- TRUE
sim.every.hours <- 24 / n.hours.per.day

min.lon <- sim.extent[1] ; max.lon <- sim.extent[2] ; min.lat <- sim.extent[3] ; max.lat <- sim.extent[4]
ratio <- abs(sim.extent[1]) +  abs(sim.extent[2]) : abs(sim.extent[4]) - abs(sim.extent[3])

source.sink.id <- 1:nrow(source.sink.xy)
particles.reference <- data.table( id = 1:(n.particles.per.cell * nrow(source.sink.xy) ) )
particles.reference[ , start.cell := as.numeric( sapply( source.sink.id ,function(x) { rep(x,n.particles.per.cell) })) ]
particles.to.sql.id <- particles.reference[ start.cell %in% movie.sites.id , id ]

# ---------------------------------------------------------------------------------------------------------

land.polygon <- shapefile("Data/Shapefiles/Global Landmass.shp")
crs(land.polygon) <- "+proj=longlat +ellps=WGS84"
land.polygon <- gBuffer(land.polygon, byid=TRUE, width=0)

land.polygon <- crop(land.polygon, extent(sim.extent + c(-2,+2,-2,+2)) )
plot(land.polygon, col="grey")

# ------------------

if( ! "Video" %in% list.files(paste0(project.folder,"/Results")) ) { dir.create(file.path(paste0(project.folder,"/Results/Video"))) }

# ------------------------------------------------------

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  while (format(date, format="%m") == m) { date <- date + 1 }
  return(as.integer(format(date - 1, format="%d")))
}

# --------------------------------------

days.months <- data.frame(
  day = unlist( sapply(sort(months),function(x) { 1:numberOfDays( as.Date(  paste(unique(movie.year),"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) } )) ,
  month = unlist(sapply(sort(months),function(x) { rep(x[1],( numberOfDays( as.Date(  paste(unique(movie.year),"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) ) ) } ))
)

# ------------------------------------------------------------------------------------------------------

t.steps <- nrow(particles.lat)
change.day.vect <- rep(1:unique(n.hours.per.day),length.out=t.steps+unique(sim.every.hours))[1:t.steps]
change.day <- rep(FALSE,length(change.day.vect))
change.day[change.day.vect == 1] <- TRUE

# ---------------------------------

distinctColors <- function(n) {
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(sample(col_vector, n))
}

cells.colors <- movie.sites.id
#cells.colors <- unique(particles.reference[ id %in% movie.sites.id , start.cell ])
cells.colors <- data.frame(cell=cells.colors,color=distinctColors(length(cells.colors)), stringsAsFactors = FALSE)

# ---------------------------------
# Aggregate colors

dist <- spDists( as.matrix( source.sink.xy[source.sink.xy$cells.id %in% cells.colors$cell , c(2,3) ] ) , as.matrix( source.sink.xy[source.sink.xy$cells.id %in% cells.colors$cell , c(2,3)  ] ) )
dist <- as.dist(dist)
mds.coor <- cmdscale(dist)
plot(mds.coor)
plot(hclust(dist(1-dist), method="single"))

k <- length(unique(cells.colors$cell))
hc <- hclust(dist, "single")

cells.colors.i <- distinctColors( k )
cells.colors <- data.frame(cell=cells.colors$cell,color=sapply(cutree(hc, k = k),function(x) cells.colors.i[x]))

# ---------------------------------

particles.lon.t <- as.matrix(particles.lon)
particles.lat.t <- as.matrix(particles.lat)

polygon.region.interest.xx <-  c( sim.extent[1] , sim.extent[1] , sim.extent[2] , sim.extent[2] )
polygon.region.interest.yy <-  c( sim.extent[3] , sim.extent[4] , sim.extent[4] , sim.extent[3] )

# ---------------------------------

if( extent(land.polygon)[1] < min.lon | extent(land.polygon)[2] < max.lon | extent(land.polygon)[3] < min.lat | extent(land.polygon)[4] < max.lat ) {
  
  fakePoint <- as.matrix(data.frame( Lon=c(min.lon-2,min.lon-2.0001,min.lon-2,min.lon-2.0001),Lat= c(min.lat-2,min.lat-2,min.lat-2.0001,min.lat-2.0001)))
  fakePoint <- spPolygons(as.matrix(fakePoint))
  crs(fakePoint) <- dt.projection 
  
  land.polygon <- as(land.polygon,"SpatialPolygons")
  crs(land.polygon) <- dt.projection 
  land.polygon <- gUnion(fakePoint,land.polygon)
  plot(land.polygon, col="grey")
  
}

# ---------------------------------

show.polygon.region.interest <- TRUE
show.additional.landmass.shp <- FALSE
print.day <- 0

if(show.additional.landmass.shp) {
  
  additional.landmass.shp <- shapefile(paste0("../",additional.landmass.shp))
  
}

png(file=paste0(project.folder,"/Results/Video/Video map_%02d.png"), width=1280, height=720)
par( mar=c(0,0,3,0) , bg="#ffffff")

for( t in 1:t.steps) {
  
  if( change.day[t] ) { print.day <- print.day + 1 }
  
  print.date.sim <- format(as.Date(  paste(movie.year,"-",days.months[print.day,2],"-",days.months[print.day,1],sep="")  , "%Y-%m-%d"), "%d %b %Y")

  moving.ids <- particles.to.sql.id[which(particles.lon[,t] != 0)]
  moving.cell.ids <- particles.reference[ id %in% moving.ids , start.cell]
  moving.colors <- as.character( sapply(moving.cell.ids, function(x) { cells.colors[ cells.colors$cell %in% x , "color"] } ) )
  
  moving.lons <- particles.lon[,t] ; moving.lons <- moving.lons[moving.lons != 0]
  moving.lats <- particles.lat[,t] ; moving.lats <- moving.lats[moving.lats != 0]
  
  points.plot <- data.frame(Lon = moving.lons , Lat = moving.lats , color=moving.colors )
  points.plot <- points.plot[complete.cases(points.plot),]

  print(  plot(land.polygon , col="grey" , border="grey") )
  
  if(show.additional.landmass.shp) {
    
    print(  plot(additional.landmass.shp , col="#8D1111" , border="#8D1111",add=TRUE) )
    
  }
  
  print(  title(paste0("Simulation of Potential Connectivity: ",print.date.sim), cex = 0.5, col="black") )
  print(  points(points.plot[,1], points.plot[,2], pch=19 , col=as.character(points.plot[,3]),cex=0.8) )
  
  if(show.polygon.region.interest) {
    
    polygon(polygon.region.interest.xx,polygon.region.interest.yy, lty =2)
    
  }
  
}

dev.off()

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

paste0(project.folder,"/Results/Video/Video map_%02d.png")

system( 'ffmpeg -s 1280x720 -i "/Volumes/Jellyfish/Dropbox/Manuscripts/Halodule connectivity patterns throughout West Africa//Results/Video/Video map_%02d.png" -vcodec libx264 -r 32 -pix_fmt yuv420p Halodule.mp4 -y' )
file.remove( list.files(paste0(project.folder,"/Results/Video"),pattern="png",full.names=TRUE) )

library(magick)
library(av)

imagesFiles <- list.files(paste0(project.folder,"/Results/Video/") , full.names = T)
images <-  image_scale(image_read(imagesFiles[1]), "1280")

for(i in 2:length(imagesFiles)) {
  
  images.i <- image_read(imagesFiles[i])
  images.i <- image_scale(images.i, "1280")
  images <- c(images,images.i)
  
}

image_write_video(images, path = "test", framerate = 10)

video <- image_animate(images, fps = 1)
image_write(video, "MapIsolatedMPA.gif",quality=100)

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Video with assignments

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
cell.to.process <- 1:nrow(dbReadTable(sql, "SourceSinkSites"))
particles.reference <- as.data.table(dbReadTable(sql, "ReferenceTable"))
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell
dbDisconnect(sql)

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

Connectivity <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN")) %dopar% { # 
  
  connectivity.pairs.to.sql <- data.frame()
  connectivity.temp.m <- particles.reference[ start.cell == cell.id.ref.f , ]
  
  for(y in unique(connectivity.temp.m$start.year)) {
    
    connectivity.temp <- connectivity.temp.m[ start.year == y , ]
    
    for( cell.id.ref.t in unique(connectivity.temp[ , cell.rafted ]) ) {
      
      connectivity.pairs.to.sql <- rbind(connectivity.pairs.to.sql,
                                         
                                         data.frame(  Pair.from = cell.id.ref.f,
                                                      Pair.to = cell.id.ref.t,
                                                      Number.events = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]),
                                                      Time.mean = mean(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.min = min(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.max = max(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Time.sd = sd(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                      Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / n.particles.per.cell,
                                                      Year = simulation.year ) )
    }
    
  }
  
  connectivity.pairs.to.sql[is.na(connectivity.pairs.to.sql)] <- 0
  return( connectivity.pairs.to.sql )
}

stopCluster(cl.2) ; rm(cl.2)

# ------------------------

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "Max.Probability" , "Mean.Time" , "Max.Time" , "N.events" )
Connectivity

# ------------------------

png(file=paste0(project.folder,"/Results/Video/Video assignments_%02d.png"), width=1280, height=720)
par( mar=c(0,0,0,0) , bg="#ffffff")

for( source in unique(Connectivity[ , Pair.from ]) ) {
  
  plot(source.sink.xy[,2:3], pch=20,col="Gray")
  points(source.sink.xy[Connectivity[ Pair.from == source , Pair.to ],2:3],col="Red")
  points(source.sink.xy[source,2:3],pch=20,col="Black")
  
}

dev.off()

# ------------------------
