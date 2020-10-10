## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2020
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("0. Project Config.R")
source("Dependences.R")

## ------------------------------------------------------------------------------------------------------------------------------
##
##
## 
## ------------------------------------------------------------------------------------------------------------------------------

## Define region of interest

if( is.null(landmass.shp)) { worldmap <- getMap(resolution = "high") }
if( ! is.null(landmass.shp)) { worldmap <- shapefile(landmass.shp) }

worldmap <- gBuffer(worldmap, byid=TRUE, width=0)
worldmap <- crop(worldmap,extent(min.lon,max.lon,min.lat,max.lat))

if( ! is.null(additional.landmass.shp) ) {
  
  if( length(additional.landmass.shp) > 1 ) { stop("Code to have multiple shapefiles") }
  worldmap <- raster::union(worldmap,shapefile(additional.landmass.shp))
  
}

worldmap <- st_as_sf(worldmap)

## --------------------------------------

shape <- raster(extent(min.lon,max.lon,min.lat,max.lat))
res(shape) <- c(sim.resolution.grid,sim.resolution.grid)
crs(shape) <- st_crs(worldmap)$input
worldmap.r <- fasterize(worldmap,shape)

coordsAll <- xyFromCell(shape,1:ncell(shape))
coordsLand <- xyFromCell(worldmap.r, Which(!is.na(worldmap.r),cells=TRUE))
coordsOcean <- xyFromCell(worldmap.r, Which(is.na(worldmap.r),cells=TRUE))

hexagons.address <- unique( geo_to_h3(coordsAll[,c("y","x")], sim.resolution ) )
hexagons.address.ocean <- unique( geo_to_h3(coordsOcean[,c("y","x")], sim.resolution ) )
hexagons.address.land <- unique( hexagons.address[ ! hexagons.address %in% hexagons.address.ocean ] )

## --------------

polygons.ocean <- h3_to_geo_boundary_sf(hexagons.address.ocean)
polygons.land <- h3_to_geo_boundary_sf(hexagons.address.land)

ggplot() + geom_sf(data = h3_to_geo_boundary_sf(hexagons.address), fill=NA, colour="red", size=0.1)
ggplot() + geom_sf(data = polygons.land, fill=NA, colour="red", size=0.1)
ggplot() + geom_sf(data = polygons.ocean, fill=NA, colour="red", size=0.1)

## --------------

## Get shore hexagons

shoreTest <- function(hexagon) {
  neighbors <- k_ring(hexagon, radius = 1)
  if(sum(hexagons.address.land %in% neighbors) > 0 & sum(hexagons.address.ocean %in% neighbors) > 0 ) { return(hexagon) } else { return(NULL) }
}

cl <- makeCluster(number.cores)
clusterExport(cl, c("hexagons.address.land","shoreTest","hexagons.address.ocean"))
hexagons.address.shore <- unique(unlist(sapply(hexagons.address.land,function(x) { shoreTest(x) })))
stopCluster(cl)

hexagons.address.land <- hexagons.address.land[! hexagons.address.land %in% hexagons.address.shore]

polygons.land <- h3_to_geo_boundary_sf(hexagons.address.land)
polygons.shore <- h3_to_geo_boundary_sf(hexagons.address.shore)

## ----------

centroid.land <- st_centroid(polygons.land)
buffer.remover <- which(st_coordinates(centroid.land)[,1] <= max.lon - buffer.val &
                        st_coordinates(centroid.land)[,1] >= min.lon + buffer.val &
                        st_coordinates(centroid.land)[,2] <= max.lat - buffer.val &
                        st_coordinates(centroid.land)[,2] >= min.lat + buffer.val )

hexagons.address.land <- hexagons.address.land[buffer.remover]
polygons.land <- polygons.land[buffer.remover,]
centroid.land <- centroid.land[buffer.remover,]

centroid.shore <- st_centroid(polygons.shore)
buffer.remover <- which(st_coordinates(centroid.shore)[,1] <= max.lon - buffer.val &
                          st_coordinates(centroid.shore)[,1] >= min.lon + buffer.val &
                          st_coordinates(centroid.shore)[,2] <= max.lat - buffer.val &
                          st_coordinates(centroid.shore)[,2] >= min.lat + buffer.val )

hexagons.address.shore <- hexagons.address.shore[buffer.remover]
polygons.shore <- polygons.shore[buffer.remover,]
centroid.shore <- centroid.shore[buffer.remover,]

if( sum(hexagons.address.shore %in% hexagons.address.ocean) > 0 | sum(hexagons.address.land %in% hexagons.address.ocean) > 0 | sum(hexagons.address.ocean %in% hexagons.address.land) > 0 ) { stop("Error :: 001") }

## ----------

ggplot() + geom_sf(data = polygons.land, fill=NA, colour="red", size=0.1)
ggplot() + geom_sf(data = polygons.shore, fill="Black", colour="Black", size=0.1)

## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

## save.image(file='../Environment.RData')
## rm(list=(ls()[ls()!="v"])); gc(reset=TRUE); load('../Environment.RData')

coastline.pts <- st_coordinates(centroid.shore)
coastline.pts <- coastline.pts[!duplicated(coastline.pts[,1:2]),]

if(   unwanted.release.coastline ){ source.sink.xy <- data.frame(cells.id=1:nrow(coastline.pts),hexagons.address=hexagons.address.shore,x=coastline.pts[,1],y=coastline.pts[,2],source=0,stringsAsFactors = FALSE) }
if( ! unwanted.release.coastline ){ source.sink.xy <- data.frame(cells.id=1:nrow(coastline.pts),hexagons.address=hexagons.address.shore,x=coastline.pts[,1],y=coastline.pts[,2],source=1,stringsAsFactors = FALSE) }

## --------------

if( ! is.null(additional.source.sink.shp) ) {
  
  additional.pts <- data.frame()
  
  for(i in 1:length(additional.source.sink.shp)){
    
    additional.shp <- shapefile(additional.source.sink.shp[i])
    crs(additional.shp) <- dt.projection
    additional.shp <- gBuffer(additional.shp, byid=TRUE, width=0.001)
    additional.shp <- st_as_sf(additional.shp)
    
    additional.shp.r <- fasterize(additional.shp,shape)
    coords.additional <- xyFromCell(additional.shp.r, Which(!is.na(additional.shp.r),cells=TRUE))
    hexagons.address.additional <- unique( geo_to_h3(coords.additional[,c("y","x")], sim.resolution ) )
    
    if( additional.source.sink.shp.force.shore ) {
      
      hexagons.address.additional.i <- character(0)
      
      for( x in hexagons.address.additional ) {
        
        dist <- h3_distance(origin = x, hexagons.address.shore)
        hexagons.address.shore.i <- hexagons.address.shore
        
        if(-1 %in% dist) { 
          hexagons.address.shore.i <- hexagons.address.shore[-which(dist == -1)]
          dist <- dist[-which(dist == -1)] }

        dist <- dist[which(!hexagons.address.shore.i %in% hexagons.address.additional.i)]
        hexagons.address.shore.i <- hexagons.address.shore.i[which(!hexagons.address.shore.i %in% hexagons.address.additional.i)]
        hexagons.address.additional.i <- c(hexagons.address.additional.i,hexagons.address.shore.i[ which.min(dist) ])

      }
      
      if(length(hexagons.address.additional.i) != length(hexagons.address.additional)) { stop("Error :: 009")}
      if( sum(! hexagons.address.additional.i %in% hexagons.address.shore ) > 0) { stop("Error :: 010")}
      
      hexagons.address.additional <- hexagons.address.additional.i

      polygons.additional <- h3_to_geo_boundary_sf(hexagons.address.additional)
      polygons.pts <- st_coordinates(st_centroid(polygons.additional))
      polygons.pts <- polygons.pts[!duplicated(polygons.pts[,1:2]),]
      
      source.sink.xy <- source.sink.xy[ ! source.sink.xy$hexagons.address %in% hexagons.address.additional,]
      hexagons.address.shore <- source.sink.xy$hexagons.address
      polygons.shore <- h3_to_geo_boundary_sf(hexagons.address.shore)
    
    }
    
    if( sum(hexagons.address.land %in% hexagons.address.additional) > 0) { stop("Error :: 002")}
    if( sum(hexagons.address.shore %in% hexagons.address.additional) > 0) { stop("Error :: 003")}
    if( sum(hexagons.address.additional %in% hexagons.address.land) > 0) { stop("Error :: 004")}
    if( sum(hexagons.address.additional %in% hexagons.address.shore) > 0) { stop("Error :: 005")}
    
    ## ----------

    additional.pts <- rbind(additional.pts,data.frame(hexagons=hexagons.address.additional,polygons.pts[,c("X","Y")]))
    
  }

  additional.source.sink.xy <- data.frame(cells.id=(max(source.sink.xy$cells.id)+1):(max(source.sink.xy$cells.id)+nrow(additional.pts)),hexagons.address=additional.pts$hexagons,x=additional.pts[,"X"],y=additional.pts[,"Y"],source=1,stringsAsFactors = FALSE) ; head(additional.source.sink.xy)
  source.sink.xy <- rbind(source.sink.xy,additional.source.sink.xy)
  source.sink.xy$cells.id <- 1:nrow(source.sink.xy)

}

## -----------------------------------------------------

## Remove unwanted release sites

if( ! is.null(unwanted.release.sites.shp) ) {
  
  source.sink.xy.t <- source.sink.xy
  coordinates(source.sink.xy.t) <- c("x","y")
  crs(source.sink.xy.t) <- crs(worldmap)

  unwanted <- shapefile(paste0(project.folder,unwanted.release.sites.shp))
  unwanted <- as(unwanted,"SpatialPolygons")
  crs(unwanted) <- crs(worldmap)
  
  points.over.polygon <- as.vector(which( ! is.na( sp::over( source.sink.xy.t , unwanted , fn = NULL) )) )

  if( length(points.over.polygon) > 0 ) { source.sink.xy[points.over.polygon,"source"] <- 0 }

}
 
if(sum(duplicated(source.sink.xy[,c("x","y")])) > 0) { stop("Error :: 007") }

## ------------------

ggplot() + geom_sf(data=worldmap) + geom_point(data = source.sink.xy[source.sink.xy$source == 0,c("x","y")], aes(x = x, y = y), size = 1, shape = 1, fill = "darkred")
ggplot() + geom_sf(data=worldmap) + geom_point(data = source.sink.xy[source.sink.xy$source == 1,c("x","y")], aes(x = x, y = y), size = 1, shape = 1, fill = "darkred")

head(source.sink.xy)
tail(source.sink.xy)

## ------------------

initial.coords <- source.sink.xy[source.sink.xy$source == 1 , c("x","y") ]
source.cells.id <- source.sink.xy[source.sink.xy$source == 1,1]
source.cells.address <- source.sink.xy[source.sink.xy$source == 1,2]

## ------------------------------------------------------------------------------------------------------------------------------

## Define Available raw data and Test if data is available

raw.data.files <- list.files(data.folder,full.names = TRUE,pattern="nc")
simulation.parameters.step <- data.frame()

for( file in 1:length(raw.data.files)) {
  
  nc <- nc_open( raw.data.files[file] , verbose=FALSE )
  nc.date <- as.Date(ncvar_get( nc, "Date"), origin = as.Date("1970-01-01")) 
  nc_close( nc )
  
  simulation.parameters.step <- rbind( simulation.parameters.step, data.frame( simulation=file, file=raw.data.files[file],  year=substr(nc.date, 1, 4) , month=substr(nc.date, 6, 7) , day=substr(nc.date, 9, 10) , stringsAsFactors = FALSE) )

}

simulation.parameters.step <- simulation.parameters.step[simulation.parameters.step$year %in% as.character(from.year:to.year) & simulation.parameters.step$day %in% sapply(from.day:to.day,function(x){ ifelse(nchar(x) > 1,x,paste0("0",x))}) & simulation.parameters.step$month %in% sapply(months.all,function(x){ ifelse(nchar(x) > 1,x,paste0("0",x))}) ,  ]

if(sum( ! from.year:to.year %in% as.numeric(simulation.parameters.step[,"year"])) + sum( ! months.all %in% as.numeric(simulation.parameters.step[,"month"])) > 0 ) { stop("Data is not available for time window")}

## ------------------------------------------------------------------------------------------------------------------------------
## Prepare video (animation) points

if( ! is.null(movie.sites.xy) ) {
  
  if(class(movie.sites.xy) == "character") { movie.sites.xy <- as.data.frame(shapefile(movie.sites.xy))[,2:3] }
  
  movie.sites.xy <- movie.sites.xy[complete.cases(movie.sites.xy),]
  movie.sites.id <- unique(sort( as.vector(get.knnx( source.sink.xy[ source.sink.xy[,"source"] == 1,c("x","y") ] , movie.sites.xy , k = 1 + movie.sites.buffer , algorithm="kd_tree" )$nn.index) ))
  
  movie.sites.id <- source.sink.xy[ source.sink.xy[,"source"] == 1, "cells.id" ][movie.sites.id]
  movie.sites.xy <- source.sink.xy[ source.sink.xy[,"cells.id"] %in% movie.sites.id ,c("x","y") ]
  
}

ggplot() + geom_sf(data=worldmap) + geom_point(data = movie.sites.xy, aes(x = x, y = y), size = 2, shape = 16, fill = "darkred")

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

## Define conditions

norm.time <- data.frame()
for(n in 1:nrow(simulation.parameters.step)) {
  norm.time <- rbind(norm.time,data.frame(do.call("rbind", replicate(n.hours.per.day, simulation.parameters.step[n,c("year" , "month" , "day")] , simplify = FALSE)) , hour=1:n.hours.per.day))
}
norm.time <- data.table(norm.time)

# ------------------

release.particles.condition <- round(seq(1,n.hours.per.day,n.hours.per.day/n.new.particles.per.day))

norm.time[ hour == release.particles.condition, release.particles := TRUE ]
norm.time[ hour != release.particles.condition, release.particles := FALSE ]

if(remove.new.particles.last.days) {
  norm.time[ ( nrow(norm.time)-(remove.new.particles.last.days.n.days*n.hours.per.day)):nrow(norm.time) , release.particles := FALSE ]
}

# ------------------

n.particles.per.cell <- (nrow(simulation.parameters.step)) * n.new.particles.per.day

## --------------------------------------------------------

## Define particles

# data.table particles.reference[id,cell,state,t.start,t.finish,cell.rafted]
# 0 unborne
# 1 living
# 2 rafted 
# 3 out of space
# 4 dead by time

particles.reference <- data.table( id = 1:(n.particles.per.cell * nrow(initial.coords) ) )
particles.reference[ , start.cell := as.numeric( sapply( source.sink.xy[source.sink.xy$source == 1 , "cells.id" ] ,function(x) { rep(x,n.particles.per.cell) })) ]
particles.reference[ , start.year := 0 ]
particles.reference[ , start.month := 0 ]
particles.reference[ , start.day := 0 ]
particles.reference[ , pos.lon := as.numeric( sapply( source.sink.xy[source.sink.xy$source == 1 , "cells.id" ] ,function(c) { rep( source.sink.xy[source.sink.xy$cells.id == c , "x" ] ,n.particles.per.cell) })) ]
particles.reference[ , pos.lat := as.numeric( sapply( source.sink.xy[source.sink.xy$source == 1 , "cells.id" ] ,function(c) { rep( source.sink.xy[source.sink.xy$cells.id == c , "y" ] ,n.particles.per.cell) })) ]
particles.reference[ , pos.alt := 0 ]
particles.reference[ , state := 0 ]
particles.reference[ , t.start := 0 ]
particles.reference[ , t.finish := 0 ]
particles.reference[ , cell.rafted := 0 ]
particles.reference[ , at.sea := 0 ]

particles.reference.names <- colnames(particles.reference)

setkey(particles.reference, id )
head(particles.reference)
nrow(particles.reference)

## --------------------------------------------------------

## Out of memory objects

if(! dir.exists(paste0(project.folder,"/Results/",project.name))) { dir.create(paste0(project.folder,"/Results/",project.name)) }
if(! dir.exists(paste0(project.folder,"/Results/",project.name,"/InternalProc"))) { dir.create(paste0(project.folder,"/Results/",project.name,"/InternalProc")) }
if(! dir.exists(paste0(project.folder,"/Results/",project.name,"/SQL"))) { dir.create(paste0(project.folder,"/Results/",project.name,"/SQL")) }

clean.dump.files(clean.dump.files=TRUE,files="raw.data.",dump.folder=paste0(project.folder,"/Results/",project.name,"/InternalProc"))
clean.dump.files(clean.dump.files=TRUE,files="particles.reference.",dump.folder=paste0(project.folder,"/Results/",project.name,"/InternalProc"))
clean.dump.files(clean.dump.files=TRUE,files="particles.video.location.",dump.folder=paste0(project.folder,"/Results/",project.name,"/InternalProc"))

particles.reference.bm <- big.matrix(nrow=nrow(particles.reference),ncol=ncol(particles.reference) , backingpath=paste0(project.folder,"Results/",project.name,"/InternalProc") , backingfile = "particles.reference.bin", descriptorfile = "particles.reference.desc")
particles.reference.bm.desc <- dget( paste0(project.folder,"Results/",project.name,"/InternalProc/particles.reference.desc"))

particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
particles.reference.bm[,1] <- unlist(particles.reference[,1])
particles.reference.bm[,2] <- unlist(particles.reference[,2])
particles.reference.bm[,6] <- unlist(particles.reference[,6])
particles.reference.bm[,7] <- unlist(particles.reference[,7])

## --------------------------------------------------------

## Data table to alocate path of particles (video)

if( movie.year > 0 ) {
  
  particles.video.id <- particles.reference[ start.cell %in% movie.sites.id , id ]
  
  particles.video.location.x.bm <- filebacked.big.matrix( nrow = length(particles.video.id), 
                                                          ncol = ( sum(simulation.parameters.step[,3] == movie.year) * n.hours.per.day ), 
                                                          backingfile = "particles.video.location.x.bin",
                                                          descriptorfile = "particles.video.location.x.desc",
                                                          backingpath=paste0(project.folder,"/Results/",project.name,"/InternalProc"))
                                    
  particles.video.location.y.bm <- filebacked.big.matrix( nrow = length(particles.video.id), 
                                                          ncol = ( sum(simulation.parameters.step[,3] == movie.year) * n.hours.per.day ), 
                                                          backingfile = "particles.video.location.y.bin",
                                                          descriptorfile = "particles.video.location.y.desc",
                                                          backingpath=paste0(project.folder,"/Results/",project.name,"/InternalProc"))

  particles.video.location.x.bm.desc <- dget( paste0(project.folder,"/Results/",project.name,"/InternalProc/particles.video.location.x.desc"))
  particles.video.location.y.bm.desc <- dget( paste0(project.folder,"/Results/",project.name,"/InternalProc/particles.video.location.y.desc"))
  
}

## ------------------------------------------------------------------------------------------------------------------

if(! exists("movie.sites.id")) { movie.sites.id <- NULL}
if(! exists("particles.video.id")) { particles.video.id <- NULL}

## -----------------

global.simulation.parameters <- data.frame(   project.name = project.name,
                                              sim.years = paste(c(from.year,to.year),collapse="-"),
                                              sim.months = paste(months.all,collapse=","),
                                              n.hours.per.day = n.hours.per.day , 
                                              n.new.particles.per.day = n.new.particles.per.day , 
                                              remove.new.particles.last.days = remove.new.particles.last.days , 
                                              longevity = longevity , 
                                              particle.max.duration = particle.max.duration , 
                                              behaviour = behaviour,
                                              n.particles.per.cell = n.particles.per.cell,
                                              movie.year = movie.year, 
                                              movie.sites.id = paste(movie.sites.id,collapse=",") , 
                                              extent = paste(c(min.lon,max.lon,min.lat,max.lat),collapse=",") )       

save(source.sink.xy, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","SourceSink.RData"))
save(global.simulation.parameters, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","Parameters.RData"))

## -----------------------

rm(particles.reference.bm) ; rm(particles.reference) ; rm(coordsAll) ; rm(coordsOcean) ; rm(worldmap.r) ; rm(coordsLand) ; rm(additional.shp.r)
gc(reset=TRUE)
memory.profile()
list.memory()

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------

## save.image(file='../Environment.RData')
## rm(list=(ls()[ls()!="v"])); gc(reset=TRUE); load('../Environment.RData')

## -------------------------------------------

## Start Simulation

start_time <- Sys.time()

time.i <- character(nrow(simulation.parameters.step))

for ( simulation.step in 1:nrow(simulation.parameters.step) ) {
  
  ## --------------------------------------------------------
  ## Progress
  
  time.i[simulation.step] <- as.character(Sys.time())
  simulation.step.previous <- simulation.step - 1
  if(simulation.step.previous == 0) { simulation.step.previous <- 1 }
  
  progress.percent <- round((simulation.step / nrow(simulation.parameters.step)) * 100)
  time.take.step.all <- round(as.numeric(difftime(time.i[simulation.step], time.i[1], units = "mins")))
  time.take.step.min <- round(as.numeric(difftime(time.i[simulation.step], time.i[simulation.step.previous], units = "mins")))
  
  cat('\014')
  cat('\n')
  cat('\n Running step #',simulation.step,'| Time taken',time.take.step.min," (total: ",time.take.step.all,")",' mins.')
  cat('\n',paste0(rep("-",100),collapse = ""))
  cat('\n',paste0(rep("-",progress.percent),collapse = ""),"||",progress.percent,"%")
  cat('\n',paste0(rep("-",100),collapse = ""))
  
  ## --------------------------------------------------------
  
  simulation.year <- simulation.parameters.step[simulation.step,"year"]
  simulation.month <- simulation.parameters.step[simulation.step,"month"]
  simulation.day <- simulation.parameters.step[simulation.step,"day"]
  simulaton.raw.data.file <- simulation.parameters.step[simulation.step,"file"]
  
  # Next day data
  
  if( simulation.step < nrow(simulation.parameters.step) ) {
    
    simulation.year.n <- simulation.parameters.step[simulation.step +1,"year"]
    simulation.month.n <- simulation.parameters.step[simulation.step +1,"month"]
    simulation.day.n <- simulation.parameters.step[simulation.step +1,"day"]
    simulaton.raw.data.file.n <- simulation.parameters.step[simulation.step +1,"file"]
    
    if( ( as.numeric(simulation.month.n) - as.numeric(simulation.month) > 1 ) & ( as.numeric(simulation.year.n ) == as.numeric(simulation.year) )  ) {
      
      simulation.year.n <- simulation.year
      simulation.month.n <- simulation.month
      simulation.day.n <- simulation.day
      simulaton.raw.data.file.n <- simulaton.raw.data.file
      
      # Clean all
      particles.reference.bm.all.inject <- attach.big.matrix(particles.reference.bm.desc)
      particles.reference.moving <- mwhich(particles.reference.bm.all.inject,c(9),list(1), list('eq') )
      particles.reference.bm.all.inject[ particles.reference.moving , 9 ] <- 4
   
    }
    
    if( ( as.numeric(simulation.month.n ) - as.numeric(simulation.month) != -11 ) & (as.numeric(simulation.year.n ) > as.numeric(simulation.year) )  ) {
      
      simulation.year.n <- simulation.year
      simulation.month.n <- simulation.month
      simulation.day.n <- simulation.day
      simulaton.raw.data.file.n <- simulaton.raw.data.file
      
      # Clean all
      particles.reference.bm.all.inject <- attach.big.matrix(particles.reference.bm.desc)
      particles.reference.moving <- mwhich(particles.reference.bm.all.inject,c(9),list(1), list('eq') )
      particles.reference.bm.all.inject[ particles.reference.moving , 9 ] <- 4      
      
    }
    
  }
  
  if( simulation.step == nrow(simulation.parameters.step) ) {
    
    simulation.year.n <- simulation.year
    simulation.month.n <- simulation.month
    simulation.day.n <- simulation.day
    simulaton.raw.data.file.n <- simulaton.raw.data.file
    
  }
  
  ## --------------------------------------------------------
  
  ## Prepare environmental data (first time only)
  
  if( simulation.step == 1 ) {
    
    simulation.raw.data <- nc_open( simulaton.raw.data.file, readunlim=FALSE )
    dim.i <- ncvar_get( simulation.raw.data, "X" )
    dim.j <- ncvar_get( simulation.raw.data, "Y" )
    Longitude <- ncvar_get( simulation.raw.data, "Longitude" )
    Latitude <- ncvar_get( simulation.raw.data, "Latitude" )
    Time <- ncvar_get( simulation.raw.data, "Date" )
    Time <- as.Date(Time, origin = "1970-01-01") 
    nc_close(simulation.raw.data)
    
    if( ! as.numeric(format(Time, "%Y"))[1] %in% from.year:to.year ) { stop("Code 01: Something is wrong with raw data extraction") }
    
    norm.field <- expand.grid( i=dim.i , j=dim.j )
    raw.data.coords <- cbind( apply( norm.field , 1 , function (x) Longitude[x[1]] ) , apply( norm.field , 1 , function (x) Latitude[x[2]] ) )
    colnames(raw.data.coords) <- c("Lon","Lat")
    
  }
  
  ## --------------------------------------------------------
  
  nc.file <- nc_open( simulaton.raw.data.file, readunlim=FALSE )
  time.start.day <- ncvar_get( nc.file, "Date" )
  time.start.day <- as.Date(time.start.day, origin = "1970-01-01") 
  nc_close(nc.file)
  
  nc.file <- nc_open( simulaton.raw.data.file.n, readunlim=FALSE )
  time.next.day <- ncvar_get( nc.file, "Date" )
  time.next.day <- as.Date(time.next.day, origin = "1970-01-01") 
  nc_close(nc.file)
  
  rd.start.day <- which( format(time.start.day, "%Y") == simulation.year & format(time.start.day, "%m") == simulation.month & format(time.start.day, "%d") == simulation.day )
  rd.next.day <- which( format(time.next.day, "%Y") == simulation.year.n & format(time.next.day, "%m") == simulation.month.n & format(time.next.day, "%d") == simulation.day.n )
  
  if( length(rd.start.day) == 0 | length(rd.next.day) == 0 ) { stop("Code 01: Something is wrong with raw data extraction") }

  ## -------------------------------------------------------------------------
   
  ## U & V Components
  
  nc.file <- nc_open( simulaton.raw.data.file, readunlim=FALSE )
  velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(1,1,rd.start.day) , count=c(length(dim.i),length(dim.j),1) )
  start.day.raw.data.u <- melt(velocity.field)[,"value"]
  velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(1,1,rd.start.day) , count=c(length(dim.i),length(dim.j),1) )
  start.day.raw.data.v <-  melt(velocity.field)[,"value"]
  nc_close(nc.file)
  
  start.day.raw.data.u[start.day.raw.data.u > 100] <- NA
  start.day.raw.data.v[start.day.raw.data.v > 100] <- NA
  start.day.raw.data.u[start.day.raw.data.u < -100] <- NA
  start.day.raw.data.v[start.day.raw.data.v < -100] <- NA
  
  nc.file <- nc_open( simulaton.raw.data.file.n, readunlim=FALSE )
  velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(1,1,rd.next.day) , count=c(length(dim.i),length(dim.j),1) )
  next.day.raw.data.u <-  melt(velocity.field)[,"value"]
  velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(1,1,rd.next.day) , count=c(length(dim.i),length(dim.j),1) )
  next.day.raw.data.v <-  melt(velocity.field)[,"value"]
  nc_close(nc.file)
  
  next.day.raw.data.u[next.day.raw.data.u > 100] <- NA
  next.day.raw.data.v[next.day.raw.data.v > 100] <- NA
  next.day.raw.data.u[next.day.raw.data.u < -100] <- NA
  next.day.raw.data.v[next.day.raw.data.v < -100] <- NA
  
  raw.data.u <- cbind(raw.data.coords,start.day.raw.data.u,next.day.raw.data.u)
  raw.data.v <- cbind(raw.data.coords,start.day.raw.data.v,next.day.raw.data.v)
  
  colnames(raw.data.u) <- c("Lon","Lat","u.start","u.next")
  colnames(raw.data.v) <- c("Lon","Lat","u.start","u.next")
  
  raw.data.u <- raw.data.u[complete.cases(raw.data.u),]
  raw.data.v <- raw.data.v[complete.cases(raw.data.v),]
  
  # ggplot() + geom_sf(data=worldmap) + geom_point(data = as.data.frame(raw.data.u[,c("Lon","Lat")]), aes(x = Lon, y = Lat), size = 1, shape = 1, fill = "darkred")
  
  # --------------------------------------------------
  
  cl <- makeCluster(number.cores)
  clusterExport(cl, c("raw.data.u","n.hours.per.day","raw.data.v"))
  
  speed.u.sec <- t( parSapply(cl = cl, 1:nrow(raw.data.u), function (x) seq( from = raw.data.u[x,3] , to = raw.data.u[x,4] , length.out = n.hours.per.day + 1 ) ))
  speed.v.sec <- t( parSapply(cl = cl, 1:nrow(raw.data.v), function (x) seq( from = raw.data.v[x,3] , to = raw.data.v[x,4] , length.out = n.hours.per.day + 1 ) ))
  
  stopCluster(cl)
  
  speed.u.sec.coords <- raw.data.u[,1:2]
  speed.v.sec.coords <- raw.data.v[,1:2]
  
  ## --------------------------------------------------------
  ## --------------------------------------------------------
  
  particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
  
  # -----------------------------------------------
  
  particles.reference.moving.i <- mwhich(particles.reference.bm.all,c(9),list(1), list('eq') )
  particles.reference.moving.dt <- data.table(matrix(particles.reference.bm.all[particles.reference.moving.i, ],ncol=length(particles.reference.names)))
  colnames(particles.reference.moving.dt) <- particles.reference.names
  
  ## --------------------------------------------------------
  ## Move particles (per day)
  
  for( h in 1:n.hours.per.day ) {
    
    t.step <- ((as.numeric(simulation.step)-1) * n.hours.per.day) + h
    
    # -----------------------------------------------
    # Release new particles, if that is the case
    
    if( norm.time[t.step,release.particles] ) {   
      
      particles.reference.new.i <- mwhich(particles.reference.bm.all,c(9),list(0), list('eq') )
      particles.reference.new.i <- data.table(matrix(particles.reference.bm.all[particles.reference.new.i,c(1,2)],ncol=2))        
      particles.reference.new.i <- particles.reference.new.i[ , .SD[1] , by=V2][,V1]
      
      # particles.reference.names
      
      particles.reference.new <- data.table(matrix(particles.reference.bm.all[particles.reference.new.i, ],ncol=length(particles.reference.names)))
      colnames(particles.reference.new) <- particles.reference.names
      
      particles.reference.new[ , 9 ] <- 1
      particles.reference.new[ , 10 ] <- t.step
      particles.reference.new[ , 3 ] <- as.numeric(simulation.year)
      particles.reference.new[ , 4 ] <- as.numeric(simulation.month)
      particles.reference.new[ , 5 ] <- as.numeric(simulation.day)
      
      setkey(particles.reference.moving.dt,id)
      particles.reference.moving.dt <-  rbindlist(list(particles.reference.moving.dt, particles.reference.new))
      
    }
    
    ## --------------------------------------------------------
    
    moving.particles.condition <- nrow(particles.reference.moving.dt) > 0
    
    # --------
    
    if( ! moving.particles.condition ) { next }
    
    if(   moving.particles.condition ) {
      
      setkey(particles.reference.moving.dt,id)
      moving.particles.id <- as.vector(unlist(particles.reference.moving.dt[, "id"]))
      points.to.interp <- particles.reference.moving.dt[, .(pos.lon,pos.lat,id)]
      
      # plot(speed.u.sec.coords) ; points(points.to.interp[, .(pos.lon,pos.lat)],col="Red")
      
      source.points.to.interp.u <- data.frame(x=speed.u.sec.coords[,1],y=speed.u.sec.coords[,2],var=speed.u.sec[,h])
      source.points.to.interp.v <- data.frame(x=speed.v.sec.coords[,1],y=speed.v.sec.coords[,2],var=speed.v.sec[,h])
      
      # Interpolate speed
      
      idwPower <- 2
      idwN <- 3
      
      idw.nearest.r <- get.knnx( source.points.to.interp.v[,1:2] , points.to.interp[,.(pos.lon,pos.lat)], k=idwN , algorithm="kd_tree" )
      idw.nearest.i <- idw.nearest.r$nn.index
      idw.nearest.d <- idw.nearest.r$nn.dist
      idw.nearest.d <- idw.nearest.d * 1e9
      
      number.cores.opt <- ifelse(round(nrow(points.to.interp) / 2500) > number.cores , number.cores , round(nrow(points.to.interp) / 2500) )
      number.cores.opt <- ifelse(number.cores.opt == 0, 1, number.cores.opt)
      
      cl <- makeCluster(number.cores.opt)
      clusterExport(cl, c("source.points.to.interp.u","source.points.to.interp.v","idw.nearest.i","idw.nearest.d","idwPower"))
      speed.u <- parSapply(cl,1:nrow(points.to.interp), function(pt.i) { (sum( source.points.to.interp.u[idw.nearest.i[pt.i,],3] / idw.nearest.d[pt.i,]^idwPower , na.rm=T)) / (sum( 1 / idw.nearest.d[pt.i,]^idwPower)) })
      speed.v <- parSapply(cl,1:nrow(points.to.interp), function(pt.i) { (sum( source.points.to.interp.v[idw.nearest.i[pt.i,],3] / idw.nearest.d[pt.i,]^idwPower , na.rm=T)) / (sum( 1 / idw.nearest.d[pt.i,]^idwPower)) })
      stopCluster(cl)

      if( max(speed.u) > 100 | min(speed.u) < -100 | max(speed.v) > 100 | min(speed.v) < -100 ) { stop("Error [!]") }
      
      mov.eastward <- speed.u * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
      mov.northward <- speed.v * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
      
      # Assign temporary positions
      
      dLon <- mov.eastward / ( 6378137 * cos( pi * (  points.to.interp[,.(pos.lat)]  / 180) ) )
      dLat <- mov.northward / 6378137
      dLon <- points.to.interp[,pos.lon] + dLon * 180/pi 
      dLat <- points.to.interp[,pos.lat] + dLat * 180/pi
      
      # -----------------------------------------------
      
      setkey(particles.reference.moving.dt,id)
      particles.reference.moving.dt[,6] <- dLon
      particles.reference.moving.dt[,7] <- dLat

      # -----------------------------------------------
      # Test over Land
    
      setkey(particles.reference.moving.dt,id)
      points.to.test <- particles.reference.moving.dt[, .(pos.lat,pos.lon)]
    
      cl <- makeCluster(number.cores)
      clusterExport(cl, c("points.to.test","geo_to_h3","sim.resolution"))
      hexagons.address.test <- parSapply(cl,1:nrow(points.to.test), function(x) { unique( geo_to_h3(points.to.test[x,], sim.resolution) ) })
      hexagons.address.test <- unlist(unique(c(hexagons.address.test)))
      stopCluster(cl)

      particles.on.sea <- which( hexagons.address.test %in% hexagons.address.ocean )
      particles.on.sea.id <- as.vector(moving.particles.id[particles.on.sea])
      setkey(particles.reference.moving.dt,id)
      particles.reference.moving.dt[ id %in% particles.on.sea.id , at.sea := 1  ]
      
      # -----------------------------------------------
      # Particles over land hexagons
      
      particles.on.land <- which( hexagons.address.test %in% hexagons.address.land )
      particles.on.land.id <- as.vector(moving.particles.id[particles.on.land])
      particles.on.land.condition <- length(particles.on.land) > 0

      if( particles.on.land.condition ) {    
        
        cells.started <- as.vector(unlist(particles.reference.moving.dt[id %in% particles.on.land.id, "start.cell"]))
        who.at.land.t.start <- as.vector(unlist(particles.reference.moving.dt[id %in% particles.on.land.id, "t.start"]))
        
        points.on.land.t <- particles.reference.moving.dt[ id %in% particles.on.land.id , .(pos.lon,pos.lat) ]
        setkey(points.to.interp,id)
        points.on.land.t.minus <- points.to.interp[ id %in% particles.on.land.id , .(pos.lon,pos.lat) ]

        cells.rafted <- get.knnx( source.sink.xy[,c("x","y")], points.on.land.t.minus, k=1 , algorithm="kd_tree" )$nn.index
        cells.rafted <- source.sink.xy[cells.rafted,"cells.id"]
        
        displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]} )
        
        # True Rafters [Distinct cell]
        
        true.rafters.id <- particles.on.land.id[ displacement != 0 ]
        true.rafters.cell <- cells.rafted[ displacement != 0 ]
        
        setkey(particles.reference.moving.dt,id)
        particles.reference.moving.dt[ id %in% true.rafters.id , state := 2 ]
        particles.reference.moving.dt[ id %in% true.rafters.id , cell.rafted := as.numeric(true.rafters.cell) ]
        particles.reference.moving.dt[ id %in% true.rafters.id , t.finish := as.numeric(t.step) ]
        
        # False Rafters [same cell of origin]
        
        setkey(particles.reference.moving.dt,id)
        non.rafters.id <- particles.on.land.id[ displacement == 0 ]
        non.rafters.cell <- cells.started[ displacement == 0 ]
        
        particles.reference.moving.dt[ id %in% non.rafters.id , state := 2 ]
        particles.reference.moving.dt[ id %in% non.rafters.id , cell.rafted := as.numeric(non.rafters.cell) ]
        particles.reference.moving.dt[ id %in% non.rafters.id , t.finish := as.numeric(t.step) ]
        
        # Allow for relocation when rafting at t.step == t.start
        
        non.rafters.t <- who.at.land.t.start[ displacement == 0 ] == t.step
        
        if( TRUE %in% non.rafters.t & allow.back.to.origin ) {

          particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , pos.lon := points.on.land.t.minus[ id %in% non.rafters.id[non.rafters.t] , pos.lon] ]
          particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , pos.lat := points.on.land.t.minus[ id %in% non.rafters.id[non.rafters.t] , pos.lat] ]

          particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , state := 1 ]
          particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , cell.rafted := as.numeric(0) ]
          particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , t.finish := as.numeric(0) ]
          
        }
        
      }
      
      # -----------------------------------------------
      # Particles over shore hexagons
      
      particles.on.shore <- which( hexagons.address.test %in% hexagons.address.shore )
      particles.on.shore.id <- as.vector(moving.particles.id[particles.on.shore])
      particles.on.shore.condition <- length(particles.on.shore) > 0
      
      if( particles.on.shore.condition ) {    
        
        cells.started <- as.vector(unlist(particles.reference.moving.dt[id %in% particles.on.shore.id, "start.cell"]))
        who.at.shore.t.start <- as.vector(unlist(particles.reference.moving.dt[id %in% particles.on.shore.id, "t.start"]))
        
        points.on.shore.t <- particles.reference.moving.dt[ id %in% particles.on.shore.id , .(pos.lon,pos.lat) ]
        setkey(points.to.interp,id)
        points.on.land.t.minus <- points.to.interp[ id %in% particles.on.shore.id , .(pos.lon,pos.lat) ]
        
        cells.rafted <- get.knnx( source.sink.xy[,c("x","y")], points.on.land.t.minus, k=1 , algorithm="kd_tree" )$nn.index
        cells.rafted <- source.sink.xy[cells.rafted,"cells.id"]
        
        displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]} )
        
        # True Rafters [Distinct cell]
        
        true.rafters.id <- particles.on.shore.id[ displacement != 0 ]
        true.rafters.cell <- cells.rafted[ displacement != 0 ]
        
        setkey(particles.reference.moving.dt,id)
        particles.reference.moving.dt[ id %in% true.rafters.id , state := 2 ]
        particles.reference.moving.dt[ id %in% true.rafters.id , cell.rafted := as.numeric(true.rafters.cell) ]
        particles.reference.moving.dt[ id %in% true.rafters.id , t.finish := as.numeric(t.step) ]
        
        # False Rafters [same cell of origin]
        
        non.rafters.id <- particles.on.shore.id[ displacement == 0 ]
        non.rafters.id.at.sea.once <- particles.reference.moving.dt[ id %in% non.rafters.id & at.sea == 1, id ]
        non.rafters.cell <- as.vector(unlist(particles.reference.moving.dt[id %in% non.rafters.id.at.sea.once, "start.cell"]))
        
        if( length(non.rafters.id.at.sea.once) > 0 ) {    
          
          particles.reference.moving.dt[ id %in% non.rafters.id.at.sea.once , state := 2 ]
          particles.reference.moving.dt[ id %in% non.rafters.id.at.sea.once, cell.rafted := as.numeric(non.rafters.cell) ]
          particles.reference.moving.dt[ id %in% non.rafters.id.at.sea.once , t.finish := as.numeric(t.step) ]
          
        }

      }
      
      # -----------------------------------------------
      # Test over additional shapefile hexagons
      
      if( ! is.null(additional.source.sink.shp) ) { 

        particles.on.additional <- which( hexagons.address.test %in% hexagons.address.additional )
        particles.on.additional.id <- as.vector(moving.particles.id[particles.on.additional])
        particles.on.additional.condition <- length(particles.on.additional) > 0
        
        if( particles.on.additional.condition ) {    
          
          cells.started <- as.vector(unlist(particles.reference.moving.dt[id %in% particles.on.additional.id, "start.cell"]))
          who.at.shore.t.start <- as.vector(unlist(particles.reference.moving.dt[id %in% particles.on.additional.id, "t.start"]))
          
          points.on.additional.t <- particles.reference.moving.dt[ id %in% particles.on.additional.id , .(pos.lon,pos.lat) ]
          setkey(points.to.interp,id)
          points.on.additional.t.minus <- points.to.interp[ id %in% particles.on.additional.id , .(pos.lon,pos.lat) ]
          
          cells.rafted <- get.knnx( source.sink.xy[,c("x","y")], points.on.additional.t.minus, k=1 , algorithm="kd_tree" )$nn.index
          cells.rafted <- source.sink.xy[cells.rafted,"cells.id"]
          
          displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]} )
          
          # True Rafters [Distinct cell]
          
          true.rafters.id <- particles.on.additional.id[ displacement != 0 ]
          true.rafters.cell <- cells.rafted[ displacement != 0 ]
          
          setkey(particles.reference.moving.dt,id)
          particles.reference.moving.dt[ id %in% true.rafters.id , state := 2 ]
          particles.reference.moving.dt[ id %in% true.rafters.id , cell.rafted := as.numeric(true.rafters.cell) ]
          particles.reference.moving.dt[ id %in% true.rafters.id , t.finish := as.numeric(t.step) ]
          
          # False Rafters [same cell of origin]
          
          non.rafters.id <- particles.on.additional.id[ displacement == 0 ]
          non.rafters.id.at.sea.once <- particles.reference.moving.dt[ id %in% non.rafters.id & at.sea == 1, id ]
          non.rafters.cell <- as.vector(unlist(particles.reference.moving.dt[id %in% non.rafters.id.at.sea.once, "start.cell"]))
          
          if( length(non.rafters.id.at.sea.once) > 0 ) {    
            
            particles.reference.moving.dt[ id %in% non.rafters.id.at.sea.once , state := 2 ]
            particles.reference.moving.dt[ id %in% non.rafters.id.at.sea.once, cell.rafted := as.numeric(non.rafters.cell) ]
            particles.reference.moving.dt[ id %in% non.rafters.id.at.sea.once , t.finish := as.numeric(t.step) ]
            
          }
        }
      }
      
      # -----------------------------------------------
      # Out of space (study region), if TRUE, place particles on hold
      
      out.of.space <- dLon > max.lon | dLon < min.lon | dLat > max.lat | dLat < min.lat
      out.of.space.ids <- moving.particles.id[out.of.space]
      
      setkey(particles.reference.moving.dt,id)
      particles.reference.moving.dt[ id %in% out.of.space.ids , state := 3  ]
      
    }
    
    # -----------------------------------------------
    # End of day, Kill by Longevity
    
    if ( longevity ) {   
      
      max.duration.id <- particles.reference.moving.dt[ ( t.step - t.start ) > ( n.hours.per.day * particle.max.duration ) , id ]
      setkey(particles.reference.moving.dt, id )
      particles.reference.moving.dt[ id %in% max.duration.id , state := 4 ]
      
    }
    
    ## ---------------------------------------------------------------
    ## Save positions to Video matrix (if condition matched) 
    
    if( ! is.null(movie.sites.xy) & movie.year == as.numeric(simulation.year) ) {
      
      ## --------------------------------------------------------
      
      particles.video.location.x.bm.i <- attach.big.matrix(particles.video.location.x.bm.desc)
      particles.video.location.y.bm.i <- attach.big.matrix(particles.video.location.y.bm.desc)
      
      ## --------------------------------------------------------
      
      setkey(particles.reference.moving.dt,id)
      particles.video.id.moving <- particles.reference.moving.dt[ state == 1 & id %in% particles.video.id,id]
      particles.video.id.moving.condition <- length(particles.video.id.moving) > 0
      
      if( particles.video.id.moving.condition ) { 
        
        t.step.movie <- ((as.numeric(which( which(as.numeric(simulation.parameters.step[,3]) == movie.year) == simulation.step))-1) * n.hours.per.day) + h
        
        setkey(particles.reference.moving.dt,id)
        particles.video.location.x.bm.i[ which(particles.video.id %in% particles.video.id.moving) , t.step.movie ] <- unlist(particles.reference.moving.dt[ id %in% particles.video.id.moving,pos.lon])
        particles.video.location.y.bm.i[ which(particles.video.id %in% particles.video.id.moving) , t.step.movie ] <- unlist(particles.reference.moving.dt[ id %in% particles.video.id.moving,pos.lat])
        
      }
    }
    
  }
  
  ## ---------------------------------------------------
  ## Inject particles to final object [ particles.reference.bm.all ] 
  
  setkey(particles.reference.moving.dt,id)
  particles.reference.bm.all.inject <- attach.big.matrix(particles.reference.bm.desc)
  sect.loop <- particles.reference.moving.dt[ state != 0 , ]

  particles.reference.bm.all.inject[ sect.loop$id , 3 ] <- as.numeric(unlist(sect.loop[  , "start.year"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 4 ] <- as.numeric(unlist(sect.loop[  , "start.month"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 5 ] <- as.numeric(unlist(sect.loop[  , "start.day"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 6 ] <- as.numeric(unlist(sect.loop[  , "pos.lon"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 7 ] <- as.numeric(unlist(sect.loop[  , "pos.lat"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 9 ] <- as.numeric(unlist(sect.loop[  , "state"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 10 ] <- as.numeric(unlist(sect.loop[  , "t.start"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 11 ] <- as.numeric(unlist(sect.loop[  , "t.finish"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 12 ] <- as.numeric(unlist(sect.loop[  , "cell.rafted"] ))
  particles.reference.bm.all.inject[ sect.loop$id , 13 ] <- as.numeric(unlist(sect.loop[  , "at.sea"] ))
  
  rm(particles.reference.moving.dt)
  rm(sect.loop)
  
  gc(reset=TRUE)
  
}

end_time <- Sys.time()
end_time - start_time

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
# Save Reference Table in SQL

# Erase those that did not acomplish 
# 0 unborne 
# 1 living 
# 2 rafted *** 
# 3 on hold out of space 
# 4 dead by time 

particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
if( length(mwhich(particles.reference.bm,c(9),list(0), list('eq'))) > 0 ) { print("An error may have occurred [!]")}

zeros <- mwhich(particles.reference.bm,c(9),list(0), list('eq'))
unique(particles.reference.bm[zeros,2])
plot(particles.reference.bm[zeros,6:7])
nrow(particles.reference.bm)

particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
particles.reference.bm.i <- mwhich(particles.reference.bm,c(9),list(0), list('neq'))
particles.reference.bm <- particles.reference.bm[particles.reference.bm.i,]
colnames(particles.reference.bm) <- particles.reference.names
nrow(particles.reference.bm)

particles.reference.bm <- particles.reference.bm[particles.reference.bm[,"cell.rafted"] != 0,]
ReferenceTable <- data.frame( particles.reference.bm , travel.time= ( 1 + particles.reference.bm[,11] - particles.reference.bm[,10] ) / n.hours.per.day )
head(ReferenceTable)
nrow(ReferenceTable)

save(ReferenceTable, file = paste0(project.folder,"/Results/",project.name,"/InternalProc/","ReferenceTable.RData"))

## ------------------------------------------------------------------------------------------------------------------
# Time taken

seq.t <- numeric(length(time.i)) ; seq.t[1] <- 0
for( t.seq in 2:length(seq.t)) {
  seq.t[t.seq] <- as.numeric(difftime(time.i[t.seq], time.i[t.seq-1], units = "mins"))
}
plot(1:length(seq.t),seq.t)

## ---------------------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------------------
## End of Code