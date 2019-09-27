## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("0. Project Config.R")

## ------------------------------------------------------------------------------------------------------------------------------
##
##
## 
## ------------------------------------------------------------------------------------------------------------------------------

## Define region of interest

options(warn=-1)

landmass <- shapefile(landmass.shp)
crs(landmass) <- dt.projection
landmass <- gBuffer(landmass, byid=TRUE, width=0)

coastline <- shapefile(coastline.shp)
crs(coastline) <- dt.projection

clipper <- as(extent(min.lon,max.lon,min.lat,max.lat), "SpatialPolygons")
crs(clipper) <- dt.projection

landmass <- gIntersection(landmass, clipper, byid=TRUE)
coastline <- gIntersection(coastline, clipper, byid=TRUE)

## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

## Define source and sink locations

coastline.pts <- as(coastline, "SpatialPointsDataFrame")
coastline.pts <- as.data.frame(coastline.pts)[,c("x","y")]

coastline.pts <- trim.by.distance(coastline.pts,source.sink.dist,TRUE)

if(   unwanted.release.coastline ){ source.sink.xy <- data.frame(cells.id=1:nrow(coastline.pts),x=coastline.pts[,1],y=coastline.pts[,2],source=0,stringsAsFactors = FALSE) }
if( ! unwanted.release.coastline ){ source.sink.xy <- data.frame(cells.id=1:nrow(coastline.pts),x=coastline.pts[,1],y=coastline.pts[,2],source=1,stringsAsFactors = FALSE) }

## --------------

if( ! is.null(additional.sourcesink.shp) ) {
  
  additional.pts <- data.frame()
  
  for(i in 1:length(additional.sourcesink.shp)){
    
    if(exists("additional.shp")) { rm(additional.shp)} 
    
    additional.shp <- shapefile(paste0(project.folder,additional.sourcesink.shp[i]))

    if(class(additional.shp) == "SpatialPolygons") { additional.shp <- as(additional.shp, "SpatialLines") }
    if(class(additional.shp) == "SpatialPolygonsDataFrame") { additional.shp <- as(additional.shp, "SpatialLinesDataFrame") }
    
    if(class(additional.shp) == "SpatialPointsDataFrame") { 
      
      crs(additional.shp) <- crs(landmass)
      over.land <- over(additional.shp,landmass)
      additional.shp <- additional.shp[which(is.na(over.land)),]
      
      }
    
    additional.shp <- as.data.frame(as(additional.shp, "SpatialPointsDataFrame"))
    additional.shp <- data.frame(x=additional.shp[,ncol(additional.shp)-1],y=additional.shp[,ncol(additional.shp)])
    additional.pts <- rbind(additional.pts,additional.shp)

  }

  additional.pts <- trim.by.distance(additional.pts,source.sink.dist,TRUE)

  additional.source.sink.xy <- data.frame(cells.id=(nrow(coastline.pts)+1):(nrow(additional.pts)+nrow(coastline.pts)),x=additional.pts[,1],y=additional.pts[,2],source=1,stringsAsFactors = FALSE) ; head(additional.source.sink.xy)
  source.sink.xy <- rbind(source.sink.xy,additional.source.sink.xy)

}

## --------------

source.sink.xy <- source.sink.xy[source.sink.xy$x >= extent(clipper)[1] & source.sink.xy$x <= extent(clipper)[2] & source.sink.xy$y >= extent(clipper)[3] & source.sink.xy$y <= extent(clipper)[4],]
head(source.sink.xy)

## -----------------------------------------------------

## Remove unwanted release sites

if( ! is.null(unwanted.release.sites.shp) ) {
  
  source.sink.xy.t <- source.sink.xy[,2:3]
  coordinates(source.sink.xy.t) <- c("x","y")
  crs(source.sink.xy.t) <- dt.projection

  unwanted <- shapefile(paste0(project.folder,unwanted.release.sites.shp))
  unwanted <- as(unwanted,"SpatialPolygons")
  
  points.over.polygon <- as.vector(which( ! is.na( sp::over( source.sink.xy.t , unwanted , fn = NULL) )) )
  source.points <- as.vector(which( is.na( sp::over( source.sink.xy.t , unwanted , fn = NULL) )))
  
  if( length(points.over.polygon) > 0 ) {
    
    source.sink.xy <- rbind( data.frame(cells.id=1:length(source.points),x=source.sink.xy[source.points,2],y=source.sink.xy[source.points,3],source=1) ,
                             data.frame(cells.id=(length(source.points)+1):(length(source.points)+length(points.over.polygon)),x=source.sink.xy[points.over.polygon,2],y=source.sink.xy[points.over.polygon,3],source=0) )
  }
 
}
 
options(warn=0)

## ------------------

plot(landmass,box=FALSE,legend=FALSE,col=c("grey"))
lines(coastline,col=c("yellow"))

points(source.sink.xy[source.sink.xy$source == 1,2:3],col=c("black"),pch=16)
points(source.sink.xy[source.sink.xy$source == 0,2:3],col=c("yellow"))

## ------------------

initial.coords <- source.sink.xy[source.sink.xy$source == 1 , 2:3 ]
source.cells.id <- source.sink.xy[source.sink.xy$source == 1,1]

## ------------------------------------------------------------------------------------------------------------------------------

## Define Available raw data and Test if data is available

raw.data.files <- list.files(paste0(project.folder,"/Data"),full.names = TRUE,pattern="nc")
simulation.parameters.step <- data.frame()

for( file in 1:length(raw.data.files)) {
  
  nc <- nc_open( raw.data.files[file] , verbose=FALSE )
  nc.date <- as.Date(ncvar_get( nc, "Time"), origin = "1970-01-01")
  nc_close( nc )
  
  simulation.parameters.step <- rbind( simulation.parameters.step, data.frame( simulation=file, file=raw.data.files[file],  year=substr(nc.date, 1, 4) , month=substr(nc.date, 6, 7) , day=substr(nc.date, 9, 10) , stringsAsFactors = FALSE) )

}

simulation.parameters.step <- simulation.parameters.step[simulation.parameters.step$year %in% as.character(from.year:to.year) & simulation.parameters.step$day %in% sapply(from.day:to.day,function(x){ ifelse(nchar(x) > 1,x,paste0("0",x))}) & simulation.parameters.step$month %in% sapply(months.all,function(x){ ifelse(nchar(x) > 1,x,paste0("0",x))}) ,  ]

if(sum( ! from.year:to.year %in% as.numeric(simulation.parameters.step[,"year"])) + sum( ! months.all %in% as.numeric(simulation.parameters.step[,"month"])) > 0 ) { stop("Data is not available for time window")}

## ------------------------------------------------------------------------------------------------------------------------------
## Prepare video (animation) points

if( !is.null(movie.sites.xy) ) {
  
  if(class(movie.sites.xy) =="character") { movie.sites.xy <- as.data.frame(shapefile(paste0(project.folder,movie.sites.xy)))[,2:3] }
  
  movie.sites.xy <- movie.sites.xy[complete.cases(movie.sites.xy),]
  movie.sites.xy <- sort( as.vector(get.knnx( source.sink.xy[,c("x","y") ] , movie.sites.xy , k = 1 + movie.sites.buffer , algorithm="kd_tree" )$nn.index) )
  movie.sites.id <- unique(movie.sites.xy)
  movie.sites.xy <- source.sink.xy[ unique(movie.sites.xy) ,c("x","y") ]
  
  points(movie.sites.xy,col=c("green"),pch=16)
  
}

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
particles.reference[ , start.cell := as.numeric( sapply( source.sink.xy[source.sink.xy$source == 1 , 1 ] ,function(x) { rep(x,n.particles.per.cell) })) ]
particles.reference[ , start.year := 0 ]
particles.reference[ , start.month := 0 ]
particles.reference[ , start.day := 0 ]
particles.reference[ , pos.lon := 0 ]
particles.reference[ , pos.lat := 0 ]
particles.reference[ , pos.alt := 0 ]
particles.reference[ , state := 0 ]
particles.reference[ , t.start := 0 ]
particles.reference[ , t.finish := 0 ]
particles.reference[ , cell.rafted := 0 ]

particles.reference.names <- colnames(particles.reference)

setkey(particles.reference, id )

for( c in source.sink.xy[source.sink.xy$source == 1 , 1 ] ) {
  
  particles.reference[ start.cell == c, pos.lon := source.sink.xy[source.sink.xy$cells.id == c , 2 ] ]
  particles.reference[ start.cell == c, pos.lat := source.sink.xy[source.sink.xy$cells.id == c , 3 ] ]
  
}

head(particles.reference)
nrow(particles.reference)

## --------------------------------------------------------

## Out of memory objects

clean.dump.files(clean.dump.files=TRUE,files="raw.data.",dump.folder=paste0(project.folder,"/InternalProc"))
clean.dump.files(clean.dump.files=TRUE,files="particles.reference.",dump.folder=paste0(project.folder,"/InternalProc"))
clean.dump.files(clean.dump.files=TRUE,files="particles.video.location.",dump.folder=paste0(project.folder,"/InternalProc"))

particles.reference.bm <- big.matrix(nrow=nrow(particles.reference),ncol=ncol(particles.reference) , backingpath=paste0(project.folder,"InternalProc") , backingfile = "particles.reference.bin", descriptorfile = "particles.reference.desc")
particles.reference.bm.desc <- dget( paste0(project.folder,"InternalProc/particles.reference.desc"))

particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
particles.reference.bm[,1] <- unlist(particles.reference[,1])
particles.reference.bm[,2] <- unlist(particles.reference[,2])
particles.reference.bm[,6] <- unlist(particles.reference[,6])
particles.reference.bm[,7] <- unlist(particles.reference[,7])

## --------------------------------------------------------

## Data table to alocate path of particles (video)

if( movie.year > 0 ) {
  
  particles.to.sql.id <- particles.reference[ start.cell %in% movie.sites.id , id ]
  
  particles.video.location.x.bm <- filebacked.big.matrix( nrow = length(particles.to.sql.id), 
                                                          ncol = ( sum(simulation.parameters.step[,3] == movie.year) * n.hours.per.day ), 
                                                          backingfile = "particles.video.location.x.bin",
                                                          descriptorfile = "particles.video.location.x.desc",
                                                          backingpath=paste0(project.folder,"/InternalProc"))
                                    
  particles.video.location.y.bm <- filebacked.big.matrix( nrow = length(particles.to.sql.id), 
                                                          ncol = ( sum(simulation.parameters.step[,3] == movie.year) * n.hours.per.day ), 
                                                          backingfile = "particles.video.location.y.bin",
                                                          descriptorfile = "particles.video.location.y.desc",
                                                          backingpath=paste0(project.folder,"/InternalProc"))

  particles.video.location.x.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.video.location.x.desc"))
  particles.video.location.y.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.video.location.y.desc"))
  
}

## ------------------------------------------------------------------------------------------------------------------

## Generate regions for simulation
## Parallel.computational.sections : latitudinal section
 
sections.lat <- data.frame( sect.from = seq(min.lat,max.lat,length.out = parallel.computational.sections+1)[-(parallel.computational.sections+1)] , 
                            sect.to = seq(min.lat,max.lat,length.out = parallel.computational.sections+1)[-1] )

## Generate polygons defining land regions for simulation

gClip <- function(shp, bb){
       if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
       else b_poly <- as(extent(bb), "SpatialPolygons")
       gIntersection(shp, b_poly, byid = T)
 }

list.of.polygons <- character()

## -------------------

if( ! is.null(landmass.shp.2) ) { 
  
  landmass <- shapefile(paste0(project.folder,landmass.shp.2)) 
  clipper <- as(extent(min.lon,max.lon,min.lat,max.lat), "SpatialPolygons")
  crs(clipper) <- dt.projection
  landmass <- gIntersection(landmass, clipper, byid=TRUE)
  
}

## -------------------

for(i in 1:parallel.computational.sections){
  
  cat('\014')
  cat('\n')
  cat('\n')
  cat('\n Progress:', 100 - round((length(parallel.computational.sections) / i) * 100 , digits = 2) , '%' )
  
  clipper <- as(extent(min.lon,max.lon, sections.lat[i,1] - parallel.computational.buffer , sections.lat[i,2] + parallel.computational.buffer ), "SpatialPolygons")
  crs(clipper) <- dt.projection 

  try( geometry.i <- gClip(landmass, sp::bbox(clipper)) , silent = TRUE)
  if( class(geometry.i)[1] ==  "SpatialCollections" ) { geometry.i <- geometry.i@polyobj }
  if( class(geometry.i)[1] !=  "SpatialPolygons" ) { geometry.i <- crop( landmass,clipper ) }

  if( is.null(geometry.i) ) { 
    
    fakePoint <- as.matrix(data.frame(Lon=extent(clipper)[1],Lat=extent(clipper)[3]))
    fakePoint <- mapview::coords2Polygons(fakePoint,ID=1)
    crs(fakePoint) <- dt.projection 
    
    fakeLandmass <- raster::aggregate(rbind(landmass, fakePoint))
    geometry.i <- gClip(fakeLandmass, sp::bbox(clipper))
    
    }

  assign( paste0("landmass.sect.",i) , geometry.i )
  list.of.polygons <- c(list.of.polygons,paste0("landmass.sect.",i))
  gc(reset=TRUE)

}

## ------------------------------------------------------------------------------------------------------------------

## SQL configuration

if(! exists("movie.sites.id")) { movie.sites.id <- NULL}
if(! exists("particles.to.sql.id")) { particles.to.sql.id <- NULL}

## -----------------

global.simulation.parameters <- data.frame(   project.name = project.name,
                                              sim.years = paste(c(from.year,to.year),collapse="-"),
                                              sim.months = paste(months.all,collapse=","),
                                              kill.by.raft = kill.by.raft , 
                                              n.hours.per.day = n.hours.per.day , 
                                              n.new.particles.per.day = n.new.particles.per.day , 
                                              remove.new.particles.last.days = remove.new.particles.last.days , 
                                              longevity = longevity , 
                                              particle.max.duration = particle.max.duration , 
                                              behaviour = behaviour,
                                              n.particles.per.cell = n.particles.per.cell,
                                              movie.year = movie.year, 
                                              # movie.sites.id = paste(movie.sites.id,collapse=",") , 
                                              # particles.to.sql.id = paste(particles.to.sql.id,collapse=",") , 
                                              extent = paste(c(min.lon,max.lon,min.lat,max.lat),collapse=",") )       

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
dbWriteTable(sql, "SourceSinkSites", as.data.frame(source.sink.xy)  , append=FALSE, overwrite=TRUE )
dbWriteTable(sql, "Parameters", global.simulation.parameters , append=FALSE, overwrite=TRUE )
dbDisconnect(sql)

## -----------------------

rm(landmass) ; rm(coastline) ; rm(particles.reference.bm) ; rm(particles.reference)
gc(reset=TRUE)
gc()
memory.profile()
list.memory()

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------

# save.image(file='../Environment.RData')
# rm(list=(ls()[ls()!="v"])); gc(reset=TRUE); load('../Environment.RData')

## -------------------------------------------

## Start Simulation
## 1:nrow(simulation.parameters.step)

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
  
  # Next day 
  
  if( simulation.step < nrow(simulation.parameters.step) ) {
    
    simulation.year.n <- simulation.parameters.step[simulation.step +1,"year"]
    simulation.month.n <- simulation.parameters.step[simulation.step +1,"month"]
    simulation.day.n <- simulation.parameters.step[simulation.step +1,"day"]
    simulaton.raw.data.file.n <- simulation.parameters.step[simulation.step +1,"file"]
    
  } else {
    
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
    Time <- ncvar_get( simulation.raw.data, "Time" )
    Time <- as.Date(Time, origin = "1970-01-01") 
    nc_close(simulation.raw.data)
    
    if( ! as.numeric(format(Time, "%Y"))[1] %in% from.year:to.year ) { stop("Code 01: Something is wrong with raw data extraction") }
    
  }
  
  ## --------------------------------------------------------
  
  nc.file <- nc_open( simulaton.raw.data.file, readunlim=FALSE )
  time.start.day <- ncvar_get( nc.file, "Time" )
  time.start.day <- as.Date(time.start.day, origin = "1970-01-01") 
  nc_close(nc.file)
  
  nc.file <- nc_open( simulaton.raw.data.file.n, readunlim=FALSE )
  time.next.day <- ncvar_get( nc.file, "Time" )
  time.next.day <- as.Date(time.next.day, origin = "1970-01-01") 
  nc_close(nc.file)
  
  rd.start.day <- which( format(time.start.day, "%Y") == simulation.year & format(time.start.day, "%m") == simulation.month & format(time.start.day, "%d") == simulation.day )
  rd.next.day <- which( format(time.next.day, "%Y") == simulation.year.n & format(time.next.day, "%m") == simulation.month.n & format(time.next.day, "%d") == simulation.day.n )
  
  if( length(rd.start.day) == 0 | length(rd.next.day) == 0 ) { stop("Code 01: Something is wrong with raw data extraction") }
  
  ## --------------------------------------------------------
  
  if( ! "InternalProc" %in% list.files(project.folder) ) { dir.create(file.path(paste0(project.folder,"/InternalProc"))) }
  
  ## ------------------------------------- 
  
  padlock(paste0(project.folder,"/InternalProc/"),"Unlock",-1)
  
  ## Divide computations by sections (parallel.computational.sections)
  
  gc(reset=TRUE)
  
  cl.2 <- makeCluster(number.cores)
  registerDoParallel(cl.2)
  
  sect.loop <- foreach(section=1:parallel.computational.sections, .combine=rbind, .verbose=FALSE, .export=list.of.polygons , .packages=c("ncdf4","gstat","gdata","raster","data.table","bigmemory","FNN")) %dopar% { 
    
    require(bigmemory)
    
    ## -------------------------------------------------------------------------
    
    sp.poly <- get(paste0("landmass.sect.",section))
    
    if( class(sp.poly) != "SpatialPolygons" &  class(sp.poly) != "NULL"  ) { sp.poly <- get(paste0("landmass.sect.",section))@polyobj }
    if( class(sp.poly) != "NULL"  ) { sp.poly.line <- as(sp.poly, "SpatialLines")  }
    
    ## -------------------------------------------------------------------------
    
    sections.lat.f.s <- as.numeric(sections.lat[section,1])
    sections.lat.t.s <- as.numeric(sections.lat[section,2])
    
    ## -------------------------------------------------------------------------
    
    source.sink.xy.id.s <- which(source.sink.xy$y >= sections.lat.f.s - parallel.computational.buffer & source.sink.xy$y < sections.lat.t.s + parallel.computational.buffer)
    source.sink.xy.s <- source.sink.xy[source.sink.xy.id.s,2:3]
    source.sink.xy.id.s <- source.sink.xy[source.sink.xy.id.s,1]
    
    ## -------------------------------------------------------------------------
    
    ## U & V Components
    
    velocity.field.sec <- which(Latitude >= min(c(sections.lat.t.s,sections.lat.f.s)) - parallel.computational.buffer & Latitude <= max(c(sections.lat.t.s,sections.lat.f.s)) + parallel.computational.buffer )
    velocity.field.sec <- velocity.field.sec[-length(velocity.field.sec)]
    
    norm.field <- expand.grid( i=dim.i , j=velocity.field.sec )
    raw.data.coords <- cbind( apply( norm.field , 1 , function (x) Longitude[x[1]] ) , apply( norm.field , 1 , function (x) Latitude[x[2]] ) )
    colnames(raw.data.coords) <- c("Lon","Lat")
    
    nc.file <- nc_open( simulaton.raw.data.file, readunlim=FALSE )
    velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(1,velocity.field.sec[1],rd.start.day) , count=c(length(dim.i),length(velocity.field.sec),1) )
    start.day.raw.data.u <- melt(velocity.field)[,"value"]
    velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(1,velocity.field.sec[1],rd.start.day) , count=c(length(dim.i),length(velocity.field.sec),1) )
    start.day.raw.data.v <-  melt(velocity.field)[,"value"]
    nc_close(nc.file)
    
    start.day.raw.data.u[start.day.raw.data.u > 100] <- NA
    start.day.raw.data.v[start.day.raw.data.v > 100] <- NA
    start.day.raw.data.u[start.day.raw.data.u < -100] <- NA
    start.day.raw.data.v[start.day.raw.data.v < -100] <- NA
    
    nc.file <- nc_open( simulaton.raw.data.file.n, readunlim=FALSE )
    velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(1,velocity.field.sec[1],rd.next.day) , count=c(length(dim.i),length(velocity.field.sec),1) )
    next.day.raw.data.u <-  melt(velocity.field)[,"value"]
    velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(1,velocity.field.sec[1],rd.next.day) , count=c(length(dim.i),length(velocity.field.sec),1) )
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
    
    # plot(poly)
    # plot(raw.data.u[,c("Lon","Lat")])
    
    # --------------------------------------------------
    
    speed.u.sec <- t( sapply( 1:nrow(raw.data.u) , function (x) seq( from = raw.data.u[x,3] , to = raw.data.u[x,4] , length.out = n.hours.per.day + 1 ) ))
    speed.v.sec <- t( sapply( 1:nrow(raw.data.v) , function (x) seq( from = raw.data.v[x,3] , to = raw.data.v[x,4] , length.out = n.hours.per.day + 1 ) ))
    
    speed.u.sec.coords <- raw.data.u[,1:2]
    speed.v.sec.coords <- raw.data.v[,1:2]
    
    ## --------------------------------------------------------
    ## --------------------------------------------------------
    
    particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
    
    # -----------------------------------------------
    
    particles.reference.moving.i <- mwhich(particles.reference.bm.all,c(7,7,9),list(sections.lat.f.s,sections.lat.t.s,1), list('ge', 'lt','eq') , 'AND')
    particles.reference.moving.dt <- data.table(matrix(particles.reference.bm.all[particles.reference.moving.i, ],ncol=length(particles.reference.names)))
    colnames(particles.reference.moving.dt) <- particles.reference.names
    
    ## --------------------------------------------------------
    ## Move particles (per day)
    
    for( h in 1:n.hours.per.day ) {
      
      t.step <- ((as.numeric(simulation.step)-1) * n.hours.per.day) + h
      
      # -----------------------------------------------
      # Release new particles, if that is the case
      
      if( norm.time[t.step,release.particles] ) {   
        
        particles.reference.new.i <- mwhich(particles.reference.bm.all,c(7,7,9),list(sections.lat.f.s,sections.lat.t.s,0), list('ge', 'lt','eq') , 'AND')
        particles.reference.new.i <- data.table(particles.reference.bm.all[particles.reference.new.i,c(1,2)])
        particles.reference.new.i <- particles.reference.new.i[ , .SD[1] , by=V2][,V1]
        
        # particles.reference.names
        
        particles.reference.new <- data.table(particles.reference.bm.all[particles.reference.new.i, ])
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
      
      # --------
      
      if( moving.particles.condition ) {
        
        setkey(particles.reference.moving.dt,id)
        moving.particles.id <- as.vector(unlist(particles.reference.moving.dt[, "id"]))
        
        points.to.interp <- particles.reference.moving.dt[, .(pos.lon,pos.lat)]
        
        # plot(speed.u.sec.coords) ; points(points.to.interp,col="Red")
        
        source.points.to.interp.u <- data.frame(x=speed.u.sec.coords[,1],y=speed.u.sec.coords[,2],var=speed.u.sec[,h])
        source.points.to.interp.v <- data.frame(x=speed.v.sec.coords[,1],y=speed.v.sec.coords[,2],var=speed.v.sec[,h])
        
        # Interpolate speed
        
        idwPower <- 2
        idwN <- 3
        
        idw.nearest.i <- get.knnx( source.points.to.interp.v[,1:2] , points.to.interp, k=idwN , algorithm="kd_tree" )$nn.index
        idw.nearest.d <- get.knnx( source.points.to.interp.v[,1:2] , points.to.interp, k=idwN , algorithm="kd_tree" )$nn.dist
        idw.nearest.d <- idw.nearest.d * 1e9
        
        speed.u <- numeric(nrow(points.to.interp))
        speed.v <- numeric(nrow(points.to.interp))
        
        for( pt.i in 1:nrow(points.to.interp)) {
          
          speed.u[pt.i] <- (sum( source.points.to.interp.u[idw.nearest.i[pt.i,],3] / idw.nearest.d[pt.i,]^idwPower , na.rm=T)) / (sum( 1 / idw.nearest.d[pt.i,]^idwPower))
          speed.v[pt.i] <- (sum( source.points.to.interp.v[idw.nearest.i[pt.i,],3] / idw.nearest.d[pt.i,]^idwPower , na.rm=T)) / (sum( 1 / idw.nearest.d[pt.i,]^idwPower))
          
        }
        
        if( max(speed.u) > 100 | min(speed.u) < -100 | max(speed.v) > 100 | min(speed.v) < -100 ) { stop("Error [!]") }
        
        # coordinates(points.to.interp) = ~pos.lon+pos.lat
        # coordinates(source.points.to.interp.u) = ~x+y
        # coordinates(source.points.to.interp.v) = ~x+y
        
        # speed.u <- idw(formula = var ~ 1, source.points.to.interp.u, points.to.interp, nmax = idwN, idp=idwPower)$var1.pred
        # speed.v <- idw(formula = var ~ 1, source.points.to.interp.v, points.to.interp, nmax = idwN, idp=idwPower)$var1.pred
        
        mov.eastward <- speed.u * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
        mov.northward <- speed.v * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
        
        # Assign temporary positions
        
        points.to.interp <- particles.reference.moving.dt[, .(pos.lon,pos.lat,id)]
        
        dLon <- mov.eastward / ( 6378137 * cos( pi * (  points.to.interp[,2]  / 180) ) )
        dLat <- mov.northward / 6378137
        dLon <- points.to.interp[,1] + dLon * 180/pi 
        dLat <- points.to.interp[,2] + dLat * 180/pi
        
        # ---------------------------
        
        particles.reference.moving.old.pos <- points.to.interp
        
        # ---------------------------
        
        setkey(particles.reference.moving.dt,id)
        particles.reference.moving.dt[,6] <- dLon
        particles.reference.moving.dt[,7] <- dLat
        
        # -----------------------------------------------
        # Out of space (study region), if TRUE, place particles on hold
        
        out.of.space <- dLon > max.lon | dLon < min.lon | dLat > max.lat | dLat < min.lat
        out.of.space.ids <- moving.particles.id[out.of.space]
        
        setkey(particles.reference.moving.dt,id)
        particles.reference.moving.dt[ id %in% out.of.space.ids , state := 3  ]
        
        # -----------------------------------------------
        # kill by first raft . Will eliminate particles that got to another cell - first raft event
        
        if( kill.by.raft & ! is.null(sp.poly) ) {
          
          setkey(particles.reference.moving.dt,id)
          points.to.test <- particles.reference.moving.dt[, .(pos.lon,pos.lat)]
          coordinates(points.to.test) = ~pos.lon+pos.lat
          crs(points.to.test) <- dt.projection
          
          particles.on.land <- as.vector(which(!is.na(over(points.to.test,sp.poly))))
          particles.on.land.condition <- length(particles.on.land) > 0
          
          if( particles.on.land.condition ) {    
            
            who.at.land.id <- moving.particles.id[particles.on.land]
            
            cells.started <- as.vector(unlist(particles.reference.moving.dt[id %in% who.at.land.id, "start.cell"]))
            who.at.land.t.start <- as.vector(unlist(particles.reference.moving.dt[id %in% who.at.land.id, "t.start"]))
            
            points.on.land.t <- particles.reference.moving.dt[id %in% who.at.land.id , 6:7 ]
            points.on.land.t.minus <- particles.reference.moving.old.pos[id %in% who.at.land.id , 1:2 ]
            
            points.on.land.corrected <- points.on.land.t.minus
            
            cells.rafted <- get.knnx( source.sink.xy.s, points.on.land.corrected, k=1 , algorithm="kd_tree" )$nn.index
            cells.rafted <- source.sink.xy.id.s[cells.rafted]
            
            displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]} )
            
            # For Rafters
            
            true.rafters.id <- who.at.land.id[ displacement != 0 ]
            true.rafters.cell <- cells.rafted[ displacement != 0 ]
            
            setkey(particles.reference.moving.dt,id)
            particles.reference.moving.dt[ id %in% true.rafters.id , state := 2 ]
            particles.reference.moving.dt[ id %in% true.rafters.id , cell.rafted := as.numeric(true.rafters.cell) ]
            particles.reference.moving.dt[ id %in% true.rafters.id , t.finish := as.numeric(t.step) ]
            
            # For non-Rafters New Particles
            
            non.rafters.id <- who.at.land.id[ displacement == 0 ]
            non.rafters.cell <- cells.started[ displacement == 0 ]
            non.rafters.t <- who.at.land.t.start[ displacement == 0 ] == t.step
            
            if( TRUE %in% non.rafters.t ) {    
              
              particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , pos.lon := particles.reference.moving.old.pos[ id %in% non.rafters.id[non.rafters.t] ,1] ]
              particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , pos.lat := particles.reference.moving.old.pos[ id %in% non.rafters.id[non.rafters.t]  ,2] ]
              
            }
            
            # For non-Rafters Old Particles
            
            non.rafters.t <- who.at.land.t.start[displacement == 0] < t.step
            
            if( TRUE %in% non.rafters.t ) {    
              
              particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , state := 2 ]
              particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , cell.rafted := as.numeric(non.rafters.cell[non.rafters.t] ) ]
              particles.reference.moving.dt[ id %in% non.rafters.id[non.rafters.t] , t.finish := as.numeric(t.step) ]
              
            }
          }
          
        }
        
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
      
      if( movie.year == as.numeric(simulation.year) ) {
        
        ## --------------------------------------------------------
        
        particles.video.location.x.bm.i <- attach.big.matrix(particles.video.location.x.bm.desc)
        particles.video.location.y.bm.i <- attach.big.matrix(particles.video.location.y.bm.desc)
        
        ## --------------------------------------------------------
        
        setkey(particles.reference.moving.dt,id)
        particles.to.sql.id.moving <- particles.reference.moving.dt[state == 1,id]
        particles.to.sql.id.moving <- particles.to.sql.id.moving[particles.to.sql.id.moving %in% particles.to.sql.id]
        particles.to.sql.id.moving.condition <- length(particles.to.sql.id.moving) > 0
        
        if( particles.to.sql.id.moving.condition ) { 
          
          t.step.movie <- ((as.numeric(which( which(as.numeric(simulation.parameters.step[,3]) == movie.year) == simulation.step))-1) * n.hours.per.day) + h
          
          setkey(particles.reference.moving.dt,id)
          particles.video.location.x.bm.i[ which(particles.to.sql.id %in% particles.to.sql.id.moving) , t.step.movie ] <- unlist(particles.reference.moving.dt[ id %in% particles.to.sql.id.moving,pos.lon])
          particles.video.location.y.bm.i[ which(particles.to.sql.id %in% particles.to.sql.id.moving) , t.step.movie ] <- unlist(particles.reference.moving.dt[ id %in% particles.to.sql.id.moving,pos.lat])
          
        }
      }
      
    }
    
    ## ---------------------------------------------------
    
    setkey(particles.reference.moving.dt,id)
    gc(reset=T)
    
    return(particles.reference.moving.dt)
    
    # -----------------------------------------------
    
    
  }
  
  stopCluster(cl.2) ; rm(cl.2)
  
  ## -------------------------------------
  ## Inject particles to final object [ particles.reference.bm.all ] 
  
  particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
  
  sect.loop <- sect.loop[ sect.loop$state != 0 , ]

  particles.reference.bm.all[ sect.loop$id , 3 ] <- as.numeric(unlist(sect.loop[  , "start.year"] ))
  particles.reference.bm.all[ sect.loop$id , 4 ] <- as.numeric(unlist(sect.loop[  , "start.month"] ))
  particles.reference.bm.all[ sect.loop$id , 5 ] <- as.numeric(unlist(sect.loop[  , "start.day"] ))
  particles.reference.bm.all[ sect.loop$id , 6 ] <- as.numeric(unlist(sect.loop[  , "pos.lon"] ))
  particles.reference.bm.all[ sect.loop$id , 7 ] <- as.numeric(unlist(sect.loop[  , "pos.lat"] ))
  particles.reference.bm.all[ sect.loop$id , 9 ] <- as.numeric(unlist(sect.loop[  , "state"] ))
  particles.reference.bm.all[ sect.loop$id , 10 ] <- as.numeric(unlist(sect.loop[  , "t.start"] ))
  particles.reference.bm.all[ sect.loop$id , 11 ] <- as.numeric(unlist(sect.loop[  , "t.finish"] ))
  particles.reference.bm.all[ sect.loop$id , 12 ] <- as.numeric(unlist(sect.loop[  , "cell.rafted"] ))
  
  gc(reset=TRUE)
  
}

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

## -----------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
dbWriteTable(sql, "ReferenceTable", ReferenceTable , append=FALSE, overwrite=TRUE )
dbDisconnect(sql)

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