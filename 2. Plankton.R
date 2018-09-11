## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

source("Dependences.R")

## ------------------------------------------------------------------------------------------------------------------------------
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

## -----------------------------------------------------

## Define source and sink locations

coastline.pts <- as(coastline, "SpatialPointsDataFrame")
coastline.pts <- as.data.frame(coastline.pts)[,c("x","y")]

coastline.pts.t <- coastline.pts
source.sink.xy <- data.frame()

while( nrow(coastline.pts.t) > 0 ){
  
  pt.i = coastline.pts.t[1,,drop=FALSE]
  
  source.sink.xy <- rbind(source.sink.xy,as.data.frame(pt.i))
                          
  coastline.pts.t.i <- coastline.pts.t
  coordinates(coastline.pts.t.i) <- c("x","y")
  crs(coastline.pts.t.i) <- dt.projection
  circle <- circles(pt.i, lonlat=TRUE, d=source.sink.dist*1000, dissolve=FALSE)
  circle <- geometry(circle)
  crs(circle) <- dt.projection
  
  to.extract <- which(!is.na(over(coastline.pts.t.i,circle)))
  coastline.pts.t <- coastline.pts.t[-to.extract,]
  
}

source.sink.xy <- data.frame(cells.id=1:nrow(source.sink.xy),x=source.sink.xy[,1],y=source.sink.xy[,2],source=1) ; head(source.sink.xy)

## -----------------------------------------------------

## Remove unwanted release sites

if( ! is.null(unwanted.release.sites.shp) ) {
  
  source.sink.xy.t <- source.sink.xy[,2:3]
  coordinates(source.sink.xy.t) <- c("x","y")
  crs(source.sink.xy.t) <- dt.projection

  unwanted <- shapefile(unwanted.release.sites.shp)
  unwanted <- as(unwanted,"SpatialPolygons")
  
  points.over.polygon <- as.vector(which( ! is.na( sp::over( source.sink.xy.t , unwanted , fn = NULL) )) )
  source.points <- as.vector(which( is.na( sp::over( source.sink.xy.t , unwanted , fn = NULL) )))
  
  if( length(points.over.polygon) > 0 ) {
    
    source.sink.xy <- rbind( data.frame(cells.id=1:length(source.points),x=source.sink.xy[source.points,2],y=source.sink.xy[source.points,3],source=1) ,
                             data.frame(cells.id=(length(source.points)+1):(length(source.points)+length(points.over.polygon)),x=source.sink.xy[points.over.polygon,2],y=source.sink.xy[points.over.polygon,3],source=0)
                            )
    
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

if(sum( !from.year:to.year %in% as.numeric(simulation.parameters.step[,"year"])) + sum( ! months.all %in% as.numeric(simulation.parameters.step[,"month"])) > 0 ) { stop("Data is not available for time window")}

## ------------------------------------------------------------------------------------------------------------------------------
## Prepare video (animation) points

if( !is.null(movie.sites.xy) ) {
    if( ! class(movie.sites.xy) == "matrix") {   
            movie.sites.xy <- shapefile(movie.sites.xy)
            movie.sites.xy <-  crop(movie.sites.xy, extent(landmass) )
            movie.sites.xy <- as.data.frame(movie.sites.xy)[,2:3]
    } 
    else {  movie.sites.xy <- as.data.frame(movie.sites.xy) 
    }
  
    movie.sites.xy <- sort( as.vector(get.knnx( initial.coords , movie.sites.xy , k = 1 + movie.sites.buffer , algorithm="kd_tree" )$nn.index) )
}

if( !is.null(movie.sites.xy) ) { movie.sites.id <- unique(movie.sites.xy)
                                 movie.sites.xy <- initial.coords[ unique(movie.sites.xy) , ]
 
}

points(movie.sites.xy,col=c("red"),pch=16)

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
particles.reference[ , start.cell := as.numeric( sapply( 1:nrow(initial.coords) ,function(x) { rep(x,n.particles.per.cell) })) ]
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
setkey(particles.reference, id )

for( c in source.sink.id) {
  
  particles.reference[ start.cell == c, pos.lon := initial.coords$x[c] ]
  particles.reference[ start.cell == c, pos.lat := initial.coords$y[c] ]
  
}

nrow(particles.reference)

## Out of memory objects

clean.dump.files(clean.dump.files=TRUE,files="particles.reference.",dump.folder=paste0(project.folder,"/InternalProc"))

particles.reference.bm <- as.big.matrix(as.matrix(particles.reference) , backingpath=paste0(project.folder,"/InternalProc") , backingfile = "particles.reference.bin", descriptorfile = "particles.reference.desc")
particles.reference.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.reference.desc"))

## --------------------------------------------------------

## Data table to alocate path of particles (video)

if( ! is.null(movie.year) ) {
  
  particles.to.sql.id <- particles.reference[ start.cell %in% movie.sites.id , id ]
  
  particles.video.location.x <- matrix( 0 , nrow = length(particles.to.sql.id) , ncol = (nrow(simulation.parameters.step) * n.hours.per.day ) )
  particles.video.location.y <- matrix( 0 , nrow = length(particles.to.sql.id) , ncol = (nrow(simulation.parameters.step) * n.hours.per.day ) )
  particles.video.location.z <- matrix( 0 , nrow = length(particles.to.sql.id) , ncol = (nrow(simulation.parameters.step) * n.hours.per.day ) )

  clean.dump.files(clean.dump.files=TRUE,files="particles.video.location.",dump.folder=paste0(project.folder,"/InternalProc"))
  
  particles.video.location.x.bm <- as.big.matrix(particles.video.location.x , backingpath=paste0(project.folder,"/InternalProc") , backingfile = "particles.video.location.x.bin", descriptorfile = "particles.video.location.x.desc")
  particles.video.location.x.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.video.location.x.desc"))

  particles.video.location.y.bm <- as.big.matrix(particles.video.location.y , backingpath=paste0(project.folder,"/InternalProc") , backingfile = "particles.video.location.y.bin", descriptorfile = "particles.video.location.y.desc")
  particles.video.location.y.bm.desc <- dget( paste0(project.folder,"/InternalProc/particles.video.location.y.desc"))
  
}

## ------------------------------------------------------------------------------------------------------------------

## Generate regions for simulation
## Parallel.computational.sections : latitudinal section
 
sections.lat <- data.frame( sect.from = seq(min.lat,max.lat,length.out = parallel.computational.sections+1)[-(parallel.computational.sections+1)] , 
                            sect.to = seq(min.lat,max.lat,length.out = parallel.computational.sections+1)[-1] )

## Generate polygons defining land regions for simulation

list.of.polygons <- character()

for(i in 1:parallel.computational.sections){
  
  clipper <- as(extent(min.lon,max.lon, sections.lat[i,1] - parallel.computational.buffer , sections.lat[i,2] + parallel.computational.buffer ), "SpatialPolygons")
  crs(clipper) <- dt.projection 
  
  assign( paste0("landmass.sect.",i) , gIntersection(landmass, clipper, byid=TRUE) )
  
  list.of.polygons <- c(list.of.polygons,paste0("landmass.sect.",i))

}

## ------------------------------------------------------------------------------------------------------------------

## SQL configuration

if( paste0(project.name,"SimulationResults.sql") %in% list.files(sql.directory) ) {
  
  x <- ""
  while( x != "Y" & x != "n" ) { x <- readline("SQL database already exists. Do you which to overwrite? (Y/n) ") }
  
  if (x == "Y" ) { file.remove( paste0(sql.directory,"/",project.name,"SimulationResults.sql") ) }
}

## -----------------------

if( ! paste0(project.name,"SimulationResults.sql") %in% list.files(sql.directory) ) {
  
  global.simulation.parameters <- data.frame(   project.name = project.name,
                                                sim.years = from.year:to.year,
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
                                                movie.sites.id = paste(movie.sites.id,collapse=",") , 
                                                particles.to.sql.id = paste(particles.to.sql.id,collapse=",") , 
                                                extent = paste(c(min.lon,max.lon,min.lat,max.lat),collapse=",") )       
  
  sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
  dbWriteTable(sql, "SourceSinkSites", as.data.frame(source.sink.xy)  , append=FALSE, overwrite=TRUE )
  dbWriteTable(sql, "Parameters", global.simulation.parameters , append=FALSE, overwrite=TRUE )
  dbDisconnect(sql)
  
}

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------

## Start Simulation
## 1:nrow(simulation.parameters.step)

for ( simulation.step in 1:nrow(simulation.parameters.step) ) {
                
                ## --------------------------------------------------------
    
                ## Progress
  
                time.i <- Sys.time()
                if(simulation.step == 1) { time.f <- time.i}
                progress.percent <- round((simulation.step / nrow(simulation.parameters.step)) * 100)
                time.take.step.min <- round(as.numeric(difftime(time.i, time.f, units = "mins")))
                
                cat('\014')
                cat('\n')
                cat('\n Running step #',simulation.step,'| Time taken',time.take.step.min,'mins.')
                cat('\n',paste0(rep("-",100),collapse = ""))
                cat('\n --',paste0(rep("-",progress.percent),collapse = ""),"||",progress.percent,"%")
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
                  
                  norm.field <- expand.grid( i=dim.i , j=dim.j )
                  raw.data.coords <- cbind( apply( norm.field , 1 , function (x) Longitude[x[1]] ) , apply( norm.field , 1 , function (x) Latitude[x[2]] ) )
                  colnames(raw.data.coords) <- c("Lon","Lat")
                  
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
                
                ## U & V Components
                
                nc.file <- nc_open( simulaton.raw.data.file, readunlim=FALSE )
                velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(1,1,rd.start.day) , count=c(length(dim.i),length(dim.j),1) )
                start.day.raw.data.u <- melt(velocity.field)[,"value"]
                velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(1,1,rd.start.day) , count=c(length(dim.i),length(dim.j),1) )
                start.day.raw.data.v <-  melt(velocity.field)[,"value"]
                nc_close(nc.file)

                nc.file <- nc_open( simulaton.raw.data.file.n, readunlim=FALSE )
                velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(1,1,rd.next.day) , count=c(length(dim.i),length(dim.j),1) )
                next.day.raw.data.u <-  melt(velocity.field)[,"value"]
                velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(1,1,rd.next.day) , count=c(length(dim.i),length(dim.j),1) )
                next.day.raw.data.v <-  melt(velocity.field)[,"value"]
                nc_close(nc.file)
                
                # Crop environmental data by extent (study region)
                
                raw.data.u <- data.table(raw.data.coords,u.start=start.day.raw.data.u,u.next=next.day.raw.data.u)
                raw.data.v <- data.table(raw.data.coords,v.start=start.day.raw.data.v,v.next=next.day.raw.data.v)
                
                raw.data.u <- raw.data.u[ Lon >= min.lon & Lat >= min.lat & Lon <= max.lon & Lat <= max.lat & !is.na(u.start) , ]
                raw.data.v <- raw.data.v[ Lon >= min.lon & Lat >= min.lat & Lon <= max.lon & Lat <= max.lat & !is.na(v.start) , ]
                
                # points(raw.data.u[,.(Lon,Lat)])

                ## -------------------------------------
                
                ## Out of memory objects
                
                raw.data.u <- as.matrix(raw.data.u) ; colnames(raw.data.u) <- NULL
                raw.data.v <- as.matrix(raw.data.v) ; colnames(raw.data.v) <- NULL
                
                if( ! "InternalProc" %in% list.files(project.folder) ) { dir.create(file.path(paste0(project.folder,"/InternalProc"))) }
                
                clean.dump.files(clean.dump.files=TRUE,files="raw.data.",dump.folder=paste0(project.folder,"/InternalProc"))
              
                raw.data.u.bm <- as.big.matrix( raw.data.u , backingpath=paste0(project.folder,"/InternalProc") , backingfile = "raw.data.u.bin", descriptorfile = "raw.data.u.desc")
                raw.data.u.bm.desc <- dget( paste0(project.folder,"/InternalProc/raw.data.u.desc"))
                raw.data.v.bm <- as.big.matrix(raw.data.v , backingpath=paste0(project.folder,"/InternalProc") , backingfile = "raw.data.v.bin", descriptorfile = "raw.data.v.desc")
                raw.data.v.bm.desc <- dget( paste0(project.folder,"/InternalProc/raw.data.v.desc"))

                ## -------------------------------------
                
                ## Divide computations by sections (parallel.computational.sections)
                
                cl.2 <- makeCluster(number.cores)
                registerDoParallel(cl.2)
                
                sect.loop <- foreach(section=1:parallel.computational.sections, .verbose=FALSE, .export=list.of.polygons , .packages=c("gstat","gdata","raster","data.table","bigmemory")) %dopar% { 
                  
                      sections.lat.f.s <- as.numeric(sections.lat[section,1])
                      sections.lat.t.s <- as.numeric(sections.lat[section,2])
                  
                      ## --------------------------------------------------------
                      
                      particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
                      particles.reference.bm.i <- mwhich(particles.reference.bm.all,c(7,7),list(sections.lat.f.s,sections.lat.t.s), list('ge', 'lt') , 'AND')
                      particles.reference.bm.sec <- as.data.table(particles.reference.bm.all[particles.reference.bm.i,])
                      
                      ## --------------------------------------------------------
                      
                      raw.data.u.bm.sec <- attach.big.matrix(raw.data.u.bm.desc)
                      raw.data.v.bm.sec <- attach.big.matrix(raw.data.v.bm.desc)
                    
                      raw.data.u.i <- mwhich(raw.data.u.bm.sec,c(2,2),list(sections.lat.f.s-parallel.computational.buffer,sections.lat.t.s+parallel.computational.buffer), list('ge', 'le') , 'AND')
                      raw.data.v.i <- mwhich(raw.data.v.bm.sec,c(2,2),list(sections.lat.f.s-parallel.computational.buffer,sections.lat.t.s+parallel.computational.buffer), list('ge', 'le') , 'AND')
                      
                      speed.u.sec <- t( sapply( raw.data.u.i , function (x) seq( from = raw.data.u.bm.sec[x,3] , to = raw.data.u.bm.sec[x,4] , length.out = n.hours.per.day + 1 ) ))
                      speed.v.sec <- t( sapply( raw.data.v.i , function (x) seq( from = raw.data.v.bm.sec[x,3] , to = raw.data.v.bm.sec[x,4] , length.out = n.hours.per.day + 1 ) ))
                      
                      speed.u.sec.coords <- raw.data.u.bm.sec[raw.data.u.i,1:2]
                      speed.v.sec.coords <- raw.data.v.bm.sec[raw.data.v.i,1:2]
                    
                      ## --------------------------------------------------------
                      ## Move particles (per day)
                          
                      for( h in 1:n.hours.per.day ) {
    
                            t.step <- ((as.numeric(simulation.step)-1) * n.hours.per.day) + h
                            
                            # -----------------------------------------------
                            # Release new particles, if that is the case
                            
                            if( norm.time[t.step,release.particles] ) {   
                              
                                  new.particles.id <- particles.reference.bm.sec[ state == 0 , .SD[1] , by=start.cell][,id]
                                  new.particles.cells <- particles.reference.bm.sec[ state == 0 , .SD[1] , by=start.cell][,start.cell]
                                  particles.reference.bm.sec[ id %in% new.particles.id , state := 1 ]
                                  particles.reference.bm.sec[ id %in% new.particles.id , t.start := t.step ]
                                  particles.reference.bm.sec[ id %in% new.particles.id , start.year := as.numeric(simulation.year) ]
                                  particles.reference.bm.sec[ id %in% new.particles.id , start.month := as.numeric(simulation.month) ]
                                  particles.reference.bm.sec[ id %in% new.particles.id , start.day := as.numeric(simulation.day) ]
                                  
                                  
                            }
    
                            # -----------------------------------------------
                            # Which to move and speed
                            
                            moving.particles.xy <- particles.reference.bm.sec[ state == 1 , ]
                            moving.particles.condition <- nrow(moving.particles.xy) > 0
    
                            if( moving.particles.condition ) {
    
                                      moving.particles.ids <- moving.particles.xy[,id]
                                      moving.particles.start.cell <- moving.particles.xy[,start.cell]
                                      moving.particles.t.start <- moving.particles.xy[,t.start]
                                      
                                      points.to.interp <- moving.particles.xy[, .(pos.lon,pos.lat) ]
                                      coordinates(points.to.interp) = ~pos.lon+pos.lat
                                              
                                      source.points.to.interp.u <- data.frame(x=speed.u.sec.coords[,1],y=speed.u.sec.coords[,2],var=speed.u.sec[,h])
                                      coordinates(source.points.to.interp.u) = ~x+y
                                      
                                      source.points.to.interp.v <- data.frame(x=speed.v.sec.coords[,1],y=speed.v.sec.coords[,2],var=speed.v.sec[,h])
                                      coordinates(source.points.to.interp.v) = ~x+y
                          
                                      # Interpolate speed
    
                                      invisible( speed.u <- idw(formula = var ~ 1, source.points.to.interp.u, points.to.interp, nmax=3)$var1.pred )
                                      invisible( speed.v <- idw(formula = var ~ 1, source.points.to.interp.v, points.to.interp, nmax=3)$var1.pred )
                                      mov.eastward <- speed.u * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
                                      mov.northward <- speed.v * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
                              
                                      # Assign temporary positions
    
                                      points.to.interp <- as.data.frame(points.to.interp)
                                      
                                      dLon <- mov.eastward / ( 6378137 * cos( pi * (  points.to.interp[,2]  / 180) ) )
                                      dLat <- mov.northward / 6378137
                                      dLon <- points.to.interp[,1] + dLon * 180/pi 
                                      dLat <- points.to.interp[,2] + dLat * 180/pi
    
                                      setkey(moving.particles.xy,id)
                                      
                                      moving.particles.xy[ , pos.lon := dLon ]
                                      moving.particles.xy[ , pos.lat := dLat ]
                                      
                                      moving.particles.xy[ , old.pos.lon := points.to.interp[,1] ]
                                      moving.particles.xy[ , old.pos.lat := points.to.interp[,2] ]
                                      
                                      # -----------------------------------------------
                                      # Out of space (study region), if TRUE, place particles on hold
    
                                      out.of.space <- dLon > max.lon | dLon < min.lon | dLat > max.lat | dLat < min.lat
                                      out.of.space.ids <- moving.particles.ids[out.of.space]
                                      moving.particles.xy[ id %in% out.of.space.ids , state := 3 ] 
                              
                                      # -----------------------------------------------
                                      # kill by first raft . Will eliminate particles that got to another cell - first raft event
    
                                      if( kill.by.raft ) {
                                        
                                        points.to.test <- moving.particles.xy[ ,.(pos.lon,pos.lat) ]
                                        coordinates(points.to.test) <- c("pos.lon","pos.lat")
                                        crs(points.to.test) <- dt.projection
    
                                        particles.on.land <- as.vector(which(!is.na(over(points.to.test,get(paste0("landmass.sect.",section))))))
                                        particles.on.land.condition <- length(particles.on.land) > 0
    
                                        if( particles.on.land.condition ) {    
                                
                                                  who.at.land.id <- moving.particles.xy[ particles.on.land  , id ]
                                                  who.at.land.t.start <- moving.particles.xy[ particles.on.land , t.start ]
                                                  cells.started <- moving.particles.xy[ particles.on.land ,start.cell] 
    
                                                  dist.to.nearest.cell <- spDists(as.matrix(moving.particles.xy[ particles.on.land ,.(pos.lon,pos.lat) ]) , as.matrix(initial.coords) , longlat = TRUE)      
                                                  cells.rafted <- apply(dist.to.nearest.cell,1,which.min)
                                  
                                                  displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]})
                                                  
                                                  # For Rafters
                                                  
                                                  true.rafters.id <- who.at.land.id[ displacement != 0 ]
                                                  true.rafters.cell <- cells.rafted[ displacement != 0 ]
                                                  moving.particles.xy[ id %in% true.rafters.id , c("cell.rafted","state","t.finish") := list(true.rafters.cell,2,t.step) ] 
                                  
                                                  # For non-Rafters New Particles
                                                  
                                                  non.rafters.id <- who.at.land.id[ displacement == 0 ]
                                                  non.rafters.cell <- cells.started[ displacement == 0 ]
                                                  non.rafters.t <- who.at.land.t.start[ displacement == 0 ] == t.step
    
                                                  if( TRUE %in% non.rafters.t ) {    
    
                                                    moving.particles.xy[id %in% non.rafters.id[non.rafters.t], c("pos.lon","pos.lat") := .(old.pos.lon,old.pos.lat) ]
    
                                                  }
                                                                  
                                                  # For non-Rafters Old Particles
    
                                                  non.rafters.t <- who.at.land.t.start[displacement == 0] < t.step
                                                  
                                                  if( TRUE %in% non.rafters.t ) {    
    
                                                    moving.particles.xy[ id %in% non.rafters.id[non.rafters.t] , c("cell.rafted","state","t.finish") := list( non.rafters.cell[non.rafters.t] , 2 , t.step ) ] 
                                                    
                                                  }
    
                                        }
                                        
                                      }
                                      
                            }
    
                            # -----------------------------------------------
                            # End of day, Kill by Longevity
                            
                            if ( longevity ) {   
                              
                              max.duration.id <- moving.particles.xy[ ( t.step - t.start ) > ( n.hours.per.day * particle.max.duration ) , id ]
                              moving.particles.xy[ id %in% max.duration.id , c("state") := 4 ]
                              
                            }
                            
                            ## ---------------------------------------------------
                            # Inject particles to temporary object [ particles.reference.bm.sec ] 
                            
                            setkey(particles.reference.bm.sec, id )
                            
                            particles.reference.bm.sec[ id %in% moving.particles.xy[,id] , pos.lon := moving.particles.xy[,pos.lon] ]
                            particles.reference.bm.sec[ id %in% moving.particles.xy[,id] , pos.lat := moving.particles.xy[,pos.lat] ]
                            particles.reference.bm.sec[ id %in% moving.particles.xy[,id] , state := moving.particles.xy[,state] ]
                            particles.reference.bm.sec[ id %in% moving.particles.xy[,id] , t.start := moving.particles.xy[,t.start] ]
                            particles.reference.bm.sec[ id %in% moving.particles.xy[,id] , t.finish := moving.particles.xy[,t.finish] ]
                            particles.reference.bm.sec[ id %in% moving.particles.xy[,id] , cell.rafted := moving.particles.xy[,cell.rafted] ]
                            
                            ## ---------------------------------------------------------------
                            ## Save positions to Video matrix (if condition matched)

                            if( movie.year == simulation.year ) {

                                ## --------------------------------------------------------
                                
                                particles.video.location.x.bm.i <- attach.big.matrix(particles.video.location.x.bm.desc)
                                particles.video.location.y.bm.i <- attach.big.matrix(particles.video.location.y.bm.desc)
                                
                                ## --------------------------------------------------------
                                
                                particles.to.sql.id.moving <- particles.reference.bm.sec[state==1,id]
                                particles.to.sql.id.moving <- particles.to.sql.id.moving[particles.to.sql.id.moving %in% particles.to.sql.id]
                                particles.to.sql.id.moving.condition <- length(particles.to.sql.id.moving) > 0
                                
                                if( particles.to.sql.id.moving.condition ) { 
  
                                  particles.video.location.x.bm.i[ which(particles.to.sql.id %in% particles.to.sql.id.moving) , t.step ] <- unlist(particles.reference.bm.sec[ id %in% particles.to.sql.id.moving,6])
                                  particles.video.location.y.bm.i[ which(particles.to.sql.id %in% particles.to.sql.id.moving) , t.step ] <- unlist(particles.reference.bm.sec[ id %in% particles.to.sql.id.moving,7])
                                  
                                }
                            }
                                    
                      }
                      
                      ## ---------------------------------------------------
                      ## Inject particles to final object [ particles.reference.bm.all ] 

                      particles.reference.bm.i <- numeric()
                      
                      for( bm.i in particles.reference.bm.sec[state!=0,id] ) {
                        
                        particles.reference.bm.i <- c( particles.reference.bm.i , mwhich(particles.reference.bm.all,1,list(bm.i),list('eq')) )
                        
                      }

                      particles.reference.bm.all[particles.reference.bm.i , 3 ] <- particles.reference.bm.sec[state!=0,start.year]
                      particles.reference.bm.all[particles.reference.bm.i , 4 ] <- particles.reference.bm.sec[state!=0,start.month]
                      particles.reference.bm.all[particles.reference.bm.i , 5 ] <- particles.reference.bm.sec[state!=0,start.day]
                      
                      particles.reference.bm.all[particles.reference.bm.i , 6 ] <- particles.reference.bm.sec[state!=0,pos.lon]
                      particles.reference.bm.all[particles.reference.bm.i , 7 ] <- particles.reference.bm.sec[state!=0,pos.lat]
                      particles.reference.bm.all[particles.reference.bm.i , 9 ] <- particles.reference.bm.sec[state!=0,state]
                      particles.reference.bm.all[particles.reference.bm.i , 10 ] <- particles.reference.bm.sec[state!=0,t.start]
                      particles.reference.bm.all[particles.reference.bm.i , 11 ] <- particles.reference.bm.sec[state!=0,t.finish]
                      particles.reference.bm.all[particles.reference.bm.i , 12 ] <- particles.reference.bm.sec[state!=0,cell.rafted]
                    
                      # -----------------------------------------------

                      return(NULL)
                }
                
                stopCluster(cl.2) ; rm(cl.2)
                
                ## -------------------------------------
                
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

particles.reference.bm.i <- mwhich(particles.reference.bm,c(9),list(2), list('eq'))

ReferenceTable <- data.frame( particles.reference.bm[particles.reference.bm.i,] , 
                              travel.time= ( 1 + particles.reference.bm[particles.reference.bm.i,11] - particles.reference.bm[particles.reference.bm.i,10] ) / n.hours.per.day )

## -----------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))

dbWriteTable(sql, "ReferenceTable", ReferenceTable , append=FALSE, overwrite=TRUE )
dbWriteTable(sql, "MovieLon", as.data.frame(t(as.matrix(particles.video.location.x.bm))) , append = FALSE)
dbWriteTable(sql, "MovieLat", as.data.frame(t(as.matrix(particles.video.location.y.bm))) , append = FALSE)

# dbWriteTable(sql, "MovieAlt", as.data.frame(t(as.matrix(particles.video.location.z.bm))) , append = FALSE)
  
dbDisconnect(sql)

##  ---------------------------------------------------------------------------------------------------------------------------------
##  ---------------------------------------------------------------------------------------------------------------------------------
## End of Code
