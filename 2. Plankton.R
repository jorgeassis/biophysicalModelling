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

options(warn=0)

plot(landmass,box=FALSE,legend=FALSE,col=c("grey"))
lines(coastline,box=FALSE,legend=FALSE,col=c("red"))

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
  
  # dist = spDists( as.matrix(pt.i),as.matrix(coastline.pts.t),longlat=TRUE)
  # to.get <- which(dist <= source.sink.dist)
  # to.get.xy <- coastline.pts.t[to.get,]
  # source.sink.xy <- rbind(source.sink.xy,data.frame(x=sort(to.get.xy[,1])[length(to.get.xy[,1])/2],y=sort(to.get.xy[,2])[length(to.get.xy[,2])/2]))
  # coastline.pts.t <- coastline.pts.t[-to.get,]

}

coordinates(source.sink.xy) <- c("x","y")
crs(source.sink.xy) <- dt.projection
source.sink.id <- 1:length(source.sink.xy)

## -----------------------------------------------------

## Remove unwanted release sites

if( !is.null(unwanted.release.sites.shp) ) {
  
  unwanted <- shapefile(unwanted.release.sites.shp)
  unwanted <- as(unwanted,"SpatialPolygons")
  
  points.over.polygon <- as.vector(which(sp::over( source.sink.xy , unwanted , fn = NULL) == 1))
  source.sink.xy <- source.sink.xy[-points.over.polygon]
  
}

points(source.sink.xy,box=FALSE,legend=FALSE,col=c("black"),pch=16)

## ------------------

initial.coords <- as.data.frame(source.sink.xy)

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
  
    movie.sites.xy <- sort( as.vector(get.knnx( as.data.frame(source.sink.xy) , movie.sites.xy , k = 1 + movie.sites.buffer , algorithm="kd_tree" )$nn.index) )
}

if( !is.null(movie.sites.xy) ) { movie.sites.id <- unique(movie.sites.xy)
                                 movie.sites.xy <- as.data.frame(source.sink.xy[ unique(movie.sites.xy) , ] ) 
 
}

points(movie.sites.xy,box=FALSE,legend=FALSE,col=c("red"),pch=16)

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

## SQL configuration

if( paste0(project.name,"SimulationResults.sql") %in% list.files(sql.directory) ) {
  
  x <- ""
  while( x != "Y" & x != "n" ) { x <- readline("SQL database already exists. Do you which to overwrite? (Y/n) ") }
  
  if (x == "Y" ) { file.remove( paste0(sql.directory,"/",project.name,"SimulationResults.sql") ) }
}
  
## -----------------------

if( ! paste0(project.name,"SimulationResults.sql") %in% list.files(sql.directory) ) {
    
    global.simulation.parameters <- data.frame(   kill.by.raft = kill.by.raft , 
                                                  n.hours.per.day = n.hours.per.day , 
                                                  n.new.particles.per.day = n.new.particles.per.day , 
                                                  remove.new.particles.last.days = remove.new.particles.last.days , 
                                                  longevity = longevity , 
                                                  particle.max.duration = particle.max.duration , 
                                                  behaviour = behaviour   )       
    
    sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))
    dbWriteTable(sql, "ReleaseSites", as.data.frame(source.sink.xy)  , append=FALSE, overwrite=TRUE )
    dbWriteTable(sql, "Parameters", global.simulation.parameters , append=FALSE, overwrite=TRUE )
    dbDisconnect(sql)
    
}


## ------------------------------------------------------------------------------------------------------------------

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

particles.reference <- data.table( id = 1:(n.particles.per.cell * length(source.sink.xy) ) )
particles.reference[ , start.cell := as.numeric( sapply( source.sink.id ,function(x) { rep(x,n.particles.per.cell) })) ]
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

for(i in 1:parallel.computational.sections){
  
  clipper <- as(extent(min.lon,max.lon, sections.lat[i,1] - parallel.computational.buffer , sections.lat[i,2] + parallel.computational.buffer ), "SpatialPolygons")
  crs(clipper) <- dt.projection 
  
  assign( paste0("landmass.sect.",i) , gIntersection(landmass, clipper, byid=TRUE) )
  
}

## ------------------------------------------------------------------------------------------------------------------

## Start Simulation

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
                cat('\n||',paste0(rep("-",progress.percent),collapse = ""),"(",progress.percent,"% )")
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
                
                ptm <- proc.time()
                cl.2 <- makeCluster(number.cores)
                registerDoParallel(cl.2)
                
                for( section in 1:parallel.computational.sections ) {
                  
                      sections.lat.f.s <- as.numeric(sections.lat[section,1])
                      sections.lat.t.s <- as.numeric(sections.lat[section,2])
                  
                      ## --------------------------------------------------------
                      
                      particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
                      particles.reference.bm.i <- mwhich(particles.reference.bm.all,c(4,4),list(sections.lat.f.s,sections.lat.t.s), list('ge', 'lt') , 'AND')
                      particles.reference.bm.sec <- as.data.table(particles.reference.bm.all[particles.reference.bm.i,])
                      
                      ## --------------------------------------------------------
                      
                      particles.video.location.x.bm.i <- attach.big.matrix(particles.video.location.x.bm.desc)
                      particles.video.location.y.bm.i <- attach.big.matrix(particles.video.location.y.bm.desc)
                      
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
                                                  moving.particles.xy[ id %in% true.rafters.id , c("cell.rafted","state") := list(true.rafters.cell,2) ] 
                                  
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
    
                                                    moving.particles.xy[ id %in% non.rafters.id[non.rafters.t] , c("cell.rafted","state") := list( non.rafters.cell[non.rafters.t] , 2 ) ] 
                                                    
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
                            # Move Particles
                            
                            particles.reference.bm.i <- mwhich(particles.reference.bm.all,rep(1,nrow(moving.particles.xy)),list(moving.particles.xy[,id]),list(rep('eq',nrow(moving.particles.xy))),'OR')
                            
                            particles.reference.bm.all[particles.reference.bm.i , 3 ] <- moving.particles.xy$pos.lon
                            particles.reference.bm.all[particles.reference.bm.i , 4 ] <- moving.particles.xy$pos.lat
                            particles.reference.bm.all[particles.reference.bm.i , 6 ] <- moving.particles.xy$state
                            particles.reference.bm.all[particles.reference.bm.i , 7 ] <- moving.particles.xy$t.start
                            particles.reference.bm.all[particles.reference.bm.i , 8 ] <- moving.particles.xy$t.finish
                            particles.reference.bm.all[particles.reference.bm.i , 9 ] <- moving.particles.xy$cell.rafted
                            
                            ## ---------------------------------------------------------------
                            ## Save positions to Video matrix (if condition matched)

                            if( movie.year == simulation.year ) {
                              
                                particles.reference.bm.i <- numeric()
                                
                                for( bm.i in particles.to.sql.id) {
                                  
                                  particles.reference.bm.i <- c(particles.reference.bm.i,mwhich(particles.reference.bm.all,c(1,6),list(bm.i,1),list('eq','eq'),'AND'))
                                  
                                }
    
                              particles.to.sql.id.moving <- particles.reference.bm.all[particles.reference.bm.i,1]
                              particles.to.sql.id.moving.condition <- length(particles.to.sql.id.moving) > 0
                              
                                    if( particles.to.sql.id.moving.condition ) { 
    
                                      particles.video.location.x.bm.i[ which(particles.to.sql.id %in% particles.to.sql.id.moving) , t.step ] <- particles.reference.bm.all[particles.reference.bm.i,3]
                                      particles.video.location.y.bm.i[ which(particles.to.sql.id %in% particles.to.sql.id.moving) , t.step ] <- particles.reference.bm.all[particles.reference.bm.i,4]
                                      
                                    }
                            }
                                    
                      }

                      # -----------------------------------------------

                }
                
                stopCluster(cl.2) ; rm(cl.2)
                
                ## -------------------------------------
                
}
                
                
                







                
                
                # ------------------------------------------------------------------------------------
                # Save Reference Table in SQL
                
                particles.reference[ , travel.time := ( ( 1 + t.finish-t.start ) / n.hours.per.day ) ]
                particles.reference[ , "Year" := simulation.year ]
                
                sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",results.files,".reference.particles.sql"))
                dbWriteTable(sql, "Particles_reference_table", particles.reference , append=TRUE, overwrite=FALSE )
                
                if( particles.to.sql.years == simulation.year ) {
                  
                dbWriteTable(sql, "lon", data.frame(particles.video.location.x) , append = TRUE)
                dbWriteTable(sql, "lat", data.frame(particles.video.location.y) , append = TRUE)
                dbWriteTable(sql, "alt", data.frame(particles.video.location.z) , append = TRUE)
                
                }
                
                dbDisconnect(sql)
                
                # -------------------------------------------------------------------------------------------------------------------------
                # Resolve connectivity

                # Erase those that did not acomplish
                # 0 unborne
                # 1 living
                # 2 rafted 
                # 3 on hold out of space
                # 4 dead by time
                
                particles.reference <- particles.reference[ state == 2 & cell.rafted != -9 ]
                setkey(particles.reference, cell )

                n.particles.alldays <- length(simulation.parameters.step$day) * n.new.particles.per.day
                progress.full <- length(unique(particles.reference[,cell]))
                cell.to.process <- unique(particles.reference[,cell])
                
                cat('\014')
                cat('\n')
                cat('\n')
                cat('\n Saving results to SQL')
                cat('\n')
                cat('\n')
                
                cl.2 <- makeCluster(number.cores)
                registerDoParallel(cl.2)
                
                all.connectivity.pairs.to.sql <- foreach(cell.id.ref.f=cell.to.process, .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN")) %dopar% { # 
                  
                          connectivity.pairs.to.sql <- data.frame()
                          connectivity.temp <- particles.reference[ cell == cell.id.ref.f , ]
        
                          for( cell.id.ref.t in unique(connectivity.temp[ , cell.rafted ]) ) {
                            
                                  connectivity.pairs.to.sql <- rbind(connectivity.pairs.to.sql,
                                                                     
                                                                     data.frame(  Pair.from = cell.id.ref.f,
                                                                                  Pair.to = cell.id.ref.t,
                                                                                  Number.events = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]),
                                                                                  Time.mean = mean(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                                  Time.min = min(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                                  Time.max = max(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                                  Time.sd = sd(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                                  Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / n.particles.alldays,
                                                                                  Year = simulation.year ) 
                                                                     
                                                                     )
                          }
                          
                          connectivity.pairs.to.sql[is.na(connectivity.pairs.to.sql)] <- 0
                          return( connectivity.pairs.to.sql )
                }
                
                stopCluster(cl.2) ; rm(cl.2)
                
                # Save pairs to SQL
                
                sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",results.files,".reference.particles.sql"))
                dbWriteTable(sql, "Connectivity", all.connectivity.pairs.to.sql , overwrite=FALSE, append=TRUE)
                dbDisconnect(sql)
                
                # ------------------------------------------------------------------
                # Remove heavy objects from memory
                
                objects.to.rm <- c("particles.video.location.x" , "particles.video.location.y" , "all.connectivity.pairs.to.sql" , "velocity.field" , "positions.fate" ,
                                   "particles.reference" , "particles.video.location.z" , "raw.data.u" , "raw.data.v" , "speed.u.rk" , "speed.v.rk" , "raw.data.u.t" , "raw.data.v.t" )
                
                for( i in 1:length(objects.to.rm) ) { if( exists(objects.to.rm[i]) ) { rm(list=objects.to.rm[i]) }  }
         
                
                time.f <- Sys.time()
                

}

                stopCluster(cl.2) ; rm(cl.2); gc(reset=TRUE)
                
                ## -------------------------------------
                
                clean.dump.files(clean.dump.files=TRUE,files="raw.data.",dump.folder=paste0(project.folder,"/InternalProc"))
                clean.dump.files(clean.dump.files=TRUE,files="particles.reference.",dump.folder=paste0(project.folder,"/InternalProc"))
                file.remove( paste0(project.folder,"/InternalProc"))
                
}
                
##  ---------------------------------------------------------------------------------------------------------------------------------
##  ---------------------------------------------------------------------------------------------------------------------------------
## End of Code
