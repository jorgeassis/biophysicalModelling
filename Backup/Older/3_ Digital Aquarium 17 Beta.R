## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
##
## Digital Aquarium 2016 Versio 3.0
## Simulation of Dispersal using ocean fields
##

## Change Log:
## Parallel save results to SQL (Try open Db before parallelization)

## ------------------------------------------------------------------------------------------------------------------
## Main Configuration

main.directory <- "/Volumes/Albacora/Dropbox/Gist/One Aquarium V2.0" # Albacora Jellyfish
raw.data.directory <- "Data/Currents"
results.directory <- "Results"
results.files <- "NAtlantic"

## ------------------------

ocean.shape.file <- "Data/ocean.tif"
coast.line.file <- "Data/coast_line.tif"
unwated.regions.file <- "Data/Shapefiles/regionsexcluded.shp" # NULL

new.extent <- NULL # extent( -15 , 10 , 30 , 47.5 )

## ------------------------

number.cores <- 10

kill.by.raft <- TRUE                              # Will eliminate particles that got to another cell - first raft event # May need a new particle every day
n.hours.per.day <- 12                             # Uncoded for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- TRUE            # If last days (particle.max.duration) are not to deliver new particles  

longevity <- TRUE
particle.max.duration <- 30                       # Days allowed to travel
behaviour <- NULL                                 # Uncoded

dump.unresolved.to.sql <- TRUE                   # Default; Save all paired events to SQL

# Ilustration (movie)

parcticles.to.sql.xy <- "Data/Shapefiles/points_video.shp" # NULL # matrix( c(  -2.20 , 36.73 , -5.42 , 36.17 , -15.92 , 23.76 ) , ncol=2 , byrow=TRUE)
parcticles.to.sql.years <- 2003
parcticles.to.sql.buffer <- 0 # Nearby cells to include, 0 for xy only

## ----------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------

setwd(main.directory)
source("Dependences.R")

# Job listing

raw.data.currents.files <- list.files(raw.data.directory,full.names = TRUE,pattern="nc")
available.raw.data <- data.frame()

for( file in 1:length(raw.data.currents.files)) {
  nc <- nc_open( raw.data.currents.files[file] , verbose=FALSE )
  nc.date <- as.Date(ncvar_get( nc, "Time"), origin = "1970-01-01")
  extent.raw.data <- c(  min(ncvar_get( nc, "Longitude")) , max(ncvar_get( nc, "Longitude")) , min(ncvar_get( nc, "Latitude")) , max(ncvar_get( nc, "Latitude")) )
  nc_close( nc )
  available.raw.data <- rbind( available.raw.data, data.frame( simulation=file, year=substr(nc.date, 1, 4) , month=substr(nc.date, 6, 7) , day=substr(nc.date, 9, 10)) )
}

## -----------------------------------------------

coast.line <- raster(coast.line.file)
ocean <- raster(ocean.shape.file)

if( is.null(new.extent) ) { new.extent <- extent( ocean ) }

coast.line <- crop(coast.line, new.extent ); plot(coast.line,col="black")
ocean <- crop(ocean, new.extent ) ; plot(ocean,col=c("black"))
simulation.resolution <- res(ocean)[1]
  
## ------------------------------------------------------------------
## Prepare data and Analysis

parcticles.to.sql.xy <- shapefile(parcticles.to.sql.xy)
parcticles.to.sql.xy <-  crop(parcticles.to.sql.xy, new.extent )
parcticles.to.sql.xy <- as.data.frame(parcticles.to.sql.xy)[,2:3]

if( ! "Results" %in% list.files() ) { dir.create(file.path("Results")) }
if( ! "SQLite" %in% list.files(results.directory) ) { dir.create(file.path(results.directory,"SQLite")) }

## ------------------------------------------------------------------
## Prepare Sites and Surfaces

if( !is.null(unwated.regions.file) ) {
  
  unwanted <- shapefile(unwated.regions.file)
  unwanted <- rasterize(unwanted, ocean)

  unwanted.1 <- unwanted
  unwanted.1[is.na(unwanted.1)] <- -9 ; unwanted.1[unwanted.1 != -9 ] <- NA ; unwanted.1[unwanted.1 == -9 ] <- 1
  study.region.used <- mask(coast.line,unwanted.1)
  study.region.used[which(!is.na(getValues(study.region.used)))] <- 1:length(which(!is.na(getValues(study.region.used))))
  
  unwanted.2 <- unwanted
  unwanted.2[!is.na(unwanted.2)] <- 1
  study.region.unused <- mask(coast.line,unwanted.2)
  study.region.unused[which(!is.na(getValues(study.region.unused)))] <- -9
  
  coast.line <- merge(study.region.used,study.region.unused)
  rm( unwanted ) ; rm( unwanted.1 ) ; rm( unwanted.2 ) ; rm( study.region.used ) ; rm( study.region.unused )

}

coast.line.val <- getValues(coast.line)
coast.line.val <- which( !is.na(coast.line.val) & coast.line.val != -9 )
coast.line[coast.line.val] <- 1:(length(coast.line.val)) # writeRaster(coast.line,filename="Data/temp_coatline",format="GTiff",overwrite=T)

## ------------------------------------------------------------------
## Ocean region
## -1 for ocean, 0 for land, -9 for unwanted cells and numbered cells

ocean.region.dt <- ocean
ocean.region.dt[ocean.region.dt == 1] <- (-1)
ocean.region.dt <- calc(stack(coast.line,ocean.region.dt) , function(x) { ifelse( !is.na(x[1]) , x[1] , x[2] )  } ) 
ocean.region.dt[is.na(ocean.region.dt)] <- 0 # writeRaster(ocean.region.dt,filename="Data/temp_ocean",format="GTiff",overwrite=T)
plot(ocean.region.dt)

## ------------------------------------------------------------------

cells <- as.data.frame(ocean.region.dt,xy=TRUE)
cells <- cells[! cells$layer %in% c(-1,0,-9) ,1:2]

if( min(cells[,1]) < extent.raw.data[1] | max(cells[,1]) > extent.raw.data[2] | min(cells[,2]) < extent.raw.data[3] | max(cells[,2]) > extent.raw.data[4] ) { print("ERROR: Extent wider than Raw Data") } else { print("OK: Extent within Raw Data")}

cells.id <- as.data.frame(ocean.region.dt,xy=TRUE)
cells.id <- cells.id[! cells.id$layer %in% c(-1,0,-9) ,3]
number.cells <- nrow(cells)

## ------------------------------------------------------------------------------------------------------------------
## Video configuration

if ( !is.null(parcticles.to.sql.xy) ) { particles.to.sql.cells <- sort( as.vector(get.knnx( cells , parcticles.to.sql.xy , k = 1 + parcticles.to.sql.buffer , algorithm="kd_tree" )$nn.index) ) }

plot(cells,col="grey")
points(parcticles.to.sql.xy,col="black")

## Save baseline information

if( paste0(results.files,".reference.particles.sql") %in% list.files(paste0(results.directory,"/SQLite/")) ) {
  
  x <- ""
  
  while( x == "" ) { x <- readline("SQL database already exists. Do you which to overwrite? (y/n) ") 
  
  if ( x == "y" | x == "n"  ) { break ; } else { x <- "" }
  
  }
  
  if (x == "y" ) {
    file.remove( paste0(results.directory,"/SQLite/",results.files,".reference.particles.sql") )
    
    if( paste0(results.files,".particles.video.sql") %in% list.files(paste0(results.directory,"/SQLite/")) ) {
      file.remove( paste0(results.directory,"/SQLite/",results.files,".particles.video.sql") )
    }
    
  }
}
  
## -------------------------------------------------------------------
## Remove Years from Database

sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/SQLite/",results.files,".reference.particles.sql"))
dbListFields(sql, "Connectivity")
dbListFields(sql, "Particles_reference_table")
all.connectivity <- dbReadTable(sql, "Connectivity") ; unique(all.connectivity$Year)
all.particles.reference <- dbReadTable(sql, "Particles_reference_table") ; unique(all.particles.reference$Year)
dbExecute(sql, 'DELETE FROM Connectivity WHERE "Year" == 2006')
dbExecute(sql, 'DELETE FROM Particles_reference_table WHERE "Year" == 2006')
dbDisconnect(sql)

## -------------------------------------------------------------------
## Generate Database

if( ! paste0(results.files,".reference.particles.sql") %in% list.files(paste0(results.directory,"/SQLite/")) ) {
    sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/SQLite/",results.files,".reference.particles.sql"))
    dbWriteTable(sql, "Cells_coordinates", data.table(cell=cells.id,x=cells[,1],y=cells[,2]) , append=TRUE, overwrite=FALSE )
    dbDisconnect(sql)
}

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------

## Start Simulation
## unique(available.raw.data$simulation)

for ( simulation.step in 7:10 ) {
                
                cat('\014') ; cat('\n')
                cat('\n Preparing Environmental data')
  
                ## --------------------------------------------------------
                
                simulation.parameters.step <- available.raw.data[which( available.raw.data$simulation == simulation.step ),]
                simulation.year <- as.numeric(as.character(unique(simulation.parameters.step$year)))

                ## --------------------------------------------------------
                ## Prepare environmental data

                simulation.raw.data <- nc_open( raw.data.currents.files[1], readunlim=FALSE )
                dim.i <- ncvar_get( simulation.raw.data, "X" )
                dim.j <- ncvar_get( simulation.raw.data, "Y" )
                Longitude <- ncvar_get( simulation.raw.data, "Longitude" )
                Latitude <- ncvar_get( simulation.raw.data, "Latitude" )
                
                nc_close(simulation.raw.data)
                
                norm.field <- expand.grid( i=dim.i , j=dim.j )
                raw.data.coords <- cbind( apply( norm.field , 1 , function (x) Longitude[x[1]] ) , apply( norm.field , 1 , function (x) Latitude[x[2]] ) )
                colnames(raw.data.coords) <- c("x","y")
                
                raw.data.u <- matrix(NA , nrow=nrow(raw.data.coords) , ncol=nrow(simulation.parameters.step) )
                raw.data.v <- matrix(NA , nrow=nrow(raw.data.coords) , ncol=nrow(simulation.parameters.step) )

                for( rd.i in 1:nrow(simulation.parameters.step) ) {
                  
                          rd.file <- simulation.parameters.step[rd.i,1]
                          nc.opened.file <- nc_open( raw.data.currents.files[rd.file], readunlim=FALSE )
                          
                          # u component
                          
                          velocity.field <- ncvar_get(  nc.opened.file , "UComponent", start=c(1,1,rd.i) , count=c(length(dim.i),length(dim.j),1) )
                          raw.data.u[,rd.i] <- apply( norm.field , 1 , function (y) velocity.field[y[1],y[2]] )
                          
                          # v component
                          
                          velocity.field <- ncvar_get(  nc.opened.file , "VComponent", start=c(1,1,rd.i) , count=c(length(dim.i),length(dim.j),1) )
                          raw.data.v[,rd.i] <- apply( norm.field , 1 , function (y) velocity.field[y[1],y[2]] )
                          
                          nc_close(nc.opened.file)

                }

                raw.data.u <- data.table(raw.data.coords,raw.data.u)
                setnames(raw.data.u,names(raw.data.u),c("Lon","Lat",sapply(1:(ncol(raw.data.u)-2),function(x) { paste("day.",x,sep="") })))
                
                raw.data.v <- data.table(raw.data.coords,raw.data.v)
                setnames(raw.data.v,names(raw.data.v),c("Lon","Lat",sapply(1:(ncol(raw.data.v)-2),function(x) { paste("day.",x,sep="") })))
                
                coords.cells <- data.table(id=cells.id,cells)
                setnames(coords.cells,names(coords.cells),c("Cell","Lon","Lat"))
                
                # Crop environmental data by extent (study region)
                
                raw.data.u <- raw.data.u[ Lon >= new.extent[1] & Lat >= new.extent[3] & Lon <= new.extent[2] & Lat <= new.extent[4]  , ]
                raw.data.v <- raw.data.v[ Lon >= new.extent[1] & Lat >= new.extent[3] & Lon <= new.extent[2] & Lat <= new.extent[4]  , ]

                ## --------------------------------------------------------
                ## Define conditions

                norm.time <- expand.grid( hour=1:n.hours.per.day , day=1:(nrow(simulation.parameters.step)) )
                norm.time <- data.table(norm.time)
                
                new.day <- norm.time[,hour == 1]
                end.of.day <- norm.time[,hour == n.hours.per.day]
                
                # Remove the last day from new.day because there is no day #367 in the data for currents
                
                last.new.day <- which(new.day)
                new.day[last.new.day[length(last.new.day)]] <- FALSE
                
                # ------
                
                release.particles.t <- seq(from=1,to=n.hours.per.day,by=(n.hours.per.day/n.new.particles.per.day))
                release.particles.condition <- norm.time[,hour %in% release.particles.t] 
                
                if( remove.new.particles.last.days ) { release.particles.condition[(length(release.particles.condition)-(n.hours.per.day*particle.max.duration)):length(release.particles.condition)] <- FALSE }

                n.particles.per.cell <- (nrow(simulation.parameters.step)) * length(release.particles.t)
                n.simulation.steps <- nrow(norm.time)
                  
                runge.kutta.sequence <- c(sapply( 1:(n.hours.per.day), function (x) rep( x ,  n.hours.per.day / ( n.hours.per.day ) ) ))
                runge.kutta.sorter <- 1:n.hours.per.day

                ## ---------------------------
                
                initial.coords <- data.table(cell=cells.id,x=cells[,1],y=cells[,2])
                
                all.but.first.day <- rep(TRUE,n.simulation.steps)
                all.but.first.day[1] <- FALSE
                
                ## --------------------------------------------------------
                ## Define particles
                
                # data.table particles.reference[id,cell,state,t.start,t.finish,cell.rafted]
                # 0 unborne
                # 1 living
                # 2 rafted 
                # 3 on hold out of space
                # 4 dead by time
                
                particles.reference <- data.table(id = 1:(n.particles.per.cell * number.cells))
                particles.reference[ , cell := as.numeric(sapply( cells.id ,function(x) { rep(x,n.particles.per.cell) })) ]
                particles.reference[ , state := 0 ]
                particles.reference[ , t.start := 0 ]
                particles.reference[ , t.finish := 0 ]
                particles.reference[ , cell.rafted := 0 ]
                setkey(particles.reference, id )
                
                # data.table particles.location[id,time]

                particles.location.x <- matrix( c(0,rep(NA,n.simulation.steps)) , nrow=1 , ncol=n.simulation.steps+1 )
                colnames(particles.location.x) <- c("id", sapply(1:n.simulation.steps , function(x) { paste0("t.",x ) }   ) )
                particles.location.x <- as.data.table( particles.location.x[-1,] )
                setkey(particles.location.x, id )
                
                particles.location.y <- matrix( c(0,rep(NA,n.simulation.steps)) , nrow=1 , ncol=n.simulation.steps+1 )
                colnames(particles.location.y) <- c("id", sapply(1:n.simulation.steps , function(x) { paste0("t.",x ) }   ) )
                particles.location.y <- as.data.table( particles.location.y[-1,] )
                setkey(particles.location.y, id )

                particles.location.z <- matrix( c(0,rep(NA,n.simulation.steps)) , nrow=1 , ncol=n.simulation.steps+1 )
                colnames(particles.location.z) <- c("id", sapply(1:n.simulation.steps , function(x) { paste0("t.",x ) }   ) )
                particles.location.z <- as.data.table( particles.location.z[-1,] )
                setkey(particles.location.z, id )
                
                particles.location.temp <- matrix( c(NA,rep(NA,n.simulation.steps)) , nrow=number.cells , ncol=n.simulation.steps+1 )
                colnames(particles.location.temp) <- c("id", sapply(1:n.simulation.steps , function(x) { paste0("t.",x ) }   ) )
                
                ## --------------------------------------------------------
                ## Move particles
   
                ptm <- proc.time()
                cl.2 <- makeCluster(number.cores)
                registerDoParallel(cl.2)
                
                # 1:n.simulation.steps
                
                for( t.step in 1:n.simulation.steps ) {
                  
                        ## -----------------------

                        hour <- norm.time[t.step,hour]
                        day <- norm.time[t.step,day]
                        
                        ## -----------------------
                        
                        # Environmental data between days
                        
                        if( new.day[t.step] ) {   
                                                  raw.data.u.t <- subset( raw.data.u[ , .(Lon , Lat , get(paste0("day.",day)) , get(paste0("day.",day+1)) ) ] , !is.na(V3) )
                                                  raw.data.v.t <- subset( raw.data.v[ , .(Lon , Lat , get(paste0("day.",day)) , get(paste0("day.",day+1)) ) ] , !is.na(V3) )

                                                  comb.to <- round ( seq(  nrow(raw.data.u.t) / number.cores , nrow(raw.data.u.t) , length.out=number.cores) )
                                                  comb.to <- c(comb.to[-length(comb.to)] , nrow(raw.data.u.t))
                                                  comb.from <- c(1,comb.to + 1)       
                                                  combinations <- data.frame(from = comb.from[-length(comb.from)] , to = comb.to )

                                                  speed.u.rk <- foreach(s=1:nrow(combinations), .combine='rbind', .verbose=FALSE, .packages=c("gstat","raster","data.table")) %dopar% {
                                                    
                                                                speed.u.rk <- t( sapply( combinations$from[s]:combinations$to[s] , function (x) seq( from = raw.data.u.t[x, V3 ] , to = raw.data.u.t[x, V4 ] , length.out = n.hours.per.day + 1 ) ))
                                                                speed.u.rk <- data.table(speed.u.rk)
                                                                setnames(speed.u.rk,sapply( 1:(n.hours.per.day + 1) , function(x) { paste0( "t",x) } ))
                                                                speed.u.rk[, c("Lon","Lat") := list(raw.data.u.t[ combinations$from[s]:combinations$to[s] , Lon ] , raw.data.u.t[combinations$from[s]:combinations$to[s],Lat])]
                                                                setcolorder(speed.u.rk,c("Lon","Lat","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13"))
                                                                return(speed.u.rk)
                                                                
                                                  }
                                                  
                                                  speed.v.rk <- foreach(s=1:nrow(combinations), .combine='rbind', .verbose=FALSE, .packages=c("gstat","raster","data.table")) %dopar% {
                                                    
                                                                speed.v.rk <- t( sapply( combinations$from[s]:combinations$to[s] , function (x) seq( from = raw.data.v.t[x, V3 ] , to = raw.data.v.t[x, V4 ] , length.out = n.hours.per.day + 1 )))
                                                                speed.v.rk <- data.table(speed.v.rk)
                                                                setnames(speed.v.rk,sapply( 1:(n.hours.per.day + 1) , function(x) { paste0( "t",x) } ))
                                                                speed.v.rk[, c("Lon","Lat") := list(raw.data.v.t[ combinations$from[s]:combinations$to[s] , Lon ] , raw.data.v.t[combinations$from[s]:combinations$to[s],Lat])]
                                                                setcolorder(speed.v.rk,c("Lon","Lat","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13"))
                                                                return(speed.v.rk)
                                                                
                                                  }
                                                  
                        }
                  
                        # -----------------------------------------------                        
                        # Get positions onto t.step
                        
                        if( all.but.first.day[t.step] ) {   particles.location.x[ , paste0("t.",t.step) := get( paste0("t.",t.step-1) ) ]
                                                            particles.location.y[ , paste0("t.",t.step) := get( paste0("t.",t.step-1) ) ]
                                                        }
                        
                        # -----------------------------------------------
                        # Release new particles, if that is the case
                       
                        if( release.particles.condition[t.step] ) {   
                          
                                                new.particles.id <- particles.reference[ state == 0 , .SD[1] , by=cell][,id]
                                                new.particles.cells <- particles.reference[ state == 0 , .SD[1] , by=cell][,cell]
                                                
                                                particles.reference[ id %in% new.particles.id , state := 1 ]
                                                particles.reference[ id %in% new.particles.id , t.start := t.step ]

                                                particles.location.temp.x <- particles.location.temp
                                                particles.location.temp.x[,1] <- new.particles.id
                                                particles.location.temp.x <- as.data.table( particles.location.temp.x )
                                                particles.location.temp.x[ , paste0("t.",t.step) := initial.coords$x[new.particles.cells] ]
                                                
                                                particles.location.temp.y <- particles.location.temp
                                                particles.location.temp.y[,1] <- new.particles.id
                                                particles.location.temp.y <- as.data.table( particles.location.temp.y )
                                                particles.location.temp.y[ , paste0("t.",t.step) := initial.coords$y[new.particles.cells] ]

                                                particles.location.x <- rbindlist( list( particles.location.x , particles.location.temp.x )  , use.names=TRUE )
                                                setkey(particles.location.x, id )
                                                particles.location.y <- rbindlist( list( particles.location.y , particles.location.temp.y )  , use.names=TRUE )
                                                setkey(particles.location.y, id )   
                        }

                        # -----------------------------------------------
                        # Which to move and speed
                        
                        moving.particles.condition <- nrow(particles.location.x) != 0

                        if( moving.particles.condition ) {

                                  # Sort Particles By Latitude
                                   
                                  moving.particles.xy <- merge( particles.location.x[ , .(id,get( paste0("t.",t.step))) ] ,
                                                                particles.location.y[ , .(id,get( paste0("t.",t.step))) ] , 
                                                                by="id" , all=TRUE )
                                  
                                  setnames(moving.particles.xy,c("id","x","y"))
                                  moving.particles.xy <- moving.particles.xy[order(y, decreasing = TRUE), ]
                                  
                                  moving.particles.ids <- moving.particles.xy[,id]
                                  moving.particles.cells <- particles.reference[ id %in% moving.particles.ids , cell ]
                                  moving.particles.tstart <- particles.reference[ id %in% moving.particles.ids , t.start ]
                                  
                                  ## -----------------------
                                  
                                  comb.to <- round ( seq(  nrow(moving.particles.xy) / number.cores , nrow(moving.particles.xy) , length.out=number.cores) )
                                  comb.to <- c(comb.to[-length(comb.to)] , nrow(moving.particles.xy))
                                  comb.from <- c(1,comb.to + 1)       
                                  combinations <- data.frame(from = comb.from[-length(comb.from)] , to = comb.to )

                                  positions.fate <- foreach(s=1:nrow(combinations), .verbose=FALSE, .combine = rbind ,  .packages=c("gstat","raster","data.table","FNN")) %dopar% { # 
                                    
                                          ids.para <- moving.particles.ids[combinations$from[s]:combinations$to[s]]
                                          cells.para <- moving.particles.cells[combinations$from[s]:combinations$to[s]]
                                          tstart.para <- moving.particles.tstart[combinations$from[s]:combinations$to[s]]
                                          
                                          points.to.interp <- moving.particles.xy[combinations$from[s]:combinations$to[s], .(x,y) ]
                                          coordinates(points.to.interp) = ~x+y
                                          
                                          extent.para <- c( extent(points.to.interp)[1] - ( simulation.resolution * 10 ) ,
                                                            extent(points.to.interp)[2] + ( simulation.resolution * 10 ) ,
                                                            extent(points.to.interp)[3] - ( simulation.resolution * 10 ) ,
                                                            extent(points.to.interp)[4] + ( simulation.resolution * 10 ) 
                                                            )

                                          ocean.region.dt.t <- crop(ocean.region.dt,extent.para)
                                          
                                          source.points.to.interp.u <- speed.u.rk[ Lon >= extent.para[1] & Lon <= extent.para[2] & Lat >= extent.para[3] & Lat <= extent.para[4] , .(Lon,Lat,get(paste0("t",hour))) ]
                                          setnames(source.points.to.interp.u,c("x","y","var"))
                                          coordinates(source.points.to.interp.u) = ~x+y
                                   
                                          source.points.to.interp.v <- speed.v.rk[ Lon >= extent.para[1] & Lon <= extent.para[2] & Lat >= extent.para[3] & Lat <= extent.para[4] , .(Lon,Lat,get(paste0("t",hour))) ]
                                          setnames(source.points.to.interp.v,c("x","y","var"))
                                          coordinates(source.points.to.interp.v) = ~x+y

                                          sink("/dev/null"); speed.u <- idw(formula = var ~ 1, source.points.to.interp.u, points.to.interp, nmax=3)$var1.pred ; sink()
                                          sink("/dev/null"); speed.v <- idw(formula = var ~ 1, source.points.to.interp.v, points.to.interp, nmax=3)$var1.pred ; sink()
                                          
                                          mov.eastward <- speed.u * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
                                          mov.northward <- speed.v * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
                                          
                                          # Assign temporary positions

                                          points.to.interp <- as.data.frame(points.to.interp)
                                          
                                          dLon <- mov.eastward / ( 6378137 * cos( pi * (  points.to.interp[,2]  / 180) ) )
                                          dLat <- mov.northward / 6378137
                                          dLon <- points.to.interp[,1] + dLon * 180/pi 
                                          dLat <- points.to.interp[,2] + dLat * 180/pi

                                          moving.particles.xy.t <- data.table( id=ids.para , old.x = points.to.interp[,1] , old.y = points.to.interp[,2] , x=dLon , y=dLat , init.cell=cells.para , t.started = tstart.para, rafted.cell=0 , fate=0)
                                          setkey(moving.particles.xy.t,id)
                                          
                                          # Test Positions
                                          # Test 1. Out of space (study region), if TRUE, place particles on hold

                                          out.of.space <- dLon > new.extent[2] | dLon < new.extent[1] | dLat > new.extent[4] | dLat < new.extent[3]
                                          out.of.space.ids <- ids.para[out.of.space]
                                          moving.particles.xy.t[ id %in% out.of.space.ids , fate := 3 ] 
                                          
                                          # Test 2. kill by first raft . Will eliminate particles that got to another cell - first raft event

                                          if( kill.by.raft ) {
                                            
                                                    particles.on.land <- extract(ocean.region.dt.t,moving.particles.xy.t[ ,.(x,y) ])
                                                    particles.on.land.condition <- sum(particles.on.land == 0,na.rm=T) > 0

                                                    if( particles.on.land.condition ) {    
                                            
                                                              who.at.land.id <- moving.particles.xy.t[ particles.on.land == 0 , id ]
                                                              who.at.land.tstart <- moving.particles.xy.t[ particles.on.land == 0 , t.started ]
                                                              cells.started <- moving.particles.xy.t[ particles.on.land == 0 ,init.cell] 
                                                              
                                                              cells.rafted <- extract(ocean.region.dt.t,moving.particles.xy.t[ particles.on.land == 0 ,.(old.x,old.y) ])
                                                              displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]})
                                                              
                                                              # For Rafters
                                                              
                                                              true.rafters.id <- who.at.land.id[ displacement != 0 ]
                                                              true.rafters.cell <- cells.rafted[ displacement != 0 ]
                                                              moving.particles.xy.t[ id %in% true.rafters.id , c("rafted.cell","fate") := list(true.rafters.cell,2) ] 
                                                              
                                                              # For non-Rafters New Particles
                                                              
                                                              non.rafters.id <- who.at.land.id[ displacement == 0 ]
                                                              non.rafters.cell <- cells.started[ displacement == 0 ]
                                                              non.rafters.t <- who.at.land.tstart[ displacement == 0 ] == t.step

                                                              if( TRUE %in% non.rafters.t ) {    
                                                                
                                                                moving.particles.xy.t[id %in% non.rafters.id[non.rafters.t], c("x","y") := list( initial.coords[cell %in% non.rafters.cell, x ] , initial.coords[cell %in% non.rafters.cell, y ] ) ]
  
                                                              }
                                                              
                                                              # For non-Rafters Old Particles
                                                              
                                                              non.rafters.id <- who.at.land.id[displacement == 0]
                                                              non.rafters.t <- who.at.land.tstart[displacement == 0] < t.step
                                                              
                                                              if( TRUE %in% non.rafters.t ) {    
                                                                
                                                                moving.particles.xy.t[id %in% non.rafters.id[non.rafters.t], c("x","y") := list( old.x , old.y ) ]
                                                                
                                                              }
                                                    }
                                          }
                                          
                                          return( moving.particles.xy.t[,.(id,x,y,rafted.cell,fate)] )
                                          
                                  }

                                  ## ---------------------------------------------------
                                  # Move Particles
                                  
                                  setkey(positions.fate, id)
                                  setkey(particles.location.x, id)
                                  setkey(particles.reference, id)
                                  
                                  particles.location.x[ id %in% positions.fate[,id] , paste0("t.",t.step) := positions.fate[,x] ]
                                  particles.location.y[ id %in% positions.fate[,id] , paste0("t.",t.step) := positions.fate[,y] ]

                                  # -----------------------------------------------
                                  # Kill by Longevity
                                  
                                  if ( longevity ) {   
                                    
                                                max.duration.id <- particles.reference[ state == 1 & ( t.step - t.start ) > ( n.hours.per.day * particle.max.duration ) , id ]
                                                particles.reference[ id %in% max.duration.id , c("state","t.finish") := list(4, t.step) ]
                                  }

                                  ## ---------------------------------------------------------------
                                  # Inject to SQL and clean objects

                                  rafters <- positions.fate[ fate == 2 ,id]
                                  rafters.position <- positions.fate[fate==2,rafted.cell]

                                  particles.reference[ id %in% rafters , c("state","cell.rafted","t.finish") := list(2,rafters.position,t.step) ]

                                  # ---------
                                  
                                  holders <- positions.fate[ fate == 3 , id ]
                                  particles.reference[ id %in% holders , c("state","t.finish") := list(3,t.step) ]
                                  
                                  # ---------
                                  
                                  inject.and.remove <- unique(c( rafters , 
                                                                 holders , 
                                                                 max.duration.id))
                                  
                                  # Save to SQL if that is the condition

                                  inject.and.remove.cells <- particles.reference[id %in% inject.and.remove,cell]
                                  
                                  if( parcticles.to.sql.years == simulation.year & TRUE %in% ( particles.to.sql.cells %in% inject.and.remove.cells ) ) {

                                            # Get Maximum of number of days * number hours per day
                                            # Save sql
                                            # In video get starting day to put particle on the loop
                                    
                                            parcticles.to.sql.video.id <- particles.reference[ id %in% inject.and.remove & cell %in% particles.to.sql.cells , id ]
                                    
                                            parcticles.to.sql.video <- data.frame(particles.location.x[ id %in% parcticles.to.sql.video.id , ])
                                            parcticles.to.sql.video.x <- parcticles.to.sql.video[,-1]
                                            
                                            parcticles.to.sql.video <- data.frame(particles.location.y[ id %in% parcticles.to.sql.video.id , ])
                                            parcticles.to.sql.video.y <- parcticles.to.sql.video[,-1]

                                            parcticles.to.sql.video.x.sql <- matrix(NA,nrow=length(parcticles.to.sql.video.id),ncol=(particle.max.duration*n.hours.per.day)+n.hours.per.day)
                                            parcticles.to.sql.video.y.sql <- matrix(NA,nrow=length(parcticles.to.sql.video.id),ncol=(particle.max.duration*n.hours.per.day)+n.hours.per.day)
        
                                            for( particles.to.sql.cells.i in 1:length(parcticles.to.sql.video.id) ) {
                                              
                                                  x.data <- parcticles.to.sql.video.x[particles.to.sql.cells.i,][!is.na(parcticles.to.sql.video.x[particles.to.sql.cells.i,])]
                                                  y.data <- parcticles.to.sql.video.y[particles.to.sql.cells.i,][!is.na(parcticles.to.sql.video.y[particles.to.sql.cells.i,])]
                                                  
                                                  parcticles.to.sql.video.x.sql[particles.to.sql.cells.i,1:length(x.data)] <- x.data
                                                  parcticles.to.sql.video.y.sql[particles.to.sql.cells.i,1:length(y.data)] <- y.data
                                                  
                                            }

                                            sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/SQLite/",results.files,".particles.video.sql"))
                                            dbWriteTable(sql, "lon", data.frame(parcticles.to.sql.video.x.sql) , append = TRUE)
                                            dbWriteTable(sql, "lat", data.frame(parcticles.to.sql.video.y.sql) , append = TRUE)
                                            dbWriteTable(sql, "reference", data.frame(parcticles.to.sql.video.id) , append = TRUE)
                                            
                                            temp.ids <- dbReadTable(sql, "reference")
                                            temp.ids <- particles.reference[ id %in% c(unlist(temp.ids)) , cell ]
                                            
                                            if( FALSE %in% (temp.ids %in% particles.to.sql.cells) ) { dbDisconnect(sql) ; stop("\nMovie Cells not included in first set") }
                                            
                                            dbDisconnect(sql)
                                        
                                  }

                                  # Remove from main objects
                                  
                                  particles.location.x <- particles.location.x[ ! id %in% inject.and.remove , ]
                                  setkey(particles.location.x, id)
                                  particles.location.y <- particles.location.y[ ! id %in% inject.and.remove , ]
                                  setkey(particles.location.y, id)
                                  
                         }

                        # -----------------------------------------------
                        # Test Errors, print Log

                        if( end.of.day[t.step] ) {

                              t.step.time.log <- proc.time() - ptm

                              cat('\014')
                              cat('\n')
                              cat('\n')
                              cat(' Simulation of Dispersal (running)')
                              cat(paste0(" | Year " , simulation.year , ", day: " , day , " | Time Taken: " , round(t.step.time.log[3]) , " (memory: ", round(sum(list.memory()$Size)),"kb)" ))
                              cat('\n')
                              cat('\n')
                              progress.partial <- floor( round( day / nrow(simulation.parameters.step) , digits=2) * 100 )
                              cat(paste0(" | ",paste(rep("-", progress.partial ),collapse="")," | ",progress.partial,"%"))

                              ptm <- proc.time()
                              
                              if( sum(is.na( particles.location.y[  , get( paste0("t.",t.step) ) ])) > 0 ) { stop("\nNAs is positions") }
                              if( ! setequal( (particles.reference[state==1,id]) , (particles.location.y[,id]) ) ) { stop("\nParameters Bigger than  positions") }
                              
                        }
                }
                
                stopCluster(cl.2) ; rm(cl.2)
                
                # End of t.step cycles
                
                particles.reference[ , "Year" := simulation.year ]
                
                # ------------------------------------------------------------------------------------
                # Save Reference Table to SQL
                
                if( dump.unresolved.to.sql ) {
                  sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/SQLite/",results.files,".reference.particles.sql"))
                  dbWriteTable(sql, "Particles_reference_table", particles.reference , append=TRUE, overwrite=FALSE )
                  dbDisconnect(sql)
                }

                # -------------------------------------------------------------------------------------------------------------------------
                # Resolve connectivity

                # Erase those that did not acomplish
                # 0 unborne
                # 1 living
                # 2 rafted 
                # 3 on hold out of space
                # 4 dead by time
                
                particles.reference <- particles.reference[ state == 2 & cell.rafted != -9 ]
                particles.reference[ , travel.time := ( ( 1 + t.finish-t.start ) / n.hours.per.day ) ]
                
                setkey(particles.reference, cell )

                n.particles.alldays <- length(simulation.parameters.step$day) * n.new.particles.per.day
                progress.full <- length(unique(particles.reference[,cell]))
                cell.to.process <- unique(particles.reference[,cell])

                # ---------
                
                for( cell.id.ref.f in cell.to.process ) {

                      cat('\014')
                      cat('\n')
                      cat('\n')
                      cat('\n Saving results to SQL')
                      cat('\n')
                      cat('\n')
                      
                      progress.partial <- floor( which( cell.to.process == cell.id.ref.f ) / progress.full )
                      cat(paste0(" | ",paste(rep("-", progress.partial ),collapse="")," | ",progress.partial,"%"))
                  
                      # Generate temporary data table and Remove cell from reference
                      
                      connectivity.temp <- particles.reference[ cell == cell.id.ref.f , ]
                      particles.reference <- particles.reference[ cell != cell.id.ref.f , ]
                      
                      for( cell.id.ref.t in unique(connectivity.temp[ , cell.rafted ]) ) {

                              connectivity.pairs.to.sql <- data.frame(  Pair.from = cell.id.ref.f,
                                                                        Pair.to = cell.id.ref.t,
                                                                        Time.mean = mean(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                        Time.min = min(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                        Time.max = max(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                        Time.sd = sd(connectivity.temp[ cell.rafted == cell.id.ref.t,]$travel.time),
                                                                        Probability = nrow(connectivity.temp[ cell.rafted == cell.id.ref.t,]) / n.particles.alldays,
                                                                        Year = simulation.year )
                              
                              # Save pairs to SQL
                              
                              sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/SQLite/",results.files,".reference.particles.sql"))
                              dbWriteTable(sql, "Connectivity", connectivity.pairs.to.sql , overwrite=FALSE, append=TRUE)
                              dbDisconnect(sql)
                              
                      } }

                # ------------------------------------------------------------------
                # Remove heavy objects from memory
                
                rm(connectivity.pairs.to.sql) ; rm(particles.reference) ; rm(particles.location.temp) ; rm(particles.location.temp.x) ; rm(particles.location.temp.y) ; rm(particles.location.x) ; rm(particles.location.y) ; rm(raw.data.u) ; rm(raw.data.v) ; rm(speed.u.rk) ; rm(speed.v.rk) ; rm(raw.data.u.t) ; rm(raw.data.v.t)
                
}

##  ---------------------------------------------------------------------------------------------------------------------------------
##  ---------------------------------------------------------------------------------------------------------------------------------
## End of Code