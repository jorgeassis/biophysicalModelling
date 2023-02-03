## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

for( simulation.step in 1:nrow(simulation.parameters.step) ) {
  
  ## --------------------------------------------------------
  ## Progress
  
  if( simulation.step== 1) { 
    start_time <- Sys.time()
    time.i <- character(nrow(simulation.parameters.step))
  }
  
  time.i[simulation.step] <- as.character(Sys.time())
  simulation.step.previous <- simulation.step - 1
  if(simulation.step.previous == 0) { simulation.step.previous <- 1 }
  
  progress.percent <- round((simulation.step / nrow(simulation.parameters.step)) * 100)
  time.take.step.all <- round(as.numeric(difftime(time.i[simulation.step], time.i[1], units = "mins")))
  time.take.step.min <- round(as.numeric(difftime(time.i[simulation.step], time.i[simulation.step.previous], units = "mins")))
  
  cat('\014')
  cat('\n')
  cat('\n Running biophysical model step #',simulation.step,'| Time taken:',time.take.step.min,"mins. | Total time taken: ",time.take.step.all,"mins.")
  cat('\n',paste0(rep("-",100),collapse = ""))
  cat('\n',paste0(rep("-",progress.percent),collapse = ""),"||",progress.percent,"%")
  cat('\n',paste0(rep("-",100),collapse = ""))
  
  ## ---------------------
  
  ## Log
  
  if( simulation.step == 1 & file.exists(paste0(results.folder,"/","process.log")) ) { file.remove(paste0(results.folder,"/","process.log"))  }
  lineToWrite <- paste0("Log ",Sys.time()," | Running step: #",simulation.step," out of #", nrow(simulation.parameters.step), " | Time taken: ",time.take.step.min," min. | Total time taken: ",time.take.step.all," mins.")
  write(lineToWrite,file=paste0(results.folder,"/","process.log"),append=TRUE)
  
  ## --------------------------------------------------------
  
  simulation.year <- simulation.parameters.step[simulation.step,"year"]
  simulation.month <- simulation.parameters.step[simulation.step,"month"]
  simulation.day <- simulation.parameters.step[simulation.step,"day"]
  simulation.raw.data.file <- simulation.parameters.step[simulation.step,"file"]
  
  # Next day data
  
  simulation.year.n <- simulation.parameters.step[simulation.step + ifelse(simulation.step < nrow(simulation.parameters.step) , 1 , 0) ,"year"]
  simulation.month.n <- simulation.parameters.step[simulation.step + ifelse(simulation.step < nrow(simulation.parameters.step) , 1 , 0) ,"month"]
  simulation.day.n <- simulation.parameters.step[simulation.step + ifelse(simulation.step < nrow(simulation.parameters.step) , 1 , 0) ,"day"]
  simulation.raw.data.file.n <- simulation.parameters.step[simulation.step + ifelse(simulation.step < nrow(simulation.parameters.step) , 1 , 0) ,"file"]
  
  ## --------------------------------------------------------
  
  if( ( as.numeric(simulation.month.n) - as.numeric(simulation.month) > 1 ) & ( as.numeric(simulation.year.n ) == as.numeric(simulation.year) )  ) {
    
    stop("Review :: 991")
    
  }
  
  if( ( as.numeric(simulation.month.n ) - as.numeric(simulation.month) != -11 ) & (as.numeric(simulation.year.n ) > as.numeric(simulation.year) )  ) {
    
    simulation.year.n <- simulation.year
    simulation.month.n <- simulation.month
    simulation.day.n <- simulation.day
    simulation.raw.data.file.n <- simulation.raw.data.file
    
    particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
    particles.reference.moving <- mwhich(particles.reference.bm.all,c(9),list(1), list('eq') )
    particles.reference.bm.all[ particles.reference.moving , 9 ] <- 4
    rm(particles.reference.bm.all)
    gc(reset=TRUE)
    
  }
  
  ## --------------------------------------------------------
  
  ## Test data coherence
  
  simulation.raw.data <- nc_open( simulation.raw.data.file, readunlim=FALSE )
  Longitude <- ncvar_get( simulation.raw.data, "Longitude" )
  Latitude <- ncvar_get( simulation.raw.data, "Latitude" )
  nc_close(simulation.raw.data)
  
  simulation.raw.data <- nc_open( simulation.raw.data.file.n, readunlim=FALSE )
  Longitude.n <- ncvar_get( simulation.raw.data, "Longitude" )
  Latitude.n <- ncvar_get( simulation.raw.data, "Latitude" )
  nc_close(simulation.raw.data)
  
  if( ! all.equal(Longitude,Longitude.n) | ! all.equal(Latitude,Latitude.n) ) { stop(" Error :: Code 104") }
  
  ## --------------------------------------------------------
  
  nc.file <- nc_open( simulation.raw.data.file, readunlim=FALSE )
  time.start.day <- ncvar_get( nc.file, "Date" )
  time.start.day <- as.Date(time.start.day, origin = "1970-01-01") 
  nc_close(nc.file)
  
  nc.file <- nc_open( simulation.raw.data.file.n, readunlim=FALSE )
  time.next.day <- ncvar_get( nc.file, "Date" )
  time.next.day <- as.Date(time.next.day, origin = "1970-01-01") 
  nc_close(nc.file)
  
  if( ! as.numeric(format(time.start.day, "%Y"))[1] %in% from.year:to.year ) { stop(" Error :: Code 101") }
  if( ! as.numeric(format(time.next.day, "%Y"))[1] %in% from.year:to.year ) { stop(" Error :: Code 102") }
  
  rd.start.day <- which( format(time.start.day, "%Y") == simulation.year & format(time.start.day, "%m") == simulation.month & format(time.start.day, "%d") == simulation.day )
  rd.next.day <- which( format(time.next.day, "%Y") == simulation.year.n & format(time.next.day, "%m") == simulation.month.n & format(time.next.day, "%d") == simulation.day.n )
  
  if( length(rd.start.day) == 0 | length(rd.next.day) == 0 ) { stop(" Error :: Code 105") }
  
  ## --------------------------------------------------------------------------------
  ## --------------------------------------------------------------------------------
  
  if( simulation.step == 1 ) {
    
    options(warn=-1)
    xy <- source.sink.xy[source.sink.xy$source == 1 , c("x","y","cells.id") ]
    xy.sp <- st_as_sf(xy, coords = c("x", "y"), crs = 4326)
    grd <- sf::st_make_grid(xy.sp, n=c(round(sqrt(number.cores) + 0.5),round(sqrt(number.cores) - 0.5))  ) # cellsize = ( round( abs( diff(c(max.lon,min.lon)) / number.cores ) + 1))
    grd <- as_Spatial(grd)
    grd$id <- 1:length(grd)
    grdCentroid <- st_centroid(st_as_sf(grd))
    grdCentroid <- do.call(rbind, st_geometry(grdCentroid))
    options(warn=0)
    
    cat("\nParallel structure defined with", length(grd),"dimentions\n")
    
    # plot(as_Spatial(worldmap), main="")
    # plot(grd,add=T)
    # points(source.sink.xy[,c("x","y")], col="red")

  }
  
  parallelChunk <- data.frame(cells.id=xy[,"cells.id"],chunk=get.knnx( grdCentroid , xy[, c("x","y") ], k=1 , algorithm="kd_tree" )$nn.index)
  if( sum(is.na(parallelChunk$chunk)) > 0 | FALSE %in% (source.cells.id %in% parallelChunk$cells.id ) ) { stop("Error :: Code 1999") }
  
  # points(source.sink.xy[source.sink.xy$cells.id %in% parallelChunk[parallelChunk$chunk == 4, 1],2:3], col="red")
  # points(source.sink.xy[source.sink.xy$cells.id %in% parallelChunk[parallelChunk$chunk == 1, 1],2:3], col="green")
  
  ## -------
  ## -------
  
  particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
  
  # Get particle reference of moving cells
  
  particles.reference.bm.moving.id <- mwhich(particles.reference.bm,9,list(1),"eq")
  particles.reference.bm.moving <- particles.reference.bm[particles.reference.bm.moving.id,]
  particles.reference.bm.moving <- data.table(particles.reference.bm.moving)
  names(particles.reference.bm.moving) <- particles.reference.names
  setkey(particles.reference.bm.moving, id )
  
  # Get ids of (unborne) available cells
  
  particles.reference.bm.unborne.id <- mwhich(particles.reference.bm,9,list(0),"eq")
  particles.reference.bm.unborne <- particles.reference.bm[particles.reference.bm.unborne.id,]
  particles.reference.bm.unborne <- data.table(particles.reference.bm.unborne)
  names(particles.reference.bm.unborne) <- particles.reference.names
  setkey(particles.reference.bm.unborne, id )
  
  particles.reference.bm.unborne <- particles.reference.bm.unborne[ , .SD[1] , by=start.cell][,.(id,start.cell)]
  
  if(n.new.particles.per.day > 1) {
    stop("Review :: 032")
    for(j in 2:n.new.particles.per.day) {
      particles.reference.bm.unborne <- rbindlist(list(particles.reference.bm.unborne, particles.reference.bm.unborne[ , .SD[j] , by=start.cell][,.(id,start.cell)] ))
    }
  }

  ## -------
  ## -------
  
  Cluster <- makeCluster( number.cores )
  registerDoParallel( Cluster )
  
  parallelProcess <- foreach(chunk=sort(unique(parallelChunk$chunk)), .verbose=FALSE, .packages=c("FNN","ncdf4","h3","data.table","zoo","bigmemory","rgeos","dismo","h3js")) %dopar% {
    
    start.cell.i <- parallelChunk[parallelChunk$chunk == chunk,"cells.id"]
    
    if( length(start.cell.i) == 0) { stop("Error :: 917") }
    
    ## -------
    
    particles.reference.dt.video <- data.table()
    particles.reference.moving <- particles.reference.bm.moving[ start.cell %in% start.cell.i, ]
    setkey(particles.reference.moving,id)
    
    ## --------------------------------------------------------
    ## Move particles (per day)
    
    for( h in 1:n.hours.per.day ) {
      
      t.step <- ((as.numeric(simulation.step)-1) * n.hours.per.day) + h
      
      # -----------------------------------------------
      # Release new particles, if that is the case
      
      if( norm.time[t.step,release.particles] ) {   
        
        if( h != 1) { stop("Review :: 033") }
        
        particles.reference.unborne <- particles.reference.bm.unborne[ start.cell %in% start.cell.i, ]
        particles.reference.new <- particles.reference.template[ start.cell %in% start.cell.i, ]
        particles.reference.new$id <- particles.reference.unborne$id
        
        particles.reference.new[ , "state" ] <- 1
        particles.reference.new[ , "t.start" ] <- t.step
        particles.reference.new[ , "start.year" ] <- as.numeric(simulation.year)
        particles.reference.new[ , "start.month" ] <- as.numeric(simulation.month)
        particles.reference.new[ , "start.day" ] <- as.numeric(simulation.day)
      
        particles.reference.moving <-  rbindlist(list(particles.reference.moving, particles.reference.new))

      }
      
      # --------
      
      if( nrow(particles.reference.moving) == 0 ) { stop("Error :: Code 982") }
      
      # -----------------------------------------------
      
      if( h == 1 ) {
        
        dim.i.subseter <- which( Longitude >= min(particles.reference.moving$pos.lon) - 5 &
                                 Longitude <= max(particles.reference.moving$pos.lon) + 5 )
        
        dim.j.subseter <- which( Latitude >= min(particles.reference.moving$pos.lat) - 2.5 &
                                 Latitude <= max(particles.reference.moving$pos.lat) + 2.5 )
        
        norm.field <- expand.grid( i=dim.i.subseter , j=dim.j.subseter )
        raw.data.coords <- cbind( apply( norm.field , 1 , function (x) Longitude[x[1]] ) , apply( norm.field , 1 , function (x) Latitude[x[2]] ) )
        colnames(raw.data.coords) <- c("Lon","Lat")
        
        dim.i.subseter <- c(min(dim.i.subseter),max(dim.i.subseter))
        dim.j.subseter <- c(min(dim.j.subseter),max(dim.j.subseter))
        
        nc.file <- nc_open( simulation.raw.data.file, readunlim=FALSE )
        velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(dim.i.subseter[1],dim.j.subseter[1],rd.start.day) , count=c(dim.i.subseter[2] - dim.i.subseter[1] + 1,dim.j.subseter[2] - dim.j.subseter[1] + 1,1) )
        velocity.field[velocity.field == -999] <- NA
        start.day.raw.data.u <- reshape2::melt(velocity.field)[,"value"]

        velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(dim.i.subseter[1],dim.j.subseter[1],rd.start.day) , count=c(dim.i.subseter[2] - dim.i.subseter[1] + 1,dim.j.subseter[2] - dim.j.subseter[1] + 1,1) )
        velocity.field[velocity.field == -999] <- NA
        start.day.raw.data.v <-  reshape2::melt(velocity.field)[,"value"]
        nc_close(nc.file)
        
        nc.file <- nc_open( simulation.raw.data.file, readunlim=FALSE )
        velocity.field <- ncvar_get(  nc.file , "UComponent", start=c(dim.i.subseter[1],dim.j.subseter[1],rd.next.day) , count=c(dim.i.subseter[2] - dim.i.subseter[1] + 1,dim.j.subseter[2] - dim.j.subseter[1] + 1,1) )
        velocity.field[velocity.field == -999] <- NA
        next.day.raw.data.u <- reshape2::melt(velocity.field)[,"value"]
        
        velocity.field <- ncvar_get(  nc.file , "VComponent", start=c(dim.i.subseter[1],dim.j.subseter[1],rd.next.day) , count=c(dim.i.subseter[2] - dim.i.subseter[1] + 1,dim.j.subseter[2] - dim.j.subseter[1] + 1,1) )
        velocity.field[velocity.field == -999] <- NA
        next.day.raw.data.v <-  reshape2::melt(velocity.field)[,"value"]
        nc_close(nc.file)
        
        raw.data.u <- cbind(start.day.raw.data.u,next.day.raw.data.u)
        raw.data.v <- cbind(start.day.raw.data.v,next.day.raw.data.v)
        colnames(raw.data.u) <- c("u.start","u.next")
        colnames(raw.data.v) <- c("v.start","v.next")
        
        raw.data.coords <- raw.data.coords[which(!is.na(raw.data.u[,1])),]
        raw.data.u <- raw.data.u[which(!is.na(raw.data.u[,1])),]
        raw.data.v <- raw.data.v[which(!is.na(raw.data.v[,1])),]
        
        if( sum(is.na(raw.data.u)) > 0 | sum(is.na(raw.data.v)) > 0 ) { 
          
          raw.data.u[,1] <- na.approx(raw.data.u[,1], rule = 2)
          raw.data.u[,2] <- na.approx(raw.data.u[,2], rule = 2) 
          raw.data.v[,1] <- na.approx(raw.data.v[,1], rule = 2) 
          raw.data.v[,2] <- na.approx(raw.data.v[,2], rule = 2) 
          
        }
        
        if( nrow(raw.data.coords) != nrow(raw.data.u) | nrow(raw.data.u) != nrow(raw.data.v) ) { stop(" Error :: Code 111") }
        
        rm(start.day.raw.data.u); rm(next.day.raw.data.u); rm(start.day.raw.data.v); rm(next.day.raw.data.v)
        
      }
      
      # -----------------------------------------------
      
      setkey(particles.reference.moving,id)
      moving.particles.id <- as.vector(unlist(particles.reference.moving[, id]))
      points.to.interp <- particles.reference.moving[, .(pos.lon,pos.lat,id)]
      
      # --------
      
      # Interpolate speed
      
      idwPower <- 2
      idwN <- 4
      
      idw.nearest.r <- get.knnx( raw.data.coords , points.to.interp[,.(pos.lon,pos.lat)], k=idwN , algorithm="kd_tree" )
      idw.nearest.i <- idw.nearest.r$nn.index
      idw.nearest.d <- idw.nearest.r$nn.dist
      idw.nearest.d <- idw.nearest.d * 1e9
      
      speed.u.sec <- apply( data.frame(from=raw.data.u[,1][idw.nearest.i],to=raw.data.u[,2][idw.nearest.i]) , 1, function(x) { seq( from = x[[1]] , to = x[[2]] , length.out = n.hours.per.day + 1 ) })
      speed.v.sec <- apply( data.frame(from=raw.data.v[,1][idw.nearest.i],to=raw.data.v[,2][idw.nearest.i]) , 1, function(x) { seq( from = x[[1]] , to = x[[2]] , length.out = n.hours.per.day + 1 ) })
      
      speedData.u <- apply( matrix(speed.u.sec[h,],ncol=idwN) / idw.nearest.d^idwPower,1,sum) / apply(1 / idw.nearest.d^idwPower,1,sum)
      speedData.v <- apply( matrix(speed.v.sec[h,],ncol=idwN) / idw.nearest.d^idwPower,1,sum) / apply(1 / idw.nearest.d^idwPower,1,sum)
      
      if( max(speedData.u) > 100 | min(speedData.u) < -100 | max(speedData.v) > 100 | min(speedData.v) < -100 ) { 
        
        speedData.u[speedData.u > 100] <- NA
        speedData.u[speedData.u < -100] <- NA
        speedData.v[speedData.v > 100] <- NA
        speedData.v[speedData.v < -100] <- NA
        
        speedData.u <- na.approx(speedData.u, rule = 2)
        speedData.v <- na.approx(speedData.v, rule = 2)
        
      }
      
      mov.eastward <- speedData.u * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
      mov.northward <- speedData.v * 60 * 60 * ( 24 / n.hours.per.day ) # Was as m/s
      
      # Assign temporary positions
      
      particles.position.lat <- unlist(points.to.interp[,.(pos.lat)] )
      particles.position.lat[particles.position.lat > 89.99] <- 89.99
      
      dLon <- mov.eastward / ( 6378137 * cos( pi * (  particles.position.lat  / 180) ) )
      names(dLon) <- NULL
      dLat <- mov.northward / 6378137
      
      dLon <- points.to.interp[,pos.lon] + dLon * 180/pi 
      dLat <- points.to.interp[,pos.lat] + dLat * 180/pi
      
      # Condition for regional transition
      
      dLon[ dLon > 180 ] <- dLon[dLon > 180 ] - 360
      dLon[ dLon < -180 ] <- dLon[dLon < -180 ] + 360
      dLat[ dLat > 90 ] <- 90 + (90 - dLat[ dLat > 90 ])
      
      if( max(dLon) > 180 | min(dLon) < -180 ) { stop(" Error :: Code 421") }
      if( max(dLat) > 90 | min(dLat) < -90 ) { stop(" Error :: Code 431") }
      
      # -----------------------------------------------
      
      setkey(particles.reference.moving,id)
      particles.reference.moving[,"pos.lon"] <- dLon
      particles.reference.moving[,"pos.lat"] <- dLat
      
      # -----------------------------------------------
      # -----------------------------------------------
      # Test over sea hexagons
      
      setkey(particles.reference.moving,id)
      points.to.test <- particles.reference.moving[, .(pos.lat,pos.lon)]

      hexagons.address.test <- sapply(1:nrow(points.to.test), function(x) { unique( geo_to_h3(c(points.to.test[x,pos.lat],points.to.test[x,pos.lon]), sim.resolution  )  ) })

      particles.on.sea <- which( hexagons.address.test %in% hexagons.address.ocean )
      particles.on.sea.id <- as.vector(moving.particles.id[particles.on.sea])
      setkey(particles.reference.moving,id)
      particles.reference.moving[ id %in% particles.on.sea.id , at.sea := 1  ]
      
      # -----------------------------------------------
      # -----------------------------------------------
      # Particles over land hexagons
      
      particles.on.land <- which( hexagons.address.test %in% hexagons.address.land )
      particles.on.land.id <- as.vector(moving.particles.id[particles.on.land])

      if( length(particles.on.land) > 0 ) {    
        
        setkey(particles.reference.moving,id)
        setkey(points.to.interp,id)
        
        cells.started <- as.vector(unlist(particles.reference.moving[id %in% particles.on.land.id, "start.cell"]))
        who.at.land.t.start <- as.vector(unlist(particles.reference.moving[id %in% particles.on.land.id, "t.start"]))
        
        points.on.land.t <- particles.reference.moving[ id %in% particles.on.land.id , .(pos.lon,pos.lat) ]
        points.on.land.t.minus <- points.to.interp[ id %in% particles.on.land.id , .(pos.lon,pos.lat) ]
        
        cells.rafted <- get.knnx( source.sink.xy[,c("x","y")], points.on.land.t.minus, k=1 , algorithm="kd_tree" )$nn.index
        cells.rafted <- source.sink.xy[cells.rafted,"cells.id"]
        
        displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]} )
        
        # True Rafters [Distinct cell]
        
        true.rafters.id <- particles.on.land.id[ displacement != 0 ]
        true.rafters.cell <- cells.rafted[ displacement != 0 ]
        
        setkey(particles.reference.moving,id)
        particles.reference.moving[ id %in% true.rafters.id , state := 2 ]
        particles.reference.moving[ id %in% true.rafters.id , cell.rafted := as.numeric(true.rafters.cell) ]
        particles.reference.moving[ id %in% true.rafters.id , t.finish := as.numeric(t.step) ]
        
        # False Rafters [same cell of origin]
        
        non.rafters.id <- particles.on.land.id[ displacement == 0 ]
        non.rafters.cell <- cells.started[ displacement == 0 ]
        
        particles.reference.moving[ id %in% non.rafters.id , state := 2 ]
        particles.reference.moving[ id %in% non.rafters.id , cell.rafted := as.numeric(non.rafters.cell) ]
        particles.reference.moving[ id %in% non.rafters.id , t.finish := as.numeric(t.step) ]
        
        # Allow for relocation when rafting at t.step == t.start
        
        non.rafters.t <- who.at.land.t.start[ displacement == 0 ] == t.step
        
        if( TRUE %in% non.rafters.t & allow.back.to.origin ) {
          
          stop("Revise code :: 001")
          
          particles.reference.moving[ id %in% non.rafters.id[non.rafters.t] , pos.lon := points.on.land.t.minus[ id %in% non.rafters.id[non.rafters.t] , pos.lon] ]
          particles.reference.moving[ id %in% non.rafters.id[non.rafters.t] , pos.lat := points.on.land.t.minus[ id %in% non.rafters.id[non.rafters.t] , pos.lat] ]
          
          particles.reference.moving[ id %in% non.rafters.id[non.rafters.t] , state := 1 ]
          particles.reference.moving[ id %in% non.rafters.id[non.rafters.t] , cell.rafted := as.numeric(0) ]
          particles.reference.moving[ id %in% non.rafters.id[non.rafters.t] , t.finish := as.numeric(0) ]
          
        }
        
      }
      
      # -----------------------------------------------
      # -----------------------------------------------
      # Particles over shore hexagons [ includes additional shapefile hexagons, if the case]
      
      particles.on.shore <- which( hexagons.address.test %in% hexagons.address.sourcesink )
      particles.on.shore.id <- as.vector(moving.particles.id[particles.on.shore])

      if( length(particles.on.shore) > 0 ) {    
        
        setkey(particles.reference.moving,id)
        setkey(points.to.interp,id)
        
        cells.started <- as.vector(unlist(particles.reference.moving[id %in% particles.on.shore.id, "start.cell"]))
        who.at.shore.t.start <- as.vector(unlist(particles.reference.moving[id %in% particles.on.shore.id, "t.start"]))
        
        points.on.shore.t <- particles.reference.moving[ id %in% particles.on.shore.id , .(pos.lon,pos.lat) ]
        points.on.shore.t.minus <- points.to.interp[ id %in% particles.on.shore.id , .(pos.lon,pos.lat) ]
        
        cells.rafted <- get.knnx( source.sink.xy[,c("x","y")], points.on.shore.t.minus, k=1 , algorithm="kd_tree" )$nn.index
        cells.rafted <- source.sink.xy[cells.rafted,"cells.id"]
        
        displacement <- apply( cbind( cells.started, cells.rafted) , 1 , function(x) { x[2] - x[1]} )
        
        # True Rafters [Distinct cell]
        
        true.rafters.id <- particles.on.shore.id[ displacement != 0 ]
        true.rafters.cell <- cells.rafted[ displacement != 0 ]
        
        particles.reference.moving[ id %in% true.rafters.id , state := 2 ]
        particles.reference.moving[ id %in% true.rafters.id , cell.rafted := as.numeric(true.rafters.cell) ]
        particles.reference.moving[ id %in% true.rafters.id , t.finish := as.numeric(t.step) ]
        
        # True Rafters by Retention [same cell of origin, but off in the ocean]
        
        non.rafters.id <- particles.on.shore.id[ displacement == 0 ]
        non.rafters.id.at.sea.once <- particles.reference.moving[ id %in% non.rafters.id & at.sea == 1, id ]
        non.rafters.cell <- as.vector(unlist(particles.reference.moving[id %in% non.rafters.id.at.sea.once, "start.cell"]))
        
        if( length(non.rafters.id.at.sea.once) > 0 ) {    
          
          particles.reference.moving[ id %in% non.rafters.id.at.sea.once , state := 2 ]
          particles.reference.moving[ id %in% non.rafters.id.at.sea.once, cell.rafted := as.numeric(non.rafters.cell) ]
          particles.reference.moving[ id %in% non.rafters.id.at.sea.once , t.finish := as.numeric(t.step) ]
          
        }
        
      }
      
      # -----------------------------------------------
      # Out of space (study region), if TRUE, place particles on hold
      
      out.of.space <- dLon > max.lon | dLon < min.lon | dLat > max.lat | dLat < min.lat
      out.of.space.ids.1 <- moving.particles.id[out.of.space]
      out.of.space.ids.2 <- moving.particles.id[which(! moving.particles.id %in% particles.on.land.id & ! moving.particles.id %in% particles.on.shore.id & ! moving.particles.id %in% particles.on.sea.id)]
      out.of.space.ids <- unique(c(out.of.space.ids.1,out.of.space.ids.2))
      
      setkey(particles.reference.moving,id)
      particles.reference.moving[ id %in% out.of.space.ids , state := 3  ]
      
      # -----------------------------------------------
      # -----------------------------------------------
      # End of day, Kill by Longevity
      
      if ( longevity ) {   
        
        max.duration.id <- particles.reference.moving[ ( t.step - t.start ) > ( n.hours.per.day * particle.max.duration ) , id ]
        setkey(particles.reference.moving, id )
        particles.reference.moving[ id %in% max.duration.id , state := 4 ]
        
      }
      
      ## ---------------------------------------------------------------
      ## Save positions to video DT (if condition matched) 
      
      if( ! is.null(movie.sites.xy) & movie.year == as.numeric(simulation.year) ) {
        
        setkey(particles.reference.moving,id)
        particles.video.id.moving <- particles.reference.moving[ state == 1 & id %in% particles.video.id,id]
        particles.video.id.moving.condition <- length(particles.video.id.moving) > 0
        
        if( particles.video.id.moving.condition ) { 
          
          t.step.movie <- ((as.numeric(which( which(as.numeric(simulation.parameters.step[,3]) == movie.year) == simulation.step))-1) * n.hours.per.day) + h
          setkey(particles.reference.moving,id)
          particles.reference.dt.video <- rbindlist(list(particles.reference.dt.video,data.table(t.step.movie=t.step.movie,year=simulation.year,month=simulation.month,day=simulation.day,particle.id=particles.video.id.moving,pos.lon=particles.reference.moving.dt[ id %in% particles.video.id.moving,pos.lon],pos.lat=particles.reference.moving.dt[ id %in% particles.video.id.moving,pos.lat])))
          
        }
      }

      ## ---------------------------------------------------------------
      
    }
    
    ## ---------------------------------------------------------------
    
    # inject particles movement
    
    inject.dt <- particles.reference.moving[ state != 0 , ]
    
    if( nrow(inject.dt) > 0 ) {
      
      particles.reference.bm.all <- attach.big.matrix(particles.reference.bm.desc)
      
      particles.reference.bm.all[ inject.dt$id , 3 ] <- as.numeric(unlist(inject.dt[  , "start.year"] ))
      particles.reference.bm.all[ inject.dt$id , 4 ] <- as.numeric(unlist(inject.dt[  , "start.month"] ))
      particles.reference.bm.all[ inject.dt$id , 5 ] <- as.numeric(unlist(inject.dt[  , "start.day"] ))
      
      particles.reference.bm.all[ inject.dt$id , 6 ] <- as.numeric(unlist(inject.dt[  , "pos.lon"] ))
      particles.reference.bm.all[ inject.dt$id , 7 ] <- as.numeric(unlist(inject.dt[  , "pos.lat"] ))
      particles.reference.bm.all[ inject.dt$id , 9 ] <- as.numeric(unlist(inject.dt[  , "state"] ))
      particles.reference.bm.all[ inject.dt$id , 10 ] <- as.numeric(unlist(inject.dt[  , "t.start"] ))
      particles.reference.bm.all[ inject.dt$id , 11 ] <- as.numeric(unlist(inject.dt[  , "t.finish"] ))
      particles.reference.bm.all[ inject.dt$id , 12 ] <- as.numeric(unlist(inject.dt[  , "cell.rafted"] ))
      particles.reference.bm.all[ inject.dt$id , 13 ] <- as.numeric(unlist(inject.dt[  , "at.sea"] ))
      
      # inject travel time
      
      inject.dt <- particles.reference.moving[ state == 2 , ]
      
      if( nrow(inject.dt) > 0 ) {
        
        particles.reference.bm.all[ inject.dt$id , 14 ] <- ( 1 + particles.reference.bm.all[ inject.dt$id , 11 ] - particles.reference.bm.all[ inject.dt$id , 10 ] ) / n.hours.per.day

      }
      
    }

    
    ## ---------------------------------------------------------------
    
    return( list(reference.dt.video=as.data.frame(particles.reference.dt.video) ) )
    
  }
  
  stopCluster(Cluster); rm(Cluster)
  closeAllConnections(); gc(reset=TRUE)
  
  ## ---------------------------------------------------
  ## Inject video particles
  
  for( i in 1:length(parallelProcess)) {
    
    reference.dt.video.i <- parallelProcess[[i]]$reference.dt.video
    
    if( nrow(reference.dt.video.i) > 0) {
      
      if( simulation.step != 1 | i != 1 ) {
        colnames(reference.dt.video.i) <- NULL
        write.table(reference.dt.video.i, file = paste0(results.folder,"/particles.video.location.csv"), append = TRUE, sep = ";", row.names = FALSE, col.names = FALSE)
      }

      if( simulation.step == 1 & i == 1 ) {
        write.table(reference.dt.video.i, file = paste0(results.folder,"/particles.video.location.csv"), append = FALSE, sep = ";", row.names = FALSE, col.names = TRUE)
        
      }
    }
  }
  
  # read.csv(paste0(results.folder,"/particles.video.location.csv"), sep=";")
  
  if( simulation.step == nrow(simulation.parameters.step)) { 
    end_time <- Sys.time()
    end_time - start_time
  }
  
}

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##