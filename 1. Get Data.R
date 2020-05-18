## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

source("0. Project Config.R")
source("Dependences.R")

## -----------------------------------------

files <- c("http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1" ,
           "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_90.9" ,
           "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.0" , 
           "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.1" , 
           "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2"
)

if( buffer == TRUE ) {  
  min.lon <- min.lon - buffer.val
  max.lon <- max.lon + buffer.val
  min.lat <- min.lat - buffer.val
  max.lat <- max.lat + buffer.val
}

## -----------------------------------------

fullDates <- data.frame(year=numeric(0), month=numeric(0), day=numeric(0), file=character(0),stringsAsFactors = FALSE)

for(y in from.year:to.year) {
  
  for(m in sort(months.all)) {
    
    fullDates.m <- 1:numberOfDays(as.Date(paste0(y,"-",ifelse(nchar(m) == 1 , paste0(0,m),m),"-01", "%Y-%m-%d")))
    fullDates <- rbind(fullDates,data.frame(year=y, month=m, day=fullDates.m, file="NA", stringsAsFactors = FALSE))
    
  }  
}

# ---------------

for(file in files) {
  
  time.corrected <- NULL
  nc <- NULL
  increment <- 0
  while(is.null(nc)) {
    tryCatch(  nc <- nc_open( file , verbose=FALSE ) , error = function(e) { Sys.sleep(0.25) })
    increment <- increment + 1
    if(increment == 100) { stop("100 attempts failed!")}
  }
  
  tryCatch( time.corrected <- as.Date(ncvar_get( nc, "time")/24, origin = "2000-01-01") , error = function(e) {  })
  if(is.null(time.corrected)) { tryCatch( time.corrected <- as.Date(ncvar_get( nc, "MT"), origin = "1900-12-31 00:00:00") , error = function(e) { }) }
  
  for( i in 1:nrow(fullDates)) {
    
    year <- fullDates[i,]$year
    month <- fullDates[i,]$month
    day <- fullDates[i,]$day
    
    if( as.Date(paste0(year, "-" ,ifelse(nchar(month) == 1 , paste0(0,month),month), "-" , ifelse(nchar(day) == 1 , paste0(0,day),day) )) %in% time.corrected ) { fullDates[i,"file"] <- file }
    
  }
}

# -----------------------------------------------------------------

for (y in 2014){ # unique(fullDates$year)
  
  fullDates.y <- fullDates[fullDates$year == y,]
  time.window <- apply(fullDates.y[,2:3] , 1 , function(x) { paste0(y,"-",ifelse(nchar(x[1]) == 1 , paste0(0,x[1]),x[1]),"-",ifelse(nchar(x[2]) == 1 , paste0(0,x[2]),x[2])) }  )

  for( k in 1:nrow(fullDates.y)) {
    
    day <- fullDates.y[k,]$day
    month <- fullDates.y[k,]$month
    file <- fullDates.y[k,]$file
    
    if( file == "NA" ) { 
      day <- fullDates.y[k-1,]$day
      month <- fullDates.y[k-1,]$month
      file <- fullDates.y[k-1,]$file
    }
    if( file == "NA" ) { 
      day <- fullDates.y[k+1,]$day
      month <- fullDates.y[k+1,]$month
      file <- fullDates.y[k+1,]$file
    }
    
    time.corrected <- NULL
    nc <- NULL
    increment <- 0
    while(is.null(nc)) {
      tryCatch(  nc <- nc_open( file , verbose=FALSE ) , error = function(e) { Sys.sleep(5) })
      increment <- increment + 1
      if(increment == 50) { stop("50 attempts failed!")}
    }
    
    tryCatch( time.corrected <- as.Date(ncvar_get( nc, "time")/24, origin = "2000-01-01") , error = function(e) {  })
    if(is.null(time.corrected)) { tryCatch( time.corrected <- as.Date(ncvar_get( nc, "MT"), origin = "1900-12-31 00:00:00") , error = function(e) { }) }

    t <- which(time.corrected == as.Date(paste0(y, "-" ,ifelse(nchar(month) == 1 , paste0(0,month),month), "-" , ifelse(nchar(day) == 1 , paste0(0,day),day) )) )
    
    if( length(t) == 0) { stop("Error ::  102")}
    
    depths <- ncvar_get( nc, "depth")
    depth <-  which(depths %in% depth.range)
    Latitude <- ncvar_get( nc, "lat")
    Longitude <- ncvar_get( nc, "lon")
    Longitude[Longitude > 180] <- Longitude[Longitude > 180] - 360
    
    array.region.lon <- which(Longitude >= min.lon & Longitude <= max.lon)
    array.region.lat <- which(Latitude >= min.lat & Latitude <= max.lat)
    Longitude.array <- Longitude[array.region.lon]
    Latitude.array <- Latitude[array.region.lat]
    
    i.min <- min(array.region.lon)
    j.min <- min(array.region.lat)
    i.max <- max(array.region.lon)
    j.max <- max(array.region.lat)
    
    u <- array(data = NA, dim = c( i.max-i.min+1 , j.max-j.min+1 , length(depth) ))
    v <- array(data = NA, dim = c( i.max-i.min+1 , j.max-j.min+1 , length(depth) ))

    for( d in 1:length(depth) ) {

      values.to.place.u <- NULL
      values.to.place.v <- NULL
      
      while(is.null(values.to.place.u)) {
        
        try( values.to.place.u <- ncvar_get( nc, "water_u", start=c(i.min,j.min,depth[d],t), count=c((i.max-i.min+1),(j.max-j.min+1),1,1)))
        
      }

      if( max(values.to.place.u,na.rm=T) > 100 | min(values.to.place.u,na.rm=T) < -100) { stop("Strange values") }
      
      # values.to.place.u[values.to.place.u > 100 ] <- NA
      # values.to.place.u[values.to.place.u < -100 ] <- NA
      
      while(is.null(values.to.place.v)) {
        
        try( values.to.place.v <- ncvar_get( nc, "water_v", start=c(i.min,j.min,depth[d],t), count=c((i.max-i.min+1),(j.max-j.min+1),1,1)))
        
      }
      
      if( max(values.to.place.v,na.rm=T) > 100 | min(values.to.place.v,na.rm=T) < -100) { stop("Strange values") }
      
      # values.to.place.v[values.to.place.v > 100 ] <- NA
      # values.to.place.v[values.to.place.v < -100 ] <- NA
      
      u[,,d] <- values.to.place.u
      v[,,d] <- values.to.place.v
      
    }
    
    if(dim(u)[3] == 1) { u <- u[,,1] }
    if(dim(v)[3] == 1) { v <- v[,,1] }
    
    # For 2d data
    
    if( k == 1 & is.na(dim(u)[3]) )  {  
      
      data.u <- array(data = NA, dim = c(dim(u)[1],dim(u)[2],nrow(fullDates.y)) )
      data.v <- array(data = NA, dim = c(dim(v)[1],dim(v)[2],nrow(fullDates.y)) ) 
      
    }
    
    
    if( is.na(dim(u)[3]) )  {  
      
      data.u[,,k] <- u
      data.v[,,k] <- v 
      
      }     
    
    # For 3d data
    
    if(  k == 1 & !is.na(dim(u)[3]) )  { 
      
            data.u <- array(data = NA, dim = c(dim(u)[1],dim(u)[2],dim(u)[3],nrow(fullDates.y)), dimnames = c("Lon","Lat","Depth","Time"))
            data.v <- array(data = NA, dim = c(dim(v)[1],dim(v)[2],dim(u)[3],nrow(fullDates.y)), dimnames = c("Lon","Lat","Depth","Time")) 
            
            }
    
    if( !is.na(dim(u)[3]) )  {  data.u[,,,k] <- u
                                data.v[,,,k] <- v }
    
    cat('\014')
    cat('\n')
    cat('\n')
    cat('\n Progress:', paste0(y, "-" ,ifelse(nchar(month) == 1 , paste0(0,month),month), "-" , ifelse(nchar(day) == 1 , paste0(0,day),day) ) )
    
  }
  
  #--------------------------------
  # Export to netcdf
  
  dimX <- ncdim_def( "X", "unit", 1:dim(data.u)[1] )
  dimY <- ncdim_def( "Y", "unit", 1:dim(data.u)[2] )
  dimT <- ncdim_def( "Time", "days", 1:length(time.window) )
  dimD <- ncdim_def( "Depths", "meters", depth )
  
  mv <- -999 # missing value to use
  
  var1d <- ncvar_def( "Longitude", "Degrees", dimX, mv)
  var2d <- ncvar_def( "Latitude", "Degrees", dimY, mv)
  var3d <- ncvar_def( "Date", "days", dimT, mv)
  var4d <- ncvar_def( "Depth", "meters", dimD, mv)
  
  # 2D
  
  if(final.dimensions == 2) {
    var4d <- ncvar_def( "UComponent", "m", list(dimX,dimY,dimT), mv, prec="double")
    var5d <- ncvar_def( "VComponent", "m", list(dimX,dimY,dimT), mv,  prec="double")
    nc.file <- nc_create( paste(project.folder,"/Data","/currents_2d_",project.name,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d) )
    ncvar_put( nc.file, var1d, Longitude.array )  
    ncvar_put( nc.file, var2d, Latitude.array )   
    ncvar_put( nc.file, var3d, as.numeric(as.Date(time.window)) )
    ncvar_put( nc.file, var4d, data.u )  
    ncvar_put( nc.file, var5d, data.v )
    nc_close(nc.file)
  }
  
  if(final.dimensions == 3) {
    var5d <- ncvar_def( "UComponent", "m", list(dimX,dimY,dimD,dimT), mv, prec="double")
    var6d <- ncvar_def( "VComponent", "m", list(dimX,dimY,dimD,dimT), mv,  prec="double")
    nc.file <- nc_create( paste(project.folder,"/Data","/currents_3d_",project.name,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d,var6d) )
    ncvar_put( nc.file, var1d, Longitude.array )  
    ncvar_put( nc.file, var2d, Latitude.array )   
    ncvar_put( nc.file, var3d, as.numeric(as.Date(time.window)) )
    ncvar_put( nc.file, var4d, as.numeric(depths[depth]) )
    ncvar_put( nc.file, var5d, data.u )  
    ncvar_put( nc.file, var6d, data.v )
    nc_close(nc.file)
  }
  
}
          
# -----------------------------------------------------------------
# -----------------------------------------------------------------
